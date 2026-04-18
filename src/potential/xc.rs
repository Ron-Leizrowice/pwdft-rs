//! Exchange-correlation functionals.
//!
//! LDA: Perdew-Zunger parametrization of the Ceperley-Alder correlation energy,
//! plus Slater exchange. Computed in real space from ρ(r).
//!
//! References:
//! - Exchange: Slater, Phys. Rev. 81, 385 (1951)
//! - Correlation: Perdew & Zunger, Phys. Rev. B 23, 5048 (1981)
//! - Ceperley & Alder, Phys. Rev. Lett. 45, 566 (1980)

use std::f64::consts::PI;

use rayon::prelude::*;

/// Minimum grid size at which rayon parallelism beats sequential execution
/// for `lda_xc_grid` / `lda_xc_spin_grid` on Apple M2.
///
/// Empirical calibration (`cargo bench --bench scf_benchmarks -- xc_grid`,
/// Apple M2 8-core):
///
/// | n       | sequential | parallel | net           |
/// |---------|------------|----------|---------------|
/// | 4 096   | 43 µs      | 113 µs   | +161% (worse) |
/// | 32 768  | 407 µs     | 195 µs   | −52%          |
/// | 262 144 | 3 035 µs   | 459 µs   | −85%          |
///
/// Rayon's fork/join/unzip costs ~70 µs per region on this hardware, which
/// dominates for n=4096 (43 µs sequential work) but is easily amortized by
/// n=32768.  We set the threshold at 16 384 — the safe side of the
/// crossover.  Grid sizes between 16³=4 096 and 32³=32 768 are the most
/// sensitive region; typical production FFT grids are 24³–48³ (≥13 824
/// points), so most real SCF calls take the parallel path.
const XC_PARALLEL_THRESHOLD: usize = 16_384;

/// Result of evaluating the XC functional at a single density point.
pub struct XcPoint {
    /// Exchange-correlation energy density ε_xc (eV per electron).
    pub exc: f64,
    /// Exchange-correlation potential V_xc = d(ρ·ε_xc)/dρ (eV).
    pub vxc: f64,
}

/// Evaluate the LDA exchange-correlation functional at a single density value.
///
/// `rho` is the electron density in e/ų. Must be non-negative.
///
/// Returns ε_xc (energy per electron, eV) and V_xc (potential, eV).
pub fn lda_xc(rho: f64) -> XcPoint {
    if rho < crate::consts::RHO_FLOOR {
        return XcPoint {
            exc: 0.0,
            vxc: 0.0,
        };
    }

    let (ex, vx) = slater_exchange(rho);
    let (ec, vc) = perdew_zunger_correlation(rho);

    XcPoint {
        exc: ex + ec,
        vxc: vx + vc,
    }
}

/// Compute ε_xc(r) and V_xc(r) on a real-space grid.
///
/// `rho_r`: electron density on real-space grid (e/ų).
///
/// Returns (exc_r, vxc_r): energy density and potential on the grid (eV).
///
/// Each point is an independent evaluation, so this is embarrassingly
/// parallel. We fall back to the sequential path below
/// `XC_PARALLEL_THRESHOLD` because rayon's per-region dispatch overhead
/// dominates on very small grids.
pub fn lda_xc_grid(rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    if rho_r.len() < XC_PARALLEL_THRESHOLD {
        let mut exc = Vec::with_capacity(rho_r.len());
        let mut vxc = Vec::with_capacity(rho_r.len());
        for &rho in rho_r {
            let xc = lda_xc(rho);
            exc.push(xc.exc);
            vxc.push(xc.vxc);
        }
        return (exc, vxc);
    }

    rho_r
        .par_iter()
        .map(|&rho| {
            let xc = lda_xc(rho);
            (xc.exc, xc.vxc)
        })
        .unzip()
}

/// Compute the total XC energy: E_xc = Ω/N_grid · Σ_r ρ(r) ε_xc(r)
///
/// Returns energy in eV.
pub fn lda_xc_energy(rho_r: &[f64], exc_r: &[f64], omega: f64) -> f64 {
    let n_grid = rho_r.len() as f64;
    let dvol = omega / n_grid; // volume per grid point

    rho_r
        .iter()
        .zip(exc_r.iter())
        .map(|(&rho, &exc)| rho * exc * dvol)
        .sum()
}

/// Slater exchange: ε_x = -(3/4)(3/π)^{1/3} ρ^{1/3}
///
/// V_x = (4/3) ε_x
///
/// Returns (ε_x, V_x) in Hartree. We convert to eV at the end.
fn slater_exchange(rho: f64) -> (f64, f64) {
    // ε_x in Hartree: -(3/4)(3ρ/π)^{1/3}
    // In eV: multiply by crate::consts::HA_TO_EV = 27.2114
    

    // rho [e/ų] → rho [e/Bohr³] = rho × Bohr_to_Å³ = rho × 0.529177³
    let bohr3 = crate::consts::BOHR3_TO_ANG3;
    let rho_bohr = rho * bohr3;

    // ε_x = -(3/4)(3ρ/π)^{1/3} in Hartree
    let cbrt = (3.0 * rho_bohr / PI).cbrt();
    let ex_ha = -0.75 * cbrt;
    // V_x = d(ρ·ε_x)/dρ = (4/3) ε_x
    let vx_ha = (4.0 / 3.0) * ex_ha;

    (ex_ha * crate::consts::HA_TO_EV, vx_ha * crate::consts::HA_TO_EV)
}

/// Perdew-Zunger parametrization of the Ceperley-Alder correlation energy.
///
/// Two regimes based on Wigner-Seitz radius r_s:
/// - r_s ≥ 1: ε_c = γ / (1 + β₁√r_s + β₂r_s)
/// - r_s < 1: ε_c = A ln(r_s) + B + C r_s ln(r_s) + D r_s
///
/// Returns (ε_c, V_c) in eV.
fn perdew_zunger_correlation(rho: f64) -> (f64, f64) {
    
    let bohr3 = crate::consts::BOHR3_TO_ANG3;
    let rho_bohr = rho * bohr3;

    // Wigner-Seitz radius in Bohr
    let rs = (3.0 / (4.0 * PI * rho_bohr)).cbrt();

    let (ec_ha, vc_ha);

    if rs >= 1.0 {
        // PZ parameters for r_s ≥ 1 (unpolarized)
        let gamma = -0.1423;
        let beta1 = 1.0529;
        let beta2 = 0.3334;

        let sqrt_rs = rs.sqrt();
        let denom = 1.0 + beta1 * sqrt_rs + beta2 * rs;

        ec_ha = gamma / denom;

        // V_c = ε_c - (r_s / 3) dε_c/dr_s
        // dε_c/dr_s = -γ (β₁/(2√r_s) + β₂) / denom²
        let d_ec = -gamma * (beta1 / (2.0 * sqrt_rs) + beta2) / (denom * denom);
        vc_ha = (rs / 3.0).mul_add(-d_ec, ec_ha);
    } else {
        // PZ parameters for r_s < 1 (unpolarized)
        let a: f64 = 0.0311;
        let b: f64 = -0.048;
        let c: f64 = 0.0020;
        let d: f64 = -0.0116;

        let ln_rs = rs.ln();
        ec_ha = d.mul_add(rs, (c * rs).mul_add(ln_rs, a.mul_add(ln_rs, b)));

        // V_c = ε_c - (r_s / 3) dε_c/dr_s
        // dε_c/dr_s = a/r_s + c(ln(r_s) + 1) + d
        let d_ec = c.mul_add(ln_rs + 1.0, a / rs) + d;
        vc_ha = (rs / 3.0).mul_add(-d_ec, ec_ha);
    }

    (ec_ha * crate::consts::HA_TO_EV, vc_ha * crate::consts::HA_TO_EV)
}

// ---------------------------------------------------------------------------
// Spin-polarized LSDA (collinear)
// ---------------------------------------------------------------------------

/// Result of spin-polarized XC evaluation at a single point.
pub struct XcSpinPoint {
    /// Exchange-correlation energy density ε_xc (eV per electron).
    pub exc: f64,
    /// XC potential for spin up (eV).
    pub vxc_up: f64,
    /// XC potential for spin down (eV).
    pub vxc_down: f64,
}

/// Spin-polarized LDA XC at a single point.
///
/// `rho_up`, `rho_down` in e/ų. Returns energy density and potentials in eV.
/// Evaluate spin-polarized LDA XC at a single point.
///
/// Exchange: ε_x^σ = -(3/4)(6ρ_σ/π)^{1/3} (fully polarized gas per channel).
/// Correlation: interpolated between unpolarized (ζ=0) and fully polarized (ζ=1)
/// using the von Barth-Hedin interpolation function:
///   f(ζ) = [(1+ζ)^{4/3} + (1-ζ)^{4/3} - 2] / [2^{4/3} - 2]
///   ε_c(r_s, ζ) = ε_c^unpol + f(ζ) [ε_c^pol - ε_c^unpol]
///
/// `rho_up`, `rho_down` in e/ų. Returns (ε_xc, V_xc↑, V_xc↓) in eV.
pub fn lda_xc_spin(rho_up: f64, rho_down: f64) -> XcSpinPoint {
    let rho = rho_up + rho_down;
    if rho < crate::consts::RHO_FLOOR {
        return XcSpinPoint { exc: 0.0, vxc_up: 0.0, vxc_down: 0.0 };
    }

    let (ex, vx_up, vx_down) = slater_exchange_spin(rho_up, rho_down);
    let (ec, vc_up, vc_down) = pz_correlation_spin(rho_up, rho_down);

    XcSpinPoint {
        exc: ex + ec,
        vxc_up: vx_up + vc_up,
        vxc_down: vx_down + vc_down,
    }
}

/// Spin-polarized XC on a real-space grid.
///
/// Returns (exc_r, vxc_up_r, vxc_down_r) in eV.
///
/// Parallelized for grids at or above `XC_PARALLEL_THRESHOLD` points.
/// Rayon's `unzip` handles only 2-tuples, so we unzip into an
/// ((exc, vxc_up), vxc_down) shape and then flatten.
pub fn lda_xc_spin_grid(
    rho_up_r: &[f64],
    rho_down_r: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = rho_up_r.len();
    debug_assert_eq!(n, rho_down_r.len(), "spin channels must share grid size");

    if n < XC_PARALLEL_THRESHOLD {
        let mut exc = Vec::with_capacity(n);
        let mut vxc_up = Vec::with_capacity(n);
        let mut vxc_down = Vec::with_capacity(n);
        for (&ru, &rd) in rho_up_r.iter().zip(rho_down_r.iter()) {
            let xc = lda_xc_spin(ru, rd);
            exc.push(xc.exc);
            vxc_up.push(xc.vxc_up);
            vxc_down.push(xc.vxc_down);
        }
        return (exc, vxc_up, vxc_down);
    }

    let ((exc, vxc_up), vxc_down): ((Vec<f64>, Vec<f64>), Vec<f64>) = rho_up_r
        .par_iter()
        .zip(rho_down_r.par_iter())
        .map(|(&ru, &rd)| {
            let xc = lda_xc_spin(ru, rd);
            ((xc.exc, xc.vxc_up), xc.vxc_down)
        })
        .unzip();

    (exc, vxc_up, vxc_down)
}

/// Spin-polarized Slater exchange.
///
/// ε_x = (1/2)[(1+ζ)·ε_x(2ρ_up) + (1-ζ)·ε_x(2ρ_down)]
/// where ε_x(ρ) = -(3/4)(3ρ/π)^{1/3} is the unpolarized exchange per electron.
///
/// V_x_σ = (4/3)·ε_x(2ρ_σ)·2^{1/3}  [derivative of the spin-scaled exchange]
///
/// Returns (ε_x, V_x_up, V_x_down) in eV.
fn slater_exchange_spin(rho_up: f64, rho_down: f64) -> (f64, f64, f64) {
    
    let bohr3 = crate::consts::BOHR3_TO_ANG3;

    let rho = rho_up + rho_down;
    if rho < crate::consts::RHO_FLOOR {
        return (0.0, 0.0, 0.0);
    }

    let rho_up_bohr = rho_up * bohr3;
    let rho_down_bohr = rho_down * bohr3;

    // Exchange energy per electron for each spin channel (fully polarized formula)
    // ε_x(ρ_σ) for a single spin channel = -(3/4)(6ρ_σ/π)^{1/3}
    // This is the exchange of a fully-polarized gas with density ρ_σ
    let ex_up_ha = if rho_up_bohr > crate::consts::RHO_FLOOR {
        -0.75 * (6.0 * rho_up_bohr / PI).cbrt()
    } else {
        0.0
    };
    let ex_down_ha = if rho_down_bohr > crate::consts::RHO_FLOOR {
        -0.75 * (6.0 * rho_down_bohr / PI).cbrt()
    } else {
        0.0
    };

    // Total exchange energy density: ε_x = (ρ_up·ε_x_up + ρ_down·ε_x_down) / ρ
    let rho_bohr = rho * bohr3;
    let ex_ha = rho_up_bohr.mul_add(ex_up_ha, rho_down_bohr * ex_down_ha) / rho_bohr;

    // Potentials: V_x_σ = d(ρ·ε_x)/dρ_σ = (4/3)·ε_x(ρ_σ)
    let vx_up_ha = (4.0 / 3.0) * ex_up_ha;
    let vx_down_ha = (4.0 / 3.0) * ex_down_ha;

    (ex_ha * crate::consts::HA_TO_EV, vx_up_ha * crate::consts::HA_TO_EV, vx_down_ha * crate::consts::HA_TO_EV)
}

/// Spin-polarized Perdew-Zunger correlation via spin interpolation.
///
/// ε_c(rs, ζ) = ε_c^unpol(rs) + f(ζ)·[ε_c^pol(rs) - ε_c^unpol(rs)]
/// where f(ζ) = [(1+ζ)^{4/3} + (1-ζ)^{4/3} - 2] / [2^{4/3} - 2]
///
/// Returns (ε_c, V_c_up, V_c_down) in eV.
fn pz_correlation_spin(rho_up: f64, rho_down: f64) -> (f64, f64, f64) {
    
    let bohr3 = crate::consts::BOHR3_TO_ANG3;

    let rho = rho_up + rho_down;
    if rho < crate::consts::RHO_FLOOR {
        return (0.0, 0.0, 0.0);
    }

    let rho_bohr = rho * bohr3;
    let rs = (3.0 / (4.0 * PI * rho_bohr)).cbrt();
    let zeta = ((rho_up - rho_down) / rho).clamp(-1.0, 1.0);

    // Unpolarized correlation
    let (ec_unpol, vc_unpol) = pz_correlation_rs(rs, false);
    // Fully polarized correlation
    let (ec_pol, vc_pol) = pz_correlation_rs(rs, true);

    // Spin interpolation function f(ζ) and its derivative
    let f_denom = 2.0_f64.cbrt().mul_add(2.0, -2.0); // 2^{4/3} - 2

    let op = (1.0 + zeta).max(0.0);
    let om = (1.0 - zeta).max(0.0);
    let cbrt_op = op.cbrt();
    let cbrt_om = om.cbrt();
    let fz = (cbrt_op * op + cbrt_om * om - 2.0) / f_denom; // (1+ζ)^{4/3} + (1-ζ)^{4/3} - 2

    // df/dζ = (4/3) [(1+ζ)^{1/3} - (1-ζ)^{1/3}] / (2^{4/3} - 2)
    let dfz = (4.0 / 3.0) * (cbrt_op - cbrt_om) / f_denom;

    // Energy density
    let ec_ha = ec_unpol + fz * (ec_pol - ec_unpol);

    // Potentials (chain rule via rs and ζ)
    // V_c_σ = ε_c - (rs/3)·dε_c/drs + (±1 - ζ)·dε_c/dζ
    // where dε_c/drs = dε_c^u/drs + fz·(dε_c^p/drs - dε_c^u/drs)
    // and dε_c/dζ = dfz·(ε_c^p - ε_c^u)
    let vc_rs_ha = vc_unpol + fz * (vc_pol - vc_unpol); // this is ε_c - (rs/3)·dε_c/drs
    let dec_dzeta = dfz * (ec_pol - ec_unpol);

    let vc_up_ha = (1.0 - zeta).mul_add(dec_dzeta, vc_rs_ha);
    let vc_down_ha = (1.0 + zeta).mul_add(-dec_dzeta, vc_rs_ha);

    (ec_ha * crate::consts::HA_TO_EV, vc_up_ha * crate::consts::HA_TO_EV, vc_down_ha * crate::consts::HA_TO_EV)
}

/// PZ correlation at given rs for unpolarized (polarized=false) or
/// fully polarized (polarized=true) electron gas.
///
/// Returns (ε_c, V_c) where V_c = ε_c - (rs/3)·dε_c/drs, in Hartree.
fn pz_correlation_rs(rs: f64, polarized: bool) -> (f64, f64) {
    let (ec_ha, vc_ha);

    if rs >= 1.0 {
        let (gamma, beta1, beta2): (f64, f64, f64) = if polarized {
            (-0.0843, 1.3981, 0.2611) // PZ fully polarized parameters
        } else {
            (-0.1423, 1.0529, 0.3334) // PZ unpolarized parameters
        };

        let sqrt_rs = rs.sqrt();
        let denom = beta2.mul_add(rs, beta1.mul_add(sqrt_rs, 1.0));
        ec_ha = gamma / denom;
        let d_ec = -gamma * (beta1 / (2.0 * sqrt_rs) + beta2) / (denom * denom);
        vc_ha = (rs / 3.0).mul_add(-d_ec, ec_ha);
    } else {
        let (a, b, c, d): (f64, f64, f64, f64) = if polarized {
            (0.01555, -0.0269, 0.0007, -0.0048) // PZ fully polarized
        } else {
            (0.0311, -0.048, 0.0020, -0.0116) // PZ unpolarized
        };

        let ln_rs = rs.ln();
        ec_ha = d.mul_add(rs, (c * rs).mul_add(ln_rs, a.mul_add(ln_rs, b)));
        let d_ec = c.mul_add(ln_rs + 1.0, a / rs) + d;
        vc_ha = (rs / 3.0).mul_add(-d_ec, ec_ha);
    }

    (ec_ha, vc_ha)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::float_cmp)] // Zero-density XC returns exact 0.0
    fn test_lda_zero_density() {
        let xc = lda_xc(0.0);
        assert_eq!(xc.exc, 0.0);
        assert_eq!(xc.vxc, 0.0);
    }

    #[test]
    fn test_lda_exchange_sign() {
        // Exchange energy should be negative
        let xc = lda_xc(0.05);
        assert!(xc.exc < 0.0, "ε_xc should be negative: {}", xc.exc);
        assert!(xc.vxc < 0.0, "V_xc should be negative: {}", xc.vxc);
    }

    #[test]
    fn test_lda_monotonic_in_density() {
        // |V_xc| should increase with density
        let v1 = lda_xc(0.01).vxc.abs();
        let v2 = lda_xc(0.1).vxc.abs();
        assert!(
            v2 > v1,
            "|V_xc| should increase with ρ: {v1} vs {v2}"
        );
    }

    #[test]
    fn test_lda_known_values() {
        // At ρ = 0.044 e/ų (≈ 0.00652 e/Bohr³, r_s ≈ 3.33 Bohr):
        // ε_x ≈ -3.76 eV, ε_c ≈ -0.96 eV, ε_xc ≈ -4.72 eV
        let rho = 0.044;
        let xc = lda_xc(rho);
        assert!(
            xc.exc > -6.0 && xc.exc < -3.0,
            "ε_xc at ρ=0.044: expected ~ -4.7 eV, got {}",
            xc.exc
        );
        // Average Si valence density ~0.20 e/ų (8 electrons / 40 ų):
        // r_s ≈ 2.0 Bohr, ε_xc ≈ -7.5 eV
        let rho_avg = 0.20;
        let xc_avg = lda_xc(rho_avg);
        assert!(
            xc_avg.exc > -10.0 && xc_avg.exc < -5.0,
            "ε_xc at avg Si density: expected ~ -7.5 eV, got {}",
            xc_avg.exc
        );
    }

    #[test]
    fn test_lsda_unpolarized_limit() {
        // lda_xc_spin(ρ/2, ρ/2) must equal lda_xc(ρ) — the ζ=0 limit
        for &rho in &[0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0] {
            let unpol = lda_xc(rho);
            let spin = lda_xc_spin(rho / 2.0, rho / 2.0);

            assert!(
                (unpol.exc - spin.exc).abs() < 1e-10,
                "ε_xc mismatch at ρ={rho}: unpol={}, spin={}",
                unpol.exc, spin.exc
            );
            // V_xc should be equal for both spins and match unpolarized
            assert!(
                (unpol.vxc - spin.vxc_up).abs() < 1e-10,
                "V_xc_up mismatch at ρ={rho}: unpol={}, spin={}",
                unpol.vxc, spin.vxc_up
            );
            assert!(
                (spin.vxc_up - spin.vxc_down).abs() < 1e-10,
                "V_xc_up != V_xc_down at ρ={rho}: up={}, down={}",
                spin.vxc_up, spin.vxc_down
            );
        }
    }

    #[test]
    fn test_lsda_fully_polarized() {
        // Fully polarized: all spin up (ζ=1)
        let rho = 0.1;
        let xc = lda_xc_spin(rho, 0.0);
        assert!(xc.exc < 0.0, "ε_xc should be negative");
        // Exchange should be more negative for polarized than unpolarized
        // (Pauli exclusion reduces exchange hole)
        let unpol = lda_xc(rho);
        assert!(
            xc.exc < unpol.exc,
            "Polarized ε_xc ({}) should be more negative than unpolarized ({})",
            xc.exc, unpol.exc
        );
    }

    #[test]
    fn test_lsda_symmetry() {
        // Swapping up/down should swap potentials but keep exc the same
        let rho_up = 0.15;
        let rho_down = 0.05;
        let xc1 = lda_xc_spin(rho_up, rho_down);
        let xc2 = lda_xc_spin(rho_down, rho_up);

        assert!(
            (xc1.exc - xc2.exc).abs() < 1e-12,
            "ε_xc not symmetric: {} vs {}", xc1.exc, xc2.exc
        );
        assert!(
            (xc1.vxc_up - xc2.vxc_down).abs() < 1e-12,
            "V_xc swap failed: up1={}, down2={}", xc1.vxc_up, xc2.vxc_down
        );
        assert!(
            (xc1.vxc_down - xc2.vxc_up).abs() < 1e-12,
            "V_xc swap failed: down1={}, up2={}", xc1.vxc_down, xc2.vxc_up
        );
    }

    #[test]
    fn test_lda_both_regimes() {
        // Low density: r_s > 1
        let xc_low = lda_xc(0.001);
        assert!(xc_low.exc < 0.0);

        // High density: r_s < 1 (very high density, rare in practice)
        // r_s < 1 Bohr requires ρ > 3/(4π) e/Bohr³ ≈ 0.239 e/Bohr³ ≈ 1.61 e/ų
        let xc_high = lda_xc(2.0);
        assert!(xc_high.exc < 0.0);
    }
}
