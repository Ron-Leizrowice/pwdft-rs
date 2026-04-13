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
    if rho < 1e-30 {
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
pub fn lda_xc_grid(rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut exc = Vec::with_capacity(rho_r.len());
    let mut vxc = Vec::with_capacity(rho_r.len());

    for &rho in rho_r {
        let xc = lda_xc(rho);
        exc.push(xc.exc);
        vxc.push(xc.vxc);
    }

    (exc, vxc)
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
    // In eV: multiply by HA_TO_EV = 27.2114
    const HA_TO_EV: f64 = 27.211386245988;

    // rho [e/ų] → rho [e/Bohr³] = rho × Bohr_to_Å³ = rho × 0.529177³
    let bohr3 = 0.529177210903_f64.powi(3);
    let rho_bohr = rho * bohr3;

    // ε_x = -(3/4)(3ρ/π)^{1/3} in Hartree
    let cbrt = (3.0 * rho_bohr / PI).powf(1.0 / 3.0);
    let ex_ha = -0.75 * cbrt;
    // V_x = d(ρ·ε_x)/dρ = (4/3) ε_x
    let vx_ha = (4.0 / 3.0) * ex_ha;

    (ex_ha * HA_TO_EV, vx_ha * HA_TO_EV)
}

/// Perdew-Zunger parametrization of the Ceperley-Alder correlation energy.
///
/// Two regimes based on Wigner-Seitz radius r_s:
/// - r_s ≥ 1: ε_c = γ / (1 + β₁√r_s + β₂r_s)
/// - r_s < 1: ε_c = A ln(r_s) + B + C r_s ln(r_s) + D r_s
///
/// Returns (ε_c, V_c) in eV.
fn perdew_zunger_correlation(rho: f64) -> (f64, f64) {
    const HA_TO_EV: f64 = 27.211386245988;
    let bohr3 = 0.529177210903_f64.powi(3);
    let rho_bohr = rho * bohr3;

    // Wigner-Seitz radius in Bohr
    let rs = (3.0 / (4.0 * PI * rho_bohr)).powf(1.0 / 3.0);

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
        vc_ha = ec_ha - rs / 3.0 * d_ec;
    } else {
        // PZ parameters for r_s < 1 (unpolarized)
        let a = 0.0311;
        let b = -0.048;
        let c = 0.0020;
        let d = -0.0116;

        let ln_rs = rs.ln();
        ec_ha = a * ln_rs + b + c * rs * ln_rs + d * rs;

        // V_c = ε_c - (r_s / 3) dε_c/dr_s
        // dε_c/dr_s = a/r_s + c(ln(r_s) + 1) + d
        let d_ec = a / rs + c * (ln_rs + 1.0) + d;
        vc_ha = ec_ha - rs / 3.0 * d_ec;
    }

    (ec_ha * HA_TO_EV, vc_ha * HA_TO_EV)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
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
