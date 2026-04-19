//! Exchange-correlation functionals.
//!
//! LDA: Perdew-Zunger parametrization of the Ceperley-Alder correlation energy,
//! plus Slater exchange. Computed in real space from ρ(r).
//!
//! # Dispatch
//!
//! The SCF driver does not call [`lda_xc_grid`] / [`lda_xc_spin_grid`]
//! directly. It holds an [`XcEvaluator`] value — a *data* enum whose
//! variants name the functional — and dispatches via a single `match`
//! inside [`XcEvaluator::eval`] / [`XcEvaluator::eval_spin`]. Variants
//! carry only plain data (no closures, no `Box<dyn Fn>`, no trait objects);
//! this keeps each variant independently implementable and leaves the
//! Hamiltonian assembly free to see ψ, which hybrid functionals require.
//!
//! [`XcEvaluator::Pz`] (LDA) is fully implemented; it calls
//! [`lda_xc_grid`] / [`lda_xc_spin_grid`] verbatim.
//! [`XcEvaluator::Pbe`] is partially implemented: the exchange half
//! (`pbe_exchange` in this module) is ported from QE's `pbex` with
//! `iflag=1` and is unit-tested. [`XcEvaluator::eval`] /
//! [`XcEvaluator::eval_spin`] still return
//! [`crate::error::PwdftError::NotImplemented`] with
//! `what = "pbe_correlation"` because PW92-based PBE correlation has
//! not yet landed; Phase C wires it in.
//!
//! References:
//! - Exchange: Slater, Phys. Rev. 81, 385 (1951)
//! - Correlation: Perdew & Zunger, Phys. Rev. B 23, 5048 (1981)
//! - Ceperley & Alder, Phys. Rev. Lett. 45, 566 (1980)

use std::f64::consts::PI;

use rayon::prelude::*;

use crate::{
    error::{PwdftError, Result},
    settings::XcFunctional,
};

/// Result of evaluating the (non-spin) LDA exchange-correlation
/// functional at a single real-space density point.
///
/// The two scalars implement the pointwise decomposition
///
/// ```text
///     E_xc[ρ] = ∫ ρ(r) · ε_xc(ρ(r)) d³r
///     V_xc(r) = δE_xc / δρ(r) = d[ρ · ε_xc(ρ)] / dρ|_{ρ(r)}
///             = ε_xc(ρ) + ρ · dε_xc/dρ
/// ```
///
/// with `ε_xc = ε_x + ε_c` built from Slater exchange
/// `ε_x(ρ) = −(3/4)(3ρ/π)^{1/3}` and the Perdew-Zunger parametrization
/// of the Ceperley-Alder correlation energy (see the module header for
/// the full piecewise formulas and citations).
///
/// Units: `exc` in eV per electron; `vxc` in eV. Both are scalar
/// real-space quantities (LDA is local, so the functional derivative is
/// itself a pointwise function of ρ).
///
/// Reference: Slater, *Phys. Rev.* **81**, 385 (1951) for exchange;
/// Perdew & Zunger, *Phys. Rev. B* **23**, 5048 (1981) Eq. (C1) for
/// the correlation parametrization.
pub struct XcPoint {
    /// Exchange-correlation energy per electron `ε_xc(ρ)` in eV.
    pub exc: f64,
    /// Exchange-correlation potential
    /// `V_xc = d[ρ·ε_xc(ρ)]/dρ = ε_xc + ρ·dε_xc/dρ` in eV.
    pub vxc: f64,
}

/// Evaluate the non-spin LDA exchange-correlation functional at a single
/// real-space density `rho`.
///
/// Returns
///
/// ```text
///     ε_xc(ρ) = ε_x(ρ) + ε_c(ρ)
///     V_xc(ρ) = (4/3) · ε_x(ρ) + [ε_c(ρ) − (r_s/3) · dε_c/dr_s]
/// ```
///
/// where
/// - `ε_x(ρ) = −(3/4)(3ρ/π)^{1/3}` is the Slater exchange energy per
///   electron (Hartree atomic units, converted to eV internally);
/// - `ε_c(ρ)` is the Perdew-Zunger parametrization of the Ceperley-Alder
///   correlation energy, with the two `r_s` regimes written out in
///   [`XcPoint`]'s module header (see references);
/// - `r_s = (3/(4π ρ))^{1/3}` is the Wigner-Seitz radius in Bohr and
///   `ρ` is converted from e/Å³ to e/Bohr³ before entering the formulas
///   to keep the dimensionless `(3/π)^{1/3}` constant correct.
///
/// Inputs:
/// - `rho`: electron density in e/Å³. Must satisfy `rho >= 0`; values
///   below [`crate::consts::RHO_FLOOR`] short-circuit to zero to avoid a
///   cube-root singularity and a spurious −∞ log in the `r_s < 1`
///   branch.
///
/// Returns [`XcPoint`] with `exc` in eV per electron and `vxc` in eV.
///
/// Reference: Slater, *Phys. Rev.* **81**, 385 (1951);
/// Perdew & Zunger, *Phys. Rev. B* **23**, 5048 (1981) Eq. (C1)
/// and surrounding discussion;
/// Ceperley & Alder, *Phys. Rev. Lett.* **45**, 566 (1980) for the
/// quantum-Monte-Carlo correlation data fit by PZ.
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

/// Evaluate the non-spin LDA functional at every point of a real-space
/// density grid.
///
/// This is the pointwise lift of [`lda_xc`] to an array:
///
/// ```text
///     exc_r[i] = ε_xc(ρ(r_i))          (eV per electron)
///     vxc_r[i] = d[ρ · ε_xc] / dρ      (eV)
/// ```
///
/// so that the total LDA energy assembled by [`lda_xc_energy`] is
///
/// ```text
///     E_xc = (Ω / N_grid) · Σ_i ρ(r_i) · ε_xc(ρ(r_i)).
/// ```
///
/// Inputs:
/// - `rho_r`: electron density on the real-space grid (e/Å³); each
///   entry must be non-negative (values below
///   [`crate::consts::RHO_FLOOR`] are short-circuited).
///
/// Returns `(exc_r, vxc_r)`, each a fresh `Vec<f64>` of length
/// `rho_r.len()`, in eV. The grid indexing matches `rho_r` 1:1.
///
/// Parallelization: pointwise-independent, so the loop is parallelized
/// unconditionally via rayon — consistent with every other grid kernel
/// in the engine.
///
/// Reference: as for [`lda_xc`] — Perdew & Zunger, *Phys. Rev. B*
/// **23**, 5048 (1981).
pub fn lda_xc_grid(rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    rho_r
        .par_iter()
        .map(|&rho| {
            let xc = lda_xc(rho);
            (xc.exc, xc.vxc)
        })
        .unzip()
}

/// Real-space quadrature of the LDA exchange-correlation energy.
///
/// Returns
///
/// ```text
///     E_xc = ∫ ρ(r) · ε_xc(ρ(r)) d³r
///          ≈ (Ω / N_grid) · Σ_{i=0}^{N_grid−1} ρ(r_i) · ε_xc(ρ(r_i))
/// ```
///
/// where `Ω / N_grid` is the real-space volume element (uniform FFT grid,
/// simple rectangle rule — exact for a band-limited integrand up to the
/// grid's Nyquist, which is the standard PW-DFT convention).
///
/// Inputs:
/// - `rho_r`: electron density on the FFT grid in e/Å³;
/// - `exc_r`: `ε_xc(ρ(r_i))` in eV per electron, typically the first
///   return value of [`lda_xc_grid`]; must have the same length as
///   `rho_r`;
/// - `omega`: cell volume Ω in Å³.
///
/// Returns `E_xc` in eV.
///
/// Reference: any LDA reference implementation, e.g. Martin,
/// *Electronic Structure*, §8.3; the rectangle-rule FFT-grid quadrature
/// is the canonical PW-DFT choice since Ihm, Zunger & Cohen,
/// *J. Phys. C* **12**, 4409 (1979).
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

/// Result of evaluating the collinear spin-polarized LSDA
/// exchange-correlation functional at a single real-space point.
///
/// The pointwise decomposition is
///
/// ```text
///     E_xc[ρ↑, ρ↓] = ∫ (ρ↑ + ρ↓) · ε_xc(ρ↑, ρ↓) d³r
///     V_xc^σ(r)    = δE_xc / δρ_σ(r)
///                  = ε_xc + (ρ↑ + ρ↓) · ∂ε_xc/∂ρ_σ
/// ```
///
/// with `σ ∈ {↑, ↓}`. The functional `ε_xc(ρ↑, ρ↓)` is built from the
/// fully-polarized Slater exchange per channel and a von Barth-Hedin
/// interpolation of PZ correlation between the `ζ = 0` (unpolarized)
/// and `ζ = 1` (fully polarized) gases; see [`lda_xc_spin`] for the
/// formulas.
///
/// Units: `exc` in eV per electron; `vxc_up`, `vxc_down` in eV.
///
/// Reference: von Barth & Hedin, *J. Phys. C* **5**, 1629 (1972) for
/// the collinear LSDA framework and the spin-interpolation function;
/// Perdew & Zunger, *Phys. Rev. B* **23**, 5048 (1981) for the
/// `ζ = 0` and `ζ = 1` endpoint parametrizations.
pub struct XcSpinPoint {
    /// Exchange-correlation energy per electron `ε_xc(ρ↑, ρ↓)` in eV.
    pub exc: f64,
    /// Spin-up channel potential
    /// `V_xc↑ = δE_xc/δρ↑ = ε_xc + (ρ↑ + ρ↓)·∂ε_xc/∂ρ↑` in eV.
    pub vxc_up: f64,
    /// Spin-down channel potential (same formula, with `↑↔↓`), in eV.
    pub vxc_down: f64,
}

/// Evaluate the collinear LSDA exchange-correlation functional at a
/// single real-space point.
///
/// Exchange. By the spin-scaling relation
/// `E_x[ρ↑, ρ↓] = ½(E_x[2ρ↑] + E_x[2ρ↓])` (Oliver & Perdew, *Phys. Rev.
/// A* **20**, 397 (1979)), evaluating the unpolarized Slater exchange
/// per-channel at `2ρ_σ` gives the fully-polarized exchange per
/// electron of that channel:
///
/// ```text
///     ε_x(2ρ_σ) = −(3/4) · (6 ρ_σ / π)^{1/3}            (eV per electron)
///     ε_x(ρ↑, ρ↓) = [ρ↑ · ε_x(2ρ↑) + ρ↓ · ε_x(2ρ↓)] / (ρ↑ + ρ↓)
///     V_x^σ      = δE_x/δρ_σ = (4/3) · ε_x(2ρ_σ)        (eV)
/// ```
///
/// Correlation. Interpolate Perdew-Zunger between the paramagnetic
/// (`ζ = 0`) and ferromagnetic (`ζ = 1`) parametrizations using the
/// von Barth-Hedin form:
///
/// ```text
///     ζ        = (ρ↑ − ρ↓) / (ρ↑ + ρ↓)                  ∈ [−1, 1]
///     f(ζ)     = [(1+ζ)^{4/3} + (1−ζ)^{4/3} − 2] / [2^{4/3} − 2]
///     ε_c(r_s, ζ) = ε_c^unpol(r_s) + f(ζ) · [ε_c^pol(r_s) − ε_c^unpol(r_s)]
/// ```
///
/// Correlation potential (chain rule through `r_s` and `ζ`):
///
/// ```text
///     V_c^σ = ε_c(r_s, ζ) − (r_s/3) · dε_c/dr_s
///                         + ( δ_{σ↑} · (1−ζ) − δ_{σ↓} · (1+ζ) ) · dε_c/dζ
/// ```
///
/// with `dε_c/dζ = f'(ζ) · [ε_c^pol − ε_c^unpol]` and
/// `f'(ζ) = (4/3)[(1+ζ)^{1/3} − (1−ζ)^{1/3}] / (2^{4/3} − 2)`.
///
/// Inputs:
/// - `rho_up`, `rho_down`: spin-channel densities in e/Å³; both must be
///   non-negative. When `ρ↑ + ρ↓ < ` [`crate::consts::RHO_FLOOR`] the
///   routine short-circuits to zero in every component.
///
/// Returns [`XcSpinPoint`] with `exc` in eV per electron and
/// `vxc_up` / `vxc_down` in eV.
///
/// Reference: von Barth & Hedin, *J. Phys. C* **5**, 1629 (1972)
/// Eq. (5.9) for the spin-interpolation function; Oliver & Perdew,
/// *Phys. Rev. A* **20**, 397 (1979) for the spin-scaling relation;
/// Perdew & Zunger, *Phys. Rev. B* **23**, 5048 (1981) §III for the
/// `ζ = 0, 1` endpoint parametrizations.
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

/// Evaluate the collinear LSDA functional at every point of a spin-
/// polarized density grid.
///
/// This is the pointwise lift of [`lda_xc_spin`] to arrays:
///
/// ```text
///     exc_r[i]      = ε_xc(ρ↑(r_i), ρ↓(r_i))          (eV per electron)
///     vxc_up_r[i]   = δE_xc/δρ↑ at r_i                 (eV)
///     vxc_down_r[i] = δE_xc/δρ↓ at r_i                 (eV)
/// ```
///
/// so that the total LSDA energy is
///
/// ```text
///     E_xc = (Ω / N_grid) · Σ_i (ρ↑(r_i) + ρ↓(r_i)) · ε_xc(ρ↑, ρ↓)(r_i).
/// ```
///
/// Inputs:
/// - `rho_up_r`, `rho_down_r`: per-channel densities on the FFT grid in
///   e/Å³; both slices must have the same length (debug-asserted).
///   Each entry must be non-negative (values summing below
///   [`crate::consts::RHO_FLOOR`] short-circuit to zero).
///
/// Returns `(exc_r, vxc_up_r, vxc_down_r)`, three fresh `Vec<f64>`s of
/// length `rho_up_r.len()`, in eV.
///
/// Parallelization: pointwise-independent, so the loop is parallelized
/// unconditionally via rayon — consistent with every other grid kernel
/// in the engine. Rayon's `unzip` only handles 2-tuples, so the
/// implementation unzips to `((exc, vxc_up), vxc_down)` and re-binds
/// the pieces.
///
/// Reference: von Barth & Hedin, *J. Phys. C* **5**, 1629 (1972);
/// Perdew & Zunger, *Phys. Rev. B* **23**, 5048 (1981) §III.
pub fn lda_xc_spin_grid(
    rho_up_r: &[f64],
    rho_down_r: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    debug_assert_eq!(
        rho_up_r.len(),
        rho_down_r.len(),
        "spin channels must share grid size",
    );

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
/// Energy density per electron:
/// ```text
/// ε_x = (1/2)[(1+ζ)·ε_x(2ρ_up) + (1-ζ)·ε_x(2ρ_down)]
/// ```
/// where `ε_x(ρ) = -(3/4)(3ρ/π)^{1/3}` is the unpolarized exchange per
/// electron (evaluated at the scaled argument `2ρ_σ`).
///
/// Potential for spin channel σ (derivation):
/// ```text
/// E_x^σ[ρ_σ] = ρ_σ · ε_x(2ρ_σ)
/// V_x_σ      = δE_x/δρ_σ
///            = ε_x(2ρ_σ) + ρ_σ · 2 · (dε_x/du)|_{u=2ρ_σ}
///            = ε_x(2ρ_σ) + u · (dε_x/du)|_{u=2ρ_σ}
///            = ε_x(2ρ_σ) + ε_x(2ρ_σ)/3           [since ε_x ∝ u^{1/3} ⇒ u·dε_x/du = ε_x/3]
///            = (4/3) · ε_x(2ρ_σ)
/// ```
/// i.e. the spin-channel exchange potential is `(4/3)·ε_x(2ρ_σ)` — no
/// `2^{1/3}` factor; that factor would appear only if the derivation were
/// expressed in terms of `ε_x(ρ_σ)` (unscaled argument) via the chain rule.
///
/// Returns `(ε_x, V_x_up, V_x_down)` in eV.
fn slater_exchange_spin(rho_up: f64, rho_down: f64) -> (f64, f64, f64) {
    
    let bohr3 = crate::consts::BOHR3_TO_ANG3;

    let rho = rho_up + rho_down;
    if rho < crate::consts::RHO_FLOOR {
        return (0.0, 0.0, 0.0);
    }

    let rho_up_bohr = rho_up * bohr3;
    let rho_down_bohr = rho_down * bohr3;

    // Exchange energy per electron for each spin channel (fully polarized formula).
    // We store ε_x(2ρ_σ) = -(3/4)(6ρ_σ/π)^{1/3}, i.e. the unpolarized
    // ε_x(u) = -(3/4)(3u/π)^{1/3} evaluated at the scaled argument u = 2ρ_σ.
    // This is the exchange per electron of a fully-polarized gas of density ρ_σ.
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

    // Potentials: V_x_σ = δ(ρ·ε_x)/δρ_σ = (4/3)·ε_x(2ρ_σ).
    // `ex_up_ha` / `ex_down_ha` already hold ε_x(2ρ_σ) (see derivation in
    // the function docstring); multiply by 4/3 to get V_x_σ.
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

// ---------------------------------------------------------------------------
// GGAP Phase B: PBE exchange (non-spin)
// ---------------------------------------------------------------------------
//
// `pbe_exchange` is the non-spin PBE exchange port. It is unit-tested
// against the PBE 1996 paper's enhancement-factor form directly (see
// `pbe_exchange_*` below). The full grid evaluation lives in
// `XcEvaluator::Pbe::eval`, which currently returns
// `NotImplemented { what: "pbe_correlation" }` — the matching
// correlation half lands in Phase C, at which point the Phase-B
// function will be folded into the per-grid-point loop. The Pbe arm
// already exercises `pbe_exchange` on the first grid point as a
// defensive smoke call so refactors can't silently bitrot the Phase-B
// port.

/// Density floor below which the PBE exchange integrand is clamped to zero
/// (QE `rho_threshold_gga`, in e/Bohr³ after conversion).
///
/// QE uses `rho_threshold_gga = 1.E-6` e/Bohr³ as its default (see
/// `qe-7.5/XClib/dft_setting_params.f90:85`). The check is applied in
/// atomic units after the input-unit conversion so the short-circuit
/// threshold tracks QE's exactly.
const PBE_RHO_THRESHOLD_AU: f64 = 1.0e-6;

/// Gradient-magnitude-squared floor (same convention as QE's
/// `grho_threshold_gga = 1.E-10` in (e/Bohr⁴)²). Below this, the
/// gradient correction is suppressed and PBE exchange reduces exactly
/// to the Slater/LDA limit.
const PBE_GRHO2_THRESHOLD_AU: f64 = 1.0e-10;

/// PBE exchange constant κ ("LO" Lieb-Oxford bound). PBE eq. 14.
const PBE_KAPPA: f64 = 0.804;

/// PBE exchange constant μ = β π² / 3. PBE eq. 12 + §III.
const PBE_MU: f64 = 0.219_514_972_764_517_1;

/// Perdew-Burke-Ernzerhof (1996) exchange, non-spin, with the canonical
/// `iflag = 1` parameter choice (the original PBE, *not* revPBE or
/// PBEsol).
///
/// Given a single grid point with density `ρ` (e/Å³) and gradient
/// magnitude `|∇ρ|` (e/Å⁴), returns the triple
///
/// ```text
///     ε_x(r)  = ρ · ε_x^PBE(ρ, s)           (eV / Å³ — energy density)
///     v1_x(r) = ∂(ρ · ε_x^PBE) / ∂ρ         (eV)
///     v2_x(r) = "QE's v2x convention" — the scalar that contracts with
///               ∇ρ into the semilocal h-vector by h(r) = v2_x(r) · ∇ρ(r).
/// ```
///
/// PBE eq. 14 uses the enhancement-factor form
///
/// ```text
///     ε_x^PBE(ρ, s) = ε_x^LDA(ρ) · F_x(s)
///     F_x(s)        = 1 + κ − κ / (1 + μ s² / κ)
///     s             = |∇ρ| / (2 k_F ρ)
///     k_F           = (3 π² ρ)^{1/3}
///     ε_x^LDA(ρ)    = −(3 / (4π)) · k_F
/// ```
///
/// with constants `κ = 0.804` and `μ = 0.21951492776…` pinned against
/// the `k(1)` / `mu(1)` entries of QE 7.5's `pbex` subroutine.
///
/// The v2 return follows QE's convention exactly: the semilocal h-vector
/// used by the driver is `h(r) = v2_x(r) · ∇ρ(r)` (no extra factor of
/// two), which dimensionally means `v2_x = 2 · ∂(ρ·ε_x) / ∂(|∇ρ|²)` — QE
/// has already absorbed the `2 ∇ρ` from the chain rule into `v2x`.
/// Staying faithful to this convention means the Phase-D driver can
/// read QE's `gcxc_spin` and `v_of_rho.f90:306,343-344` as canonical
/// without any re-scaling.
///
/// For `ρ < 1e-6 e/Bohr³` (QE's `rho_threshold_gga`) or
/// `|∇ρ|² < 1e-10 (e/Bohr⁴)²` (QE's `grho_threshold_gga`) the gradient
/// enhancement is skipped. For very low density the return is
/// `(0, 0, 0)`; for low gradient we return the Slater/LDA limit
/// (F_x = 1, v2 = 0). The thresholds are applied in atomic units to
/// track QE line-for-line; the function is internally pure AU.
///
/// Returns `(eps_x, v1_x, v2_x)` in the units given above.
//
// Source: qe-7.5/XClib/qe_funct_exch_gga.f90::pbex lines 111-331, CASE
// DEFAULT branch (iflag = 1 matches the default because case arms 4-9
// handle revPBE / PBEsol / PBEQ2D / optB88 / optB86b / EV / RPBE /
// W31X). The internal AU→eV/Å conversion uses HA_TO_EV and
// BOHR3_TO_ANG3 from `crate::consts`.
#[inline]
fn pbe_exchange(rho: f64, grad_rho_mag: f64) -> (f64, f64, f64) {
    // Convert inputs to atomic units so the port is line-for-line
    // identical to QE's pbex. QE: rho in e/Bohr³, |∇ρ| in e/Bohr⁴
    // (grho = |∇ρ|²).
    let bohr3 = crate::consts::BOHR3_TO_ANG3;
    let bohr = crate::consts::BOHR_TO_ANG;

    let rho_au = rho * bohr3;
    // 1 e/Å⁴ = (BOHR_TO_ANG)⁴ e/Bohr⁴, so the AU-side magnitude is
    // `grad_rho_mag · BOHR_TO_ANG⁴`. Squared for QE's `grho`.
    let bohr4 = bohr3 * bohr;
    let agrho_au = grad_rho_mag * bohr4;
    let grho_au = agrho_au * agrho_au;

    // Low-density short-circuit (QE `rho_threshold_gga`).
    if rho_au <= PBE_RHO_THRESHOLD_AU {
        return (0.0, 0.0, 0.0);
    }

    // QE-equivalent locals. `c1 = 3/(4π)` (Slater exchange prefactor in
    // Hartree), `c2 = (3π²)^{1/3}` (the k_F prefactor).
    let c1 = 0.75 / PI;
    // (3π²)^(1/3) computed at compile time from the true constant π.
    // QE uses the literal 3.093667726280136 at line 153; we recompute
    // from PI to avoid a literal drift if a future toolchain widens π.
    let c2: f64 = (3.0 * PI * PI).cbrt();
    let c5 = 4.0 / 3.0;

    // QE: kf = c2 * rho^(1/3).
    let kf = c2 * rho_au.cbrt();
    // exunif = ε_x^LDA (Hartree, per electron).
    let exunif = -c1 * kf;

    // Low-gradient short-circuit: drop the gradient enhancement and
    // return bare LDA. This matches QE's `grho_threshold_gga` guard and
    // also protects the `/agrho_au` inside the v2 assembly below.
    if grho_au <= PBE_GRHO2_THRESHOLD_AU {
        // sx = rho · ε_x^LDA (Hartree · e/Bohr³).
        let sx_ha = rho_au * exunif;
        // v1 = d(ρ ε_x^LDA)/dρ = (4/3) ε_x^LDA.
        let v1_ha = c5 * exunif;
        // Convert: Ha · e/Bohr³ → eV/Å³  (multiply by HA_TO_EV / BOHR3_TO_ANG3);
        // Ha → eV for v1.
        let ha = crate::consts::HA_TO_EV;
        return (sx_ha * ha / bohr3, v1_ha * ha, 0.0);
    }

    // QE CASE DEFAULT (iflag = 1), lines 301-324.
    //
    // dsg = 0.5 / kf,  s = |∇ρ| · dsg / ρ,  s² = s·s.
    let dsg = 0.5 / kf;
    let s1 = agrho_au * dsg / rho_au;
    let s2 = s1 * s1;

    // QE's `fx` for the default arm is the *gradient correction* to
    // the enhancement factor, `fx = F_x^{PBE}(s) − 1 = κ − κ/(1 + μs²/κ)`.
    // The task's enhancement factor F_x(s) from PBE eq. 14 is therefore
    // `1 + fx` — identical value, different labelling.
    let f1 = s2 * PBE_MU / PBE_KAPPA;
    let f2 = 1.0 + f1;
    let f3 = PBE_KAPPA / f2;
    let fx = PBE_KAPPA - f3;

    // Full PBE exchange energy density ρ·ε_x^PBE = ρ·ε_x^LDA·F_x(s).
    // QE only stores the gradient-only part (ρ·exunif·fx) in `sx`; we
    // add the LDA piece explicitly so downstream callers get the full
    // ε_x in one shot.
    let sx_full_ha = rho_au * exunif * (1.0 + fx);

    // Derivatives — keep QE's assembly verbatim for the gradient piece
    // and add the LDA ∂/∂ρ externally. QE's v1x equals d(ρ·exunif·fx)/dρ
    // (the gradient-only partial), via the chain-rule decomposition
    //   d(ρ·exunif·fx)/dρ = exunif·fx + (exunif/3)·fx + exunif·(dfx/ds)·(ρ·ds/dρ)
    //                     = sx_s    + dxunif·fx      + exunif·dfx·ds
    // where ρ·ds/dρ = −(4/3)·s = `ds` in QE's locals. The LDA v1
    // contribution `(4/3)·exunif` is added once below.
    let dxunif = exunif / 3.0;
    // dfx/ds for the default branch: dfx1 = (1 + μs²/κ)², dfx = 2μs/dfx1.
    let dfx1 = f2 * f2;
    let dfx = 2.0 * PBE_MU * s1 / dfx1;
    // ρ·ds/dρ = −(4/3)·s. (Used as ds in QE's naming.)
    let ds = -c5 * s1;

    // QE's v1x for the gradient piece equals d(sx_qe)/dρ where
    // sx_qe = rho · exunif · fx. We extend this to the full PBE
    // partial d(rho · exunif · F_x)/dρ by adding the LDA
    // derivative d(rho · exunif)/dρ = (4/3) · exunif, since
    // `d(exunif)/dρ = exunif / (3ρ)` and `rho · d(exunif)/dρ = exunif/3`,
    // so total d(rho · exunif)/dρ = exunif + exunif/3 = (4/3) exunif.
    let sx_s = exunif * fx;
    let v1_grad_ha = sx_s + dxunif * fx + exunif * dfx * ds;
    let v1_lda_ha = c5 * exunif;
    let v1_ha = v1_lda_ha + v1_grad_ha;

    // v2 follows QE's convention (h = v2 · ∇ρ). The LDA exchange is
    // ∇ρ-independent so v2 has no LDA contribution.
    let v2_ha = exunif * dfx * dsg / agrho_au;

    // Unit conversion back to pwdft-rs native units (eV, Å).
    //   ε_x (Ha · e/Bohr³)   → ε_x (eV/Å³):  × HA_TO_EV / BOHR3_TO_ANG3
    //   v1  (Ha)             → v1  (eV):     × HA_TO_EV
    //   v2  (Ha · Bohr⁵ / e) → v2 (eV·Å⁵/e): × HA_TO_EV · BOHR_TO_ANG⁵
    //
    // The target v2 unit is fixed by the driver contract h = v2·∇ρ:
    //   h[eV·Å]  = v2[eV·Å⁵/e] · ∇ρ[e/Å⁴]
    //   ∇·h[eV]  = (1/Å) · h[eV·Å]  →  adds cleanly to v1 giving V_xc [eV].
    let ha = crate::consts::HA_TO_EV;
    let bohr5 = bohr4 * bohr;
    let eps_x = sx_full_ha * ha / bohr3;
    let v1 = v1_ha * ha;
    let v2 = v2_ha * ha * bohr5;

    (eps_x, v1, v2)
}

// ---------------------------------------------------------------------------
// GGAP Phase A: data-enum dispatch
// ---------------------------------------------------------------------------

/// Runtime-selected exchange-correlation functional.
///
/// This is a **data** enum: each variant carries parameters only (no
/// closures, no `Box<dyn Fn>`, no trait objects). Dispatch lives in
/// [`XcEvaluator::eval`] and [`XcEvaluator::eval_spin`] as a single
/// `match`. The shape is deliberately open-coded rather than polymorphic
/// so that hybrid functionals (PBE0, HSE06) can extend this enum cleanly
/// — they need `(ρ, ψ)` access inside the SCF loop, and a closure-based
/// shape would paint the dispatch into a corner.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XcEvaluator {
    /// Perdew-Zunger 81 LDA + Ceperley-Alder correlation + Slater exchange.
    /// The only fully-implemented variant today.
    Pz,
    /// Perdew-Burke-Ernzerhof GGA (1996). The exchange half
    /// (`pbe_exchange`) is implemented and unit-tested in this module.
    /// PW92-based PBE correlation still has to land before the full
    /// functional is usable; until then [`XcEvaluator::eval`] /
    /// [`XcEvaluator::eval_spin`] return
    /// [`PwdftError::NotImplemented`] with
    /// `what = "pbe_correlation"`.
    Pbe,
}

impl XcEvaluator {
    /// Construct an evaluator from the YAML-level [`XcFunctional`] choice.
    ///
    /// This is the single site that concentrates the "is this functional
    /// implemented yet?" check. `Pbe0` and `Hse06` map to
    /// [`PwdftError::NotImplemented`]; `Pbe` constructs successfully but
    /// its `eval` / `eval_spin` methods currently return
    /// [`PwdftError::NotImplemented`] with `what = "pbe_correlation"`
    /// (exchange has landed, correlation is Phase C).
    ///
    /// Returning an error at this construction site (rather than at first
    /// evaluation) lets `scf::run_scf` fail fast before any compute work.
    ///
    /// # Errors
    ///
    /// Returns [`PwdftError::NotImplemented`] for `Pbe0` and `Hse06`;
    /// hybrid functionals are not yet implemented. `Pz` and `Pbe`
    /// always succeed.
    pub fn from_settings(xc: XcFunctional) -> Result<Self> {
        match xc {
            XcFunctional::Pz => Ok(Self::Pz),
            XcFunctional::Pbe => Ok(Self::Pbe),
            XcFunctional::Pbe0 => Err(PwdftError::NotImplemented { what: "xc_functional 'pbe0'".into() }),
            XcFunctional::Hse06 => Err(PwdftError::NotImplemented { what: "xc_functional 'hse06'".into() }),
        }
    }

    /// Whether this functional requires ∇ρ on the real-space grid.
    ///
    /// Drivers call this once per SCF iteration to decide whether to spend
    /// the gradient FFTs. `false` for LDA (zero FFT work beyond the LDA
    /// pipeline); `true` for any GGA.
    #[must_use]
    pub fn needs_gradient(&self) -> bool {
        match self {
            Self::Pz => false,
            Self::Pbe => true,
        }
    }

    /// Evaluate the (non-spin) exchange-correlation energy density and
    /// potential on a real-space density grid.
    ///
    /// For the LDA variant ([`XcEvaluator::Pz`]) this returns
    ///
    /// ```text
    ///     exc_r[i] = ε_xc(ρ_xc(r_i))           (eV per electron)
    ///     v1_r[i]  = d[ρ · ε_xc(ρ)]/dρ|_{ρ_xc(r_i)}  (eV)
    ///     v2_r     = None                      (no GGA channel)
    /// ```
    ///
    /// For GGAs the return shape additionally populates `v2_r` with the
    /// contracted semilocal-∇ρ partial derivative; the PBE variant is
    /// declared here so drivers compile without `match` arms changing,
    /// but `eval` returns [`PwdftError::NotImplemented`] until the PBE
    /// integrator lands.
    ///
    /// Inputs:
    /// - `rho_r`: electron density on the FFT grid in e/Å³. Under NLCC
    ///   the caller passes `ρ_val + ρ_core` here (see
    ///   `scf::potentials::compute_core_density`); without NLCC it is
    ///   `ρ_val` alone.
    /// - `rho_grad_r`: the three Cartesian components of ∇ρ at each
    ///   grid point (Å⁻¹ · e/Å³ = e/Å⁴). Ignored for [`Self::Pz`];
    ///   required (not yet used) for [`Self::Pbe`].
    ///
    /// Returns [`XcGridResult`] (see that struct for the full shape
    /// contract). `v2_r.is_none()` on the LDA path short-circuits the
    /// caller's semilocal V_xc assembly, so the LDA pipeline does zero
    /// gradient-FFT work.
    ///
    /// Reference: Slater, *Phys. Rev.* **81**, 385 (1951);
    /// Perdew & Zunger, *Phys. Rev. B* **23**, 5048 (1981).
    ///
    /// # Errors
    /// Returns [`PwdftError::NotImplemented`] for any variant whose
    /// functional has not yet been ported (currently `Pbe`).
    pub fn eval(
        &self,
        rho_r: &[f64],
        rho_grad_r: Option<&[[f64; 3]]>,
    ) -> Result<XcGridResult> {
        match self {
            Self::Pz => {
                // `v2_r = None` makes the caller's semilocal V_xc
                // assembly short-circuit to `v1_r`.
                let _ = rho_grad_r; // LDA ignores the gradient.
                let (exc_r, v1_r) = lda_xc_grid(rho_r);
                Ok(XcGridResult { exc_r, v1_r, v2_r: None })
            }
            Self::Pbe => {
                // GGAP Phase B: `pbe_exchange` (the non-spin exchange half)
                // is implemented and unit-tested in this module, but the
                // matching correlation is still pending Phase C. Without
                // PW92-based PBE correlation the functional is not
                // well-defined. If the caller supplies the density
                // gradient we exercise `pbe_exchange` on the first grid
                // point (defensive call so a future refactor cannot
                // silently bitrot the Phase-B routine), then fail fast
                // with a scoped NotImplemented until Phase C lands.
                if let Some(grad) = rho_grad_r
                    && let Some(&rho0) = rho_r.first()
                    && let Some(&g0) = grad.first()
                {
                    let grad_mag0 = (g0[0] * g0[0] + g0[1] * g0[1] + g0[2] * g0[2]).sqrt();
                    let _ = pbe_exchange(rho0, grad_mag0);
                }
                Err(PwdftError::NotImplemented { what: "pbe_correlation".into() })
            }
        }
    }

    /// Evaluate the collinear spin-polarized exchange-correlation
    /// energy density and per-channel potential on a real-space density
    /// grid.
    ///
    /// For the LDA variant ([`XcEvaluator::Pz`]) this returns
    ///
    /// ```text
    ///     exc_r[i]      = ε_xc(ρ↑(r_i), ρ↓(r_i))    (eV per electron)
    ///     v1_up_r[i]    = δE_xc / δρ↑ at r_i         (eV)
    ///     v1_down_r[i]  = δE_xc / δρ↓ at r_i         (eV)
    ///     v2_up_r       = v2_down_r = None           (no GGA channels)
    /// ```
    ///
    /// using the von Barth-Hedin spin-interpolated PZ parametrization
    /// of [`lda_xc_spin_grid`]. The PBE branch is declared here but
    /// returns [`PwdftError::NotImplemented`] until the semilocal
    /// integrator lands.
    ///
    /// Inputs:
    /// - `rho_up_r`, `rho_down_r`: per-channel densities on the FFT
    ///   grid in e/Å³. Under NLCC the caller adds `ρ_core/2` to each
    ///   channel (the core is assumed spin-unpolarized).
    /// - `rho_grad_up_r`, `rho_grad_down_r`: per-channel ∇ρ on the
    ///   same grid (e/Å⁴). Ignored for [`Self::Pz`]; required for
    ///   [`Self::Pbe`] once that variant is implemented.
    ///
    /// Returns [`XcSpinGridResult`]; `v2_*_r` is `None` on the LDA
    /// path, so the caller's semilocal V_xc assembly degenerates to
    /// `V_xc^σ = v1_σ` with no gradient-FFT work.
    ///
    /// Reference: von Barth & Hedin, *J. Phys. C* **5**, 1629 (1972);
    /// Perdew & Zunger, *Phys. Rev. B* **23**, 5048 (1981) §III.
    ///
    /// # Errors
    /// Returns [`PwdftError::NotImplemented`] for any variant whose
    /// functional has not yet been ported (currently `Pbe`).
    pub fn eval_spin(
        &self,
        rho_up_r: &[f64],
        rho_down_r: &[f64],
        rho_grad_up_r: Option<&[[f64; 3]]>,
        rho_grad_down_r: Option<&[[f64; 3]]>,
    ) -> Result<XcSpinGridResult> {
        match self {
            Self::Pz => {
                let _ = (rho_grad_up_r, rho_grad_down_r); // LDA ignores gradients.
                let (exc_r, v1_up_r, v1_down_r) = lda_xc_spin_grid(rho_up_r, rho_down_r);
                Ok(XcSpinGridResult {
                    exc_r,
                    v1_up_r,
                    v1_down_r,
                    v2_up_r: None,
                    v2_down_r: None,
                })
            }
            Self::Pbe => {
                // GGAP Phase B: exchange half implemented (see
                // `pbe_exchange` in this module) but PBE correlation is
                // still pending Phase C. Spin scaling will be wired up
                // in Phase D once Phase C lands; for now the spin path
                // bails with the same scoped NotImplemented marker.
                let _ = (rho_grad_up_r, rho_grad_down_r);
                Err(PwdftError::NotImplemented { what: "pbe_correlation".into() })
            }
        }
    }
}

/// Result of a non-spin XC evaluation on the FFT grid.
///
/// All grids have length `n_grid`. Energies in eV.
///
/// `v2_r` is the GGA semilocal-∇ρ partial derivative, already contracted
/// with ∇ρ into the vector field `h(r) = 2 · (∂ρ·ε_xc/∂σ) · ∇ρ(r)` that
/// appears inside the divergence term of V_xc. LDA leaves it `None`; the
/// caller's V_xc assembly then degenerates to the LDA expression
/// `V_xc = v1_r` with no FFT work.
#[derive(Debug, Clone)]
pub struct XcGridResult {
    /// Exchange-correlation energy density ε_xc(r) (eV per electron).
    pub exc_r: Vec<f64>,
    /// V_xc contribution `∂(ρ·ε_xc)/∂ρ` on the FFT grid (eV).
    pub v1_r: Vec<f64>,
    /// GGA-only contribution: `h(r) = 2 · (∂(ρ·ε_xc)/∂σ) · ∇ρ(r)` (3-vector
    /// per grid point). The caller feeds this into `∇·h` to complete the
    /// semilocal V_xc. `None` for LDA.
    pub v2_r: Option<Vec<[f64; 3]>>,
}

/// Result of a spin-polarized XC evaluation on the FFT grid.
///
/// Mirrors [`XcGridResult`] but with per-channel V_xc^σ grids.
#[derive(Debug, Clone)]
pub struct XcSpinGridResult {
    /// Shared energy density ε_xc(r) on the FFT grid (eV).
    pub exc_r: Vec<f64>,
    /// V_xc contribution for the ↑ channel: `∂(ρ·ε_xc)/∂ρ↑` (eV).
    pub v1_up_r: Vec<f64>,
    /// V_xc contribution for the ↓ channel.
    pub v1_down_r: Vec<f64>,
    /// GGA-only per-channel `h_↑(r)`. `None` for LDA.
    pub v2_up_r: Option<Vec<[f64; 3]>>,
    /// GGA-only per-channel `h_↓(r)`. `None` for LDA.
    pub v2_down_r: Option<Vec<[f64; 3]>>,
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

    // -----------------------------------------------------------------------
    // XcEvaluator (GGAP Phase A) dispatcher tests.
    //
    // These pin the contract that the enum stays *data-driven* — variants
    // name functionals, dispatch is a `match`, no trait objects, no closures.
    // The LDA path must produce bit-for-bit the same arrays as the pre-Phase-A
    // direct-`lda_xc_grid` call site; any drift here would be a regression
    // against every QE-validated LDA integration test.
    // -----------------------------------------------------------------------

    #[test]
    fn xc_evaluator_from_settings_maps_variants_correctly() {
        assert_eq!(
            XcEvaluator::from_settings(XcFunctional::Pz).unwrap(),
            XcEvaluator::Pz,
        );
        assert_eq!(
            XcEvaluator::from_settings(XcFunctional::Pbe).unwrap(),
            XcEvaluator::Pbe,
        );

        // Hybrids fail fast at construction time with a NotImplemented error.
        let err = XcEvaluator::from_settings(XcFunctional::Pbe0).unwrap_err();
        match err {
            PwdftError::NotImplemented { what } => assert_eq!(what, "xc_functional 'pbe0'"),
            other => panic!("expected NotImplemented, got {other:?}"),
        }
        let err = XcEvaluator::from_settings(XcFunctional::Hse06).unwrap_err();
        match err {
            PwdftError::NotImplemented { what } => assert_eq!(what, "xc_functional 'hse06'"),
            other => panic!("expected NotImplemented, got {other:?}"),
        }
    }

    #[test]
    fn xc_evaluator_needs_gradient_flags() {
        assert!(!XcEvaluator::Pz.needs_gradient(), "LDA must not need ∇ρ");
        assert!(XcEvaluator::Pbe.needs_gradient(), "PBE needs ∇ρ");
    }

    #[test]
    fn xc_evaluator_pz_eval_matches_direct_lda_grid() {
        // Densities spanning both PZ regimes (rs >= 1 and rs < 1) with the
        // same shape the SCF driver produces after add_core_density.
        let rho_r: Vec<f64> = (1..=2000).map(|i| 0.001 + f64::from(i) * 0.0005).collect();
        let (direct_exc, direct_v1) = lda_xc_grid(&rho_r);

        let eval = XcEvaluator::Pz;
        let result = eval.eval(&rho_r, None).expect("PZ path must not error");

        assert_eq!(result.exc_r.len(), direct_exc.len());
        assert_eq!(result.v1_r.len(), direct_v1.len());
        assert!(
            result.v2_r.is_none(),
            "LDA must leave v2_r absent so the semilocal assembly short-circuits"
        );

        // Bit-for-bit: `eval` is a simple pass-through on `Pz`, so nothing
        // above the tolerance of `lda_xc_grid`'s own determinism is allowed.
        for (i, (&direct, &routed)) in direct_exc.iter().zip(result.exc_r.iter()).enumerate() {
            assert!(
                (direct - routed).abs() < 1e-15,
                "ε_xc mismatch at {i}: direct={direct}, eval={routed}"
            );
        }
        for (i, (&direct, &routed)) in direct_v1.iter().zip(result.v1_r.iter()).enumerate() {
            assert!(
                (direct - routed).abs() < 1e-15,
                "v1_r mismatch at {i}: direct={direct}, eval={routed}"
            );
        }
    }

    #[test]
    fn xc_evaluator_pz_eval_ignores_gradient_argument() {
        // Regression pin: the LDA path must not observe rho_grad_r. If a
        // future refactor accidentally threads ∇ρ into the LDA branch, the
        // result should be indistinguishable from passing None.
        let rho_r: Vec<f64> = (1..=256).map(|i| 0.01 + f64::from(i) * 0.001).collect();
        let fake_grad: Vec<[f64; 3]> = rho_r.iter().map(|&r| [r * 0.5, -r, 2.0 * r]).collect();

        let no_grad = XcEvaluator::Pz.eval(&rho_r, None).unwrap();
        let with_grad = XcEvaluator::Pz.eval(&rho_r, Some(&fake_grad)).unwrap();

        for (i, (&a, &b)) in no_grad.exc_r.iter().zip(with_grad.exc_r.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-15,
                "LDA eval must be ∇ρ-invariant at {i}: {a} vs {b}"
            );
        }
        for (i, (&a, &b)) in no_grad.v1_r.iter().zip(with_grad.v1_r.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-15,
                "LDA v1_r must be ∇ρ-invariant at {i}: {a} vs {b}"
            );
        }
    }

    #[test]
    fn xc_evaluator_pz_eval_spin_matches_direct_lda_spin_grid() {
        let n = 1024;
        let rho_up: Vec<f64> = (1..=n).map(|i| 0.01 + f64::from(i) * 0.001).collect();
        let rho_down: Vec<f64> = (1..=n).map(|i| 0.005 + f64::from(i) * 0.0007).collect();

        let (direct_exc, direct_up, direct_down) = lda_xc_spin_grid(&rho_up, &rho_down);
        let result = XcEvaluator::Pz
            .eval_spin(&rho_up, &rho_down, None, None)
            .expect("PZ spin path must not error");

        assert!(result.v2_up_r.is_none() && result.v2_down_r.is_none());
        for (i, (&d, &r)) in direct_exc.iter().zip(result.exc_r.iter()).enumerate() {
            assert!((d - r).abs() < 1e-15, "ε_xc spin mismatch at {i}");
        }
        for (i, (&d, &r)) in direct_up.iter().zip(result.v1_up_r.iter()).enumerate() {
            assert!((d - r).abs() < 1e-15, "v1_up mismatch at {i}");
        }
        for (i, (&d, &r)) in direct_down.iter().zip(result.v1_down_r.iter()).enumerate() {
            assert!((d - r).abs() < 1e-15, "v1_down mismatch at {i}");
        }
    }

    #[test]
    fn xc_evaluator_pbe_eval_returns_not_implemented() {
        // GGAP Phase B landed the exchange half (see the `pbe_exchange_*`
        // unit tests below); correlation + PW92 is Phase C. The grid
        // evaluator therefore bails with a scoped NotImplemented marker
        // pointing specifically at the correlation half so the next-phase
        // author can search for it and knows where to drop the fix in.
        let rho_r = vec![0.1; 32];
        let err = XcEvaluator::Pbe
            .eval(&rho_r, None)
            .expect_err("PBE eval must be NotImplemented until Phase C");
        match err {
            PwdftError::NotImplemented { what } => {
                assert_eq!(what, "pbe_correlation");
            }
            other => panic!("expected NotImplemented, got {other:?}"),
        }

        // Same shape on the spin path.
        let rho_down = vec![0.05; 32];
        let err = XcEvaluator::Pbe
            .eval_spin(&rho_r, &rho_down, None, None)
            .expect_err("PBE spin eval must be NotImplemented until Phase C");
        match err {
            PwdftError::NotImplemented { what } => {
                assert_eq!(what, "pbe_correlation");
            }
            other => panic!("expected NotImplemented, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------------
    // GGAP Phase B: PBE exchange (non-spin) unit tests.
    //
    // Four shapes of test:
    //   1. F_x(s) pinned against analytic PBE eq. 14 at s ∈ {0, 0.1, 1, 5, 10}
    //      to 1e-14 (pure analytic expression, no unit conversion drift).
    //   2. s=0 reduces exactly to the LDA/Slater exchange energy density
    //      (this is the Slater limit of the PBE enhancement factor).
    //   3. One-point cross-check at (ρ = 0.1 e/Å³, |∇ρ| = 0.05 e/Å⁴)
    //      against a hand-derived value from the PBE 1996 paper formula.
    //   4. Low-density short-circuit (ρ ≪ QE threshold).
    // -----------------------------------------------------------------------

    /// PBE enhancement factor, eq. (14) of Perdew-Burke-Ernzerhof
    /// *PRL* **77**, 3865 (1996). Used by the unit tests to pin the
    /// hard-coded F_x(s) values we claim — kept here so the tests pull
    /// from the paper's formula, not from the function under test.
    fn fx_pbe_analytic(s: f64) -> f64 {
        let kappa = 0.804;
        let mu = 0.219_514_972_764_517_1;
        1.0 + kappa - kappa / (1.0 + mu * s * s / kappa)
    }

    #[test]
    fn pbe_exchange_fx_matches_pbe_paper_formula() {
        // At five values of s, the PBE F_x should match the analytic
        // PBE 1996 paper formula to machine precision. We extract F_x
        // from the function under test by evaluating at a fixed ρ and
        // inverting eps_x = ρ · ε_x^LDA · F_x.
        let rho = 0.10_f64; // e/Å³
        let bohr3 = crate::consts::BOHR3_TO_ANG3;
        let rho_bohr = rho * bohr3;
        let kf_bohr = (3.0 * PI * PI * rho_bohr).cbrt(); // Bohr⁻¹
        let ex_lda_ha = -(0.75 / PI) * kf_bohr; // Ha per electron
        // Energy density ρ·ε_x^LDA in eV/Å³. Per-electron exchange is
        // unit-agnostic (Ha/e = HA_TO_EV eV/e), so no Bohr↔Å factor is
        // needed — just multiply ρ[e/Å³] by ε_x^LDA[eV/e].
        let eps_lda = rho * ex_lda_ha * crate::consts::HA_TO_EV; // eV/Å³

        for &s in &[0.0_f64, 0.1, 1.0, 5.0, 10.0] {
            // s = |∇ρ|/(2 k_F ρ)  ⇒  |∇ρ| = s · 2 k_F · ρ with k_F and
            // |∇ρ| in matched units. Convert k_F to Å⁻¹ to get |∇ρ| in
            // e/Å⁴.
            let kf_inv_ang = kf_bohr / crate::consts::BOHR_TO_ANG;
            let grad_mag = s * 2.0 * kf_inv_ang * rho; // e/Å⁴

            let (eps_x, _v1, _v2) = pbe_exchange(rho, grad_mag);
            let fx_numeric = eps_x / eps_lda;
            let fx_expected = fx_pbe_analytic(s);

            assert!(
                (fx_numeric - fx_expected).abs() < 1e-14,
                "F_x(s={s}) numeric={fx_numeric}, analytic={fx_expected}"
            );
        }
    }

    #[test]
    fn pbe_exchange_zero_gradient_reduces_to_lda() {
        // |∇ρ| = 0 is the deep Slater limit: eps_x = ρ·ε_x^LDA, v1 =
        // (4/3)·ε_x^LDA, v2 = 0. Matches the LDA path bit-for-bit
        // modulo the f64·multiply reordering between `lda_xc_grid` and
        // the explicit AU→eV/Å conversion inside `pbe_exchange`.
        for &rho in &[0.01_f64, 0.05, 0.1, 0.5, 1.0, 2.0] {
            let (eps_pbe, v1_pbe, v2_pbe) = pbe_exchange(rho, 0.0);
            // LDA reference: ε_x^LDA computed in Hartree at ρ_au, then
            // converted to per-electron eV and multiplied by ρ[e/Å³]
            // to get energy density eV/Å³.
            let bohr3 = crate::consts::BOHR3_TO_ANG3;
            let kf = (3.0 * PI * PI * rho * bohr3).cbrt();
            let exunif_ha = -(0.75 / PI) * kf;
            let eps_lda = rho * exunif_ha * crate::consts::HA_TO_EV;
            let v1_lda = (4.0 / 3.0) * exunif_ha * crate::consts::HA_TO_EV;

            assert!(
                (eps_pbe - eps_lda).abs() < 1e-12 * eps_lda.abs().max(1.0),
                "eps_x(∇ρ=0) mismatch at ρ={rho}: pbe={eps_pbe}, lda={eps_lda}"
            );
            assert!(
                (v1_pbe - v1_lda).abs() < 1e-12 * v1_lda.abs().max(1.0),
                "v1_x(∇ρ=0) mismatch at ρ={rho}: pbe={v1_pbe}, lda={v1_lda}"
            );
            assert!(
                v2_pbe.abs() < 1e-30,
                "v2_x(∇ρ=0) must be exactly 0 at ρ={rho}: got {v2_pbe}"
            );
        }
    }

    #[test]
    fn pbe_exchange_single_point_reference() {
        // One-point reference at (ρ = 0.1 e/Å³, |∇ρ| = 0.05 e/Å⁴).
        //
        // Reference is hand-derived from the PBE 1996 paper formula
        // (eq. 14 + ε_x^LDA expression) with constants κ=0.804,
        // μ=0.2195149727645171. We recompute the chain ρ → ρ_au → k_F
        // → exunif → s → F_x → eps_x inline in the test so a future
        // unit-conversion refactor can't silently desync from the paper
        // formula.
        let rho = 0.10_f64; // e/Å³
        let grad = 0.05_f64; // e/Å⁴

        let bohr3 = crate::consts::BOHR3_TO_ANG3;
        let bohr = crate::consts::BOHR_TO_ANG;
        let ha = crate::consts::HA_TO_EV;

        let rho_au = rho * bohr3;
        let grho_mag_au = grad * bohr3 * bohr;
        let kf_au = (3.0 * PI * PI * rho_au).cbrt();
        let exunif_ha = -(0.75 / PI) * kf_au;
        let s = grho_mag_au / (2.0 * kf_au * rho_au);
        let fx = fx_pbe_analytic(s);

        // Expected energy density (eV/Å³): full PBE.
        let eps_expected = rho_au * exunif_ha * fx * ha / bohr3;

        // Expected v1 = d(ρ·ε_x^LDA)/dρ + d(ρ·ε_x^LDA·(F_x−1))/dρ.
        // The gradient part equals QE's v1x (see `pbe_exchange` source
        // comment), so re-derive from that decomposition.
        let kappa = 0.804_f64;
        let mu = 0.219_514_972_764_517_1_f64;
        let s2 = s * s;
        let f2 = 1.0 + mu * s2 / kappa;
        let fx_qe = kappa - kappa / f2; // QE's fx = F_x^task − 1
        let dfx = 2.0 * mu * s / (f2 * f2);
        let ds = -(4.0 / 3.0) * s;
        let dxunif = exunif_ha / 3.0;
        let v1_grad_ha = exunif_ha * fx_qe + dxunif * fx_qe + exunif_ha * dfx * ds;
        let v1_lda_ha = (4.0 / 3.0) * exunif_ha;
        let v1_expected = (v1_grad_ha + v1_lda_ha) * ha;

        // Expected v2 (QE convention, so h = v2·∇ρ):
        //   v2_ha = exunif · dfx · dsg / agrho     [Ha·Bohr⁵/e]
        //         = exunif · dfx · (0.5/kf) / agrho
        let dsg = 0.5 / kf_au;
        let v2_ha = exunif_ha * dfx * dsg / grho_mag_au;
        let v2_expected = v2_ha * ha * bohr.powi(5);

        let (eps, v1, v2) = pbe_exchange(rho, grad);

        // 1e-14 eV/Å³ is far tighter than the task's 1e-10 eV/atom
        // downstream budget; the only sources of drift here are f64
        // rounding in `.cbrt()` and the chain of multiplies.
        assert!(
            (eps - eps_expected).abs() < 1e-14,
            "eps_x mismatch at (ρ=0.1, ∇ρ=0.05): got={eps}, want={eps_expected}"
        );
        assert!(
            (v1 - v1_expected).abs() < 1e-14,
            "v1_x mismatch at (ρ=0.1, ∇ρ=0.05): got={v1}, want={v1_expected}"
        );
        assert!(
            (v2 - v2_expected).abs() < 1e-14,
            "v2_x mismatch at (ρ=0.1, ∇ρ=0.05): got={v2}, want={v2_expected}"
        );

        // Sanity: eps_x and v1_x must be negative (exchange is
        // attractive) and v2_x must be negative (the gradient correction
        // raises energy, so its σ-derivative is negative because ε_x^LDA
        // is negative — double-check by sign tracking).
        assert!(eps < 0.0, "ε_x must be negative: {eps}");
        assert!(v1 < 0.0, "V1_x must be negative: {v1}");
        assert!(v2 < 0.0, "v2_x must be negative: {v2}");
    }

    #[test]
    fn pbe_exchange_low_density_short_circuits() {
        // Below QE's rho_threshold_gga (1e-6 e/Bohr³ ≈ 6.75e-6 e/Å³)
        // the return must be exact zeros — this protects the SCF driver
        // from NaNs at density tails far from atoms.
        for &rho in &[0.0_f64, 1e-8, 1e-7, 1e-6] {
            let (eps, v1, v2) = pbe_exchange(rho, 0.01);
            assert!(
                eps.abs() < 1e-30 && v1.abs() < 1e-30 && v2.abs() < 1e-30,
                "low-density clamp failed at ρ={rho}: eps={eps}, v1={v1}, v2={v2}"
            );
        }
    }
}
