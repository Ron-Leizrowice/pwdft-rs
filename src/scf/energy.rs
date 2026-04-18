//! Total energy computation and density utilities for SCF.
//!
//! ## Nonlinear core correction (NLCC)
//!
//! Several helpers here — `add_core_density`, `xc_energy_corrected` — are
//! shared between the standard Kohn-Sham path and the NLCC path. NLCC
//! (Louie, Froyen, Cohen, *Phys. Rev. B* **26**, 1738 (1982)) restores
//! the nonlinear coupling
//! ```text
//!     E_xc[ρ_val + ρ_core] − E_xc[ρ_val]
//! ```
//! that is dropped when the core is frozen and orthogonalized out of the
//! valence problem. The implementation keeps the core *only* inside the
//! XC functional:
//!
//! - `ρ_val + ρ_core` enters `ε_xc[·]` and `v_xc[·]` (see
//!   [`xc_energy_corrected`] and QE `PW/src/v_of_rho.f90:511`).
//! - `ρ_val` alone enters the Hartree source, the electron count, and the
//!   double-counting integral `∫ ρ_val · v_xc dr`.
//! - In LSDA, `ρ_core` is spin-unpolarized and split evenly as
//!   `ρ_core/2` between the two spin channels before being added to each
//!   `ρ_σ` (see `scf::driver_spin::run_scf_spin`).
//!
//! `ρ_core` itself is built on the FFT grid by
//! [`scf::potentials::compute_core_density`](super::potentials::compute_core_density)
//! from the PP's `PP_NLCC` block (see `src/pseudopotential/upf/convert.rs`
//! for the storage-unit convention — bare ρ_core(r) in e/Å³, *not* the
//! 4πr²·ρ convention used by `PP_RHOATOM`).

use nalgebra::Vector3;
use num_complex::Complex64;

use crate::{
    basis::BasisSet,
    consts::HBAR2_OVER_2M,
    fft::FFT3D,
    potential::xc,
};

// ---------------------------------------------------------------------------
// Energy computation
// ---------------------------------------------------------------------------

/// Band energy: E_band = Σ_{n,k} f_{n,k} w_k ε_{n,k}
pub(crate) fn band_energy(
    eigenvalues: &[Vec<f64>],
    occupations: &[Vec<f64>],
    kpoint_weights: &[f64],
) -> f64 {
    eigenvalues
        .iter()
        .zip(occupations.iter())
        .zip(kpoint_weights.iter())
        .map(|((evs, occs), &w)| {
            evs.iter()
                .zip(occs.iter())
                .map(|(&e, &f)| f * w * e)
                .sum::<f64>()
        })
        .sum()
}

/// Hartree energy: E_H = (ω/2) Σ_G |ρ(G)|² 4πe²/|G|²
pub(crate) fn hartree_energy(rho_g: &[Complex64], g_squared: &[f64], omega: f64) -> f64 {
    let fourpi_e2 = 4.0 * std::f64::consts::PI * crate::consts::E2_COULOMB;
    rho_g
        .iter()
        .zip(g_squared.iter())
        .map(|(rho, &g2)| {
            if g2 > crate::consts::G2_ZERO_THRESHOLD {
                rho.norm_sqr() * fourpi_e2 / g2
            } else {
                0.0
            }
        })
        .sum::<f64>()
        * 0.5
        * omega
}

/// XC energy with double-counting correction: E_xc − E_vxc.
///
/// For the Kohn-Sham total energy with NLCC
/// (Louie, Froyen, Cohen, *Phys. Rev. B* **26**, 1738 (1982)):
/// ```text
///     E_xc − E_dc = E_xc[ρ_val + ρ_core]  −  ∫ ρ_val · v_xc[ρ_val + ρ_core] d³r
/// ```
/// where the XC potential `v_xc` is evaluated on the total density
/// (val + core) but the double-counting integrand couples it only to
/// the valence density — the core is frozen and does not appear in the
/// band sum. Without NLCC, `ρ_core = 0` and the formula collapses to
/// the usual `E_xc[ρ_val] − ∫ ρ_val · v_xc[ρ_val] d³r`.
///
/// - `rho_xc`: density for E_xc (ρ_val + ρ_core if NLCC, else ρ_val).
/// - `rho_val`: valence density only, used in E_vxc double-counting.
pub(crate) fn xc_energy_corrected(
    rho_xc: &[f64],
    rho_val: &[f64],
    exc_r: &[f64],
    vxc_r: &[f64],
    omega: f64,
) -> f64 {
    let n_grid = rho_xc.len();
    let dvol = omega / n_grid as f64;

    let e_xc = xc::lda_xc_energy(rho_xc, exc_r, omega);
    let e_vxc: f64 = rho_val
        .iter()
        .zip(vxc_r.iter())
        .map(|(&rho, &vxc)| rho * vxc * dvol)
        .sum();

    e_xc - e_vxc
}

/// Kohn-Sham total energy: E_band - E_H[rho_out] + (E_xc[rho_out] - E_vxc[rho_out]) + E_ewald.
///
/// Uses the OUTPUT density (from new wavefunctions) for double-counting corrections.
pub(crate) fn total_energy(
    e_band: f64,
    e_hartree: f64,
    e_xc_corrected: f64,
    e_ewald: f64,
) -> f64 {
    e_band - e_hartree + e_xc_corrected + e_ewald
}

/// Harris-Foulkes energy: E_band - E_H[rho_in] + (E_xc[rho_in] - E_vxc[rho_in]) + E_ewald.
///
/// Uses the INPUT density for all double-counting corrections but OUTPUT eigenvalues
/// (from diagonalizing H[rho_in]). This is stationary at self-consistency: first-order
/// density errors cancel, making E_HF converge quadratically to E_KS.
///
/// Reference: Harris, Phys. Rev. B 31, 1770 (1985);
///            Foulkes & Haydock, Phys. Rev. B 39, 12520 (1989).
pub(crate) fn harris_foulkes_energy(
    e_band: f64,
    e_hartree_in: f64,
    e_xc_corrected_in: f64,
    e_ewald: f64,
) -> f64 {
    e_band - e_hartree_in + e_xc_corrected_in + e_ewald
}

/// Add the `V_local(G=0) · N_electrons` uniform-background shift.
///
/// The G=0 component of the local pseudopotential is zeroed in
/// `v_local_fft` (see `ScfContext::new`) so the Hamiltonian matrix
/// elements remain finite. The constant background it represents is
/// then added back to the total and Harris-Foulkes energies as
/// `V_local(G=0) · N_el`. See proposals NCFX / VGCMP for the
/// derivation; this helper collapses the expression duplicated across
/// the spin and non-spin drivers into a single named site.
pub(crate) fn with_g0_shift(energy: f64, ctx: &super::context::ScfContext<'_>) -> f64 {
    energy + ctx.v_local_g0 * ctx.n_electrons
}

// ---------------------------------------------------------------------------
// Density utilities
// ---------------------------------------------------------------------------

/// RMS density difference between two densities.
pub(crate) fn density_diff(rho_old: &[f64], rho_new: &[f64], omega: f64, n_grid: usize) -> f64 {
    let dvol = omega / n_grid as f64;
    let sum_sq: f64 = rho_old
        .iter()
        .zip(rho_new.iter())
        .map(|(&a, &b)| (a - b).powi(2) * dvol)
        .sum();
    (sum_sq / omega).sqrt()
}

/// FFT density from real space to G-space (normalized).
pub(crate) fn density_r_to_g(fft: &mut FFT3D, rho_r: &[f64], rho_g: &mut [Complex64]) {
    for (i, &r) in rho_r.iter().enumerate() {
        rho_g[i] = Complex64::new(r, 0.0);
    }
    fft.forward(rho_g);
    let norm = 1.0 / fft.total_size() as f64;
    for v in rho_g.iter_mut() {
        *v *= norm;
    }
}

/// Convert real-space array to G-space with FFT normalization.
pub(crate) fn real_to_g_space(data_r: &[f64], fft: &mut FFT3D) -> Vec<Complex64> {
    let mut data_g = vec![Complex64::new(0.0, 0.0); data_r.len()];
    density_r_to_g(fft, data_r, &mut data_g);
    data_g
}

/// Assemble V_eff = V_local + V_H + V_xc.
pub(crate) fn assemble_v_eff(
    v_local: &[Complex64],
    v_h: &[Complex64],
    v_xc: &[Complex64],
) -> Vec<Complex64> {
    use rayon::prelude::*;
    v_local
        .par_iter()
        .zip(v_h.par_iter())
        .zip(v_xc.par_iter())
        .map(|((&vl, &vh), &vxc)| vl + vh + vxc)
        .collect()
}

/// Add NLCC core density to valence density for XC evaluation.
///
/// Implements `ρ_xc(r) = max(ρ_val(r) + ρ_core(r), 0)`. If `rho_core` is
/// empty (no NLCC), returns `rho_val` unchanged.
///
/// The clamp protects the LDA XC functional from spurious negative
/// densities that can arise from FFT-wrap round-off in ρ_core or from
/// density mixing; without it, `ρ^(1/3)` in the exchange term would
/// produce NaN. Mirrors QE `PW/src/v_of_rho.f90:511` (adds `rho_core` to
/// `rho%of_r(ir,1)` before calling `xc_lda`).
///
/// Reference: Louie, Froyen, Cohen, *Phys. Rev. B* **26**, 1738 (1982).
pub(crate) fn add_core_density(rho_val: &[f64], rho_core: &[f64]) -> Vec<f64> {
    if rho_core.is_empty() {
        rho_val.to_vec()
    } else {
        rho_val
            .iter()
            .zip(rho_core.iter())
            .map(|(&v, &c)| (v + c).max(0.0))
            .collect()
    }
}

/// Compute Hartree potential on the full FFT grid.
pub(crate) fn hartree_on_fft_grid(rho_g: &[Complex64], g_squared: &[f64]) -> Vec<Complex64> {
    use rayon::prelude::*;
    let fourpi_e2 = 4.0 * std::f64::consts::PI * crate::consts::E2_COULOMB;

    rho_g
        .par_iter()
        .zip(g_squared.par_iter())
        .map(|(&rho, &g2)| {
            if g2 > crate::consts::G2_ZERO_THRESHOLD {
                rho * fourpi_e2 / g2
            } else {
                Complex64::new(0.0, 0.0)
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Per-component diagnostics (VGC5)
// ---------------------------------------------------------------------------

/// Kinetic expectation value: Σ_{n,k} f·w·⟨ψ_{n,k}|T|ψ_{n,k}⟩ (eV).
///
/// For plane-wave coefficients c_{n,k}(G) (columns of `wavefunctions[ik]`):
/// ⟨ψ|T|ψ⟩ = Σ_G |c(G)|² · (ℏ²/2m) · |k+G|²
///
/// Coefficients are assumed orthonormal: Σ_G |c(G)|² = 1.
pub(crate) fn kinetic_expectation(
    basis: &BasisSet,
    k_points: &[Vector3<f64>],
    kpoint_weights: &[f64],
    wavefunctions: &[faer::Mat<Complex64>],
    occupations: &[Vec<f64>],
) -> f64 {
    use rayon::prelude::*;
    let g_vecs = basis.g_vectors();
    (0..k_points.len())
        .into_par_iter()
        .map(|ik| {
            let k = &k_points[ik];
            let wfn = &wavefunctions[ik];
            let occ = &occupations[ik];
            let w = kpoint_weights[ik];
            let n_pw = wfn.nrows();
            let n_bands = wfn.ncols();
            let mut e = 0.0f64;
            for nb in 0..n_bands {
                let f = occ[nb];
                if f == 0.0 {
                    continue;
                }
                let mut t_band = 0.0f64;
                for ig in 0..n_pw {
                    let c = wfn[(ig, nb)];
                    let ke = HBAR2_OVER_2M * (k + g_vecs[ig]).norm_squared();
                    t_band += c.norm_sqr() * ke;
                }
                e += f * w * t_band;
            }
            e
        })
        .sum()
}

/// Local-PP expectation (G ≠ 0 piece): ∫ρ(r)·V_local(r)dr (eV).
///
/// `v_local_fft` is the FFT-grid V_local with its G=0 component already
/// zeroed (see `ScfContext::new`). The integral is therefore
///   ∫ρ·V_local(G≠0)dr = Σ_r ρ(r)·V_local(r)·dV.
/// The compensating `V_local(G=0)·N_el` shift is reported separately.
pub(crate) fn local_pp_energy_grid(
    rho_r: &[f64],
    v_local_fft_r: &[f64],
    omega: f64,
    n_grid: usize,
) -> f64 {
    let dvol = omega / n_grid as f64;
    rho_r
        .iter()
        .zip(v_local_fft_r.iter())
        .map(|(&r, &v)| r * v * dvol)
        .sum()
}

/// Non-local PP expectation: Σ_{n,k} f·w·⟨ψ|V_NL|ψ⟩ (eV).
///
/// Uses the cached `NonlocalPotential` for each k-point to build an
/// ephemeral H_NL matrix, then computes ⟨ψ|H_NL|ψ⟩ for each band.
pub(crate) fn nonlocal_expectation(
    basis: &BasisSet,
    crystal: &crate::crystal::Crystal,
    k_points: &[Vector3<f64>],
    kpoint_weights: &[f64],
    wavefunctions: &[faer::Mat<Complex64>],
    occupations: &[Vec<f64>],
    vnl_cache: &[crate::potential::nonlocal::NonlocalPotential],
) -> f64 {
    use rayon::prelude::*;
    (0..k_points.len())
        .into_par_iter()
        .map(|ik| {
            let k = &k_points[ik];
            let wfn = &wavefunctions[ik];
            let occ = &occupations[ik];
            let w = kpoint_weights[ik];
            let n_pw = wfn.nrows();
            let n_bands = wfn.ncols();

            let mut h_nl = faer::Mat::<Complex64>::zeros(n_pw, n_pw);
            vnl_cache[ik].add_to_hamiltonian(&mut h_nl, crystal, basis, k);

            let mut e = 0.0f64;
            for nb in 0..n_bands {
                let f = occ[nb];
                if f == 0.0 {
                    continue;
                }
                // ⟨ψ|H_NL|ψ⟩ = Σ_{G,G'} c*(G) H_NL[G,G'] c(G')
                let mut acc = Complex64::new(0.0, 0.0);
                for ig in 0..n_pw {
                    let mut row_sum = Complex64::new(0.0, 0.0);
                    for jg in 0..n_pw {
                        row_sum += h_nl[(ig, jg)] * wfn[(jg, nb)];
                    }
                    acc += wfn[(ig, nb)].conj() * row_sum;
                }
                e += f * w * acc.re;
            }
            e
        })
        .sum()
}

/// Bare XC energy: ∫ρ(r)·ε_xc(r)dr (eV). Same sign as QE's "xc contribution".
/// `rho_xc` = ρ_val + ρ_core (for NLCC) or ρ_val otherwise.
pub(crate) fn xc_energy_bare(rho_xc: &[f64], exc_r: &[f64], omega: f64) -> f64 {
    xc::lda_xc_energy(rho_xc, exc_r, omega)
}

// ---------------------------------------------------------------------------
// Per-component decomposition (VGC5 diagnostic)
// ---------------------------------------------------------------------------

/// Per-component energy decomposition of a converged SCF total energy.
///
/// All values in eV. Identity (at convergence):
/// ```text
/// E_total = e_kinetic
///         + e_local
///         + e_local_g0_shift   (= V_local(G=0) · N_el)
///         + e_nonlocal
///         + e_hartree
///         + e_xc
///         + e_ewald
/// ```
/// and by the Kohn-Sham double-counting identity:
/// ```text
/// e_band = e_kinetic + e_local + e_nonlocal + 2·e_hartree + e_vxc
/// ```
/// where `e_vxc = ∫ρ(r)·V_xc(r)dr`. The `V_local(G=0)·N_el` background shift
/// is the compensating term for zeroing the G=0 component of the local
/// pseudopotential in the Hamiltonian; see `scf::context::ScfContext::new`.
///
/// Mirrors QE's `pw.x` standard-output decomposition:
/// ```text
///   one-electron contribution = e_kinetic + e_local + e_nonlocal + e_local_g0_shift
///   hartree    contribution = e_hartree
///   xc         contribution = e_xc
///   ewald      contribution = e_ewald
/// ```
/// Intended for validation (see proposal VGC5) rather than routine SCF use.
/// Computed on the final iteration by one extra pass over wavefunctions,
/// V_local on the FFT grid, and the V_NL operator.
#[derive(Debug, Clone)]
pub struct EnergyComponents {
    /// Band energy: Σ_{n,k} f_{n,k} w_k ε_{n,k}.
    pub e_band: f64,
    /// Kinetic: Σ_{n,k} f·w·⟨ψ|T|ψ⟩ = Σ_{n,k} f·w·Σ_G |c_G|² · ℏ²/(2m)·|k+G|².
    pub e_kinetic: f64,
    /// Local PP (G ≠ 0): ∫ρ(r)·V_local(r)dr on the FFT grid (G=0 excluded).
    pub e_local: f64,
    /// Local PP G=0 compensating shift: V_local(G=0)·N_el.
    /// Constant background subtracted from `v_local_fft` at setup to keep the
    /// Hamiltonian diagonal finite.
    pub e_local_g0_shift: f64,
    /// Non-local (KB separable): Σ_{n,k} f·w·⟨ψ|V_NL|ψ⟩.
    pub e_nonlocal: f64,
    /// Hartree: (Ω/2) Σ_G |ρ(G)|² · 4πe²/|G|² (from OUTPUT density).
    pub e_hartree: f64,
    /// XC energy: ∫ρ(r)·ε_xc(r)dr (from OUTPUT density; same sign as QE's
    /// "xc contribution"). NLCC: ρ here is ρ_val + ρ_core.
    pub e_xc: f64,
    /// Ewald ion-ion energy (spin- and density-independent).
    pub e_ewald: f64,
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "ERR2 § Phase 0: in-src test modules are allowed to panic")]
mod tests {
    use super::*;

    use approx::relative_eq;
    use crate::fft::FFT3D;

    #[test]
    fn test_real_to_g_space_dc_component() {
        // A constant real-space function f(r) = C should give
        // F(G=0) = C and F(G≠0) = 0.
        let mut fft = FFT3D::new(8, 8, 8);
        let n = fft.total_size();
        let c = 3.5;
        let data_r = vec![c; n];
        let data_g = real_to_g_space(&data_r, &mut fft);

        // G=0 component (index 0) should be C
        assert!(
            relative_eq!(data_g[0].re, c, epsilon = 1e-10),
            "DC component: expected {c}, got {}", data_g[0].re
        );
        assert!(data_g[0].im.abs() < 1e-10);

        // All other G-components should be ~0
        for (i, &v) in data_g.iter().enumerate().skip(1) {
            assert!(
                v.norm() < 1e-10,
                "G≠0 component at {i}: expected ~0, got {v}"
            );
        }
    }

    #[test]
    fn test_real_to_g_space_roundtrip() {
        let mut fft = FFT3D::new(8, 8, 8);
        let n = fft.total_size();
        let data_r: Vec<f64> = (0..n).map(|i| (i as f64 * 0.1).sin()).collect();
        let data_g = real_to_g_space(&data_r, &mut fft);

        // Inverse FFT should recover original (unnormalized → need N factor)
        let mut data_back = data_g;
        fft.inverse(&mut data_back);
        // real_to_g_space divides by N, inverse multiplies by N → should recover original
        for (i, (&orig, &back)) in data_r.iter().zip(data_back.iter()).enumerate() {
            assert!(
                relative_eq!(orig, back.re, epsilon = 1e-10),
                "Roundtrip failed at {i}: original={orig}, recovered={}", back.re
            );
            assert!(back.im.abs() < 1e-10, "Imaginary part at {i}: {}", back.im);
        }
    }

    #[test]
    fn test_assemble_v_eff_adds_correctly() {
        let n = 100;
        let v1: Vec<Complex64> = (0..n).map(|i| Complex64::new(i as f64, 0.0)).collect();
        let v2: Vec<Complex64> = (0..n).map(|i| Complex64::new(0.0, i as f64 * 0.1)).collect();
        let v3: Vec<Complex64> = (0..n).map(|i| Complex64::new(-(i as f64) * 0.5, 0.0)).collect();

        let result = assemble_v_eff(&v1, &v2, &v3);

        for i in 0..n {
            let expected = v1[i] + v2[i] + v3[i];
            assert!(
                (result[i] - expected).norm() < 1e-14,
                "V_eff mismatch at {i}: expected {expected}, got {}", result[i]
            );
        }
    }

    #[test]
    fn test_hartree_on_fft_grid_g0_zero() {
        // V_H(G=0) should be zero (no divergence)
        let rho_g = vec![Complex64::new(1.0, 0.0); 10];
        let g_squared = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let v_h = hartree_on_fft_grid(&rho_g, &g_squared);
        assert!(v_h[0].norm() < 1e-15, "V_H(G=0) should be zero, got {}", v_h[0]);
        // V_H(G≠0) should be finite and positive real for positive ρ
        for &v in &v_h[1..] {
            assert!(v.re > 0.0, "V_H should be positive for positive ρ: {v}");
        }
    }

    #[test]
    fn test_density_diff_identical() {
        let rho = vec![1.0; 100];
        let diff = density_diff(&rho, &rho, 40.0, 100);
        assert!(diff < 1e-15, "Identical densities should give zero diff: {diff}");
    }

    #[test]
    fn test_density_diff_known() {
        let omega = 40.0;
        let n = 100;
        let rho_a = vec![1.0; n];
        let rho_b = vec![2.0; n];
        // diff = sqrt(Σ(1.0)² × dvol / omega) = sqrt(n × dvol / omega) = sqrt(dvol × n / omega)
        // dvol = omega / n = 0.4
        // diff = sqrt(0.4 * 100 / 40) = sqrt(1.0) = 1.0
        let diff = density_diff(&rho_a, &rho_b, omega, n);
        assert!(
            relative_eq!(diff, 1.0, epsilon = 1e-10),
            "Expected diff=1.0, got {diff}"
        );
    }
}
