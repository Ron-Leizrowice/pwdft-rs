//! Total energy computation and density utilities for SCF.

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

/// XC energy with double-counting correction: E_xc - E_vxc.
///
/// `rho_xc`: density for E_xc (ρ_val + ρ_core if NLCC).
/// `rho_val`: valence density only for E_vxc double-counting.
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
/// Clamps to non-negative to avoid NaN in XC.
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
