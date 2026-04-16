//! Total energy computation and density utilities for SCF.

use num_complex::Complex64;

use crate::{
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
