//! Charge density construction from wavefunctions.
//!
//! ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²

use nalgebra::DMatrix;
use num_complex::Complex64;

use crate::{basis::BasisSet, fft::FFT3D, kpoints::KPoint};

/// Compute the charge density on the real-space FFT grid from wavefunctions.
///
/// For each (k-point, band):
/// 1. Place PW coefficients onto the FFT grid
/// 2. Inverse FFT to get ψ(r)
/// 3. Accumulate f_{n,k} × w_k × |ψ(r)|²
///
/// The result is normalized so that ∫ρ(r)dr = N_electrons.
pub fn compute_density(
    basis: &BasisSet,
    kpoints: &[KPoint],
    wavefunctions: &[DMatrix<Complex64>],
    occupations: &[Vec<f64>],
    g_to_fft: &[usize],
    fft: &FFT3D,
    n_electrons: f64,
    omega: f64,
) -> Vec<f64> {
    let n_grid = fft.total_size();
    let n_pw = basis.len();
    let mut rho_r = vec![0.0; n_grid];

    for (ik, kp) in kpoints.iter().enumerate() {
        let wfn = &wavefunctions[ik];
        let occ = &occupations[ik];
        let n_bands = occ.len();

        for ib in 0..n_bands {
            let f = occ[ib] * kp.weight;
            if f < 1e-15 {
                continue;
            }

            // Place PW coefficients on FFT grid
            let mut psi_g = vec![Complex64::new(0.0, 0.0); n_grid];
            for ig in 0..n_pw {
                psi_g[g_to_fft[ig]] = wfn[(ig, ib)];
            }

            // Inverse FFT: ψ(G) → ψ(r)
            fft.inverse(&mut psi_g);
            // Note: no normalization needed for |ψ|² because it cancels

            // Accumulate |ψ(r)|²
            for (i, &psi) in psi_g.iter().enumerate() {
                rho_r[i] += f * psi.norm_sqr();
            }
        }
    }

    // Normalize: the integral ∫ρ(r)dr should equal N_electrons
    // ∫ρ(r)dr = (Ω/N_grid) Σ_r ρ(r) = N_electrons
    // The FFT convention introduces a factor of N_grid from the inverse FFT
    let dvol = omega / n_grid as f64;
    let integral: f64 = rho_r.iter().sum::<f64>() * dvol;

    if integral.abs() > 1e-15 {
        let scale = n_electrons / integral;
        for v in &mut rho_r {
            *v *= scale;
        }
    }

    rho_r
}
