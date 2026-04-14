//! Charge density construction from wavefunctions.
//!
//! ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²

use nalgebra::DMatrix;
use num_complex::Complex64;
use rayon::prelude::*;

use crate::{basis::BasisSet, fft::FFT3D, kpoints::KPoint};

/// Compute the charge density on the real-space FFT grid from wavefunctions.
///
/// For each (k-point, band):
/// 1. Place PW coefficients onto the FFT grid
/// 2. Inverse FFT to get ψ(r)
/// 3. Accumulate f_{n,k} × w_k × |ψ(r)|²
///
/// The result is normalized so that ∫ρ(r)dr = N_electrons.
///
/// K-point contributions are computed in parallel and reduced.
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

    // Each k-point computes its contribution independently, then we reduce
    let rho_per_k: Vec<Vec<f64>> = (0..kpoints.len())
        .into_par_iter()
        .map(|ik| {
            let kp = &kpoints[ik];
            let wfn = &wavefunctions[ik];
            let occ = &occupations[ik];
            let n_bands = occ.len();
            let mut rho_k = vec![0.0; n_grid];

            // Each k-point needs its own FFT instance (plans are Arc-shared internally)
            let fft_local = FFT3D::new(fft.dims()[0], fft.dims()[1], fft.dims()[2]);

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
                fft_local.inverse(&mut psi_g);

                // Accumulate |ψ(r)|²
                for (i, &psi) in psi_g.iter().enumerate() {
                    rho_k[i] += f * psi.norm_sqr();
                }
            }
            rho_k
        })
        .collect();

    // Reduce: sum contributions from all k-points
    let mut rho_r = vec![0.0; n_grid];
    for rho_k in &rho_per_k {
        for (i, &v) in rho_k.iter().enumerate() {
            rho_r[i] += v;
        }
    }

    // Normalize: the integral ∫ρ(r)dr should equal N_electrons
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
