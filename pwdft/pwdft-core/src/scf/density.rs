//! Charge density construction from wavefunctions.
//!
//! ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²

use num_complex::Complex64;
use rayon::prelude::*;

use crate::{basis::BasisSet, fft::FFT3D, kpoints::KPoint};

/// Bands with effective weight `occ × k-weight` below this threshold are skipped
/// during density reconstruction — their contribution is below eigensolver
/// round-off, and the forward FFT plus |ψ(r)|² accumulation would only add noise.
const OCCUPATION_SKIP_THRESHOLD: f64 = 1e-15;

/// Lower bound on ∫ρ(r)dr before the normalization rescale is applied. Guards
/// against dividing by a vanishingly small integral when the density on the
/// grid has underflowed to zero (e.g. all bands skipped as unoccupied).
const NORMALIZATION_INTEGRAL_FLOOR: f64 = 1e-15;

/// Grid-level parameters needed for density construction.
pub struct DensityGrid<'a> {
    pub basis: &'a BasisSet,
    pub g_to_fft: &'a [usize],
    pub fft: &'a mut FFT3D,
    pub n_electrons: f64,
    pub omega: f64,
}

/// Compute the charge density on the real-space FFT grid from wavefunctions.
///
/// For each (k-point, band):
/// 1. Place PW coefficients onto the FFT grid
/// 2. Inverse FFT to get ψ(r)
/// 3. Accumulate f_{n,k} × w_k × |ψ(r)|²
///
/// The result is normalized so that ∫ρ(r)dr = N_electrons.
/// Output density is in e/ų.
///
/// K-point contributions are computed in parallel and reduced.
pub fn compute_density(
    grid: &mut DensityGrid<'_>,
    kpoints: &[KPoint],
    wavefunctions: &[faer::Mat<Complex64>],
    occupations: &[Vec<f64>],
) -> Vec<f64> {
    let basis = grid.basis;
    let g_to_fft = grid.g_to_fft;
    let fft = &mut grid.fft;
    let n_electrons = grid.n_electrons;
    let omega = grid.omega;
    let n_grid = fft.total_size();
    let n_pw = basis.len();

    // Parallel map-reduce: each k-point accumulates into a thread-local buffer,
    // then rayon reduces by summing the buffers. Avoids n_kpoints intermediate
    // allocations — only allocates one buffer per rayon worker thread.
    let [dnx, dny, dnz] = fft.dims();
    let mut rho_r: Vec<f64> = (0..kpoints.len())
        .into_par_iter()
        .fold(
            || vec![0.0; n_grid],
            |mut rho_acc, ik| {
                let kp = &kpoints[ik];
                let wfn = &wavefunctions[ik];
                let occ = &occupations[ik];
                let n_bands = occ.len();
                let mut fft_local = FFT3D::new(dnx, dny, dnz);

                for ib in 0..n_bands {
                    let f = occ[ib] * kp.weight;
                    if f < OCCUPATION_SKIP_THRESHOLD {
                        continue;
                    }

                    let mut psi_g = vec![Complex64::new(0.0, 0.0); n_grid];
                    for ig in 0..n_pw {
                        psi_g[g_to_fft[ig]] = wfn[(ig, ib)];
                    }
                    fft_local.inverse(&mut psi_g);

                    for (i, &psi) in psi_g.iter().enumerate() {
                        rho_acc[i] += f * psi.norm_sqr();
                    }
                }
                rho_acc
            },
        )
        .reduce(
            || vec![0.0; n_grid],
            |mut a, b| {
                for (ai, &bi) in a.iter_mut().zip(b.iter()) {
                    *ai += bi;
                }
                a
            },
        );

    // Normalize: the integral ∫ρ(r)dr should equal N_electrons
    let dvol = omega / n_grid as f64;
    let integral: f64 = rho_r.iter().sum::<f64>() * dvol;

    if integral.abs() > NORMALIZATION_INTEGRAL_FLOOR {
        let scale = n_electrons / integral;
        for v in &mut rho_r {
            *v *= scale;
        }
    }

    rho_r
}
