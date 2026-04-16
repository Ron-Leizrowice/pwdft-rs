//! Potential computation on the FFT grid for SCF.
//!
//! V_local (pseudopotential), NLCC core density, and Hamiltonian construction.

use nalgebra::Vector3;
use num_complex::Complex64;
use rayon::prelude::*;

use crate::{
    basis::BasisSet,
    consts::HBAR2_OVER_2M,
    crystal::Crystal,
    pseudopotential::PseudopotentialData,
};

use super::grid::{FftGrid, g_vector_at_dims, miller_to_idx};

/// Compute local pseudopotential V_local(G) on the full FFT grid.
pub(crate) fn compute_v_local(
    crystal: &Crystal,
    grid: &FftGrid,
    pseudopotentials: &[&PseudopotentialData],
    omega: f64,
) -> Vec<Complex64> {
    let atom_data: Vec<(Vector3<f64>, &PseudopotentialData)> = crystal
        .atoms
        .iter()
        .map(|atom| {
            let pp = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials);
            (atom.cart_position(&crystal.lattice), pp)
        })
        .collect();

    let dims = grid.dims;
    let recip = grid.recip.clone();
    (0..grid.total_size())
        .into_par_iter()
        .map(|idx| {
            let g = g_vector_at_dims(idx, dims, &recip);
            let g_norm = g.norm();
            let mut v = Complex64::new(0.0, 0.0);

            for &(ref tau, pp) in &atom_data {
                let phase = -g.dot(tau);
                let sf = Complex64::cis(phase);
                let v_form = pp.v_local_of_g(g_norm, omega);
                v += sf * v_form;
            }
            v
        })
        .collect()
}

/// Compute NLCC core density on the real-space FFT grid.
///
/// Returns empty vec if no PP has NLCC.
pub(crate) fn compute_core_density(
    crystal: &Crystal,
    grid: &mut FftGrid,
    pseudopotentials: &[&PseudopotentialData],
) -> Vec<f64> {
    let any_nlcc = pseudopotentials.iter().any(|pp| pp.has_nlcc());
    if !any_nlcc {
        return vec![];
    }

    let n_grid = grid.total_size();
    let omega = crystal.lattice.volume();
    let mut rho_core_g = vec![Complex64::new(0.0, 0.0); n_grid];

    for atom in &crystal.atoms {
        let pp = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials);
        if !pp.has_nlcc() || pp.core_charge.is_empty() {
            continue;
        }

        let tau = atom.cart_position(&crystal.lattice);

        for (idx, rho_g_val) in rho_core_g.iter_mut().enumerate() {
            let g = grid.g_vector_at(idx);
            let g_norm = g.norm();

            let mut integral = 0.0;
            for ((&rho_c, &r), &dr) in pp
                .core_charge
                .iter()
                .zip(pp.r_grid.iter())
                .zip(pp.rab.iter())
            {
                let gr = g_norm * r;
                let j0 = if gr < 1e-10 {
                    1.0 - gr * gr / 6.0
                } else {
                    gr.sin() / gr
                };
                integral += rho_c * j0 * dr;
            }

            let phase = -g.dot(&tau);
            let sf = Complex64::cis(phase);
            *rho_g_val += sf * (integral / omega);
        }
    }

    grid.fft.inverse(&mut rho_core_g);
    rho_core_g.iter().map(|c| c.re).collect()
}

/// Build Hamiltonian matrix: kinetic + V_eff(G-G') from FFT grid.
pub(crate) fn build_hamiltonian_with_v_eff(
    basis: &BasisSet,
    k: &Vector3<f64>,
    v_eff_fft: &[Complex64],
    grid_dims: [usize; 3],
) -> faer::Mat<Complex64> {
    let n = basis.len();
    let mut h = faer::Mat::<Complex64>::zeros(n, n);

    for (i, g) in basis.g_vectors().iter().enumerate() {
        let ke = HBAR2_OVER_2M * (k + g).norm_squared();
        h[(i, i)] = Complex64::new(ke, 0.0);
    }

    let miller_idx = basis.miller_indices();
    for i in 0..n {
        for j in 0..n {
            let dn1 = miller_idx[i][0] - miller_idx[j][0];
            let dn2 = miller_idx[i][1] - miller_idx[j][1];
            let dn3 = miller_idx[i][2] - miller_idx[j][2];
            let fft_idx = miller_to_idx(grid_dims, dn1, dn2, dn3);
            h[(i, j)] += v_eff_fft[fft_idx];
        }
    }

    h
}
