//! Potential computation on the FFT grid for SCF.
//!
//! V_local (pseudopotential), NLCC core density, and Hamiltonian construction.

use std::collections::HashMap;

use elements_rs::Element;
use nalgebra::Vector3;
use num_complex::Complex64;
use rayon::prelude::*;

use super::grid::{FftGrid, g_vector_at_dims, miller_to_idx};
use crate::{
    basis::BasisSet, consts::HBAR2_OVER_2M, crystal::Crystal, error::Result, pseudopotential::UpfPseudoPotential,
};

/// Compute local pseudopotential V_local(G) on the full FFT grid.
pub(crate) fn compute_v_local(
    crystal: &Crystal,
    grid: &FftGrid,
    pseudopotentials: &HashMap<Element, UpfPseudoPotential>,
    omega: f64,
) -> Result<Vec<Complex64>> {
    let atom_data: Vec<_> = crystal
        .atoms
        .iter()
        .map(|atom| {
            let pp = pseudopotentials.get(&atom.symbol).expect("BUG: PP missing");
            (atom.cart_position(&crystal.lattice), pp)
        })
        .collect();

    let recip = &grid.recip;
    let dims = grid.dims;

    let v_local = (0..grid.total_size())
        .into_par_iter()
        .map(|idx| {
            let g = g_vector_at_dims(idx, dims, recip);
            let g_norm = g.norm();

            // Use fold for a more functional summation
            atom_data.iter().fold(Complex64::default(), |acc, (tau, pp)| {
                let phase = -g.dot(tau);
                acc + Complex64::cis(phase) * pp.v_local_of_g(g_norm, omega)
            })
        })
        .collect();

    Ok(v_local)
}

/// Compute NLCC core density on the real-space FFT grid.
///
/// For a spherically symmetric radial density `ρ_core(r)`, the Fourier
/// transform per unit cell is
///
/// ```text
///     ρ_core(G) = (4π / Ω) · ∫₀^∞ ρ_core(r) · j₀(|G| r) · r² dr · S(G),
/// ```
///
/// where `S(G) = exp(−i G·τ)` is the atomic structure factor.
/// `pp.core_charge` stores the bare `ρ_core(r)` in `e/Å³`
/// (see `PseudopotentialData::core_charge`); the `r²` weight and `4π`
/// prefactor are supplied here.
///
/// Units: `r` and `rab` in Å, `G` in 1/Å, `ρ_core` in `e/Å³`; the
/// Simpson integral has units `e/Å³ · Å² · Å = e`, so `4π·I/Ω` is
/// `e/Å³` (real-space density after the inverse FFT).
///
/// Returns empty vec if no PP has NLCC.
pub(crate) fn compute_core_density(
    crystal: &Crystal,
    grid: &mut FftGrid,
    pseudopotentials: &HashMap<Element, UpfPseudoPotential>,
) -> Vec<f64> {
    if !pseudopotentials.values().any(|pp| pp.has_nlcc()) {
        return vec![];
    }

    let omega = crystal.lattice.volume();
    let four_pi = 4.0 * std::f64::consts::PI;

    // Filter atoms that actually have NLCC to avoid wasted math in the inner loop
    let nlcc_atoms: Vec<_> = crystal
        .atoms
        .iter()
        .filter_map(|atom| {
            let pp = pseudopotentials.get(&atom.symbol)?;
            if pp.has_nlcc() {
                Some((atom.cart_position(&crystal.lattice), pp))
            } else {
                None
            }
        })
        .collect();

    let mut rho_g: Vec<Complex64> = (0..grid.total_size())
        .into_par_iter() // Parallelize the heavy radial integrations
        .map(|idx| {
            let g = grid.g_vector_at(idx);
            let g_norm = g.norm();

            nlcc_atoms.iter().fold(Complex64::default(), |acc, (tau, pp)| {
                // Inline the Simpson integration for clarity or keep in helper
                let integral = integrate_core_charge(pp, g_norm);
                let phase = -g.dot(tau);
                acc + Complex64::cis(phase) * (four_pi * integral / omega)
            })
        })
        .collect();

    grid.fft.inverse(&mut rho_g);
    rho_g.into_iter().map(|c| c.re).collect()
}

// Helper to keep the main loop clean
fn integrate_core_charge(pp: &UpfPseudoPotential, g_norm: f64) -> f64 {
    let integrand: Vec<f64> = pp
        .core_charge
        .iter()
        .zip(&pp.r_grid)
        .map(|(&rho, &r)| {
            let gr = g_norm * r;
            let j0 = if gr < 1e-10 { 1.0 - gr * gr / 6.0 } else { gr.sin() / gr };
            rho * r * r * j0
        })
        .collect();
    crate::numerics::simpson_integrate(&integrand, &pp.rab)
}

/// Assemble the kinetic + local part of the Kohn-Sham Hamiltonian at
/// k-point `k`, fully overwriting a caller-supplied `faer::Mat`.
///
/// Every entry `(i, j)` is written exactly once:
/// ```text
///     H[i, j] = V_eff(G_i − G_j)                       (off-diagonal)
///     H[i, i] = (ℏ²/2m) · |k + G_i|² + V_eff(0)        (diagonal)
/// ```
/// so `h` does **not** need to be pre-zeroed. This matters because the
/// non-local KB term
/// ([`crate::potential::nonlocal::NonlocalPotential::add_to_hamiltonian`])
/// is layered on top via an accumulating `matmul(..., Accum::Add, ...)`:
/// the contract is "fill first, then accumulate". Under the
/// QE-compatible gauge (see `ScfContext::new`), `V_local(G=0)` is
/// kept on the Hamiltonian diagonal — and since `V_H(G=0) = 0` and
/// `V_xc(G=0)` is a real scalar, `V_eff(0) = V_local(G=0) + V_xc(G=0)`
/// is a finite uniform shift that enters every KS eigenvalue.
///
/// `h` must already be sized `basis.len() × basis.len()`; the shape is
/// checked with `debug_assert!`. Kept separate from the public
/// `crate::hamiltonian::build_kinetic` (which constructs the kinetic-only
/// Hamiltonian used by the free-electron band-structure path and does not
/// share this tight inner loop).
///
/// **Allocation contract.** This routine performs no heap allocations.
/// The backing `faer::Mat` is owned by the caller — typically one slot
/// per `(ispin, ik)` in `ScfContext::h_scratch` — so the SCF loop does
/// not allocate an `n_pw × n_pw` `Mat::<Complex64>` each iteration.
/// Callers that want an owned matrix (free-electron band structure,
/// tests) should combine this with `faer::Mat::<Complex64>::zeros(n, n)`
/// inline.
pub(crate) fn fill_hamiltonian_with_v_eff(
    h: &mut faer::Mat<Complex64>,
    basis: &BasisSet,
    k: &Vector3<f64>,
    v_eff_fft: &[Complex64],
    grid_dims: [usize; 3],
) {
    let n = basis.len();
    debug_assert_eq!(h.nrows(), n, "fill_hamiltonian_with_v_eff: row count mismatch");
    debug_assert_eq!(h.ncols(), n, "fill_hamiltonian_with_v_eff: col count mismatch");

    let miller_idx = basis.miller_indices();
    let g_vectors = basis.g_vectors();

    for i in 0..n {
        let ke_i = HBAR2_OVER_2M * (k + g_vectors[i]).norm_squared();
        let mi = miller_idx[i];
        for j in 0..n {
            let mj = miller_idx[j];
            let dn1 = mi[0] - mj[0];
            let dn2 = mi[1] - mj[1];
            let dn3 = mi[2] - mj[2];
            let fft_idx = miller_to_idx(grid_dims, dn1, dn2, dn3);
            let v = v_eff_fft[fft_idx];
            h[(i, j)] = if i == j { Complex64::new(ke_i, 0.0) + v } else { v };
        }
    }
}
