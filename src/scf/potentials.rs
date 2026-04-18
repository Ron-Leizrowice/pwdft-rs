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
    error::{PwdftError, Result},
    pseudopotential::PseudopotentialData,
};

use super::grid::{FftGrid, g_vector_at_dims, miller_to_idx};

/// Compute local pseudopotential V_local(G) on the full FFT grid.
pub(crate) fn compute_v_local(
    crystal: &Crystal,
    grid: &FftGrid,
    pseudopotentials: &[&PseudopotentialData],
    omega: f64,
) -> Result<Vec<Complex64>> {
    let atom_data: Vec<(Vector3<f64>, &PseudopotentialData)> = crystal
        .atoms
        .iter()
        .map(|atom| {
            let pp = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials)
                .ok_or_else(|| PwdftError::MissingPseudopotential(
                    format!("Z={} not found in loaded pseudopotentials", atom.z)
                ))?;
            Ok((atom.cart_position(&crystal.lattice), pp))
        })
        .collect::<Result<Vec<_>>>()?;

    let dims = grid.dims;
    let recip = grid.recip.clone();
    Ok((0..grid.total_size())
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
        .collect())
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
/// prefactor are supplied here. This mirrors QE's `init_tab_rhc` at
/// `qe-7.5/upflib/rhoc_mod.f90:107-115`:
///
/// ```text
///     aux(ir)     = upf%rho_atc(ir) * rgrid%r2(ir) * sin(qr)/(qr)
///     tab_rhc(iq) = fpi * simpson(aux, rab) / omega
/// ```
///
/// Units: `r` and `rab` in Å, `G` in 1/Å, `ρ_core` in `e/Å³`; the
/// Simpson integral has units `e/Å³ · Å² · Å = e`, so `4π·I/Ω` is
/// `e/Å³` (real-space density after the inverse FFT).
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
    let four_pi = 4.0 * std::f64::consts::PI;
    let mut rho_core_g = vec![Complex64::new(0.0, 0.0); n_grid];

    for atom in &crystal.atoms {
        // SAFETY: if we reached this point, ScfContext::new already validated
        // that all atoms have matching pseudopotentials. A missing PP here
        // would be a programming error, not a user input error.
        let Some(pp) = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials) else {
            continue;
        };
        if !pp.has_nlcc() || pp.core_charge.is_empty() {
            continue;
        }

        let tau = atom.cart_position(&crystal.lattice);

        for (idx, rho_g_val) in rho_core_g.iter_mut().enumerate() {
            let g = grid.g_vector_at(idx);
            let g_norm = g.norm();

            // Integrand: ρ_core(r) · r² · j₀(|G| r)
            let integrand: Vec<f64> = pp.core_charge.iter().zip(pp.r_grid.iter())
                .map(|(&rho_c, &r)| {
                    let gr = g_norm * r;
                    let j0 = if gr < 1e-10 {
                        1.0 - gr * gr / 6.0
                    } else {
                        gr.sin() / gr
                    };
                    rho_c * r * r * j0
                })
                .collect();
            let integral = crate::numerics::simpson_integrate(&integrand, &pp.rab);

            let phase = -g.dot(&tau);
            let sf = Complex64::cis(phase);
            *rho_g_val += sf * (four_pi * integral / omega);
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
            // Miller entries are `i16` (TYPE-A). Widen to `i32` before
            // subtraction so differences cannot overflow even at the
            // i16 boundary; `miller_to_idx` takes `i32` natively.
            let dn1 = i32::from(miller_idx[i][0]) - i32::from(miller_idx[j][0]);
            let dn2 = i32::from(miller_idx[i][1]) - i32::from(miller_idx[j][1]);
            let dn3 = i32::from(miller_idx[i][2]) - i32::from(miller_idx[j][2]);
            let fft_idx = miller_to_idx(grid_dims, dn1, dn2, dn3);
            h[(i, j)] += v_eff_fft[fft_idx];
        }
    }

    h
}
