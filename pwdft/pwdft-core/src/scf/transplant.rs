//! VGCH-2 Part B diagnostic — single-iteration SCF replica seeded from
//! a pre-computed G-space density.
//!
//! Rationale: the heavy-atom VGCH-2 residual (Cu 16.6 eV, Fe 11.5 eV, etc.)
//! shows a linear-response Δone-e / −ΔE_H signature across the per-term
//! decomposition — pwdft-core and QE converge to different self-consistent
//! densities within each code's numerical basin. VGCH-2 Part B tests the
//! "different basin" hypothesis (H3) by transplanting QE's converged
//! density ρ_QE(G) into pwdft-core, running **exactly one** SCF iteration
//! (build V_eff from ρ_QE, assemble and diagonalize H, reconstruct ρ_out,
//! compute per-term energies), and comparing each term against QE's
//! converged per-term output.
//!
//! - If per-term values match QE at iter-1 within ~1 meV/term, the driver
//!   reproduces QE's decomposition at QE's fixed point. The 16.6 eV
//!   residual at pwdft-core self-consistency is then mixer-basin (H3).
//! - If per-term values do NOT match at iter-1 even with the transplanted
//!   density, H3 is cleared and the bug is finer-grained (H_{G,G'}
//!   assembly, ψ_nk reconstruction, occupation smearing, or a subtle v_xc
//!   / v_H term only visible on heavy-atom density structure).
//!
//! This module is a `#[doc(hidden)]` public API: it is not part of
//! pwdft-core's engineering surface. The routine mirrors the first
//! iteration of the non-spin SCF driver line-for-line, so the
//! diagnostic measures the same code path production SCF executes.
//! Any behavioral drift between this module and the driver at iter 0
//! would be a bug to fix; both sites intentionally share the same
//! helpers (`add_core_density`, `hartree_on_fft_grid`, `assemble_v_eff`,
//! `fill_hamiltonian_with_v_eff`, etc.) to keep them in sync.

use num_complex::Complex64;
use rayon::prelude::*;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    eigensolver::EigensolverKind,
    error::{PwdftError, Result},
    fft::compute_density_gradient,
    kpoints::KPoint,
    potential::xc::{XcEvaluator, assemble_semilocal_vxc},
    pseudopotential::PseudopotentialData,
};

use super::context::ScfContext;
use super::driver::{compute_occupations, diagonalize_dispatch};
use super::energy::{
    EnergyComponents, add_core_density, assemble_v_eff, band_energy, density_r_to_g,
    harris_foulkes_energy, hartree_energy, hartree_on_fft_grid, kinetic_expectation,
    local_pp_energy_grid, nonlocal_expectation, real_to_g_space, total_energy,
    xc_energy_bare, xc_energy_corrected,
};
use super::potentials::fill_hamiltonian_with_v_eff;
use super::{density, smearing, ScfParams};

/// Result of a single transplanted iteration. Field units:
/// all energies in eV, densities in e/Å³ on the FFT grid.
#[derive(Debug)]
#[doc(hidden)]
pub struct TransplantIter1Result {
    /// Total Mermin free energy `F = E − TS` evaluated with the OUTPUT
    /// density (same accounting as `ScfResult::total_energy`).
    pub total_energy: f64,
    /// Harris-Foulkes estimator (paired with INPUT density = ρ_transplant).
    pub harris_foulkes_energy: f64,
    /// Fermi level from this iteration's occupations.
    pub fermi_energy: f64,
    /// Entropy `TS` in eV (positive).
    pub entropy_ts: f64,
    /// Eigenvalues per-k: `[ik][nb]`.
    pub eigenvalues: Vec<Vec<f64>>,
    /// VGC5 per-term breakdown (see [`EnergyComponents`] for exact
    /// definitions). The OUTPUT density is used throughout.
    pub components: EnergyComponents,
    /// Output density on the real-space FFT grid (after PCFX).
    pub rho_r_out: Vec<f64>,
    /// Output density on the FFT-grid G-space ordering.
    pub rho_g_out: Vec<Complex64>,
    /// Input density re-materialised for caller-side sanity checks. Same
    /// FFT-grid real-space layout as `rho_r_out`.
    pub rho_r_in: Vec<f64>,
    /// `Δρ_rms` between input and output densities (e/Å³).
    pub delta_rho: f64,
}

/// VGCH-2 Part B — run one SCF iteration seeded from a transplanted
/// G-space density.
///
/// `rho_g_fft_in` must be a length-`n_grid` complex array on the dense
/// FFT grid with the pwdft-core ordering (see `scf::grid::g_vector_at_dims`).
/// Units are e/Å³ and the Fourier convention is
/// `ρ(G) = (1/N) Σ_r ρ(r) exp(−iG·r)` — the same convention QE uses via
/// `fwfft('Rho', …, dfftp)` (verified in `qe-7.5/PW/src/v_of_rho.f90`)
/// and the same convention pwdft-core' `density_r_to_g` applies. Callers
/// transplanting QE's densities must convert units (multiply by
/// `1/BOHR_TO_ANG³`) and build the dense-grid layout via Miller-index
/// lookup through `miller_to_idx`.
///
/// The implementation copy-pastes the first iteration of the non-spin
/// SCF driver with one change — the initial density comes from
/// `rho_g_fft_in` instead of SAD. All shared helpers are called
/// identically so the iteration is bit-identical to what production
/// SCF would do on its first step from the same density.
///
/// # Errors
/// - `PwdftError::InvalidInput` if `rho_g_fft_in.len() != ctx.n_grid`.
/// - Any error propagated from the eigensolver, XC evaluator, or the
///   per-calculation context constructor.
#[doc(hidden)]
pub fn run_scf_iter1_from_rho_g_fft(
    crystal: &Crystal,
    basis: &BasisSet,
    kpoints: &[KPoint],
    pseudopotentials: &[&PseudopotentialData],
    params: &ScfParams,
    symmetry: &crate::symmetry::SymmetryInfo,
    rho_g_fft_in: &[Complex64],
) -> Result<TransplantIter1Result> {
    params.validate()?;
    if params.nspin != 1 {
        return Err(PwdftError::InvalidInput(
            "VGCH-2B transplant diagnostic only supports nspin=1".into(),
        ));
    }
    let xc_evaluator = XcEvaluator::from_settings(params.xc_functional)?;

    let mut ctx = ScfContext::new(crystal, basis, kpoints, pseudopotentials, params, symmetry)?;

    if rho_g_fft_in.len() != ctx.n_grid {
        return Err(PwdftError::InvalidInput(format!(
            "transplanted rho_g_fft has {} entries; expected n_grid={}",
            rho_g_fft_in.len(),
            ctx.n_grid
        )));
    }

    // Materialise ρ_in on the real-space grid. pwdft-core' inverse FFT is
    // unnormalised (matches `density_r_to_g`'s 1/N forward), so
    // `ρ_r(r) = Σ_G ρ_g(G) exp(+i G·r) = N · IFFT(ρ_g)` produces the
    // expected real-space density when ρ_g is in the `1/N`-forward
    // convention.
    let mut rho_g: Vec<Complex64> = rho_g_fft_in.to_vec();
    let mut rho_r_complex: Vec<Complex64> = rho_g.clone();
    ctx.grid.fft.inverse(&mut rho_r_complex);
    let rho_r: Vec<f64> = rho_r_complex.iter().map(|c| c.re).collect();

    // Re-forward-FFT to get a consistent rho_g that pairs bit-exactly
    // with rho_r under the pwdft-core 1/N convention. The caller's
    // rho_g_fft_in might have tiny numerical drift from an exact
    // (forward ∘ inverse) trip on our FFT; this re-round brings the pair
    // into the same state as the production driver's iter-N state.
    density_r_to_g(&mut ctx.grid.fft, &rho_r, &mut rho_g);

    // Keep a snapshot of the input density for caller diagnostics.
    let rho_r_in = rho_r.clone();

    // --- Begin replica of `driver::run_scf_unpolarized` iter 0 ---

    // 1. Hartree potential from INPUT density.
    let v_h_fft = hartree_on_fft_grid(&rho_g, &ctx.g_squared);

    // 2. XC from INPUT density (with NLCC if present).
    let rho_for_xc = add_core_density(&rho_r, &ctx.rho_core_r);
    let rho_grad_for_xc: Option<Vec<[f64; 3]>> = xc_evaluator.needs_gradient().then(|| {
        let mut grad_val = compute_density_gradient(&rho_r, &mut ctx.grid.fft, &ctx.g_vectors);
        if let Some(ref grad_core) = ctx.rho_core_grad_r {
            for (g, gc) in grad_val.iter_mut().zip(grad_core.iter()) {
                g[0] += gc[0];
                g[1] += gc[1];
                g[2] += gc[2];
            }
        }
        grad_val
    });
    let xc_in = xc_evaluator.eval(&rho_for_xc, rho_grad_for_xc.as_deref())?;
    let exc_r_in = xc_in.exc_r;
    let vxc_r = if let Some(ref h_r) = xc_in.v2_r {
        assemble_semilocal_vxc(&xc_in.v1_r, h_r, &mut ctx.grid.fft, &ctx.g_vectors)
    } else {
        xc_in.v1_r
    };
    let vxc_g = real_to_g_space(&vxc_r, &mut ctx.grid.fft);

    // 3. V_eff = V_local + V_H + V_xc
    let v_eff_fft = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_g);

    // 4. Diagonalize per-k. Mirrors driver's `par_iter_mut()` loop.
    let eigensolver_kind = ctx.params.eigensolver;
    let n_bands = ctx.params.n_bands;
    let basis_ref = ctx.basis;
    let grid_dims = ctx.grid.dims;
    let crystal_ref = ctx.crystal;
    let kpoints_ref = ctx.kpoints;
    let vnl_cache = &ctx.vnl_cache;
    let n_k = kpoints_ref.len();
    let h_scratch_up = &mut ctx.h_scratch[..n_k];
    // Iter-0 from transplant has no prior wavefunctions; dense subspace
    // path with `v_prev == None` degrades to a cold full diag, and the
    // iterative path uses its deterministic seed.
    let kpoint_results: Result<Vec<_>> = h_scratch_up
        .par_iter_mut()
        .zip(kpoints_ref.par_iter())
        .zip(vnl_cache.par_iter())
        .map(|((h, kp), vnl)| {
            fill_hamiltonian_with_v_eff(h, basis_ref, &kp.k, &v_eff_fft, grid_dims);
            vnl.add_to_hamiltonian(h, crystal_ref, basis_ref, &kp.k);
            diagonalize_dispatch(
                h,
                n_bands,
                eigensolver_kind,
                matches!(eigensolver_kind, EigensolverKind::Dense) && ctx.params.wfrx_subspace,
                None,
            )
        })
        .collect();
    let kpoint_results = kpoint_results?;
    let eigenvalues_all: Vec<Vec<f64>> =
        kpoint_results.iter().map(|r| r.eigenvalues.clone()).collect();
    let all_kpoint_wavefns: Vec<_> =
        kpoint_results.into_iter().map(|r| r.eigenvectors).collect();

    // 5. Occupations + Fermi level
    let fermi_energy = smearing::find_fermi_energy(
        &eigenvalues_all,
        &ctx.kpt_weights,
        ctx.n_electrons,
        ctx.params.smearing_sigma,
        ctx.params.smearing_scheme,
        ctx.spin_factor,
    );
    let occupations = compute_occupations(
        &eigenvalues_all,
        ctx.params.smearing_scheme,
        fermi_energy,
        ctx.params.smearing_sigma,
        ctx.spin_factor,
    );

    // 6. Output density + PCFX
    let mut rho_r_new = density::compute_density(
        &mut density::DensityGrid {
            basis: ctx.basis,
            g_to_fft: &ctx.g_to_fft,
            fft: &mut ctx.grid.fft,
            n_electrons: ctx.n_electrons,
            omega: ctx.omega,
        },
        ctx.kpoints,
        &all_kpoint_wavefns,
        &occupations,
    );
    crate::symmetry::density::symmetrize_density_g(
        &mut rho_r_new,
        ctx.grid.dims,
        &mut ctx.grid.fft,
        ctx.symmetry,
    );

    // 7. ρ-diff, E_band, entropy
    let dvol_total: f64 = ctx.omega / ctx.n_grid as f64;
    let mut sum_sq = 0.0f64;
    for (a, b) in rho_r.iter().zip(rho_r_new.iter()) {
        let d = a - b;
        sum_sq += d * d * dvol_total;
    }
    let delta_rho = (sum_sq / ctx.omega).sqrt();

    let mut rho_g_new = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
    density_r_to_g(&mut ctx.grid.fft, &rho_r_new, &mut rho_g_new);

    let rho_new_for_xc = add_core_density(&rho_r_new, &ctx.rho_core_r);
    let rho_grad_out: Option<Vec<[f64; 3]>> = xc_evaluator.needs_gradient().then(|| {
        let mut grad_val =
            compute_density_gradient(&rho_r_new, &mut ctx.grid.fft, &ctx.g_vectors);
        if let Some(ref grad_core) = ctx.rho_core_grad_r {
            for (g, gc) in grad_val.iter_mut().zip(grad_core.iter()) {
                g[0] += gc[0];
                g[1] += gc[1];
                g[2] += gc[2];
            }
        }
        grad_val
    });
    let xc_out = xc_evaluator.eval(&rho_new_for_xc, rho_grad_out.as_deref())?;
    let exc_r = xc_out.exc_r;
    let vxc_r_energy = if let Some(ref h_r) = xc_out.v2_r {
        assemble_semilocal_vxc(&xc_out.v1_r, h_r, &mut ctx.grid.fft, &ctx.g_vectors)
    } else {
        xc_out.v1_r
    };

    let e_band = band_energy(&eigenvalues_all, &occupations, &ctx.kpt_weights);
    let ts = smearing::entropy_ts(
        &eigenvalues_all,
        &ctx.kpt_weights,
        fermi_energy,
        ctx.params.smearing_sigma,
        ctx.params.smearing_scheme,
        ctx.spin_factor,
    );
    let e_smearing = -ts;

    // Post-VGCH-SiEF-B1: V_loc(G=0) is kept on the Hamiltonian diagonal
    // (QE gauge), so the band sum already picks up its contribution and
    // no external `with_g0_shift` is applied here.
    let e_total = total_energy(
        e_band,
        hartree_energy(&rho_g_new, &ctx.g_squared, ctx.omega),
        xc_energy_corrected(&rho_new_for_xc, &rho_r_new, &exc_r, &vxc_r_energy, ctx.omega),
        ctx.e_ewald,
        e_smearing,
    );
    let e_harris = harris_foulkes_energy(
        e_band,
        hartree_energy(&rho_g, &ctx.g_squared, ctx.omega),
        xc_energy_corrected(&rho_for_xc, &rho_r, &exc_r_in, &vxc_r, ctx.omega),
        ctx.e_ewald,
        e_smearing,
    );

    // VGC5 per-term decomposition on the OUTPUT density.
    let k_vecs: Vec<nalgebra::Vector3<f64>> = ctx.kpoints.iter().map(|kp| kp.k).collect();
    let e_kinetic = kinetic_expectation(
        ctx.basis,
        &k_vecs,
        &ctx.kpt_weights,
        &all_kpoint_wavefns,
        &occupations,
    );
    let mut v_local_cplx = ctx.v_local_fft.clone();
    ctx.grid.fft.inverse(&mut v_local_cplx);
    let v_local_r: Vec<f64> = v_local_cplx.iter().map(|c| c.re).collect();
    let e_local = local_pp_energy_grid(&rho_r_new, &v_local_r, ctx.omega, ctx.n_grid);
    // Post-VGCH-SiEF-B1: V_loc(G=0) lives on the Hamiltonian diagonal,
    // so ∫ρ·V_local already includes the uniform-background piece and
    // `e_local_g0_shift` is zero by construction.
    let e_local_g0_shift = 0.0;
    let e_nonlocal = nonlocal_expectation(
        ctx.basis,
        ctx.crystal,
        &k_vecs,
        &ctx.kpt_weights,
        &all_kpoint_wavefns,
        &occupations,
        &ctx.vnl_cache,
    );
    let e_hartree_term = hartree_energy(&rho_g_new, &ctx.g_squared, ctx.omega);
    let e_xc_term = xc_energy_bare(&rho_new_for_xc, &exc_r, ctx.omega);
    let dvol = ctx.omega / ctx.n_grid as f64;
    let e_vxc_term: f64 = rho_r_new
        .iter()
        .zip(vxc_r_energy.iter())
        .map(|(&rho, &vxc)| rho * vxc * dvol)
        .sum();

    let components = EnergyComponents {
        e_band,
        e_kinetic,
        e_local,
        e_local_g0_shift,
        e_nonlocal,
        e_hartree: e_hartree_term,
        e_xc: e_xc_term,
        e_vxc: e_vxc_term,
        e_ewald: ctx.e_ewald,
        e_smearing,
    };

    // Move rho_r (the INPUT density, consumed only in construction;
    // rho_r_in already holds a clone we hand back) out of scope so the
    // return value only references the OUTPUT density.
    let _ = rho_r;

    Ok(TransplantIter1Result {
        total_energy: e_total,
        harris_foulkes_energy: e_harris,
        fermi_energy,
        entropy_ts: ts,
        eigenvalues: eigenvalues_all,
        components,
        rho_r_out: rho_r_new,
        rho_g_out: rho_g_new,
        rho_r_in,
        delta_rho,
    })
}
