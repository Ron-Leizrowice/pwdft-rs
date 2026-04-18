//! Non-spin-polarized (nspin=1) self-consistent field driver.
//!
//! Single hot loop over SCF iterations:
//! Hartree → LDA XC → V_eff → diagonalize at each k-point → Fermi
//! occupations → reconstruct density → symmetrize → check convergence →
//! mix. On convergence, the final pass computes the per-component
//! energy decomposition (VGC5) and returns a populated [`ScfResult`].
//!
//! Shared helpers (`diagonalize_dispatch`, `compute_occupations`,
//! `scf_progress_bar`) are `pub(super)` for reuse by `driver_spin`.

use log::info;
use num_complex::Complex64;
use rayon::prelude::*;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    eigensolver::{EigenResult, EigensolverKind, dense, iterative},
    error::{PwdftError, Result},
    kpoints::KPoint,
    potential::xc::XcEvaluator,
    pseudopotential::PseudopotentialData,
};
#[cfg(feature = "gpu")]
use crate::potential::xc::XcGridResult;

use super::energy::{
    EnergyComponents, add_core_density, assemble_v_eff, band_energy, density_diff,
    density_r_to_g, harris_foulkes_energy, hartree_energy, hartree_on_fft_grid,
    kinetic_expectation, local_pp_energy_grid, nonlocal_expectation,
    real_to_g_space, total_energy, with_g0_shift, xc_energy_bare, xc_energy_corrected,
};
use super::potentials::build_hamiltonian_with_v_eff;
use super::report::{log_components, log_convergence_summary, log_entropy, log_iteration, IterationReport};
use super::{ScfParams, ScfResult, context, density, initial_density, mixing, smearing};

/// Dispatch a per-k-point diagonalization to the configured backend.
///
/// Transparent fallback: if the iterative solver fails to converge
/// `n_bands` eigenpairs within its restart budget, this helper emits a
/// `log::warn!` and retries on the dense path. The SCF loop never sees a
/// convergence-style failure from ITEV — only a genuine panic would
/// escape, which faer's upstream tests exercise heavily.
pub(super) fn diagonalize_dispatch(
    h: &faer::Mat<Complex64>,
    n_bands: usize,
    kind: EigensolverKind,
) -> Result<EigenResult> {
    match kind {
        EigensolverKind::Dense => dense::diagonalize_lowest(h, n_bands),
        EigensolverKind::Iterative => {
            match iterative::diagonalize_lowest_iterative(
                h,
                n_bands,
                None,
                iterative::DEFAULT_TOL,
            ) {
                Ok(r) => Ok(r),
                Err(e) => {
                    // Fires per-k-point per-SCF-iteration on persistent failure;
                    // warn once, downgrade the rest to debug to avoid log spam.
                    static ITERATIVE_FALLBACK_WARNED: std::sync::atomic::AtomicBool =
                        std::sync::atomic::AtomicBool::new(false);
                    if !ITERATIVE_FALLBACK_WARNED.swap(true, std::sync::atomic::Ordering::Relaxed) {
                        log::warn!(
                            "iterative eigensolver failed ({e}); falling back to dense \
                             (subsequent failures on this run will be logged at debug level)"
                        );
                    } else {
                        log::debug!("iterative eigensolver failed ({e}); falling back to dense");
                    }
                    dense::diagonalize_lowest(h, n_bands)
                }
            }
        }
    }
}

/// Compute occupation numbers for all k-points from eigenvalues and Fermi energy.
pub(super) fn compute_occupations(
    eigenvalues: &[Vec<f64>],
    scheme: smearing::SmearingScheme,
    fermi_energy: f64,
    sigma: f64,
    spin_factor: f64,
) -> Vec<Vec<f64>> {
    eigenvalues
        .iter()
        .map(|evs| {
            evs.iter()
                .map(|&e| smearing::occupation(scheme, e, fermi_energy, sigma, spin_factor))
                .collect()
        })
        .collect()
}

/// Evaluate the non-spin XC functional, preferring the GPU fast path for LDA.
///
/// GGAP Phase A wiring: the GPU kernel in `gpu::GpuAccelerator::lda_xc` is
/// LDA-specific (`src/gpu/shaders/lda_xc.wgsl`). When the active functional
/// is `XcEvaluator::Pz` and a GPU is available, we keep the direct f32
/// kernel call — bit-identical to pre-Phase-A output on that code path.
/// For any other functional, we fall through to the CPU evaluator so the
/// dispatch stays a single `match` on the data enum (no hidden GPU-only
/// override). Phase E adds a GPU PBE shader; this helper is the seam where
/// that branch slots in.
#[cfg(feature = "gpu")]
fn eval_xc_with_gpu(
    xc_evaluator: &XcEvaluator,
    rho_r: &[f64],
    rho_grad_r: Option<&[[f64; 3]]>,
    gpu: Option<&crate::gpu::GpuAccelerator>,
) -> Result<XcGridResult> {
    if let (XcEvaluator::Pz, Some(gpu)) = (xc_evaluator, gpu) {
        let (exc_r, v1_r) = gpu.lda_xc(rho_r);
        return Ok(XcGridResult { exc_r, v1_r, v2_r: None });
    }
    xc_evaluator.eval(rho_r, rho_grad_r)
}

pub(super) fn scf_progress_bar(max_iter: usize) -> indicatif::ProgressBar {
    let pb = indicatif::ProgressBar::new(max_iter as u64);
    pb.set_style(
        indicatif::ProgressStyle::with_template(
            "SCF [{bar:30}] {pos}/{len}  {msg}  [{elapsed_precise} elapsed]",
        )
        // SAFETY: This is a static, valid template string -- with_template cannot fail.
        .expect("BUG: invalid progress bar template")
        .progress_chars("##-"),
    );
    pb
}

/// Run the non-spin-polarized (nspin=1) SCF loop.
///
/// Precondition: `params.nspin == 1` (caller dispatches on nspin).
/// `xc_evaluator` is pre-constructed by `run_scf` and held by value here;
/// the driver dispatches XC kernels through it on each iteration.
pub(crate) fn run_scf_unpolarized(
    crystal: &Crystal,
    basis: &BasisSet,
    kpoints: &[KPoint],
    pseudopotentials: &[&PseudopotentialData],
    params: &ScfParams,
    symmetry: &crate::symmetry::SymmetryInfo,
    xc_evaluator: XcEvaluator,
) -> Result<ScfResult> {
    let mut ctx = context::ScfContext::new(
        crystal, basis, kpoints, pseudopotentials, params, symmetry,
    )?;

    // Try to initialize GPU if compiled with gpu feature
    #[cfg(feature = "gpu")]
    let mut gpu = crate::gpu::GpuAccelerator::try_new();
    #[cfg(feature = "gpu")]
    if gpu.is_some() {
        info!("GPU acceleration enabled for grid operations");
    }
    // Prepare persistent GPU buffers if GPU is available
    #[cfg(feature = "gpu")]
    if let Some(ref mut g) = gpu {
        g.prepare_buffers(ctx.n_grid, &ctx.g_squared);
    }

    // Initial density: superposition of atomic densities (SAD)
    let init_config = initial_density::InitialDensityConfig::non_magnetic(ctx.crystal.atoms.len());
    let mut rho_r = initial_density::generate_initial_density(
        ctx.crystal, &mut ctx.grid, ctx.pseudopotentials, ctx.n_electrons, &init_config,
    );
    let mut rho_g = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
    density_r_to_g(&mut ctx.grid.fft, &rho_r, &mut rho_g);
    info!("Initial density: superposition of atomic densities (Gaussian model)");

    let mut mixer = mixing::Mixer::new(
        ctx.params.mixing_beta,
        ctx.params.mixing_ndim,
        &ctx.params.mixing_mode,
        Some(&ctx.g_squared),
        ctx.n_electrons,
        ctx.omega,
        ctx.params.adaptive_beta,
    );
    let mut eigenvalues_all: Vec<Vec<f64>>;
    let mut fermi_energy;
    let mut e_prev: Option<f64> = None;
    let mut last_delta = f64::INFINITY;
    let pb = scf_progress_bar(ctx.params.max_iter);

    for iter in 0..ctx.params.max_iter {
        // Steps 1-3: Hartree, XC, V_eff assembly.
        // GPU path uses f32 for Hartree, XC, and V_eff; CPU path uses f64 + rayon.
        // XC FFT normalization is shared between both paths.

        // 1. Hartree potential
        #[cfg(feature = "gpu")]
        let v_h_fft = if let Some(ref gpu) = gpu {
            let fourpi_e2 = 4.0 * std::f64::consts::PI * crate::consts::E2_COULOMB;
            gpu.hartree_potential(&rho_g, &ctx.g_squared, fourpi_e2)
        } else {
            hartree_on_fft_grid(&rho_g, &ctx.g_squared)
        };
        #[cfg(not(feature = "gpu"))]
        let v_h_fft = hartree_on_fft_grid(&rho_g, &ctx.g_squared);

        // 2. XC potential: compute in real space, FFT to G-space.
        // NLCC: add core density to valence density for XC evaluation.
        //
        // GGAP Phase A: route LDA through `xc_evaluator.eval` so the
        // dispatcher is the single entry point for every functional. The
        // GPU fast path stays wired to the `Pz` variant via a direct
        // `gpu.lda_xc` call — same numbers as pre-Phase-A — because the
        // GPU kernel is LDA-specific. Once Phase E adds a GPU PBE shader
        // this branch extends; the CPU side of the enum carries the shape.
        let rho_for_xc = add_core_density(&rho_r, &ctx.rho_core_r);
        #[cfg(feature = "gpu")]
        let xc_in = eval_xc_with_gpu(
            &xc_evaluator, &rho_for_xc, None, gpu.as_ref(),
        )?;
        #[cfg(not(feature = "gpu"))]
        let xc_in = xc_evaluator.eval(&rho_for_xc, None)?;
        let exc_r_in = xc_in.exc_r;
        let vxc_r = xc_in.v1_r;
        // v2_r stays None for LDA — Phase B wires it into the semilocal
        // divergence assembly. No change in behaviour in Phase A.
        debug_assert!(
            xc_in.v2_r.is_none(),
            "GGAP Phase A: LDA eval must produce v2_r = None"
        );

        let vxc_g = real_to_g_space(&vxc_r, &mut ctx.grid.fft);

        // 3. V_eff = V_local + V_H + V_xc
        #[cfg(feature = "gpu")]
        let v_eff_fft = if let Some(ref gpu) = gpu {
            gpu.v_eff_assembly(&ctx.v_local_fft, &v_h_fft, &vxc_g)
        } else {
            assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_g)
        };
        #[cfg(not(feature = "gpu"))]
        let v_eff_fft = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_g);

        // 4. Solve eigenvalue problem at each k-point (parallel over k-points)
        let eigensolver_kind = ctx.params.eigensolver;
        let kpoint_results: Result<Vec<_>> = ctx.kpoints
            .par_iter()
            .enumerate()
            .map(|(ik, kp)| {
                let mut h = build_hamiltonian_with_v_eff(ctx.basis, &kp.k, &v_eff_fft, ctx.grid.dims);
                ctx.vnl_cache[ik].add_to_hamiltonian(&mut h, ctx.crystal, ctx.basis, &kp.k);
                diagonalize_dispatch(&h, ctx.params.n_bands, eigensolver_kind)
            })
            .collect();
        let kpoint_results = kpoint_results?;

        eigenvalues_all = kpoint_results.iter().map(|r| r.eigenvalues.clone()).collect();
        let all_kpoint_wavefns: Vec<_> = kpoint_results.into_iter().map(|r| r.eigenvectors).collect();

        // 5. Occupations (configurable smearing scheme)
        fermi_energy = smearing::find_fermi_energy(
            &eigenvalues_all,
            &ctx.kpt_weights,
            ctx.n_electrons,
            ctx.params.smearing_sigma,
            ctx.params.smearing_scheme,
            ctx.spin_factor,
        );
        let occupations = compute_occupations(
            &eigenvalues_all, ctx.params.smearing_scheme, fermi_energy,
            ctx.params.smearing_sigma, ctx.spin_factor,
        );

        // 6. New density
        let mut rho_r_new = density::compute_density(
            &mut density::DensityGrid {
                basis: ctx.basis, g_to_fft: &ctx.g_to_fft, fft: &mut ctx.grid.fft,
                n_electrons: ctx.n_electrons, omega: ctx.omega,
            },
            ctx.kpoints, &all_kpoint_wavefns, &occupations,
        );

        // 6b. Symmetrize density in G-space (PCFX). The G-space form applies
        //     each space-group fractional translation τ_S as an analytic
        //     phase factor `exp(-i G · τ_S)`, so it is exact for non-
        //     symmorphic groups on arbitrary FFT grids. The real-space
        //     form used pre-PCFX rounded τ to the nearest grid point and
        //     smeared density across neighbours whenever `N_i · τ_i` was
        //     not integer (e.g. Fd-3m τ=(¼,¼,¼) on an 18³ grid),
        //     producing a ~1.2 eV per-component residual on Si.
        //     `symmetrize_density_g` short-circuits for identity-only
        //     groups (bit-identical to the pre-PCFX skip).
        crate::symmetry::density::symmetrize_density_g(
            &mut rho_r_new,
            ctx.grid.dims,
            &mut ctx.grid.fft,
            ctx.symmetry,
        );

        // 7. Convergence check (dual criterion: density AND energy)
        let delta = density_diff(&rho_r, &rho_r_new, ctx.omega, ctx.n_grid);
        last_delta = delta;

        // Compute energy every iteration for convergence monitoring
        let mut rho_g_new = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
        density_r_to_g(&mut ctx.grid.fft, &rho_r_new, &mut rho_g_new);

        let rho_new_for_xc = add_core_density(&rho_r_new, &ctx.rho_core_r);
        // Recompute XC from the OUTPUT density for the Kohn-Sham total
        // energy (see SPXC/VGC5 rationale in `scf::energy`). Routed through
        // the same evaluator as the input-density path so a future PBE
        // swap cannot leave one site on LDA.
        let xc_out = xc_evaluator.eval(&rho_new_for_xc, None)?;
        let exc_r = xc_out.exc_r;
        let vxc_r_energy = xc_out.v1_r;

        let e_band = band_energy(&eigenvalues_all, &occupations, &ctx.kpt_weights);

        // Kohn-Sham energy: double-counting from OUTPUT density
        let e_total = with_g0_shift(total_energy(
            e_band,
            hartree_energy(&rho_g_new, &ctx.g_squared, ctx.omega),
            xc_energy_corrected(&rho_new_for_xc, &rho_r_new, &exc_r, &vxc_r_energy, ctx.omega),
            ctx.e_ewald,
        ), &ctx);

        // Harris-Foulkes energy: double-counting from INPUT density
        // rho_g, rho_for_xc, exc_r_in, vxc_r are all from the input density
        let e_harris = with_g0_shift(harris_foulkes_energy(
            e_band,
            hartree_energy(&rho_g, &ctx.g_squared, ctx.omega),
            xc_energy_corrected(&rho_for_xc, &rho_r, &exc_r_in, &vxc_r, ctx.omega),
            ctx.e_ewald,
        ), &ctx);

        let hf_diff = (e_harris - e_total).abs();

        let de = e_prev.map(|ep| (e_total - ep).abs());
        e_prev = Some(e_total);

        let rho_converged = delta < ctx.params.conv_threshold;
        let energy_converged = de.is_some_and(|de| de < ctx.params.energy_threshold);

        log_iteration(&pb, &IterationReport {
            iter,
            e_total,
            e_harris,
            hf_diff,
            de,
            delta,
            beta: ctx.params.adaptive_beta.then(|| mixer.current_beta()),
            spin: None,
        });

        // Warn if density converged but Harris-Foulkes difference is large
        if rho_converged && energy_converged && hf_diff > 0.01 {
            log::warn!(
                "Density converged but |E_HF - E_KS| = {hf_diff:.2e} eV — \
                 energy may not be reliable. Consider tightening conv_threshold."
            );
        }

        if rho_converged && energy_converged {
            pb.finish_and_clear();
            info!("SCF converged after {} iterations", iter + 1);
            rho_g = rho_g_new;
            let rho_g_basis: Vec<Complex64> = ctx.g_to_fft.iter().map(|&idx| rho_g[idx]).collect();

            // Entropy and free energy
            let ts = smearing::entropy_ts(
                &eigenvalues_all, &ctx.kpt_weights, fermi_energy,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );
            let free_energy = e_total - ts;
            let energy_sigma0 = f64::midpoint(e_total, free_energy);

            // -------------------------------------------------------------
            // VGC5 per-component decomposition (diagnostic).
            // Computed once at convergence; final wavefunctions are alive here.
            // -------------------------------------------------------------
            let k_vecs: Vec<nalgebra::Vector3<f64>> =
                ctx.kpoints.iter().map(|kp| kp.k).collect();

            let e_kinetic = kinetic_expectation(
                ctx.basis, &k_vecs, &ctx.kpt_weights,
                &all_kpoint_wavefns, &occupations,
            );

            // V_local in real space (G=0 already zeroed in v_local_fft).
            let mut v_local_cplx = ctx.v_local_fft.clone();
            ctx.grid.fft.inverse(&mut v_local_cplx);
            let v_local_r: Vec<f64> = v_local_cplx.iter().map(|c| c.re).collect();
            let e_local = local_pp_energy_grid(&rho_r_new, &v_local_r, ctx.omega, ctx.n_grid);
            let e_local_g0_shift = ctx.v_local_g0 * ctx.n_electrons;

            let e_nonlocal = nonlocal_expectation(
                ctx.basis, ctx.crystal, &k_vecs, &ctx.kpt_weights,
                &all_kpoint_wavefns, &occupations, &ctx.vnl_cache,
            );

            // NB: `rho_g` was assigned `rho_g_new` in the convergence branch
            // above; use `rho_g` here (the final-iteration density in G-space).
            let e_hartree_term = hartree_energy(&rho_g, &ctx.g_squared, ctx.omega);
            let e_xc_term = xc_energy_bare(&rho_new_for_xc, &exc_r, ctx.omega);
            let e_ewald_term = ctx.e_ewald;

            let components = EnergyComponents {
                e_band,
                e_kinetic,
                e_local,
                e_local_g0_shift,
                e_nonlocal,
                e_hartree: e_hartree_term,
                e_xc: e_xc_term,
                e_ewald: e_ewald_term,
            };

            log_convergence_summary(e_total, e_harris, hf_diff, free_energy, energy_sigma0);
            log_entropy(ts, ctx.crystal.atoms.len());
            log_components(&components, e_total, ctx.n_electrons);

            return Ok(ScfResult {
                total_energy: e_total,
                harris_foulkes_energy: e_harris,
                free_energy,
                energy_sigma0,
                entropy_ts: ts,
                eigenvalues: eigenvalues_all,
                fermi_energy,
                n_iterations: iter + 1,
                final_delta: last_delta,
                rho_g: rho_g_basis,
                magnetization: 0.0,
                nspin: 1,
                components,
            });
        }

        // 8. Mix (Kerker preconditioning applied inside if enabled)
        rho_r = mixer.mix(&rho_r, &rho_r_new, &mut ctx.grid.fft);
        density_r_to_g(&mut ctx.grid.fft, &rho_r, &mut rho_g);
    }

    pb.abandon_with_message("did not converge");
    Err(PwdftError::ConvergenceFailure {
        iterations: ctx.params.max_iter,
        delta: last_delta,
    })
}
