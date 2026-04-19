//! Non-spin-polarized (nspin=1) self-consistent field driver.
//!
//! Solves the Kohn-Sham equations
//! ```text
//!     ( −(ℏ²/2m) ∇² + V_eff[ρ](r) + V_NL ) ψ_{n,k}(r) = ε_{n,k} ψ_{n,k}(r)
//!     V_eff[ρ](r) = V_local(r) + V_H[ρ](r) + V_xc[ρ](r)
//!     ρ(r)        = Σ_{n,k} f_{n,k} · w_k · |ψ_{n,k}(r)|²
//! ```
//! by fixed-point iteration on `ρ`. Each iteration:
//!
//! 1. `V_H(G) = 4πe² · ρ(G) / |G|²` — Poisson in G-space
//!    ([`crate::scf::energy::hartree_on_fft_grid`]).
//! 2. `(ε_xc(r), v_xc(r))` from Perdew-Zunger LDA on
//!    `ρ_val(r) + ρ_core(r)` (NLCC if any PP carries `PP_NLCC`;
//!    Louie, Froyen, Cohen, *Phys. Rev. B* **26**, 1738 (1982)).
//! 3. Assemble `V_eff(G) = V_local(G) + V_H(G) + V_xc(G)` and build
//!    per-k Hamiltonians; add the separable Kleinman-Bylander non-local
//!    operator (Kleinman & Bylander, *Phys. Rev. Lett.* **48**, 1425
//!    (1982)).
//! 4. Diagonalize for the lowest `n_bands` eigenpairs (dense faer or
//!    partial Arnoldi; see [`crate::eigensolver`]).
//! 5. Fermi-Dirac / Methfessel-Paxton / cold occupations satisfying
//!    `Σ_{n,k} f_{n,k} · w_k = N_el / spin_factor`.
//! 6. Reconstruct the output density and symmetrize in G-space via
//!    analytic fractional-translation phase factors (PCFX).
//! 7. Convergence: density RMS plus energy-change threshold. The
//!    diagnostic Harris-Foulkes estimator (stationary in density error,
//!    Harris, *Phys. Rev. B* **31**, 1770 (1985)) is tracked against the
//!    KS total; a large residual at "convergence" flags a premature
//!    threshold.
//! 8. Density mixing (Anderson / Broyden / Periodic Pulay with optional
//!    Kerker preconditioning).
//!
//! On convergence the final pass assembles
//! [`crate::scf::energy::EnergyComponents`] (VGC5 per-term decomposition
//! for validation against QE's `pw.x` output). All internal quantities
//! are in eV / Å / e·Å⁻³; see [`crate::consts`] for the unit factors.
//!
//! ## Density convention for VGC5 energy terms
//!
//! All [`EnergyComponents`] Kohn-Sham terms — `e_hartree`, `e_xc`,
//! `e_kinetic`, `e_local`, `e_nonlocal`, and the assembled `e_total`
//! — are evaluated on the **post-PCFX symmetrized output density**
//! produced at step 6, *not* on the raw band-reconstructed density.
//! Step 6 (`density::compute_density`) returns the raw band sum;
//! `symmetrize_density_g` immediately rewrites it in place, and every
//! energy call that follows (`hartree_energy`, `xc_energy_corrected`,
//! `kinetic_expectation`, `local_pp_energy_grid`, `nonlocal_expectation`)
//! consumes the symmetrized density. The spin driver mirrors this order.
//!
//! A future maintainer who reorders symmetrize-vs-diagonalize, or who
//! inserts an energy evaluation between `compute_density` and
//! `symmetrize_density_g`, would silently break the direct-sum
//! identity documented on [`EnergyComponents`]:
//!
//! ```text
//! e_band = e_kinetic + e_local + e_nonlocal + 2·e_hartree + e_vxc
//! ```
//!
//! Pre-PCFX the symmetry residual on the output density was as large
//! as ~1.2 eV on Si (see proposal `PCFX`, `symmetry::density::g_space`).
//! The Harris-Foulkes estimator is the one deliberate exception: it
//! uses the *input* density (`rho_g` from the previous mix output) and
//! the corresponding input V_xc / ε_xc, because HF's O(Δρ²)
//! stationarity argument is stated about the input density.
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
use super::potentials::fill_hamiltonian_with_v_eff;
use super::report::{log_components, log_convergence_summary, log_entropy, log_iteration, IterationReport};
use super::{ScfParams, ScfResult, context, density, initial_density, mixing, smearing};

/// Dispatch a per-k-point diagonalization to the configured backend.
///
/// Transparent fallback: if the iterative solver fails to converge
/// `n_bands` eigenpairs within its restart budget, this helper emits a
/// `log::warn!` and retries on the dense path. The SCF loop never sees a
/// convergence-style failure from the iterative path — only a genuine
/// panic would escape, which faer's upstream tests exercise heavily.
///
/// ## Warm-start semantics
///
/// Both backends consume the caller-supplied `v_prev` (typically the
/// previous SCF iteration's eigenvectors at the same k-point), but on
/// different algorithmic paths:
///
/// - **Dense** with `wfrx_subspace == true` uses
///   [`dense::diagonalize_subspace`] — Rayleigh–Ritz projection onto the
///   subspace spanned by the columns of `v_prev`, with an accuracy-gate
///   fallback to [`dense::diagonalize_lowest`] if the subspace is not
///   close enough to invariant. With `wfrx_subspace == false`, `v_prev`
///   is ignored.
/// - **Iterative** passes the first column of `v_prev` as the starting
///   vector `v0` of faer's Arnoldi iteration. Without this, the cold
///   deterministic seed produces a *different* Krylov subspace on each
///   SCF iteration and drives the SCF to a different fixed point than
///   Dense (measured: 0.77 eV on Si at n_pw = 89). With the warm start,
///   the Krylov subspace is biased toward the previous iteration's
///   ground-state orbital and converges on the same fixed point as the
///   Dense path to SCF tolerance.
///
/// On the first iteration (`v_prev == None`) the iterative path uses its
/// deterministic internal seed; see
/// [`iterative::diagonalize_lowest_iterative`].
pub(super) fn diagonalize_dispatch(
    h: &faer::Mat<Complex64>,
    n_bands: usize,
    kind: EigensolverKind,
    wfrx_subspace: bool,
    v_prev: Option<&faer::Mat<Complex64>>,
) -> Result<EigenResult> {
    match kind {
        EigensolverKind::Dense => {
            if wfrx_subspace {
                dense::diagonalize_subspace(h, n_bands, v_prev)
            } else {
                dense::diagonalize_lowest(h, n_bands)
            }
        }
        EigensolverKind::Iterative => {
            // Thread the previous iteration's first-column eigenvector
            // through as the Arnoldi starting vector. Building the Vec
            // here (rather than a ColRef) keeps the iterative API's
            // slice signature clean; this is a small copy (O(n_pw))
            // next to the O(n_pw² · max_dim) Arnoldi work.
            let v0_owned: Option<Vec<Complex64>> = v_prev.and_then(|vp| {
                if vp.nrows() == h.nrows() && vp.ncols() >= 1 {
                    let n = vp.nrows();
                    let mut col = Vec::with_capacity(n);
                    for row in 0..n {
                        col.push(vp[(row, 0)]);
                    }
                    Some(col)
                } else {
                    None
                }
            });
            if v0_owned.is_some() {
                log::debug!(
                    "iterative eigensolver: using WFRX warm-start v0 (n_bands={n_bands})"
                );
            } else {
                log::debug!(
                    "iterative eigensolver: cold start (no prev_eigvecs; n_bands={n_bands})"
                );
            }
            let v0_slice = v0_owned.as_deref();
            match iterative::diagonalize_lowest_iterative(
                h,
                n_bands,
                v0_slice,
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

/// Map each eigenvalue to its occupation via the configured smearing
/// scheme.
///
/// ```text
///     f_{n,k} = spin_factor · g((ε_{n,k} − ε_F) / σ)
/// ```
/// with `g` one of Fermi-Dirac, Gaussian, Methfessel-Paxton-N, or cold
/// smearing (selected by `SmearingScheme`; see
/// [`crate::scf::smearing::occupation`]). `σ = sigma` in eV is the
/// broadening width; `ε_F` in eV is the Fermi level located by
/// [`crate::scf::smearing::find_fermi_energy`] so that
/// `Σ_{n,k} f_{n,k} · w_k = N_el / spin_factor`.
/// `spin_factor = 2` in nspin=1 (degenerate channels) and `1` in nspin=2.
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
/// The GPU kernel in `gpu::GpuAccelerator::lda_xc` is LDA-specific
/// (`src/gpu/shaders/lda_xc.wgsl`). When the active functional is
/// `XcEvaluator::Pz` and a GPU is available, we dispatch directly to the
/// f32 kernel. For any other functional, we fall through to the CPU
/// evaluator so the dispatch stays a single `match` on the data enum
/// (no hidden GPU-only override).
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

/// Run the non-spin-polarized (nspin=1) Kohn-Sham SCF loop; see the
/// module header for the equations and per-iteration pipeline.
///
/// On convergence returns a populated [`ScfResult`] with total, Harris-
/// Foulkes, free, and σ→0 energies in eV, the converged eigenvalues,
/// Fermi level, final density in G-space, and the VGC5
/// [`EnergyComponents`] decomposition. `xc_evaluator` is pre-constructed
/// by `run_scf` and held by value here; the driver dispatches XC kernels
/// through it on each iteration.
///
/// Precondition: `params.nspin == 1` (caller dispatches on nspin).
/// Errors with [`crate::error::PwdftError::ConvergenceFailure`] if the
/// density RMS criterion is not met within `max_iter` iterations.
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

    // Initial density: superposition of atomic densities (SAD).
    // CFGN Phase 1: Gaussian width threaded from `ScfParams` (which in turn
    // is threaded from `Settings::initial_density.gaussian_sigma`); defaults
    // to `initial_density::DEFAULT_GAUSSIAN_SIGMA` (1.0 Å) when YAML omits
    // the field, preserving bit-identity with the pre-CFGN hardcoded value.
    let init_config = initial_density::InitialDensityConfig {
        magnetic_moments: vec![0.0; ctx.crystal.atoms.len()],
        gaussian_sigma: Some(ctx.params.gaussian_sigma),
    };
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
    mixer.log_init(
        &ctx.params.mixing_mode,
        ctx.params.mixing_ndim,
        ctx.params.adaptive_beta,
        "",
    );
    let mut eigenvalues_all: Vec<Vec<f64>>;
    let mut fermi_energy;
    let mut e_prev: Option<f64> = None;
    let mut last_delta = f64::INFINITY;
    let pb = scf_progress_bar(ctx.params.max_iter);

    // Warm-start cache: previous iteration's eigenvectors per k-point.
    // `None` on the first iteration → cold start. Populated with length
    // `ctx.kpoints.len()` after each successful eigensolve pass; indexed
    // by `ik` inside the rayon k-point loop so each closure captures its
    // own `&faer::Mat` reference (no shared-mut aliasing).
    //
    // Two backends consume this cache on distinct algorithmic paths:
    //
    // - **Dense + `wfrx_subspace`** (the original WFRX Phase-1 path):
    //   Rayleigh–Ritz projection onto the previous iteration's subspace,
    //   with residual-gate fallback to full diag.
    // - **Iterative**: the first column of `prev_wavefunctions[ik]`
    //   becomes the Arnoldi starting vector `v0`. Without this, cold
    //   Arnoldi finds a different SCF fixed point than Dense (measured
    //   0.77 eV on Si at n_pw = 89 before this was wired).
    //
    // The cache is populated whenever *either* consumer is active. The
    // allocation is per-eigensolve clone; at n_pw = 725, n_bands = 8
    // that is ≈ 90 KiB per k-point, negligible next to the eigensolve.
    let wfrx_dense_enabled = ctx.params.wfrx_subspace
        && matches!(ctx.params.eigensolver, EigensolverKind::Dense);
    let iterative_warmstart_enabled =
        matches!(ctx.params.eigensolver, EigensolverKind::Iterative);
    let cache_prev_wavefunctions = wfrx_dense_enabled || iterative_warmstart_enabled;
    let mut prev_wavefunctions: Option<Vec<faer::Mat<Complex64>>> = None;

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

        // 4. Solve eigenvalue problem at each k-point (parallel over k-points).
        //    WFRX: when warm-start is enabled and we have a previous
        //    iteration's eigenvectors, pass them in per-k as the subspace
        //    seed. `prev_wavefunctions.as_ref()` lives outside the
        //    closure; each closure indexes by `ik` to get its own
        //    `&faer::Mat`, so there is no shared mutable aliasing.
        //
        //    ALOC F-5: `h_scratch` is an owned per-k scratch `Mat`
        //    slot on the context. `fill_hamiltonian_with_v_eff` fully
        //    overwrites every entry — no zero-fill needed — and then
        //    VNL accumulates via `matmul(Accum::Add, ...)` on top. The
        //    `par_iter_mut()` over `h_scratch[..n_k]` gives each rayon
        //    closure its own `&mut Mat` without any `Vec::split_at_mut`
        //    or cell trickery; the remaining `ctx.*` fields are
        //    captured through disjoint immutable borrows.
        let eigensolver_kind = ctx.params.eigensolver;
        let n_bands = ctx.params.n_bands;
        let prev_wfn_ref = prev_wavefunctions.as_ref();
        let basis = ctx.basis;
        let grid_dims = ctx.grid.dims;
        let crystal = ctx.crystal;
        let kpoints = ctx.kpoints;
        let vnl_cache = &ctx.vnl_cache;
        let n_k = kpoints.len();
        // Non-spin SCF uses only the first n_k slots (the spin driver
        // owns its own disjoint slice for the ↓ channel).
        let h_scratch_up = &mut ctx.h_scratch[..n_k];
        let kpoint_results: Result<Vec<_>> = h_scratch_up
            .par_iter_mut()
            .zip(kpoints.par_iter())
            .zip(vnl_cache.par_iter())
            .enumerate()
            .map(|(ik, ((h, kp), vnl))| {
                fill_hamiltonian_with_v_eff(h, basis, &kp.k, &v_eff_fft, grid_dims);
                vnl.add_to_hamiltonian(h, crystal, basis, &kp.k);
                let v_prev_k = prev_wfn_ref.map(|wfns| &wfns[ik]);
                diagonalize_dispatch(
                    h,
                    n_bands,
                    eigensolver_kind,
                    wfrx_dense_enabled,
                    v_prev_k,
                )
            })
            .collect();
        let kpoint_results = kpoint_results?;

        eigenvalues_all = kpoint_results.iter().map(|r| r.eigenvalues.clone()).collect();
        let all_kpoint_wavefns: Vec<_> = kpoint_results.into_iter().map(|r| r.eigenvectors).collect();

        // Populate the warm-start cache for the next iteration. Enabled
        // for either Dense+WFRX (Rayleigh–Ritz subspace) or Iterative
        // (Arnoldi `v0`); the `.clone()` is a deep copy of each per-k
        // faer matrix (~90 KiB at n_pw = 725, n_bands = 8 — negligible
        // next to the eigensolve). The cost is paid only when warm-start
        // is active.
        if cache_prev_wavefunctions {
            prev_wavefunctions = Some(all_kpoint_wavefns.clone());
        }

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
            // XC double-counting ∫ρ_val · V_xc d³r on the OUTPUT density —
            // same integral `xc_energy_corrected` subtracts from E_xc. The
            // valence density `rho_r_new` is integrated (core excluded),
            // mirroring the `rho_val` parameter of `xc_energy_corrected`
            // and appearing as `e_vxc` in the band-sum identity.
            let dvol = ctx.omega / ctx.n_grid as f64;
            let e_vxc_term: f64 = rho_r_new
                .iter()
                .zip(vxc_r_energy.iter())
                .map(|(&rho, &vxc)| rho * vxc * dvol)
                .sum();
            let e_ewald_term = ctx.e_ewald;

            let components = EnergyComponents {
                e_band,
                e_kinetic,
                e_local,
                e_local_g0_shift,
                e_nonlocal,
                e_hartree: e_hartree_term,
                e_xc: e_xc_term,
                e_vxc: e_vxc_term,
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
