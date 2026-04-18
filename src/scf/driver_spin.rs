//! Spin-polarized (nspin=2) self-consistent field driver (LSDA).
//!
//! Solves the spin-polarized Kohn-Sham equations
//! ```text
//!     ( −(ℏ²/2m) ∇² + V_eff^σ[ρ_↑, ρ_↓] + V_NL ) ψ_{n,k,σ} = ε_{n,k,σ} ψ_{n,k,σ}
//!     V_eff^σ[ρ_↑, ρ_↓](r) = V_local(r) + V_H[ρ_↑ + ρ_↓](r) + V_xc^σ[ρ_↑, ρ_↓](r)
//!     ρ_σ(r) = Σ_{n,k} f_{n,k,σ} · w_k · |ψ_{n,k,σ}(r)|²        (σ ∈ {↑, ↓})
//! ```
//! Hartree and `V_local` depend on the total density `ρ = ρ_↑ + ρ_↓`
//! (spin-independent operators); LSDA exchange-correlation `V_xc^σ`
//! depends on both channels through the local spin polarization
//! `ζ(r) = (ρ_↑ − ρ_↓) / ρ` (see
//! [`crate::potential::xc::lda_xc_spin_grid`]; Perdew-Zunger LSDA
//! interpolation). Under NLCC the frozen core charge is spin-unpolarized
//! and added as `ρ_core / 2` to each channel before evaluating `V_xc^σ`.
//!
//! ## CCMX: coupled-channel mixer on (ρ_total, m)
//!
//! Mixing `(ρ_↑, ρ_↓)` as two independent DIIS/Broyden streams decouples
//! two channels that are physically coupled through `V_H` and `V_xc`;
//! the independent histories let the channels oscillate in antiphase
//! without ever seeing each other's residuals, stalling convergence on
//! magnetic systems.
//!
//! CCMX (following QE's `rhoz_or_updw`,
//! `qe-7.5/PW/src/scf_mod.f90:1360-1414`) mixes in the change-of-basis
//! ```text
//!     ρ_total = ρ_↑ + ρ_↓,      m = ρ_↑ − ρ_↓
//! ```
//! the "charge + magnetization" representation that diagonalizes the
//! symmetric / antisymmetric part of the Hartree + XC response. The
//! mixed quantities are inverted back via
//! `ρ_↑ = (ρ_total + m) / 2,  ρ_↓ = (ρ_total − m) / 2`. Each basis
//! channel has its own DIIS/Broyden history, but the residuals couple
//! through a well-conditioned 2×2 rotation so the two physical channels
//! no longer drift independently.
//!
//! Kerker preconditioning
//! `K(G) = |G|² / (|G|² + q_TF²)` is the Thomas-Fermi charge-charge
//! response and has no meaning for the magnetization channel (no
//! long-wavelength charge-sloshing mode to damp). CCMX therefore applies
//! Kerker only to the `ρ_total` mixer and runs the `m` mixer with Kerker
//! disabled; the downgrade is logged once per run.
//!
//! ## Convergence (SPNC)
//!
//! Per-spin RMS difference `max(Δρ_↑, Δρ_↓)` is compared to
//! `conv_threshold` — strictly stronger than the total-density
//! difference, and essential: a spin-flip fluctuation `(+ε, −ε)` is
//! invisible to the total but leaves `ζ` inconsistent between input and
//! output, spoiling the O(Δρ²) Harris-Foulkes convergence.
//!
//! Shared primitives (`diagonalize_dispatch`, `compute_occupations`,
//! `scf_progress_bar`) are imported from `scf::driver`. All energies in
//! eV; densities in e/Å³; magnetization `m(r)` in e/Å³ (its spatial
//! integral `∫m(r)d³r` is the total magnetic moment in units of μ_B).

use log::info;
use num_complex::Complex64;
use rayon::prelude::*;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    error::{PwdftError, Result},
    kpoints::KPoint,
    potential::xc::{self, XcEvaluator},
    pseudopotential::PseudopotentialData,
};

use super::driver::{compute_occupations, diagonalize_dispatch, scf_progress_bar};
use super::energy::{
    EnergyComponents, add_core_density, assemble_v_eff, band_energy, density_diff,
    density_r_to_g, harris_foulkes_energy, hartree_energy, hartree_on_fft_grid,
    kinetic_expectation, local_pp_energy_grid, nonlocal_expectation, real_to_g_space,
    total_energy, with_g0_shift, xc_energy_bare,
};
use super::potentials::build_hamiltonian_with_v_eff;
use super::report::{log_components, log_convergence_summary, log_iteration, IterationReport, SpinIterationFields};
use super::{ScfParams, ScfResult, context, density, initial_density, mixing, smearing};

/// Run the spin-polarized (nspin=2) LSDA SCF loop; see the module header
/// for the equations, CCMX coupled-channel mixer, and SPNC per-spin
/// convergence criterion.
///
/// Two spin channels `σ ∈ {↑, ↓}` with independent densities, LSDA XC
/// potentials, and k-point Hamiltonians. Hartree and `V_local` are
/// spin-independent (functions of `ρ_↑ + ρ_↓`); NLCC core charge splits
/// evenly (`ρ_core / 2` per channel). Total magnetization
/// `M = ∫(ρ_↑ − ρ_↓) d³r` in μ_B; can be optimized freely or constrained
/// via `ScfParams::tot_magnetization`. `xc_evaluator` is the data-enum
/// XC dispatcher built by `run_scf` (currently LSDA Perdew-Zunger only;
/// GGA spin support lands in GGAP Phase C).
///
/// Precondition: `params.nspin == 2` (caller dispatches). Returns a
/// populated [`ScfResult`] with per-component [`EnergyComponents`] and
/// `nspin = 2`; errors with [`crate::error::PwdftError::ConvergenceFailure`]
/// if either channel fails to converge within `max_iter`.
pub(crate) fn run_scf_spin(
    crystal: &Crystal,
    basis: &BasisSet,
    kpoints: &[KPoint],
    pseudopotentials: &[&PseudopotentialData],
    params: &ScfParams,
    symmetry: &crate::symmetry::SymmetryInfo,
    xc_evaluator: XcEvaluator,
) -> Result<ScfResult> {
    let mut ctx = context::ScfContext::new(crystal, basis, kpoints, pseudopotentials, params, symmetry)?;

    // Determine initial spin split from starting_magnetization
    let per_atom_mag: Vec<f64> = ctx.crystal.atoms.iter().map(|a| {
        let sym = crate::atoms::from_z(a.z).map(|e| e.symbol().to_string()).unwrap_or_default();
        *ctx.params.starting_magnetization.get(&sym).unwrap_or(&0.0)
    }).collect();

    // Log the target spin populations if the user fixed the magnetization.
    // The per-spin populations aren't otherwise needed here: the CCMX
    // mixer below operates on total/magnetization and both are sized to
    // `ctx.n_electrons` (the total is conserved exactly, m can be any
    // sign). Per-spin occupations are re-computed inside the SCF loop
    // from eigenvalues.
    if let Some(tot_mag) = ctx.params.tot_magnetization {
        let n_up = f64::midpoint(ctx.n_electrons, tot_mag);
        let n_down = (ctx.n_electrons - tot_mag) / 2.0;
        info!("Fixed magnetization: n_up={n_up:.2}, n_down={n_down:.2}");
    }

    // Initial spin density via SAD with magnetic moments
    let init_config = initial_density::InitialDensityConfig {
        magnetic_moments: per_atom_mag.clone(),
        gaussian_sigma: None,
    };
    // Generate total density then split
    let rho_total = initial_density::generate_initial_density(
        ctx.crystal, &mut ctx.grid, ctx.pseudopotentials, ctx.n_electrons, &init_config,
    );
    // Split: rho_up = (1+m)/2 * rho, rho_down = (1-m)/2 * rho
    // For a uniform initial split, m=0 -> equal channels
    let avg_mag: f64 = per_atom_mag.iter().sum::<f64>() / per_atom_mag.len().max(1) as f64;
    let mut rho_up_r: Vec<f64> = rho_total.iter().map(|&r| r * (1.0 + avg_mag) / 2.0).collect();
    let mut rho_down_r: Vec<f64> = rho_total.iter().map(|&r| r * (1.0 - avg_mag) / 2.0).collect();

    // CCMX: coupled-channel mixer. Rather than mixing (ρ↑, ρ↓) independently
    // — which decouples two channels that are physically coupled by Hartree
    // and XC, leaving them free to oscillate in opposite directions without
    // ever seeing each other's residual history — we mix in the QE basis
    // (ρ_total, m) = (ρ↑ + ρ↓, ρ↑ − ρ↓). The total-density mixer uses the
    // user's full mixing_mode (Kerker/Broyden unchanged); the magnetization
    // mixer uses the same algorithm with Kerker disabled, because Kerker is
    // the charge-charge Thomas-Fermi response 4πe²/(|G|²+q_TF²) and has no
    // meaning for the spin-spin channel (no long-wavelength sloshing to
    // damp). After mixing, we invert back to (ρ↑_new, ρ↓_new). See QE
    // `rhoz_or_updw` (qe-7.5/PW/src/scf_mod.f90:1360-1414) and
    // proposals/CCMX-coupled-channel-mixer.md.
    let mag_mixing_mode = match &ctx.params.mixing_mode {
        mixing::MixingMode::Plain | mixing::MixingMode::Kerker { .. } => {
            mixing::MixingMode::Plain
        }
        mixing::MixingMode::Broyden { .. } => mixing::MixingMode::Broyden { kerker: false },
        mixing::MixingMode::PeriodicPulay { period, .. } => {
            mixing::MixingMode::PeriodicPulay { period: *period, kerker: false }
        }
    };
    // Surface the silent Kerker-off on m once — otherwise a user who configures
    // Kerker sees the charge channel preconditioned but gets no signal about
    // the magnetization channel running plain.
    if matches!(
        ctx.params.mixing_mode,
        mixing::MixingMode::Kerker { .. }
            | mixing::MixingMode::PeriodicPulay { kerker: true, .. }
            | mixing::MixingMode::Broyden { kerker: true }
    ) {
        info!(
            "nspin=2 CCMX: Kerker preconditioning applied to ρ_total channel only; \
             magnetization (m = ρ↑ − ρ↓) channel mixes plain (no charge-response kernel on spin)."
        );
    }
    let mut mixer_total = mixing::Mixer::new(
        ctx.params.mixing_beta, ctx.params.mixing_ndim, &ctx.params.mixing_mode,
        Some(&ctx.g_squared), ctx.n_electrons, ctx.omega, ctx.params.adaptive_beta,
    );
    let mut mixer_mag = mixing::Mixer::new(
        ctx.params.mixing_beta, ctx.params.mixing_ndim, &mag_mixing_mode,
        Some(&ctx.g_squared), ctx.n_electrons, ctx.omega, ctx.params.adaptive_beta,
    );

    let mut e_prev: Option<f64> = None;
    let mut last_delta = f64::INFINITY;
    let pb = scf_progress_bar(ctx.params.max_iter);

    // Precompute half core density (constant across iterations)
    let rho_core_half: Vec<f64> = ctx.rho_core_r.iter().map(|&c| c / 2.0).collect();

    for iter in 0..ctx.params.max_iter {
        // Total density for Hartree (spin-independent)
        let rho_total_r: Vec<f64> = rho_up_r.iter().zip(rho_down_r.iter())
            .map(|(&u, &d)| u + d).collect();
        let mut rho_total_g = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
        density_r_to_g(&mut ctx.grid.fft, &rho_total_r, &mut rho_total_g);

        // 1. Hartree from total density
        let v_h_fft = hartree_on_fft_grid(&rho_total_g, &ctx.g_squared);

        // 2. Spin-dependent XC (routed through the Phase A dispatcher).
        //    v2_*_r stays None for LDA — the LDA path produces zero FFT
        //    work beyond pre-Phase-A behaviour.
        let rho_up_xc = add_core_density(&rho_up_r, &rho_core_half);
        let rho_down_xc = add_core_density(&rho_down_r, &rho_core_half);
        let xc_in = xc_evaluator.eval_spin(&rho_up_xc, &rho_down_xc, None, None)?;
        debug_assert!(
            xc_in.v2_up_r.is_none() && xc_in.v2_down_r.is_none(),
            "GGAP Phase A: LDA spin eval must produce v2_*_r = None"
        );
        let exc_r = xc_in.exc_r;
        let vxc_up_r = xc_in.v1_up_r;
        let vxc_down_r = xc_in.v1_down_r;

        let vxc_up_g = real_to_g_space(&vxc_up_r, &mut ctx.grid.fft);
        let vxc_down_g = real_to_g_space(&vxc_down_r, &mut ctx.grid.fft);

        // 3. Two V_eff: V_local + V_H + V_xc_sigma
        let v_eff_up = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_up_g);
        let v_eff_down = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_down_g);

        // 4. Diagonalize both spins at each k-point. Run the two spin channels
        //    concurrently via `rayon::join`: each closure returns a `Result<Vec<_>>`
        //    from an inner `par_iter` over k-points. `ctx` is borrowed by shared
        //    reference, and every captured field is `Sync` (plain data, `Arc`, or
        //    slice references), so both closures can execute in parallel without
        //    cloning. Errors propagate after the join.
        let eigensolver_kind = ctx.params.eigensolver;
        let (kpoint_results_up, kpoint_results_down): (Result<Vec<_>>, Result<Vec<_>>) = rayon::join(
            || {
                ctx.kpoints.par_iter().enumerate().map(|(ik, kp)| {
                    let mut h = build_hamiltonian_with_v_eff(ctx.basis, &kp.k, &v_eff_up, ctx.grid.dims);
                    ctx.vnl_cache[ik].add_to_hamiltonian(&mut h, ctx.crystal, ctx.basis, &kp.k);
                    diagonalize_dispatch(&h, ctx.params.n_bands, eigensolver_kind)
                }).collect()
            },
            || {
                ctx.kpoints.par_iter().enumerate().map(|(ik, kp)| {
                    let mut h = build_hamiltonian_with_v_eff(ctx.basis, &kp.k, &v_eff_down, ctx.grid.dims);
                    ctx.vnl_cache[ik].add_to_hamiltonian(&mut h, ctx.crystal, ctx.basis, &kp.k);
                    diagonalize_dispatch(&h, ctx.params.n_bands, eigensolver_kind)
                }).collect()
            },
        );
        let kpoint_results_up = kpoint_results_up?;
        let kpoint_results_down = kpoint_results_down?;

        // Flatten eigenvalues: [up_k0, up_k1, ..., down_k0, down_k1, ...]
        let eig_up: Vec<Vec<f64>> = kpoint_results_up.iter().map(|r| r.eigenvalues.clone()).collect();
        let eig_down: Vec<Vec<f64>> = kpoint_results_down.iter().map(|r| r.eigenvalues.clone()).collect();
        let eigenvalues_all: Vec<Vec<f64>> = eig_up.iter().chain(eig_down.iter()).cloned().collect();
        let weights_all: Vec<f64> = ctx.kpt_weights.iter().chain(ctx.kpt_weights.iter()).copied().collect();

        let wfn_up: Vec<_> = kpoint_results_up.into_iter().map(|r| r.eigenvectors).collect();
        let wfn_down: Vec<_> = kpoint_results_down.into_iter().map(|r| r.eigenvectors).collect();

        // 5. Fermi energy and occupations
        let (fermi_energy, occ_up, occ_down) = if let Some(tot_mag) = ctx.params.tot_magnetization {
            // Fixed magnetization: separate Fermi energies per spin
            let n_up_target = f64::midpoint(ctx.n_electrons, tot_mag);
            let n_down_target = (ctx.n_electrons - tot_mag) / 2.0;

            let ef_up = smearing::find_fermi_energy(
                &eig_up, &ctx.kpt_weights, n_up_target,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );
            let ef_down = smearing::find_fermi_energy(
                &eig_down, &ctx.kpt_weights, n_down_target,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );

            let occ_up = compute_occupations(
                &eig_up, ctx.params.smearing_scheme, ef_up,
                ctx.params.smearing_sigma, ctx.spin_factor,
            );
            let occ_down = compute_occupations(
                &eig_down, ctx.params.smearing_scheme, ef_down,
                ctx.params.smearing_sigma, ctx.spin_factor,
            );

            // Report average Fermi energy
            (f64::midpoint(ef_up, ef_down), occ_up, occ_down)
        } else {
            // Free magnetization: single Fermi energy for both spins
            let fermi_energy = smearing::find_fermi_energy(
                &eigenvalues_all, &weights_all, ctx.n_electrons,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );

            let occ_up = compute_occupations(
                &eig_up, ctx.params.smearing_scheme, fermi_energy,
                ctx.params.smearing_sigma, ctx.spin_factor,
            );
            let occ_down = compute_occupations(
                &eig_down, ctx.params.smearing_scheme, fermi_energy,
                ctx.params.smearing_sigma, ctx.spin_factor,
            );

            (fermi_energy, occ_up, occ_down)
        };

        // 6. Reconstruct spin densities
        let n_el_up: f64 = occ_up.iter().zip(ctx.kpt_weights.iter())
            .flat_map(|(occs, &w)| occs.iter().map(move |&f| f * w))
            .sum();
        let n_el_down: f64 = occ_down.iter().zip(ctx.kpt_weights.iter())
            .flat_map(|(occs, &w)| occs.iter().map(move |&f| f * w))
            .sum();

        let rho_up_new = density::compute_density(
            &mut density::DensityGrid {
                basis: ctx.basis, g_to_fft: &ctx.g_to_fft, fft: &mut ctx.grid.fft,
                n_electrons: n_el_up, omega: ctx.omega,
            },
            ctx.kpoints, &wfn_up, &occ_up,
        );
        let rho_down_new = density::compute_density(
            &mut density::DensityGrid {
                basis: ctx.basis, g_to_fft: &ctx.g_to_fft, fft: &mut ctx.grid.fft,
                n_electrons: n_el_down, omega: ctx.omega,
            },
            ctx.kpoints, &wfn_down, &occ_down,
        );

        // Symmetrize each channel in G-space (PCFX). See the non-spin
        // run_scf for the full rationale; key point: the real-space form
        // is exact only when τ_S lands on an integer grid point, which
        // fails for Fd-3m on an 18³ grid. G-space is exact for any τ.
        let mut rho_up_sym = rho_up_new;
        let mut rho_down_sym = rho_down_new;
        crate::symmetry::density::symmetrize_density_g(
            &mut rho_up_sym,
            ctx.grid.dims,
            &mut ctx.grid.fft,
            ctx.symmetry,
        );
        crate::symmetry::density::symmetrize_density_g(
            &mut rho_down_sym,
            ctx.grid.dims,
            &mut ctx.grid.fft,
            ctx.symmetry,
        );

        // 7. Convergence check
        let rho_total_new: Vec<f64> = rho_up_sym.iter().zip(rho_down_sym.iter())
            .map(|(&u, &d)| u + d).collect();
        // SPNC: use per-spin max, not total-density diff. A spin-flip fluctuation
        // (+ε in rho_up, −ε in rho_down) is invisible to the total but keeps
        // zeta_in ≠ zeta_out, which leaves E_xc[rho, zeta] inconsistent and
        // spoils the O(Δρ²) convergence of |E_HF − E_KS|. Per-spin max is
        // strictly stronger than total and keeps the scalar `conv_threshold`
        // semantics unchanged. See proposals/SPNC-spin-per-density-convergence.md.
        let delta_up = density_diff(&rho_up_r, &rho_up_sym, ctx.omega, ctx.n_grid);
        let delta_down = density_diff(&rho_down_r, &rho_down_sym, ctx.omega, ctx.n_grid);
        let delta = delta_up.max(delta_down);
        last_delta = delta;

        // Energy from new density
        let mut rho_total_new_g = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
        density_r_to_g(&mut ctx.grid.fft, &rho_total_new, &mut rho_total_new_g);

        let rho_xc_total: Vec<f64> = add_core_density(&rho_total_new, &ctx.rho_core_r);
        let occ_all: Vec<Vec<f64>> = occ_up.iter().chain(occ_down.iter()).cloned().collect();

        // SPXC fix: recompute spin-polarized XC from OUTPUT spin densities for E_KS.
        // Previously `exc_r` (INPUT density, depends on zeta_in) was integrated against
        // `rho_xc_total` (OUTPUT total density), which is not physically meaningful and
        // introduced an O(delta_rho) error that spoiled the quadratic convergence of
        // |E_HF - E_KS|. The non-spin run_scf already recomputes XC from the output
        // density (see `let (exc_r, vxc_r_energy) = xc::lda_xc_grid(&rho_new_for_xc)`).
        // INPUT-based quantities (exc_r, vxc_up_r, vxc_down_r) remain for E_HF below.
        let rho_up_xc_out = add_core_density(&rho_up_sym, &rho_core_half);
        let rho_down_xc_out = add_core_density(&rho_down_sym, &rho_core_half);
        let xc_out = xc_evaluator.eval_spin(&rho_up_xc_out, &rho_down_xc_out, None, None)?;
        let exc_r_out = xc_out.exc_r;
        let vxc_up_r_out = xc_out.v1_up_r;
        let vxc_down_r_out = xc_out.v1_down_r;

        // Spin XC double-counting (OUTPUT density): E_vxc = integral(V_xc_up rho_up_out + V_xc_down rho_down_out) dr
        // Both V_xc and rho_sigma here come from the OUTPUT (symmetrized) density,
        // matching the non-spin run_scf convention.
        let dvol = ctx.omega / ctx.n_grid as f64;
        let e_vxc_spin_out: f64 = rho_up_sym.iter().zip(vxc_up_r_out.iter())
            .zip(rho_down_sym.iter().zip(vxc_down_r_out.iter()))
            .map(|((&ru, &vu), (&rd, &vd))| (ru * vu + rd * vd) * dvol)
            .sum();
        let e_xc_out = xc::lda_xc_energy(&rho_xc_total, &exc_r_out, ctx.omega);
        let e_xc_corrected_out = e_xc_out - e_vxc_spin_out;

        let e_band = band_energy(&eigenvalues_all, &occ_all, &weights_all);

        // Kohn-Sham energy: double-counting from OUTPUT density
        let e_total = with_g0_shift(total_energy(
            e_band,
            hartree_energy(&rho_total_new_g, &ctx.g_squared, ctx.omega),
            e_xc_corrected_out,
            ctx.e_ewald,
        ), &ctx);

        // Harris-Foulkes energy: double-counting from INPUT density
        // rho_total_g, rho_up_xc, rho_down_xc, exc_r, vxc_up_r, vxc_down_r
        // are all from the input density
        let rho_xc_total_in: Vec<f64> = add_core_density(&rho_total_r, &ctx.rho_core_r);
        let e_xc_in = xc::lda_xc_energy(&rho_xc_total_in, &exc_r, ctx.omega);
        let e_vxc_spin_in: f64 = rho_up_r.iter().zip(vxc_up_r.iter())
            .zip(rho_down_r.iter().zip(vxc_down_r.iter()))
            .map(|((&ru, &vu), (&rd, &vd))| (ru * vu + rd * vd) * dvol)
            .sum();
        let e_xc_corrected_in = e_xc_in - e_vxc_spin_in;

        let e_harris = with_g0_shift(harris_foulkes_energy(
            e_band,
            hartree_energy(&rho_total_g, &ctx.g_squared, ctx.omega),
            e_xc_corrected_in,
            ctx.e_ewald,
        ), &ctx);

        let hf_diff = (e_harris - e_total).abs();

        let de = e_prev.map(|ep| (e_total - ep).abs());
        e_prev = Some(e_total);

        let rho_converged = delta < ctx.params.conv_threshold;
        let energy_converged = de.is_some_and(|de| de < ctx.params.energy_threshold);

        let mag = (n_el_up - n_el_down).abs();
        log_iteration(&pb, &IterationReport {
            iter,
            e_total,
            e_harris,
            hf_diff,
            de,
            delta,
            beta: ctx.params.adaptive_beta.then(|| mixer_total.current_beta()),
            spin: Some(SpinIterationFields {
                delta_up,
                delta_down,
                magnetization: mag,
                mag_beta: ctx.params.adaptive_beta.then(|| mixer_mag.current_beta()),
            }),
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
            info!("Magnetization: {mag:.4} μB ({n_el_up:.4} up, {n_el_down:.4} down)");

            let rho_g_basis: Vec<Complex64> = ctx.g_to_fft.iter().map(|&idx| rho_total_new_g[idx]).collect();

            let ts = smearing::entropy_ts(
                &eigenvalues_all, &weights_all, fermi_energy,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );
            let free_energy = e_total - ts;
            let energy_sigma0 = f64::midpoint(e_total, free_energy);

            // -------------------------------------------------------------
            // VGC5 per-component decomposition (diagnostic, spin-polarized).
            // Kinetic / non-local: sum over both spin channels.
            // Local / Hartree / Ewald: built from total density (spin-
            //   independent operators).
            // XC: ∫(ρ_up+ρ_down+ρ_core)·ε_xc(ρ_up,ρ_down)dr (bare, from OUTPUT).
            // -------------------------------------------------------------
            let k_vecs: Vec<nalgebra::Vector3<f64>> =
                ctx.kpoints.iter().map(|kp| kp.k).collect();

            let e_kin_up = kinetic_expectation(
                ctx.basis, &k_vecs, &ctx.kpt_weights, &wfn_up, &occ_up,
            );
            let e_kin_down = kinetic_expectation(
                ctx.basis, &k_vecs, &ctx.kpt_weights, &wfn_down, &occ_down,
            );
            let e_kinetic = e_kin_up + e_kin_down;

            let mut v_local_cplx = ctx.v_local_fft.clone();
            ctx.grid.fft.inverse(&mut v_local_cplx);
            let v_local_r: Vec<f64> = v_local_cplx.iter().map(|c| c.re).collect();
            let e_local = local_pp_energy_grid(&rho_total_new, &v_local_r, ctx.omega, ctx.n_grid);
            let e_local_g0_shift = ctx.v_local_g0 * ctx.n_electrons;

            let e_nl_up = nonlocal_expectation(
                ctx.basis, ctx.crystal, &k_vecs, &ctx.kpt_weights,
                &wfn_up, &occ_up, &ctx.vnl_cache,
            );
            let e_nl_down = nonlocal_expectation(
                ctx.basis, ctx.crystal, &k_vecs, &ctx.kpt_weights,
                &wfn_down, &occ_down, &ctx.vnl_cache,
            );
            let e_nonlocal = e_nl_up + e_nl_down;

            let e_hartree_term = hartree_energy(&rho_total_new_g, &ctx.g_squared, ctx.omega);
            let e_xc_term = xc_energy_bare(&rho_xc_total, &exc_r_out, ctx.omega);
            // LSDA XC double-counting: `e_vxc_spin_out` is computed above
            // (line ~386) as ∫(ρ↑·V_xc↑ + ρ↓·V_xc↓) d³r on the OUTPUT
            // density and is the `e_vxc` term in the band-sum identity
            // `E_band = e_kin + e_loc + e_nl + 2·e_hartree + e_vxc`.
            let e_vxc_term = e_vxc_spin_out;
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
                magnetization: mag,
                nspin: 2,
                components,
            });
        }

        // 8. CCMX: coupled-channel mixing in the (ρ_total, m) basis.
        //    Forward basis change: ρ_total = ρ↑ + ρ↓, m = ρ↑ − ρ↓ for both
        //    the input and output (symmetrized) spin densities. Mix each
        //    mode independently with its own history. Invert with
        //    ρ↑ = (ρ_total + m) / 2, ρ↓ = (ρ_total − m) / 2 (matches QE's
        //    vi=0.5 convention in rhoz_or_updw, scf_mod.f90:1395-1396).
        let rho_total_in: Vec<f64> = rho_up_r
            .iter()
            .zip(rho_down_r.iter())
            .map(|(&u, &d)| u + d)
            .collect();
        let m_in: Vec<f64> = rho_up_r
            .iter()
            .zip(rho_down_r.iter())
            .map(|(&u, &d)| u - d)
            .collect();
        let rho_total_out: Vec<f64> = rho_up_sym
            .iter()
            .zip(rho_down_sym.iter())
            .map(|(&u, &d)| u + d)
            .collect();
        let m_out: Vec<f64> = rho_up_sym
            .iter()
            .zip(rho_down_sym.iter())
            .map(|(&u, &d)| u - d)
            .collect();

        let rho_total_mixed = mixer_total.mix(&rho_total_in, &rho_total_out, &mut ctx.grid.fft);
        let m_mixed = mixer_mag.mix(&m_in, &m_out, &mut ctx.grid.fft);

        rho_up_r = rho_total_mixed
            .iter()
            .zip(m_mixed.iter())
            .map(|(&t, &mm)| 0.5 * (t + mm))
            .collect();
        rho_down_r = rho_total_mixed
            .iter()
            .zip(m_mixed.iter())
            .map(|(&t, &mm)| 0.5 * (t - mm))
            .collect();
    }

    pb.abandon_with_message("did not converge");
    Err(PwdftError::ConvergenceFailure {
        iterations: ctx.params.max_iter,
        delta: last_delta,
    })
}
