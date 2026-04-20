//! Mixer-robustness regression pins.
//!
//! This binary collects tests that document **known mixer-algorithm failure
//! modes** on specific physical systems, pinned so that any future change
//! which silently alters the failure mode is forced to declare itself.
//!
//! The pattern follows `tests/mxba_adaptive_beta_fe.rs` — rather than
//! assert-on-convergence (positive regression), these assert-on-stall
//! (negative regression). If someone later fixes the underlying mixer
//! pathology, the assertion flips: the test fails loudly, the author has
//! to decide whether the fix is real, and the assertion is rewritten as
//! a positive convergence pin with a tolerance against the QE reference.
//!
//! ## Anderson-stall on wide-gap insulators
//!
//! Plain Anderson mixing (no preconditioning, no Broyden history of inverse
//! Jacobians) stalls on wide-gap insulators whose residual response is
//! dominated by long-wavelength modes. For C diamond at ecut = 30 Ry on a
//! 4×4×4 Γ-centered grid, the residual plateaus at Δρ ≈ 10⁻⁵ and does not
//! decay — even with 150 iterations. Switching to Kerker preconditioning,
//! Broyden (with or without Kerker), or Periodic Pulay converges cleanly
//! in 10–15 iterations (see the table in
//! `tests/qe_validation.rs::test_c_diamond_vs_qe` docstring).
//!
//! Two independent observations (Si at FFT grid 20/24 per the PCRS /
//! per-component-residual investigation; C diamond per VQEF's QE full-
//! matrix sweep) establish the pattern as a real, latent, mixer-robustness
//! limitation — not a one-off anomaly. The limitation is a *performance*
//! issue, not a correctness one: the robust-mixer alternatives already
//! exist and converge to the same physics. This test exists to document
//! the failure mode so the next person who touches the Anderson
//! implementation cannot accidentally remove or regress the pathology
//! without noticing.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use std::collections::HashMap;

use nalgebra::Vector3;
use pwdft_core::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    kpoints,
    scf::{self, ScfParams, mixing::MixingMode, smearing::SmearingScheme},
    symmetry::SymmetryInfo,
};

const RY_TO_EV: f64 = 13.605_693_122_994;

fn c_diamond() -> Crystal {
    // FCC primitive cell, lattice parameter 3.567 Å (matches qe_validation
    // reference for C diamond). Two-atom basis at (0,0,0) and (1/4,1/4,1/4).
    let a = 3.567;
    Crystal {
        lattice: Lattice::new(
            a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        ),
        atoms: vec![
            Atom::new(6, [0.00, 0.00, 0.00]),
            Atom::new(6, [0.25, 0.25, 0.25]),
        ],
    }
}

/// Pin Plain Anderson's stall on C diamond as a negative regression guard.
///
/// **Expected outcome:** SCF does NOT converge within 150 iterations — the
/// driver returns `ConvergenceFailure { delta, .. }` with `delta` *above*
/// `conv_threshold = 1e-8`. Observed residual floor on 2026-04-19 with
/// `mixing_beta = 0.3`, `mixing_ndim = 8`: Δρ ≈ 1.75e-8, stuck there for
/// the final ~50 iters. (With `β = 0.1` or `ndim = 16` the floor is
/// higher at ≈ 1.4e-5 / 2.0e-5; see `tests/qe_validation.rs::
/// test_c_diamond_vs_qe` docstring table.) Any mixer with Kerker
/// preconditioning, Broyden, or Periodic Pulay converges cleanly under
/// the same conditions in 10–15 iters.
///
/// If this test ever returns `Ok(result)` instead of `ConvergenceFailure`,
/// that means someone has fixed Plain Anderson's stall behaviour. At that
/// point the next step is NOT to silently weaken the test: it is to
/// rewrite the assertion as a positive convergence pin against the QE
/// reference (`E = -23.843_439_10 Ry`) and update the `MixingMode::Plain`
/// docstring accordingly.
///
/// Tier-2 because it runs a full SCF with `max_iter = 150`. Stays
/// `#[ignore]` so the default `cargo test` wall time stays green;
/// exercise it with `cargo test --test mixer_robustness -- --ignored`.
#[test]
#[ignore = "TSPL tier-2: heavy SCF, run with --ignored when touching scf/mixing/**, scf/driver*.rs, or scf/density.rs — pins Plain Anderson stall on wide-gap insulators as negative regression"]
fn test_plain_anderson_stalls_on_c_diamond() {
    let _ = env_logger::builder().is_test(true).try_init();

    let crystal = c_diamond();
    let pp_c = pwdft_core::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR"))
            .join("pseudopotentials/nc/lda/C.upf"),
    )
    .expect("C.upf must load");

    let ecut_ry = 30.0_f64;
    let ecut_ev = ecut_ry * RY_TO_EV;
    let basis = BasisSet::new(&crystal.lattice, ecut_ev);
    // Γ-centered 4×4×4 — matches QE's `K_POINTS automatic / 4 4 4 0 0 0`
    // in data/qe/c_scf.in. The stall is specific to this FFT-grid
    // / gap / mixer combination; changing the grid shift or k-density may
    // alter the failure mode.
    let kpts = kpoints::monkhorst_pack(
        4,
        4,
        4,
        kpoints::KGridShift::GammaCentered,
        &crystal.lattice,
    );

    // All other knobs match the `test_c_diamond_vs_qe` config in
    // qe_validation.rs so the only difference is the mixer choice.
    let params = ScfParams {
        n_bands: 8,
        max_iter: 150,
        conv_threshold: 1e-8,
        energy_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.01 * RY_TO_EV,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Plain,
        nspin: 1,
        starting_magnetization: HashMap::new(),
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(&crystal, &basis, &kpts, &[&pp_c], &params, &symmetry);

    // The negative regression: Plain Anderson is expected to stall on C
    // diamond, i.e. the driver should return `ConvergenceFailure` rather
    // than `Ok(..)`. Three outcomes flag a behaviour change that deserves
    // human attention:
    //
    //  1. `Ok(res)` — someone fixed Plain Anderson. Rewrite as a positive
    //     convergence pin against the QE reference.
    //  2. `Err(ConvergenceFailure { delta, .. })` with `delta <=
    //     conv_threshold` — self-contradictory, would indicate a driver
    //     bug.
    //  3. `Err(other)` — the SCF failed for a different reason (NaN,
    //     eigensolver breakdown, etc.). Investigate; don't paper over.
    match result {
        Ok(res) => {
            panic!(
                "Plain Anderson on C diamond (ecut=30 Ry, Γ-centered 4×4×4) was expected \
                 to stall (ConvergenceFailure) after 150 iters, but SCF converged with \
                 final_delta = {:.3e} in {} iters (E_KS = {:.6} eV). If you fixed this \
                 mixer pathology, replace this negative-regression assertion with a \
                 positive convergence pin against the QE reference \
                 (E = -23.843_439_10 Ry = -324.319 eV, tol 0.05 eV) and update the \
                 MixingMode::Plain docstring in src/scf/mixing/mod.rs accordingly.",
                res.final_delta, res.n_iterations, res.total_energy,
            );
        }
        Err(pwdft_core::error::PwdftError::ConvergenceFailure { iterations, delta }) => {
            // Self-consistency: if the driver reports ConvergenceFailure
            // the residual must be above conv_threshold = 1e-8, else the
            // driver itself is buggy (failing to recognize convergence).
            assert!(
                delta > 1e-8,
                "Plain Anderson on C diamond: ConvergenceFailure reported with \
                 delta = {delta:.3e} <= conv_threshold (1e-8) after {iterations} iters. \
                 The driver is reporting failure on an already-converged state — \
                 investigate the convergence check, not the mixer.",
            );
            eprintln!(
                "Plain Anderson on C diamond: documented stall hit — \
                 ConvergenceFailure(iters={iterations}, delta={delta:.3e}). Expected.",
            );
        }
        Err(other) => {
            panic!(
                "Plain Anderson on C diamond failed for an unexpected reason \
                 (not ConvergenceFailure): {other}. The documented pathology is a mixer \
                 stall, not a hard SCF failure — investigate before updating the assertion.",
            );
        }
    }
}
