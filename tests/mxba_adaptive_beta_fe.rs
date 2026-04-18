//! MXBA diagnostic: Fe BCC free-magnetization SCF with adaptive β.
//!
//! Companion to `tests/spin_polarization.rs::test_ccmx_fe_free_magnetization_converges`
//! (which runs the same system at the MXBA *default* of adaptive β = off
//! and converges in 14 iterations). This test exercises the *adaptive β on*
//! path and documents the empirical cost: on Fe BCC 4×4×4 nspin=2 Kerker,
//! adaptive β *hurts* because the residual norm plateaus at Δρ≈0.34 for
//! the first ~5 iterations while Anderson builds DIIS history. The
//! Eyert-1996 monitor interprets the flat residual as "not converging
//! well" and damps β all the way down to β_min ≈ 0.017, at which point
//! Anderson's DIIS extrapolation is starved and the SCF never escapes
//! the plateau.
//!
//! Observed trajectory (post-MXBA, adaptive ON):
//! - iter 1–9: Δρ ≈ 0.34, β = 0.3 (initial monitor warm-up)
//! - iter 10–20: β_tot steps 0.3 → 0.21 → 0.147 → 0.103 → 0.072
//! - iter 20+:  β_tot floors at ~0.017; dE drops to 1e-6 but Δρ stuck at 0.34
//! - iter 80:   ConvergenceFailure (delta = 0.34 vs threshold 1e-3)
//!
//! This is exactly the failure mode the user identified: "If adaptive β is
//! slower, the default should be disabled (user opts in)." The MXBA default
//! is therefore `adaptive_beta = false`, and this test is `#[ignore]`d to
//! preserve the documented failure as a reproducible reference without
//! blocking CI.
//!
//! Re-enable by running `cargo test test_mxba_fe_documents_adaptive_failure
//! -- --include-ignored --nocapture` and inspecting the per-iteration β log.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use pwdft_rs::{basis::BasisSet, crystal::Crystal, scf, scf::mixing::MixingMode};

fn fe_bcc() -> Crystal {
    use nalgebra::Vector3;
    use pwdft_rs::crystal::{Atom, Lattice};
    let a = 2.867;
    Crystal {
        lattice: Lattice::new(
            a * Vector3::new(1.0, 0.0, 0.0),
            a * Vector3::new(0.0, 1.0, 0.0),
            a * Vector3::new(0.0, 0.0, 1.0),
        ),
        atoms: vec![
            Atom::new(26, [0.0, 0.0, 0.0]),
            Atom::new(26, [0.5, 0.5, 0.5]),
        ],
    }
}

#[test]
#[ignore = "documents known adaptive-β failure mode on Fe BCC CCMX; \
            confirms the MXBA default-off decision. Re-run when tuning the \
            Eyert thresholds or after a mixer-topology change."]
fn test_mxba_fe_documents_adaptive_failure() {
    let _ = env_logger::builder().is_test(true).try_init();
    let crystal = fe_bcc();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("pseudopotentials/nc/lda/Fe.upf"),
    )
    .unwrap();
    let ecut = 15.0 * 13.605_693_122_994;
    let basis = BasisSet::new(&crystal.lattice, ecut);
    let kpoints = pwdft_rs::kpoints::monkhorst_pack(4, 4, 4, &crystal.lattice);

    let mut starting_mag = std::collections::HashMap::new();
    starting_mag.insert("Fe".to_string(), 0.5);

    let params = scf::ScfParams {
        n_bands: 8,
        max_iter: 80,
        conv_threshold: 1e-3,
        energy_threshold: 1e-3,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.02 * 13.605_693_122_994,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Kerker { q_tf: None },
        adaptive_beta: true, // MXBA: documented to fail on this system
        nspin: 2,
        starting_magnetization: starting_mag,
        ..Default::default()
    };

    let symmetry = pwdft_rs::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params, &symmetry);

    // Documents the failure: adaptive β at these defaults does NOT converge
    // Fe CCMX. If the Eyert tuning is reworked and this test starts passing,
    // the assertion below will flag that — at which point the `#[ignore]`
    // can be removed and this test re-purposed as a positive regression.
    match result {
        Ok(res) => {
            panic!(
                "MXBA adaptive β unexpectedly converged Fe CCMX: E_KS={:.6} eV, \
                 iters={}, M={:.3} μB. The previously-documented failure mode \
                 has changed — update proposal MXBA + unignore this test.",
                res.total_energy, res.n_iterations, res.magnetization
            );
        }
        Err(err) => {
            eprintln!("MXBA adaptive β on Fe CCMX: documented failure mode hit — {err}");
        }
    }
}
