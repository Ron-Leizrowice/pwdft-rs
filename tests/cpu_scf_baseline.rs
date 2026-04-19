//! CPU baseline SCF captures for the GPU-vs-CPU consistency gate.
//!
//! This file is built **without** the `gpu` feature so that it cannot be
//! contaminated by accidentally running through the GPU path — a defense
//! against the exact audit finding (RWHK C1) that motivated this file.
//! Under `cargo test --features gpu`, this file is skipped by the
//! `#![cfg(not(feature = "gpu"))]` gate at the top of the file; under
//! plain `cargo test`, it runs and captures the CPU-f64 baseline that
//! the GPU-side test in `tests/gpu_consistency.rs` pins against.
//!
//! The GPU side of the comparison lives in
//! `tests/gpu_consistency.rs::test_gpu_scf_matches_cpu_within_f32_tolerance`,
//! which asserts the GPU result matches the constant captured here.
//!
//! See `proposals/RWHK-reward-hacking-audit-2026-04-19.md` § C1 for the
//! original reward-hacking finding (the previous test compared GPU
//! against GPU — a tautology).

#![cfg(not(feature = "gpu"))]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use nalgebra::Vector3;

use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
};

// ---------------------------------------------------------------------------
// CPU baseline constants — mirrored in tests/gpu_consistency.rs
// ---------------------------------------------------------------------------
//
// These constants are the "Si diamond primitive, ecut=100, FFT 16³,
// Γ-only, 4 bands, LDA" SCF baseline. The GPU-vs-CPU consistency test
// in `tests/gpu_consistency.rs::test_gpu_scf_matches_cpu_within_f32_tolerance`
// pins against `SI_CPU_BASELINE_TOTAL_EV` specifically. Keep the two in
// sync; any drift here requires updating both files.
const SI_CPU_BASELINE_TOTAL_EV: f64 = -213.028_339;
const SI_CPU_BASELINE_FERMI_EV: f64 = 8.052_064;

/// Si diamond primitive cell (2 atoms) matching
/// `tests/gpu_consistency.rs::si_crystal()` exactly so the CPU baseline
/// and GPU test pin the same SCF configuration.
fn si_crystal() -> Crystal {
    let a = 5.431;
    Crystal {
        lattice: Lattice::new(
            a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        ),
        atoms: vec![
            Atom::new(14, [0.0, 0.0, 0.0]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    }
}

/// SCF parameters matching `tests/gpu_consistency.rs::si_scf_params()`.
fn si_scf_params() -> pwdft_rs::scf::ScfParams {
    pwdft_rs::scf::ScfParams {
        n_bands: 4,
        max_iter: 40,
        conv_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 4,
        smearing_sigma: 0.05,
        ecutrho_ratio: 4,
        fft_grid: Some([16, 16, 16]),
        ..Default::default()
    }
}

/// Captures the CPU-f64 SCF baseline for Si diamond at ecut=100, FFT 16³,
/// Γ-only, 4 bands — the reference against which the GPU-side
/// `test_gpu_scf_matches_cpu_within_f32_tolerance` pins. Verifies the
/// pinned baseline constants (`SI_CPU_BASELINE_TOTAL_EV`,
/// `SI_CPU_BASELINE_FERMI_EV`) in `tests/gpu_consistency.rs` still
/// reproduce bit-identically on the CPU path.
///
/// This test is **tautology-proof** for the GPU-vs-CPU comparison: the
/// whole file is gated `#![cfg(not(feature = "gpu"))]`, so
/// `pwdft_rs` is compiled without the `gpu` feature here — there is no
/// GPU init code linked in, no `GpuAccelerator::try_new()` possible. When
/// `cargo test --features gpu` runs the GPU-side matching test, it pins
/// against the constants captured from this CPU-only run.
///
/// Pinned values (captured 2026-04-19 under RWHK-FIX on `origin/main`
/// post-MODR / CAST / DWGT / FGRD / TSPL / VGCH-SiEF-B1):
///   `E_total = -213.028_339 eV`, `E_F = 8.052_064 eV`.
///
/// The E_F value is post-VGCH-SiEF-B1: every KS eigenvalue now carries
/// `V_loc(G=0)` as a DC offset (QE-compatible gauge), shifting Si E_F by
/// ≈ 1.343 eV from the pre-B1 reference. `E_total` is algebraically
/// invariant to that gauge change (see VGCH-SiEF-B1 in
/// `tests/qe_validation.rs::test_si_total_energy_bit_identity_post_siefb1`),
/// which is why it agrees with the pre-B1 -213.0283 eV pin to rounding.
///
/// Tolerance 1e-3 eV (way tighter than f32 noise) asserts the CPU path
/// is deterministic across runs — any drift caught would indicate a
/// regression in the deterministic CPU pipeline (e.g. a parallel reduction
/// became order-dependent).
#[test]
#[ignore = "TSPL Tier-2: runs Si SCF on CPU at ecut=100, 40 iters; run with cargo test -- --ignored when touching scf/, xc/, or fft paths. Paired with tests/gpu_consistency.rs::test_gpu_scf_matches_cpu_within_f32_tolerance."]
fn test_cpu_scf_baseline_convergence() {
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();
    let kpoints = vec![pwdft_rs::kpoints::KPoint {
        k: nalgebra::Vector3::zeros(),
        weight: 1.0,
        label: None,
    }];
    let params = si_scf_params();

    let result = pwdft_rs::scf::run_scf(
        &crystal,
        &basis,
        &kpoints,
        &[&pp],
        &params,
        &pwdft_rs::symmetry::SymmetryInfo::identity_only(),
    )
    .expect("CPU Si SCF must converge");

    assert!(
        result.n_iterations < params.max_iter,
        "CPU Si SCF hit max_iter={} without converging",
        params.max_iter,
    );

    eprintln!(
        "CPU baseline: {} iters, E = {:.10} eV, E_F = {:.10} eV",
        result.n_iterations, result.total_energy, result.fermi_energy,
    );

    // The baselines (module consts above) are the values the GPU
    // test pins against. If this test fails, either the CPU path has
    // drifted (physics regression — investigate before updating the pin)
    // or the pin needs to be refreshed after an intentional change
    // (update both module consts and `SI_CPU_BASELINE_TOTAL_EV` in
    // `tests/gpu_consistency.rs` together).
    let de = (result.total_energy - SI_CPU_BASELINE_TOTAL_EV).abs();
    assert!(
        de < 1e-3,
        "CPU Si baseline E_total drifted: expected {SI_CPU_BASELINE_TOTAL_EV:.6} eV, \
         got {:.6} eV, |Δ|={de:.3e} eV. If this was intentional (e.g. XC fix), \
         update both this constant and SI_CPU_BASELINE_TOTAL_EV in \
         tests/gpu_consistency.rs together.",
        result.total_energy,
    );

    let df = (result.fermi_energy - SI_CPU_BASELINE_FERMI_EV).abs();
    assert!(
        df < 1e-3,
        "CPU Si baseline E_F drifted: expected {SI_CPU_BASELINE_FERMI_EV:.6} eV, \
         got {:.6} eV, |Δ|={df:.3e} eV",
        result.fermi_energy,
    );
}
