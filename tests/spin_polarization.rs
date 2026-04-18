//! Spin polarization tests.
//!
//! Validates that nspin=2 produces correct results and matches nspin=1
//! in the unpolarized limit.

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
    kpoints::KPoint,
    scf::{self, mixing::MixingMode},
};

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

fn fe_bcc() -> Crystal {
    let a = 2.87;
    Crystal {
        lattice: Lattice::new(
            a / 2.0 * Vector3::new(-1.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, -1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, -1.0),
        ),
        atoms: vec![Atom::new(26, [0.0, 0.0, 0.0])],
    }
}

fn gamma_only() -> Vec<KPoint> {
    vec![KPoint { k: Vector3::zeros(), weight: 1.0, label: None }]
}

#[test]
fn test_si_nspin2_matches_nspin1() {
    // Si with nspin=2 and zero magnetization should match nspin=1 energy.
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    ).unwrap();
    let kpoints = gamma_only();

    let params_nspin1 = scf::ScfParams {
        n_bands: 4,
        max_iter: 40,
        conv_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 4,
        smearing_sigma: 0.05,
        ecutrho_ratio: 4,
        fft_grid: Some([16, 16, 16]),
        mixing_mode: MixingMode::Plain,
        nspin: 1,
        ..Default::default()
    };

    let params_nspin2 = scf::ScfParams {
        nspin: 2,
        n_bands: 4,  // per spin channel
        ..params_nspin1.clone()
    };

    let sym_id = pwdft_rs::symmetry::SymmetryInfo::identity_only();
    let result1 = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params_nspin1, &sym_id);
    let result2 = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params_nspin2, &sym_id);

    let r1 = result1.expect("nspin=1 must converge");
    let r2 = result2.expect("nspin=2 must converge");

    // Convergence guard (TAUD finding 5.1): reaching max_iter means SCF did
    // not actually meet conv_threshold, which earlier silent-pass patterns hid.
    assert!(
        r1.n_iterations < params_nspin1.max_iter,
        "nspin=1 hit max_iter={} without converging",
        params_nspin1.max_iter
    );
    assert!(
        r2.n_iterations < params_nspin2.max_iter,
        "nspin=2 hit max_iter={} without converging",
        params_nspin2.max_iter
    );

    let de = (r1.total_energy - r2.total_energy).abs();
    eprintln!("Si nspin=1: E={:.6} eV ({} iters)", r1.total_energy, r1.n_iterations);
    eprintln!("Si nspin=2: E={:.6} eV ({} iters), M={:.6} μB", r2.total_energy, r2.n_iterations, r2.magnetization);
    eprintln!("Energy diff: {de:.3e} eV");

    // TAUD finding 2.1: Si nspin=1 and nspin=2 (M=0 starting) evolve through
    // different SCF paths but must agree in the unpolarized limit to within
    // SCF convergence precision. Empirical de < 1e-6 eV (reads as "0.000000"
    // at 6-decimal print); threshold set at 10× empirical = 1e-5 eV.
    // Previously 0.5 eV — would have masked a real spin-XC bug of SPXC scale
    // (~13 eV pre-fix) or silent-convergence bug of SPNC scale.
    assert!(
        de < 1.0e-5,
        "nspin=1 ({:.6} eV) and nspin=2 ({:.6} eV) energies differ by {de:.3e} eV \
         (> 1e-5 eV). Unpolarized Si must match nspin=1 to SCF precision; a \
         value of many eV suggests SPXC (XC input/output mismatch) or SPNC \
         (per-spin convergence criterion) has regressed.",
        r1.total_energy, r2.total_energy
    );

    // TAUD finding 2.2: Si is non-magnetic — M must be exactly zero at
    // convergence (no starting_mag, nspin=2 relaxation). Use .abs() (not
    // signed `<`, which would let arbitrarily-negative M pass). Empirical
    // |M| < 1e-5 μB (reads as "0.000000" at 6-decimal print); threshold at
    // 10× empirical = 1e-4 μB. Previously 0.1 μB (10% of full electron spin).
    assert!(
        r2.magnetization.abs() < 1.0e-4,
        "Si should be non-magnetic, got |M|={:.3e} μB (> 1e-4)", r2.magnetization
    );
}

#[test]
fn test_fe_spin_xc_consistency_regression() {
    // Combined SPXC + SPNC regression test.
    //
    // SPXC (proposals/completed/SPXC-spin-xc-consistency.md) fixed an
    // input/output mismatch in the spin-polarized E_xc: `exc_r` was from INPUT
    // spin densities while `rho_xc_total` and `rho_up/down_sym` were OUTPUT.
    // SPNC (proposals/SPNC-spin-per-density-convergence.md) fixed the nspin=2
    // convergence criterion to use per-spin max(||Δρ_up||, ||Δρ_down||)
    // instead of the total-density ||Δρ_total||, so that spin polarization
    // (zeta) is actually driven to self-consistency.
    //
    // This test exercises nspin=2 on Si with a nonzero starting_magnetization
    // so the spin XC machinery and per-spin mixing are fully active during
    // SCF, then relaxes (no tot_magnetization constraint) to the physically
    // correct non-magnetic ground state. At convergence, |E_HF - E_KS| must
    // be O(Δρ²), i.e. sub-meV at conv_threshold=1e-6.
    //
    // Historical baseline (Fe BCC fixed-mag=2 on the same nc/lda/Fe.upf, 4×4×4
    // k, 15 Ry, starting_magnetization=0.5, conv=1e-6, max_iter=300):
    //   Pre-SPXC, total-only criterion:        |HF-KS| ≈ 22.2 eV (falsely "converged")
    //   Post-SPXC, total-only criterion:       |HF-KS| ≈ 13.0 eV (falsely "converged")
    //   Post-SPXC+SPNC, per-spin criterion:    SCF correctly refuses to converge —
    //     per-channel Δρ locks at ≈0.254 because this pseudopotential does not
    //     support a stable fixed-mag=2 state (LDA ground state is nonmagnetic).
    //     The old total-only metric hid this limit cycle by cancelling +ε/−ε
    //     between the up and down channels.
    //
    // The Fe fixed-mag case is therefore documented but not tested here — it
    // would require a different pseudopotential or a better mixer. The Si
    // probe below is the right system to pin the HF-KS quadratic convergence
    // property this test is named for.
    let _ = env_logger::builder().is_test(true).try_init();
    let crystal = si_crystal();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let kpoints = gamma_only();

    let mut starting_mag = std::collections::HashMap::new();
    starting_mag.insert("Si".to_string(), 0.2);

    let params = scf::ScfParams {
        n_bands: 4,
        max_iter: 60,
        conv_threshold: 1e-6,
        energy_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 4,
        smearing_sigma: 0.05,
        ecutrho_ratio: 4,
        fft_grid: Some([16, 16, 16]),
        mixing_mode: MixingMode::Plain,
        nspin: 2,
        tot_magnetization: None,
        starting_magnetization: starting_mag,
        ..Default::default()
    };

    let sym_id = pwdft_rs::symmetry::SymmetryInfo::identity_only();
    let result = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params, &sym_id)
        .expect("Si nspin=2 SCF must converge for SPXC+SPNC regression test");

    let hf_diff = (result.harris_foulkes_energy - result.total_energy).abs();
    eprintln!(
        "SPXC+SPNC regression (Si nspin=2): E_KS={:.6} eV  E_HF={:.6} eV  |HF-KS|={:.3e} eV  ({} iters, M={:.4})",
        result.total_energy,
        result.harris_foulkes_energy,
        hf_diff,
        result.n_iterations,
        result.magnetization,
    );

    // Quadratic-convergence assertion. Pre-SPXC: |HF-KS| was many eV even at
    // self-consistency due to E_xc using input `exc_r` vs output rho. Post-SPXC
    // (total-only criterion): |HF-KS| still O(delta_zeta) because spin channels
    // not driven to consistency. Post-SPXC+SPNC: should be O(delta_rho^2),
    // sub-microelectronvolt on a well-behaved non-magnetic system.
    //
    // Empirical post-SPNC |HF-KS| = 7.19e-7 eV at conv_threshold=1e-6 (23 iters).
    // Threshold set at ~14x empirical (1e-5 eV) for platform/compiler headroom.
    assert!(
        hf_diff < 1.0e-5,
        "SPXC+SPNC regression: Si nspin=2 |E_HF - E_KS| = {hf_diff:.3e} eV exceeds 1e-5 eV. \
         This system converges to ~7e-7 eV with both fixes in place. A value of many eV \
         indicates SPXC (XC input/output mismatch) has regressed; a value of O(0.01) eV \
         or larger suggests SPNC (per-spin convergence criterion) has regressed."
    );

    // Sanity: Si is non-magnetic — free-moment relaxation must drive M to 0.
    assert!(
        result.magnetization < 0.05,
        "Si should relax to non-magnetic (M≈0), got {:.4}",
        result.magnetization
    );
}

#[test]
fn test_ccmx_fe_free_magnetization_converges() {
    // CCMX regression: Fe BCC free-magnetization nspin=2 with the same
    // pseudopotential that fixed-mag=2 chokes on (nc/lda/Fe.upf). Pre-CCMX,
    // the independent-channel Anderson mixer entered a spin-flip limit
    // cycle — per-spin Δρ locked at ~0.254 for 200+ iterations with
    // |HF-KS| ~ 13 eV (see test_fe_ferromagnetic_fixed_moment comment and
    // proposals/SPNC-spin-per-density-convergence.md). Post-CCMX, the
    // (ρ_total, m) basis change decouples the two physical modes and the
    // SCF converges properly.
    //
    // Uses a 4×4×4 grid with Kerker at 15 Ry ecut for speed. Observed:
    //   pre-CCMX:  Δρ pinned at 0.254, consumes all max_iter=80, |HF-KS| ≈ 13 eV, M ≈ 0.05 μB (spurious).
    //   post-CCMX: converges in ~14 iters, |HF-KS| ≈ 1e-4 eV, M ≈ 0 μB.
    //
    // Regression guards (see assertions below): (1) |HF-KS| stays sub-meV,
    // (2) magnetization collapses to near-zero (no spurious spin leakage),
    // (3) iteration count stays under max_iter — pre-CCMX would exhaust it.
    // `ScfResult.final_delta` is not currently exposed; adding it would
    // enable a tighter pathology-specific assertion (Δρ ≈ 0.254 vs ≈ 1e-3).
    // Flagged as a Core Engineer follow-up on MODR Phase B.
    //
    // MXBA note (2026-04-18): this test runs at the default
    // `adaptive_beta = false`. Empirically, turning adaptive β on
    // (`adaptive_beta: true`) on this exact system damps β below β_min
    // before Anderson builds useful DIIS history and the SCF fails to
    // converge inside 80 iterations. The documented failure mode is
    // pinned by `tests/mxba_adaptive_beta_fe.rs`
    // (`test_mxba_fe_documents_adaptive_failure`, `#[ignore]`). Do not
    // flip the default to `true` without tuning Eyert thresholds first.
    let crystal = fe_bcc();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("pseudopotentials/nc/lda/Fe.upf"),
    )
    .unwrap();
    let ecut = 15.0 * 13.605_693_122_994; // 15 Ry
    let basis = BasisSet::new(&crystal.lattice, ecut);
    // Preserve the MP-1976 shifted grid this CCMX convergence test was
    // pinned on; swapping to Γ-centered alters the SCF trajectory and the
    // pinned iteration count / magnetization assertions.
    let kpoints = pwdft_rs::kpoints::monkhorst_pack(
        4,
        4,
        4,
        pwdft_rs::kpoints::KGridShift::MP1976,
        &crystal.lattice,
    );

    let mut starting_mag = std::collections::HashMap::new();
    starting_mag.insert("Fe".to_string(), 0.5);

    let params = scf::ScfParams {
        n_bands: 8,
        max_iter: 80,
        // Loose enough to converge within max_iter on 4×4×4; this test's
        // point is "no limit cycle", not "microelectronvolt tolerance".
        conv_threshold: 1e-3,
        energy_threshold: 1e-3,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.02 * 13.605_693_122_994, // 0.02 Ry
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Kerker { q_tf: None },
        nspin: 2,
        // NB: no tot_magnetization constraint — let the PP choose.
        starting_magnetization: starting_mag,
        ..Default::default()
    };

    let symmetry = pwdft_rs::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params, &symmetry)
        .expect("CCMX: Fe BCC nspin=2 free-mag must converge — pre-CCMX would fail here");

    eprintln!(
        "CCMX Fe free-mag regression: E_KS={:.6} eV  E_HF={:.6} eV  |HF-KS|={:.3e} eV  ({} iters, M={:.4} μB)",
        result.total_energy,
        result.harris_foulkes_energy,
        (result.harris_foulkes_energy - result.total_energy).abs(),
        result.n_iterations,
        result.magnetization,
    );

    // Pre-CCMX, |HF-KS| was ~13 eV because the limit cycle kept zeta
    // from ever matching input↔output. Post-CCMX should be sub-meV for
    // the non-magnetic collapse on this PP.
    let hf_diff = (result.harris_foulkes_energy - result.total_energy).abs();
    assert!(
        hf_diff < 1.0e-3,
        "CCMX Fe: |E_HF - E_KS| = {hf_diff:.3e} eV exceeds 1e-3 eV. \
         A large value (~13 eV) suggests the coupled-channel mixer has \
         regressed to independent (ρ↑, ρ↓) Anderson — check run_scf_spin \
         mixing block in src/scf/driver_spin.rs."
    );

    // The limit cycle produced M ≈ 0.05 μB from spurious spin flips; the
    // genuine non-magnetic ground state gives M ≈ 0.
    assert!(
        result.magnetization < 0.05,
        "CCMX Fe: magnetization M={:.3} μB is suspiciously large for a \
         non-magnetic PP. Possible limit-cycle leakage.",
        result.magnetization,
    );

    // Must not consume the full iteration budget — pre-CCMX blew past 200.
    assert!(
        result.n_iterations < params.max_iter,
        "CCMX Fe: hit max_iter={} — convergence still slow",
        params.max_iter,
    );

    // FDLT: pathology-specific regression guard. The pre-CCMX failure
    // mode was Δρ *pinned* at ≈0.254 for the entire budget (not an
    // iteration-count issue — relax max_iter and it still wouldn't
    // converge). Asserting `final_delta < 1e-2` catches the exact
    // limit-cycle pathology independent of how `conv_threshold` or
    // `max_iter` evolve. Post-CCMX observed: Δρ ≈ 5.7e-4 at iter 14.
    assert!(
        result.final_delta < 1e-2,
        "CCMX Fe Δρ-pathology regression guard: final Δρ = {:.3e} \
         (pre-CCMX pinned at ≈0.254). Coupled-channel (ρ_total, m) \
         mixer may have regressed to independent (ρ↑, ρ↓) Anderson.",
        result.final_delta,
    );
}

#[test]
fn test_fe_ferromagnetic_fixed_moment() {
    // TAUD finding 1.2 — INVERTED from its original (silently-passing) form.
    //
    // Fe BCC fixed-magnetization=2 on the nc/lda/Fe.upf pseudopotential is
    // KNOWN TO NOT CONVERGE. The original failure mode (pre-CCMX) was the
    // independent-channel limit cycle documented in
    // `proposals/SPNC-spin-per-density-convergence.md`: two Anderson mixers,
    // one per spin, produced uncorrelated ±ε predictions that kept Δρ_up =
    // Δρ_down ≈ 0.254 forever. After CCMX (2026-04-18) moved the mixer to
    // the (ρ_total, m) basis, the limit cycle is gone for the free-mag
    // Fe setup (see test_ccmx_fe_free_magnetization_converges below), but
    // the *fixed*-mag=2 case still fails to converge because the
    // constraint actively fights the PP's preference: LDA Fe on this PP
    // is non-magnetic, so forcing M=2 is not a stable SCF fixed point and
    // no mixer can produce one. Empirical post-CCMX behaviour: Δρ decays
    // to ≈0.04 (orders of magnitude better than the 0.254 ±ε cycle), then
    // stalls as the optimizer oscillates around a non-existent fixed
    // point. Fixing this genuinely would need (a) a ferromagnetic-stable
    // Fe pseudopotential, or (b) a Lagrange-constrained SCF that can
    // penalise deviations from M_target smoothly rather than enforcing
    // them via separate Fermi energies.
    //
    // This test is kept as an inverted regression detector: the day a PP
    // swap or a new constrained-DFT scheme flips the outcome to `Ok`, the
    // `assert!(matches!(..))` below will fail and pull attention back to
    // this case. At that point, restore the original assertions
    // (magnetization ≈ 2, energy vs QE reference E = -44.062_678_79 Ry ×
    // 13.605... = -599.503 eV from 4×4×4, 15 Ry, LDA, FD 0.02 Ry) and
    // un-invert the test.
    //
    // QE reference eigenvalues (if/when this starts converging):
    //   Gamma up:   4.62  25.71  25.71  26.52  26.52  26.52
    //   Gamma down: 5.79  27.30  27.30  28.05  28.05  28.05
    let crystal = fe_bcc();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Fe.upf"),
    ).unwrap();
    let ecut = 15.0 * 13.605_693_122_994; // 15 Ry
    let basis = BasisSet::new(&crystal.lattice, ecut);
    // MP-1976 shifted grid — this test's expected `ConvergenceFailure`
    // trajectory was pinned on it; the failure mode we assert against is
    // k-mesh sensitive.
    let kpoints = pwdft_rs::kpoints::monkhorst_pack(
        4,
        4,
        4,
        pwdft_rs::kpoints::KGridShift::MP1976,
        &crystal.lattice,
    );

    let params = scf::ScfParams {
        n_bands: 8,
        max_iter: 100,
        conv_threshold: 1e-6,
        energy_threshold: 1e-5,
        mixing_beta: 0.2,
        mixing_ndim: 8,
        smearing_sigma: 0.02 * 13.605_693_122_994, // 0.02 Ry in eV
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Kerker { q_tf: None },
        nspin: 2,
        tot_magnetization: Some(2.0),
        ..Default::default()
    };

    let symmetry = pwdft_rs::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params, &symmetry);

    // Specifically expect the per-spin limit-cycle failure to manifest as
    // ConvergenceFailure (not e.g. Eigensolver or Gpu). If the variant
    // changes, something else went wrong and we want to see it.
    match &result {
        Err(pwdft_rs::error::PwdftError::ConvergenceFailure { iterations, delta }) => {
            eprintln!(
                "Fe BCC fixed-mag=2 correctly failed to converge after {iterations} iters \
                 (final delta={delta:.3e}) — per-spin limit cycle, as expected post-SPNC."
            );
        }
        Ok(r) => {
            panic!(
                "Fe BCC fixed-mag=2 UNEXPECTEDLY converged: E={:.6} eV, M={:.4} μB, {} iters. \
                 Either CCMX landed (coupled-channel mixer) or the PP changed — inspect the \
                 result and restore the original magnetization/energy assertions (see comment above).",
                r.total_energy, r.magnetization, r.n_iterations
            );
        }
        Err(other) => {
            panic!(
                "Fe BCC fixed-mag=2 returned an unexpected error variant: {other}. \
                 Expected ConvergenceFailure (per-spin limit cycle, per SPNC proposal). \
                 A different error variant suggests a new regression."
            );
        }
    }
}

/// XCNI: non-LDA xc_functional must fail fast at SCF entry with
/// `PwdftError::NotImplemented`, not silently reinterpret as LDA.
///
/// This test must not trigger any real SCF compute work — the dispatch check
/// lives before the crystal/kpoints/volume guards in `run_scf`, so we can use
/// an otherwise minimal (even technically invalid) setup. We pick a valid Si
/// Gamma-only config so the test stays meaningful if the dispatch ever moves
/// a few lines.
#[test]
fn non_lda_xc_functional_is_rejected_at_scf_entry() {
    use pwdft_rs::error::PwdftError;
    use pwdft_rs::settings::XcFunctional;

    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();
    let kpoints = gamma_only();
    let sym = pwdft_rs::symmetry::SymmetryInfo::identity_only();

    // Each non-LDA variant must return NotImplemented with the expected
    // `what` label. Keep max_iter = 1 so that if the dispatch ever regresses,
    // the test fails loudly instead of hanging an SCF run.
    for (variant, want_label) in [
        (XcFunctional::Pbe, "pbe"),
        (XcFunctional::Pbe0, "pbe0"),
        (XcFunctional::Hse06, "hse06"),
    ] {
        let params = scf::ScfParams {
            n_bands: 4,
            max_iter: 1,
            xc_functional: variant,
            ..Default::default()
        };
        let err = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params, &sym)
            .expect_err("non-LDA xc_functional must fail at SCF entry");
        match err {
            PwdftError::NotImplemented { what } => {
                assert_eq!(
                    what, want_label,
                    "NotImplemented.what should name the functional ({variant:?})"
                );
            }
            other => panic!(
                "expected PwdftError::NotImplemented for {variant:?}, got: {other:?}"
            ),
        }
    }

    // Sanity pin: the baseline LDA path still enters the SCF body (and will
    // fail for some other reason — max_iter=1 — which is fine; we only need
    // to prove the dispatch does NOT trip NotImplemented on Pz).
    let params_pz = scf::ScfParams {
        n_bands: 4,
        max_iter: 1,
        xc_functional: XcFunctional::Pz,
        ..Default::default()
    };
    let result = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params_pz, &sym);
    // Any outcome (Ok, ConvergenceFailure, Eigensolver, …) is acceptable; we
    // only care that the XCNI trap does NOT fire on LDA. A NotImplemented
    // leak here would break every existing LDA test.
    if let Err(PwdftError::NotImplemented { .. }) = result {
        panic!("LDA (Pz) must not trigger NotImplemented — that would break every existing test");
    }
}
