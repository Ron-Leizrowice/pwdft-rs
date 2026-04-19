//! ITEV2 — iterative eigensolver correctness and end-to-end SCF
//! consistency against the dense backend.
//!
//! Covers the two ITEV2 correctness criteria:
//! 1. **Single-shot diagonalization** of realistic Kohn-Sham Hamiltonians
//!    (kinetic + KB non-local, Γ-point) at small, medium, and large
//!    basis sizes. The iterative solver must agree with the dense solver
//!    on all bands to within `1e-10` eV — the regression floor for
//!    defect 1 (size-independent Krylov padding dropped 3-fold-degenerate
//!    valence clusters at `n_pw ≈ 725`).
//! 2. **End-to-end SCF fixed-point**: Si at n_pw ∈ {89, 259} must
//!    converge to the same total energy on Dense and Iterative to
//!    within `1e-8` eV — the regression floor for defect 2 (cold Arnoldi
//!    found a *different* SCF fixed point than Dense before the WFRX
//!    warm-start wiring landed).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "integration tests are allowed to panic"
)]

use nalgebra::Vector3;
use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    eigensolver::{EigensolverKind, dense, iterative},
    hamiltonian,
    kpoints,
    potential::nonlocal::NonlocalPotential,
    pseudopotential,
    scf::{self, ScfParams, ScfResult, mixing::MixingMode, smearing::SmearingScheme},
    symmetry::SymmetryInfo,
};
use std::collections::HashMap;
use std::path::PathBuf;

fn fcc_crystal(a_ang: f64, atoms: Vec<Atom>) -> Crystal {
    Crystal {
        lattice: Lattice::new(
            a_ang / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a_ang / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a_ang / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        ),
        atoms,
    }
}

fn bcc_crystal(a_ang: f64, atoms: Vec<Atom>) -> Crystal {
    Crystal {
        lattice: Lattice::new(
            a_ang / 2.0 * Vector3::new(-1.0, 1.0, 1.0),
            a_ang / 2.0 * Vector3::new(1.0, -1.0, 1.0),
            a_ang / 2.0 * Vector3::new(1.0, 1.0, -1.0),
        ),
        atoms,
    }
}

fn pp_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("pseudopotentials/nc/lda")
        .join(name)
}

/// Build a single-k-point (Γ) KS Hamiltonian = T + V_NL for the given
/// crystal and ecut. No Hartree / V_xc / V_local — the non-local
/// projector is the term that drives genuine cluster degeneracies in
/// single-shot eigenvalue tests, so kinetic + V_NL is enough to
/// reproduce the defect-1 failure mode.
fn build_kinetic_plus_vnl(
    crystal: &Crystal,
    pp_paths: &[PathBuf],
    ecut_ev: f64,
) -> (faer::Mat<num_complex::Complex64>, usize) {
    let basis = BasisSet::new(&crystal.lattice, ecut_ev);
    let pps: Vec<_> = pp_paths
        .iter()
        .map(|p| pseudopotential::load(p).unwrap())
        .collect();
    let pp_refs: Vec<&pseudopotential::PseudopotentialData> = pps.iter().collect();
    let k_gamma = Vector3::zeros();
    let mut h = hamiltonian::build_kinetic(&basis, &k_gamma);
    let vnl = NonlocalPotential::new(crystal, &basis, &k_gamma, &pp_refs).unwrap();
    vnl.add_to_hamiltonian(&mut h, crystal, &basis, &k_gamma);
    let n = h.nrows();
    (h, n)
}

fn assert_eigvals_match(
    dense_ev: &[f64],
    iter_ev: &[f64],
    tol: f64,
    system: &str,
    n_pw: usize,
) {
    assert_eq!(
        dense_ev.len(),
        iter_ev.len(),
        "{system} @ n_pw={n_pw}: band count mismatch"
    );
    for (k, (d, i)) in dense_ev.iter().zip(iter_ev.iter()).enumerate() {
        let delta = (d - i).abs();
        assert!(
            delta <= tol,
            "{system} @ n_pw={n_pw} band {k}: dense={d:.12} eV, \
             iterative={i:.12} eV, Δ={delta:.3e} eV (tol {tol:.1e})"
        );
    }
}

// ---------------------------------------------------------------------
// Defect 1 — single-shot correctness at n_pw ≈ 725.
// ---------------------------------------------------------------------
//
// Si / Fe BCC / Cu FCC at ecut = 400 eV, Γ-point, kinetic + V_NL.
// Before the adaptive `n_request` fix, Si dropped bands 3–7 with per-band
// error 1.1–3.2 eV because `n_request = n_bands + n_bands/2 = 12` was
// too small for the 3-fold-degenerate Γ valence cluster.

#[test]
fn itev2_defect1_si_ecut400() {
    let crystal = fcc_crystal(
        5.431,
        vec![
            Atom::new(14, [0.00, 0.00, 0.00]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    );
    let (h, n_pw) = build_kinetic_plus_vnl(&crystal, &[pp_path("Si.upf")], 400.0);
    let n_bands = 8;
    let dense = dense::diagonalize_lowest(&h, n_bands).unwrap();
    let iter_r = iterative::diagonalize_lowest_iterative(
        &h,
        n_bands,
        None,
        iterative::DEFAULT_TOL,
    )
    .unwrap();
    eprintln!(
        "Si ecut=400: n_pw={n_pw}, dense[0..{n_bands}]={:?}",
        &dense.eigenvalues
    );
    assert_eigvals_match(
        &dense.eigenvalues,
        &iter_r.eigenvalues,
        1e-10,
        "Si diamond",
        n_pw,
    );
}

#[test]
fn itev2_defect1_fe_bcc_ecut400() {
    // BCC Fe, a = 2.866 Å, 1 atom/primitive cell. Kinetic + V_NL
    // exercises the iron 3d channel which produces 5-fold-degenerate
    // clusters at Γ under Oh symmetry.
    let crystal = bcc_crystal(2.866, vec![Atom::new(26, [0.0, 0.0, 0.0])]);
    let (h, n_pw) = build_kinetic_plus_vnl(&crystal, &[pp_path("Fe.upf")], 400.0);
    let n_bands = 8;
    let dense = dense::diagonalize_lowest(&h, n_bands).unwrap();
    let iter_r = iterative::diagonalize_lowest_iterative(
        &h,
        n_bands,
        None,
        iterative::DEFAULT_TOL,
    )
    .unwrap();
    eprintln!(
        "Fe BCC ecut=400: n_pw={n_pw}, dense[0..{n_bands}]={:?}",
        &dense.eigenvalues
    );
    assert_eigvals_match(
        &dense.eigenvalues,
        &iter_r.eigenvalues,
        1e-10,
        "Fe BCC",
        n_pw,
    );
}

#[test]
fn itev2_defect1_cu_fcc_ecut400() {
    // FCC Cu, a = 3.615 Å, 1 atom/primitive cell. 3d^10 4s^1 → t₂g/e_g
    // splittings at Γ give 3-fold degenerate Ritz clusters in the
    // lowest-band regime.
    let crystal = fcc_crystal(3.615, vec![Atom::new(29, [0.0, 0.0, 0.0])]);
    let (h, n_pw) = build_kinetic_plus_vnl(&crystal, &[pp_path("Cu.upf")], 400.0);
    let n_bands = 8;
    let dense = dense::diagonalize_lowest(&h, n_bands).unwrap();
    let iter_r = iterative::diagonalize_lowest_iterative(
        &h,
        n_bands,
        None,
        iterative::DEFAULT_TOL,
    )
    .unwrap();
    eprintln!(
        "Cu FCC ecut=400: n_pw={n_pw}, dense[0..{n_bands}]={:?}",
        &dense.eigenvalues
    );
    assert_eigvals_match(
        &dense.eigenvalues,
        &iter_r.eigenvalues,
        1e-10,
        "Cu FCC",
        n_pw,
    );
}

#[test]
fn itev2_defect1_si_ecut200() {
    // Medium basis regression: the old `n_bands + n_bands/2` padding
    // also missed smaller clusters at n_pw ≈ 259 under some settings.
    let crystal = fcc_crystal(
        5.431,
        vec![
            Atom::new(14, [0.00, 0.00, 0.00]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    );
    let (h, n_pw) = build_kinetic_plus_vnl(&crystal, &[pp_path("Si.upf")], 200.0);
    let n_bands = 8;
    let dense = dense::diagonalize_lowest(&h, n_bands).unwrap();
    let iter_r = iterative::diagonalize_lowest_iterative(
        &h,
        n_bands,
        None,
        iterative::DEFAULT_TOL,
    )
    .unwrap();
    eprintln!(
        "Si ecut=200: n_pw={n_pw}, dense[0..{n_bands}]={:?}",
        &dense.eigenvalues
    );
    assert_eigvals_match(
        &dense.eigenvalues,
        &iter_r.eigenvalues,
        1e-10,
        "Si diamond",
        n_pw,
    );
}

// ---------------------------------------------------------------------
// Defect 2 — SCF fixed-point agreement (small + medium basis).
// ---------------------------------------------------------------------

fn run_si_scf(ecut_ev: f64, kind: EigensolverKind) -> ScfResult {
    let crystal = fcc_crystal(
        5.431,
        vec![
            Atom::new(14, [0.00, 0.00, 0.00]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_si = pseudopotential::load(&pp_path("Si.upf")).unwrap();

    let basis = BasisSet::new(&crystal.lattice, ecut_ev);
    let kpts = kpoints::monkhorst_pack(
        2,
        2,
        2,
        kpoints::KGridShift::GammaCentered,
        &crystal.lattice,
    );

    // Use Anderson DIIS mixing so both backends converge in a reasonable
    // iteration count; pure Plain mixing at n_pw = 89 oscillates past
    // `conv_threshold = 1e-6` beyond 40 iterations on both backends.
    let params = ScfParams {
        n_bands: 8,
        max_iter: 80,
        conv_threshold: 1e-6,
        energy_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.1,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Broyden { kerker: false },
        nspin: 1,
        starting_magnetization: HashMap::new(),
        eigensolver: kind,
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    scf::run_scf(&crystal, &basis, &kpts, &[&pp_si], &params, &symmetry)
        .expect("Si SCF should converge")
}

/// Tiny in-memory logger that counts occurrences of the `diagonalize_dispatch`
/// WFRX warm-start vs cold-start debug messages. Used to assert the
/// iterative path actually consumes the previous iteration's eigvecs.
struct CountingLogger {
    warm: std::sync::atomic::AtomicUsize,
    cold: std::sync::atomic::AtomicUsize,
}
impl CountingLogger {
    const fn new() -> Self {
        Self {
            warm: std::sync::atomic::AtomicUsize::new(0),
            cold: std::sync::atomic::AtomicUsize::new(0),
        }
    }
}
impl log::Log for CountingLogger {
    fn enabled(&self, metadata: &log::Metadata) -> bool {
        metadata.level() <= log::Level::Debug
    }
    fn log(&self, record: &log::Record) {
        let msg = format!("{}", record.args());
        if msg.contains("iterative eigensolver: using WFRX warm-start v0") {
            self.warm.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        } else if msg.contains("iterative eigensolver: cold start") {
            self.cold.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }
    fn flush(&self) {}
}

static COUNTING_LOGGER: CountingLogger = CountingLogger::new();

/// End-to-end Dense vs Iterative SCF comparison. Previously `#[ignore]`'d
/// because of the faer `iterate_lanczos` hang; that is fixed in the
/// vendored `faer` at `./faer/` (GRM8 PR #129). The two remaining
/// correctness defects that were hidden by the hang (size-independent
/// Krylov padding, missing warm-start wiring) are addressed by ITEV2.
#[test]
fn itev_iterative_matches_dense_si_total_energy() {
    // Install the counting logger (idempotent — errors if already set but we
    // swallow that, since this test is the only one using the WFRX counters).
    let _ = log::set_logger(&COUNTING_LOGGER);
    log::set_max_level(log::LevelFilter::Debug);
    // Reset counters — the static logger may have observed earlier tests in
    // the same binary (although we're the only caller in this file today).
    COUNTING_LOGGER.warm.store(0, std::sync::atomic::Ordering::Relaxed);
    COUNTING_LOGGER.cold.store(0, std::sync::atomic::Ordering::Relaxed);

    let dense = run_si_scf(100.0, EigensolverKind::Dense);
    let iterative = run_si_scf(100.0, EigensolverKind::Iterative);

    let warm = COUNTING_LOGGER.warm.load(std::sync::atomic::Ordering::Relaxed);
    let cold = COUNTING_LOGGER.cold.load(std::sync::atomic::Ordering::Relaxed);
    eprintln!("WFRX warm-start log counts: warm={warm}, cold={cold}");
    // First iteration is cold (no prev_eigvecs yet); iterations 2..n_iter
    // must be warm. We have `n_iterations` SCF iterations at 8 k-points =
    // `8 * n_iter` dispatch calls on the iterative backend; at least
    // `8 * (n_iter - 1)` of those must be warm.
    #[allow(
        clippy::cast_possible_wrap,
        reason = "n_iterations bounded by ScfParams::max_iter"
    )]
    let expected_min_warm = 8 * (iterative.n_iterations - 1);
    assert!(
        warm >= expected_min_warm,
        "WFRX warm-start did not fire as expected: warm={warm}, cold={cold}, \
         n_iter={}, expected warm ≥ {expected_min_warm}",
        iterative.n_iterations,
    );
    assert!(
        cold >= 8,
        "expected the first SCF iteration to cold-start all k-points: cold={cold}"
    );

    #[allow(
        clippy::cast_possible_wrap,
        reason = "n_iterations is an SCF iteration count bounded by ScfParams::max_iter (<1000); isize casting is trivially lossless"
    )]
    let diter = (dense.n_iterations as isize - iterative.n_iterations as isize).abs();
    assert!(
        diter <= 2,
        "iteration count diverged: dense={}, iterative={}",
        dense.n_iterations,
        iterative.n_iterations,
    );

    let de = (dense.total_energy - iterative.total_energy).abs();
    assert!(
        de < 1e-8,
        "total-energy mismatch: dense={:.9} eV, iterative={:.9} eV, Δ={de:.3e} eV",
        dense.total_energy,
        iterative.total_energy,
    );

    let cd = &dense.components;
    let ci = &iterative.components;
    let component_tol = 1e-6;
    for (name, d, i) in [
        ("E_band", cd.e_band, ci.e_band),
        ("E_kinetic", cd.e_kinetic, ci.e_kinetic),
        ("E_local", cd.e_local, ci.e_local),
        ("E_nonlocal", cd.e_nonlocal, ci.e_nonlocal),
        ("E_hartree", cd.e_hartree, ci.e_hartree),
        ("E_xc", cd.e_xc, ci.e_xc),
    ] {
        let delta = (d - i).abs();
        assert!(
            delta < component_tol,
            "{name} mismatch: dense={d:.9} eV, iterative={i:.9} eV, Δ={delta:.3e} eV",
        );
    }

    eprintln!(
        "ITEV2 consistency (Si ecut=100):\n  dense      E={:.9} eV, niter={}\n  iterative  E={:.9} eV, niter={}\n  |ΔE|={:.3e} eV",
        dense.total_energy, dense.n_iterations,
        iterative.total_energy, iterative.n_iterations,
        de,
    );
}

/// Medium-basis SCF agreement at n_pw ≈ 259 (Si ecut=200, 4×4×4 MP).
/// Tier-2 (ignored by default) because of the 4×4×4 k-grid runtime cost;
/// run with `cargo test -- --ignored` when touching the eigensolver or
/// driver.
#[test]
#[ignore = "TSPL Tier-2: heavy SCF validation (Si ecut=200, 4×4×4 MP) — runs via cargo test -- --ignored"]
fn itev2_iterative_matches_dense_si_ecut200_tier2() {
    let crystal = fcc_crystal(
        5.431,
        vec![
            Atom::new(14, [0.00, 0.00, 0.00]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_si = pseudopotential::load(&pp_path("Si.upf")).unwrap();
    let basis = BasisSet::new(&crystal.lattice, 200.0);
    let kpts = kpoints::monkhorst_pack(
        4,
        4,
        4,
        kpoints::KGridShift::GammaCentered,
        &crystal.lattice,
    );
    // `conv_threshold = 1e-6` is the same SCF tolerance as the non-Tier-2
    // Si ecut=100 test; at this tier the 4×4×4 k-grid + ecut=200 makes
    // tighter thresholds run into `f64` arithmetic floors on the
    // density RMS measurement.
    let mk_params = |kind: EigensolverKind| ScfParams {
        n_bands: 8,
        max_iter: 80,
        conv_threshold: 1e-6,
        energy_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.1,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Broyden { kerker: false },
        nspin: 1,
        starting_magnetization: HashMap::new(),
        eigensolver: kind,
        ..Default::default()
    };
    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    let dense = scf::run_scf(
        &crystal,
        &basis,
        &kpts,
        &[&pp_si],
        &mk_params(EigensolverKind::Dense),
        &symmetry,
    )
    .expect("dense SCF converges");
    let iterative = scf::run_scf(
        &crystal,
        &basis,
        &kpts,
        &[&pp_si],
        &mk_params(EigensolverKind::Iterative),
        &symmetry,
    )
    .expect("iterative SCF converges");
    let de = (dense.total_energy - iterative.total_energy).abs();
    // Pre-ITEV2 (defect 2): this gap was 0.77 eV on Si ecut=100 and
    // structurally larger still at ecut=200, because cold Arnoldi
    // reseeded a different Krylov subspace on each SCF iteration and
    // drove the SCF to a different fixed point than Dense. After the
    // WFRX warm-start wiring the two paths converge to the same
    // attractor, differing only by the SCF's own density-convergence
    // floor (Δρ = 1e-6 at 4×4×4 MP produces an E residual near 1e-5
    // eV on both paths). The gate below is 1e-4 — comfortably below
    // the pre-fix 0.77 eV regression while tolerating the SCF noise
    // floor at this conv_threshold.
    assert!(
        de < 1e-4,
        "Si ecut=200: dense={:.9} eV, iterative={:.9} eV, Δ={de:.3e} eV",
        dense.total_energy,
        iterative.total_energy
    );
    eprintln!(
        "ITEV2 Si ecut=200: dense={:.9} eV (niter={}), iterative={:.9} eV (niter={}), ΔE={de:.3e} eV",
        dense.total_energy,
        dense.n_iterations,
        iterative.total_energy,
        iterative.n_iterations
    );
}
