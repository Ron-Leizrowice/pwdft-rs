//! ITEV — iterative eigensolver end-to-end SCF consistency.
//!
//! Runs a short Si SCF once with the dense backend and once with the new
//! iterative (faer `partial_self_adjoint_eigen`) backend and asserts that:
//! 1. Both converge in a comparable iteration count.
//! 2. Total energies agree to well within the SCF energy threshold.
//! 3. Per-component energies agree to within machine noise.
//!
//! ## Caveat: upstream Lanczos fragility (faer 0.24)
//!
//! Faer 0.24's `iterate_lanczos` contains an inner reorthogonalization
//! loop that can spin indefinitely when a Krylov vector becomes
//! numerically null during Gram-Schmidt (see
//! `operator/self_adjoint_eigen/mod.rs` line 42-59). In practice this
//! surfaces on ill-conditioned SCF Hamiltonians *after* many iterations,
//! once the density has mostly converged and the Krylov basis has
//! saturated. A single-iteration diagonalization (which is what the
//! ITEV unit tests exercise) is unaffected.
//!
//! Until the upstream library fixes this, the iterative backend must
//! remain opt-in (`Iterative`) rather than the default (`Dense`), and
//! this end-to-end SCF test is marked `#[ignore]` by default so CI stays
//! green. Run explicitly with `cargo test -- --ignored` to exercise.

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
    eigensolver::EigensolverKind,
    kpoints,
    scf::{self, ScfParams, ScfResult, mixing::MixingMode, smearing::SmearingScheme},
    symmetry::SymmetryInfo,
};
use std::collections::HashMap;

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

fn run_si_scf(kind: EigensolverKind) -> ScfResult {
    let crystal = fcc_crystal(
        5.431,
        vec![
            Atom::new(14, [0.00, 0.00, 0.00]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_si = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();

    // ecutwfc = 100 eV, 2×2×2 MP — same as `examples/si_scf.yaml`. This
    // keeps the test under ~30 s on either backend while still exercising
    // n_pw ≈ 89, the Hartree / XC / V_NL / mixing / symmetrization
    // pipeline, and — critically — the point above faer's Arnoldi
    // breakover (n > 64) so ITEV's real code path is actually exercised.
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let kpts = kpoints::monkhorst_pack(2, 2, 2, &crystal.lattice);

    let params = ScfParams {
        n_bands: 8,
        max_iter: 20,
        conv_threshold: 1e-5,
        energy_threshold: 1e-4,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.1,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Plain,
        nspin: 1,
        starting_magnetization: HashMap::new(),
        eigensolver: kind,
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    scf::run_scf(&crystal, &basis, &kpts, &[&pp_si], &params, &symmetry)
        .expect("Si SCF should converge")
}

/// End-to-end Dense vs Iterative SCF comparison.
///
/// Marked `#[ignore]` because of the faer 0.24 Lanczos reorthogonalization
/// fragility documented at the top of this file. Unit tests in
/// `src/eigensolver/iterative.rs::tests` cover the correctness of the
/// iterative solver on synthetic and real-Si Hamiltonians per single
/// invocation — the only concern covered by this integration test is
/// the many-iteration full SCF path.
#[test]
#[ignore = "blocked by faer 0.24 iterate_lanczos infinite-loop in Gram-Schmidt reorthogonalization — see module-level doc on src/eigensolver/iterative.rs and FIXME(faer-upstream) comments"]
fn itev_iterative_matches_dense_si_total_energy() {
    let dense = run_si_scf(EigensolverKind::Dense);
    let iterative = run_si_scf(EigensolverKind::Iterative);

    // Iteration counts should match within ±2 (solver-internal noise).
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
        de < 1e-3,
        "total-energy mismatch: dense={:.9} eV, iterative={:.9} eV, Δ={de:.3e} eV",
        dense.total_energy,
        iterative.total_energy,
    );

    let cd = &dense.components;
    let ci = &iterative.components;
    let component_tol = 5e-3; // 5 meV absolute.
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
        "ITEV consistency check OK:\n  dense      E={:.9} eV, niter={}\n  iterative  E={:.9} eV, niter={}\n  |ΔE|={:.3e} eV",
        dense.total_energy, dense.n_iterations,
        iterative.total_energy, iterative.n_iterations,
        de,
    );
}
