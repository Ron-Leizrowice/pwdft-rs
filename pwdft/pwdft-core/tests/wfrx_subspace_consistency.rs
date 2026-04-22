//! WFRX Phase 1 — subspace warm-start end-to-end SCF consistency.
//!
//! Runs a short Si SCF once with the stock dense eigensolver and once
//! with `wfrx_subspace = true` (Rayleigh-Ritz warm-start on top of the
//! same dense solver) and asserts that the two SCF outputs agree to
//! well within the SCF energy threshold. This is the non-negotiable
//! correctness gate for the WFRX proposal (§2 acceptance criterion:
//! final energy must match within 1e-8 eV on a converged Si run).
//!
//! Because the subspace path has an internal residual fallback to a
//! full diagonalization whenever the projected subspace is not invariant
//! enough, the two runs should converge to essentially the same fixed
//! point — any drift is driven solely by the ~1e-6 residual tolerance,
//! far below the conv_threshold used here.

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
    eigensolver::EigensolverKind,
    kpoints,
    pseudopotential::UpfPseudoPotential,
    scf::{self, ScfParams, ScfResult, mixing::MixingMode, smearing::SmearingScheme},
    symmetry::SymmetryInfo,
};

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

fn run_si_scf(wfrx: bool) -> ScfResult {
    let crystal = fcc_crystal(
        5.431,
        vec![Atom::new(14, [0.00, 0.00, 0.00]), Atom::new(14, [0.25, 0.25, 0.25])],
    );
    let pp_si = UpfPseudoPotential::load(
        &std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();

    // Same knobs as `itev_iterative_eigensolver.rs`: ecutwfc = 100 eV,
    // 2×2×2 MP, tight enough to converge in ~20 iters.
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let kpts = kpoints::monkhorst_pack(2, 2, 2, kpoints::KGridShift::GammaCentered, &crystal.lattice);

    let params = ScfParams {
        n_bands: 8,
        max_iter: 120,
        // Tight convergence: the WFRX proposal's 1e-8 eV total-energy
        // agreement criterion can only be met when the SCF itself is
        // converged well below that. We use Broyden mixing for fast
        // convergence on Si and push density to 1e-8 e/Å³ and energy to
        // 1e-9 eV so that the residual SCF oscillation is safely below
        // the WFRX subspace residual gate (1e-6).
        conv_threshold: 1e-8,
        energy_threshold: 1e-9,
        mixing_beta: 0.7,
        mixing_ndim: 8,
        smearing_sigma: 0.1,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Broyden { kerker: false },
        nspin: 1,
        starting_magnetization: HashMap::new(),
        eigensolver: EigensolverKind::Dense,
        wfrx_subspace: wfrx,
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    scf::run_scf(&crystal, &basis, &kpts, &[&pp_si], &params, &symmetry).expect("Si SCF should converge")
}

/// Final SCF energy must match between the reference (WFRX off) and
/// the warm-started (WFRX on) runs to within 1e-8 eV. The WFRX
/// proposal's §2 criterion is "match pre-WFRX to within 1e-8 eV on a
/// converged Si run"; we make that concrete here.
#[test]
#[ignore = "TSPL Tier-2: runs Si SCF twice (dense vs WFRX subspace) at ecut=100 conv=1e-8; run with cargo test -- --ignored when touching scf/, eigensolver/, or WFRX paths"]
fn test_wfrx_subspace_matches_dense_reference_total_energy() {
    let reference = run_si_scf(false);
    let warm = run_si_scf(true);

    let de = (reference.total_energy - warm.total_energy).abs();
    println!(
        "WFRX total_energy agreement: |ΔE| = {de:.3e} eV (reference = {:.12}, warm = {:.12}, \
         ref iters = {}, warm iters = {})",
        reference.total_energy, warm.total_energy, reference.n_iterations, warm.n_iterations,
    );
    assert!(
        de < 1e-8,
        "WFRX total_energy must match reference to 1e-8 eV; got |ΔE| = {de:.3e} eV \
         (reference = {:.12}, warm = {:.12})",
        reference.total_energy,
        warm.total_energy
    );
}

/// Free energy (Mermin) must also match; this is the quantity a user
/// reports for metals at finite T. Independent of total_energy only
/// modulo the entropy TS, which should be indistinguishable at this
/// smearing.
#[test]
#[ignore = "TSPL Tier-2: runs Si SCF twice (dense vs WFRX subspace) at ecut=100 conv=1e-8; run with cargo test -- --ignored when touching scf/, eigensolver/, or WFRX paths"]
fn test_wfrx_subspace_matches_dense_reference_free_energy() {
    let reference = run_si_scf(false);
    let warm = run_si_scf(true);

    let de = (reference.free_energy - warm.free_energy).abs();
    assert!(
        de < 1e-8,
        "WFRX free_energy must match reference to 1e-8 eV; got |ΔE| = {de:.3e} eV \
         (reference = {:.12}, warm = {:.12})",
        reference.free_energy,
        warm.free_energy
    );
}

/// Fermi energy and band-structure eigenvalues must match. The lowest
/// occupied bands carry O(1) eV weight on the total, so any drift here
/// would dwarf the 1e-8 eV total-energy gate — this is a tighter per-
/// eigenvalue sanity check.
#[test]
#[ignore = "TSPL Tier-2: runs Si SCF twice (dense vs WFRX subspace) at ecut=100 conv=1e-8; run with cargo test -- --ignored when touching scf/, eigensolver/, or WFRX paths"]
fn test_wfrx_subspace_matches_dense_reference_eigenvalues() {
    let reference = run_si_scf(false);
    let warm = run_si_scf(true);

    assert!(
        (reference.fermi_energy - warm.fermi_energy).abs() < 1e-6,
        "WFRX fermi_energy drift: reference = {:.10}, warm = {:.10}",
        reference.fermi_energy,
        warm.fermi_energy
    );
    assert_eq!(
        reference.eigenvalues.len(),
        warm.eigenvalues.len(),
        "k-point count mismatch"
    );
    for (ik, (ref_k, warm_k)) in reference.eigenvalues.iter().zip(warm.eigenvalues.iter()).enumerate() {
        assert_eq!(ref_k.len(), warm_k.len(), "band count mismatch at ik={ik}");
        for (ib, (&r, &w)) in ref_k.iter().zip(warm_k.iter()).enumerate() {
            let dev = (r - w).abs();
            assert!(
                dev < 1e-6,
                "eigenvalue drift at ik={ik}, ib={ib}: ref={r:.10}, warm={w:.10}, \
                 Δ={dev:.3e}"
            );
        }
    }
}
