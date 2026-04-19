//! VGCH Phase 1a — per-component energy accounting for heavy-atom cells.
//!
//! Follows the VGC5 methodology (see `tests/vgc5_per_component_si.rs` and
//! `proposals/completed/VGC5-per-component-energy-accounting.md`) to pin
//! the per-term breakdown of the heavy-atom SCF vs. QE 7.5. This is the
//! diagnostic infrastructure for VGCH; it does NOT assert QE agreement
//! because the residual is currently >> per-component tolerance.
//!
//! The prints below are the primary deliverable: they tell us whether
//! the residual lives in (kinetic, local, non-local, Hartree, XC, Ewald)
//! or in the one-electron sum, so the fix target is known before any
//! code change in `src/pseudopotential/`.
//!
//! Cells match `tests/qe_validation.rs` where practical. Fe uses nspin=1
//! 4×4×4 (QE ref was nspin=2 8×8×8) to avoid CCMX's extra iteration
//! budget; Cu uses nspin=1 4×4×4 (QE ref was nspin=1 8×8×8, ecut=25 Ry)
//! for the same reason. Both are informative about per-component
//! distribution of residual even without matching k-grids — the residual
//! scales smoothly with ecut/k, not with an abrupt transition.

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
    kpoints,
    pseudopotential::PseudopotentialData,
    scf::{self, ScfParams, ScfResult, mixing::MixingMode, smearing::SmearingScheme},
    symmetry::SymmetryInfo,
};
use std::collections::HashMap;

const RY_TO_EV: f64 = 13.605_693_122_994;

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

fn bcc_crystal(a_ang: f64, atom: Atom) -> Crystal {
    Crystal {
        lattice: Lattice::new(
            a_ang / 2.0 * Vector3::new(-1.0, 1.0, 1.0),
            a_ang / 2.0 * Vector3::new(1.0, -1.0, 1.0),
            a_ang / 2.0 * Vector3::new(1.0, 1.0, -1.0),
        ),
        atoms: vec![atom],
    }
}

fn load_pp(element: &str) -> PseudopotentialData {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("pseudopotentials/nc/lda")
        .join(format!("{element}.upf"));
    pwdft_rs::pseudopotential::load(&path)
        .unwrap_or_else(|e| panic!("failed to load {}: {e}", path.display()))
}

struct QeReference {
    total: f64,
    one_electron: f64,
    hartree: f64,
    xc: f64,
    ewald: f64,
}

fn print_side_by_side(label: &str, result: &ScfResult, qe: &QeReference) {
    let c = &result.components;
    let ours_one_electron = c.e_kinetic + c.e_local + c.e_local_g0_shift + c.e_nonlocal;

    eprintln!("\n===== {label} : per-component decomposition (eV) =====");
    eprintln!(
        "  {:<24}  {:>14}  {:>14}  {:>12}",
        "term", "pwdft-rs", "QE", "Δ (ours−QE)"
    );
    eprintln!(
        "  {:<24}  {:>14}  {:>14}  {:>12}",
        "----", "--------", "--", "-----------"
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14}  {:>12}",
        "E_kinetic", c.e_kinetic, "(in 1e)", ""
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14}  {:>12}",
        "E_local (G≠0)", c.e_local, "(in 1e)", ""
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14}  {:>12}",
        "E_local(G=0)*N_el", c.e_local_g0_shift, "(in 1e)", ""
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14}  {:>12}",
        "E_nonlocal", c.e_nonlocal, "(in 1e)", ""
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14.6}  {:>+12.6}",
        "one-electron (sum)",
        ours_one_electron,
        qe.one_electron,
        ours_one_electron - qe.one_electron
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14.6}  {:>+12.6}",
        "E_hartree", c.e_hartree, qe.hartree, c.e_hartree - qe.hartree
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14.6}  {:>+12.6}",
        "E_xc", c.e_xc, qe.xc, c.e_xc - qe.xc
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14.6}  {:>+12.6}",
        "E_ewald", c.e_ewald, qe.ewald, c.e_ewald - qe.ewald
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14.6}  {:>+12.6}",
        "E_total (E_KS)",
        result.total_energy,
        qe.total,
        result.total_energy - qe.total
    );

    let e_sum = c.e_kinetic
        + c.e_local
        + c.e_local_g0_shift
        + c.e_nonlocal
        + c.e_hartree
        + c.e_xc
        + c.e_ewald;
    let sum_err = e_sum - result.total_energy;
    eprintln!(
        "  [self-check] Σ(components) = {e_sum:.6} eV, E_total = {:.6} eV, Δ = {sum_err:.2e} eV",
        result.total_energy
    );
}

/// VGCH Phase 1a — Cu FCC per-component diagnostic.
///
/// Matches `qe_validation/cu_fcc_scf.in` cell (a=3.610 Å FCC, nspin=1)
/// but runs at a reduced 4×4×4 k-grid + ecut=25 Ry to keep the test
/// tractable. Prints the per-component decomposition alongside the QE
/// reference. Used by VGCH Phase 1b to localize the ~16 eV Cu residual.
///
/// Not a QE-match assertion — the residual is too large; this test
/// exists to generate the diagnostic table in the VGCH PR body.
#[test]
#[ignore = "VGCH Phase 1a diagnostic: Cu per-component run is intentionally slow; use `-- --ignored` to print the table"]
fn vgch_cu_per_component() {
    // FCC a = 6.8219 Bohr ≈ 3.610 Å (matches qe_validation/cu_fcc_scf.in).
    let a_bohr = 6.8219_f64;
    let a_ang = a_bohr * 0.529_177_210_903;
    let crystal = fcc_crystal(a_ang, vec![Atom::new(29, [0.0, 0.0, 0.0])]);
    let pp_cu = load_pp("Cu");

    let basis = BasisSet::new(&crystal.lattice, 25.0 * RY_TO_EV);
    let kpts =
        kpoints::monkhorst_pack(4, 4, 4, kpoints::KGridShift::GammaCentered, &crystal.lattice);

    let params = ScfParams {
        n_bands: 14,
        max_iter: 150,
        conv_threshold: 1e-6,
        energy_threshold: 1e-5,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.02 * RY_TO_EV,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Kerker { q_tf: None },
        nspin: 1,
        starting_magnetization: HashMap::new(),
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(&crystal, &basis, &kpts, &[&pp_cu], &params, &symmetry)
        .expect("Cu SCF should converge");

    // QE reference at 8×8×8 (from qe_validation/cu_fcc_scf.out):
    //   total = -356.73602869 Ry; one_e = -149.48760888 Ry;
    //   hartree = 76.47616676 Ry; xc = -41.09518682 Ry; ewald = -242.62085475 Ry.
    // Note: this test's 4×4×4 sampling will differ from QE's 8×8×8 by
    // k-sampling noise (~O(100 meV)), which is small against the ~16 eV
    // residual this diagnostic is meant to localize.
    let qe = QeReference {
        total: -356.736_028_69 * RY_TO_EV,
        one_electron: -149.487_608_88 * RY_TO_EV,
        hartree: 76.476_166_76 * RY_TO_EV,
        xc: -41.095_186_82 * RY_TO_EV,
        ewald: -242.620_854_75 * RY_TO_EV,
    };
    print_side_by_side("Cu FCC", &result, &qe);

    // Sanity: components sum to total (PCFX self-check).
    let c = &result.components;
    let e_sum = c.e_kinetic
        + c.e_local
        + c.e_local_g0_shift
        + c.e_nonlocal
        + c.e_hartree
        + c.e_xc
        + c.e_ewald;
    let sum_residual = (e_sum - result.total_energy).abs();
    assert!(
        sum_residual < 0.1,
        "PCFX self-check on Cu: Σ − E_total = {sum_residual:.2e} eV"
    );
    assert!(result.total_energy.is_finite());
}

/// VGCH Phase 1a — Fe BCC per-component diagnostic (8×8×8 nspin=1).
///
/// Matches `qe_validation/fe_bcc_fm_scf.in` at the full 8×8×8 k-grid
/// but with nspin=1 (the QE ref uses nspin=2, which collapses to M=0
/// anyway with PseudoDojo Fe LDA at ecut=15 Ry — see
/// `qe_validation/reference_data.toml`). Prints per-component
/// decomposition for the VGCH Phase 1 PR body.
///
/// The `vgc5_fe_per_component` test in `tests/vgc5_per_component_si.rs`
/// runs the same Fe cell at 4×4×4 for regression-guard purposes; this
/// test runs at the QE-matched 8×8×8 grid to expose how the residual
/// scales with k-sampling.
#[test]
#[ignore = "VGCH Phase 1a diagnostic: Fe 8×8×8 per-component run is intentionally slow; use `-- --ignored`"]
fn vgch_fe_per_component_8x8x8() {
    let crystal = bcc_crystal(2.87, Atom::new(26, [0.0, 0.0, 0.0]));
    let pp_fe = load_pp("Fe");

    let basis = BasisSet::new(&crystal.lattice, 15.0 * RY_TO_EV);
    let kpts =
        kpoints::monkhorst_pack(8, 8, 8, kpoints::KGridShift::GammaCentered, &crystal.lattice);

    let params = ScfParams {
        n_bands: 12,
        max_iter: 150,
        conv_threshold: 1e-6,
        energy_threshold: 1e-5,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.02 * RY_TO_EV,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Kerker { q_tf: None },
        nspin: 1,
        starting_magnetization: HashMap::new(),
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(&crystal, &basis, &kpts, &[&pp_fe], &params, &symmetry)
        .expect("Fe 8×8×8 SCF should converge");

    let qe = QeReference {
        total: -224.917_449_34 * RY_TO_EV,
        one_electron: -50.856_512_13 * RY_TO_EV,
        hartree: 26.641_158_55 * RY_TO_EV,
        xc: -28.904_021_34 * RY_TO_EV,
        ewald: -171.779_065_80 * RY_TO_EV,
    };
    print_side_by_side("Fe BCC (nspin=1, 8×8×8)", &result, &qe);

    let c = &result.components;
    let e_sum = c.e_kinetic
        + c.e_local
        + c.e_local_g0_shift
        + c.e_nonlocal
        + c.e_hartree
        + c.e_xc
        + c.e_ewald;
    let sum_residual = (e_sum - result.total_energy).abs();
    assert!(
        sum_residual < 0.1,
        "PCFX self-check on Fe 8×8×8: Σ − E_total = {sum_residual:.2e} eV"
    );
    assert!(result.total_energy.is_finite());
}
