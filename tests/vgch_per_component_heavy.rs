//! VGCH Phase 1a / VGCH-2 Part A — per-component energy accounting for
//! heavy-atom cells, extended to the full 8-system VGCH matrix.
//!
//! Follows the VGC5 methodology (see `tests/vgc5_per_component_si.rs` and
//! `proposals/completed/VGC5-per-component-energy-accounting.md`) to pin
//! the per-term breakdown of the heavy-atom SCF vs. QE 7.5. Historically
//! this file covered only Cu + Fe (Phase 1a diagnostic); VGCH-2 Part A
//! extends it to the 7 "VGCH class" systems (Fe, Cu, GaAs, NaCl, MgO, C
//! diamond, Al) plus Si as a control, and emits a joinable CSV for the
//! VGCH-2 PR comparison table.
//!
//! Tests are diagnostic, not QE-match assertions. The residual on the
//! heavy-atom cells is currently O(eV), far above per-component
//! tolerance; these tests exist to localize *which* per-term carries the
//! bulk of the residual so VGCH-2 Part B can target a fix.
//!
//! Most cells match `tests/qe_validation.rs` at the QE reference grid
//! where practical. The 2 exceptions use a reduced k-grid or nspin for
//! tractability:
//!
//! - `vgch_cu_per_component`: 4×4×4 (QE ref 8×8×8), nspin=1.
//! - `vgch_fe_per_component_8x8x8`: 8×8×8 nspin=1 (QE ref 8×8×8
//!   nspin=2 FM, collapses to NM anyway under PseudoDojo LDA Fe at
//!   ecut=15 Ry).
//!
//! Per-term CSV emission (VGCH-2 Part A): every test writes rows to
//! `<CARGO_TARGET_TMPDIR>/vgch2_per_term_trace_pwdft.csv` via the shared
//! `write_pwdft_terms` helper. The Python script
//! `scripts/validate/vgch2_per_term_trace.py` parses QE's equivalent
//! from `qe_validation/*.out`; the two CSVs are joined on
//! `(system, term_name)` for the PR-body comparison table.
//!
//! QE's printed decomposition (see `qe-7.5/PW/src/electrons.f90:1612-1621`)
//! is `total = one_electron + hartree + xc + ewald + (-TS)`. Our mapping
//! onto `EnergyComponents` is
//!
//! ```text
//!     qe.one_electron   <==>   e_kinetic + e_local + e_local_g0_shift + e_nonlocal
//!     qe.hartree        <==>   e_hartree
//!     qe.xc             <==>   e_xc              (bare, no double-counting)
//!     qe.ewald          <==>   e_ewald
//! ```
//!
//! Post-VGCH-SiEF-B1 (2026-04-19): `e_local_g0_shift` is always `0.0`
//! under the QE gauge convention — the `V_loc(G=0)·N_el` piece now
//! lives inside `e_local`. The sum on the left of the
//! `qe.one_electron` mapping above is therefore numerically unchanged,
//! but the split between `e_local` and `e_local_g0_shift` has
//! collapsed into the former.
//!
//! `e_vxc` is NOT mapped to a QE print (QE accumulates it into the
//! double-counting subtraction internally); we emit it alongside anyway
//! for the `E_band = T + V_loc + V_nl + 2·E_H + E_vxc` identity check.

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
use std::{collections::HashMap, io::Write, path::PathBuf, sync::Mutex};

const RY_TO_EV: f64 = 13.605_693_122_994;

// ---------------------------------------------------------------------------
// CSV sink (VGCH-2 Part A)
// ---------------------------------------------------------------------------
//
// Each per-component test appends rows to
// `<CARGO_TARGET_TMPDIR>/vgch2_per_term_trace_pwdft.csv` so
// `scripts/validate/vgch2_per_term_trace.py` can join against QE.
// A static mutex serializes writes if cargo runs tests in parallel.

static CSV_LOCK: Mutex<()> = Mutex::new(());

fn pwdft_csv_path() -> PathBuf {
    PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("vgch2_per_term_trace_pwdft.csv")
}

/// Append this system's per-term values to the shared pwdft CSV. Columns:
///
///     system, term_name, pwdft_value_eV
///
/// where `term_name` is one of `one_electron`, `hartree`, `xc`, `ewald`,
/// `total`, `e_band`, `e_kinetic`, `e_local`, `e_local_g0_shift`,
/// `e_nonlocal`, `e_vxc`. The header is written on the first write per
/// test invocation (we re-open the file in append mode; to regenerate a
/// fresh trace, delete the file between runs).
fn write_pwdft_terms(system: &str, result: &ScfResult) {
    let guard = CSV_LOCK.lock().unwrap();
    let path = pwdft_csv_path();
    let header_needed = !path.exists();
    let mut contents = String::new();
    if header_needed {
        contents.push_str("system,term_name,pwdft_value_eV\n");
    }
    let c = &result.components;
    let one_electron = c.e_kinetic + c.e_local + c.e_local_g0_shift + c.e_nonlocal;
    let rows: [(&str, f64); 11] = [
        ("one_electron", one_electron),
        ("hartree", c.e_hartree),
        ("xc", c.e_xc),
        ("ewald", c.e_ewald),
        ("total", result.total_energy),
        ("e_band", c.e_band),
        ("e_kinetic", c.e_kinetic),
        ("e_local", c.e_local),
        ("e_local_g0_shift", c.e_local_g0_shift),
        ("e_nonlocal", c.e_nonlocal),
        ("e_vxc", c.e_vxc),
    ];
    for (term, value) in rows {
        contents.push_str(&format!("{system},{term},{value:.6}\n"));
    }
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)
        .unwrap_or_else(|e| panic!("open {path:?}: {e}"));
    f.write_all(contents.as_bytes())
        .unwrap_or_else(|e| panic!("write {path:?}: {e}"));
    drop(guard);
}

// ---------------------------------------------------------------------------
// Crystal builders & PP loaders
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// QE reference data (parsed from `qe_validation/*.out`)
// ---------------------------------------------------------------------------

/// Per-term QE reference block in eV. Matches
/// `scripts/validate/vgch2_per_term_trace.py` column order.
struct QeReference {
    total: f64,
    one_electron: f64,
    hartree: f64,
    xc: f64,
    ewald: f64,
}

impl QeReference {
    const fn from_ry(total: f64, one_e: f64, hartree: f64, xc: f64, ewald: f64) -> Self {
        Self {
            total: total * RY_TO_EV,
            one_electron: one_e * RY_TO_EV,
            hartree: hartree * RY_TO_EV,
            xc: xc * RY_TO_EV,
            ewald: ewald * RY_TO_EV,
        }
    }
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

    // Post-TSEN: `total_energy` carries `−TS` (= `c.e_smearing`); include it
    // so the identity still closes to machine precision on metals.
    let e_sum = c.e_kinetic
        + c.e_local
        + c.e_local_g0_shift
        + c.e_nonlocal
        + c.e_hartree
        + c.e_xc
        + c.e_ewald
        + c.e_smearing;
    let sum_err = e_sum - result.total_energy;
    eprintln!(
        "  [self-check] Σ(components) = {e_sum:.6} eV, E_total = {:.6} eV, Δ = {sum_err:.2e} eV",
        result.total_energy
    );
}

fn run_and_assert_sum(
    label: &str,
    crystal: &Crystal,
    pps: &[&PseudopotentialData],
    params: &ScfParams,
    basis_ecut_ev: f64,
    nk: u32,
    qe: &QeReference,
) -> ScfResult {
    let basis = BasisSet::new(&crystal.lattice, basis_ecut_ev);
    let kpts =
        kpoints::monkhorst_pack(nk, nk, nk, kpoints::KGridShift::GammaCentered, &crystal.lattice);

    let symmetry = SymmetryInfo::from_crystal(crystal, 1e-5);
    let result = scf::run_scf(crystal, &basis, &kpts, pps, params, &symmetry)
        .unwrap_or_else(|e| panic!("{label} SCF failed: {e}"));

    print_side_by_side(label, &result, qe);

    // PCFX self-check — the direct-sum identity must hold independently
    // of how close the residual lands to QE. Post-TSEN the identity
    // includes `c.e_smearing` (= −TS).
    let c = &result.components;
    let e_sum = c.e_kinetic
        + c.e_local
        + c.e_local_g0_shift
        + c.e_nonlocal
        + c.e_hartree
        + c.e_xc
        + c.e_ewald
        + c.e_smearing;
    let sum_residual = (e_sum - result.total_energy).abs();
    assert!(
        sum_residual < 0.1,
        "PCFX self-check on {label}: Σ − E_total = {sum_residual:.2e} eV"
    );
    assert!(result.total_energy.is_finite());

    write_pwdft_terms(label, &result);
    result
}

fn default_heavy_params(n_bands: usize, mixing: MixingMode, nspin: usize) -> ScfParams {
    ScfParams {
        n_bands,
        max_iter: 150,
        conv_threshold: 1e-6,
        energy_threshold: 1e-5,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.02 * RY_TO_EV,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: mixing,
        nspin,
        starting_magnetization: HashMap::new(),
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------
// Si control — VGCH-2 Part A
// ---------------------------------------------------------------------------

/// VGCH-2 Part A — Si control cell. Si closes E_total to 33 meV (per
/// `tests/qe_validation.rs::test_si_diamond_energy_vs_qe`); this test
/// re-runs the same SCF config to capture Si's per-term breakdown in
/// the same CSV format as the VGCH-class systems, so the comparison
/// table in the PR body can tell us whether the VGCH residual
/// *structure* (which term carries it) differs from Si's.
#[test]
#[ignore = "VGCH-2 Part A diagnostic: Si per-component trace alongside VGCH systems; use `-- --ignored`"]
fn vgch2_si_per_component() {
    let crystal = fcc_crystal(
        5.431,
        vec![
            Atom::new(14, [0.00, 0.00, 0.00]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_si = load_pp("Si");

    let params = ScfParams {
        smearing_sigma: 0.01 * RY_TO_EV,
        ..default_heavy_params(8, MixingMode::Plain, 1)
    };

    // QE: qe_validation/si_scf.out — total=-17.02299344 Ry,
    // one-e=4.86744632, hartree=1.11010670, xc=-6.20301507, ewald=-16.79667313.
    let qe = QeReference::from_ry(
        -17.022_993_44,
        4.867_446_32,
        1.110_106_70,
        -6.203_015_07,
        -16.796_673_13,
    );
    run_and_assert_sum(
        "Si",
        &crystal,
        &[&pp_si],
        &params,
        15.0 * RY_TO_EV,
        4,
        &qe,
    );
}

// ---------------------------------------------------------------------------
// C diamond — VGCH light-atom 1.45 eV residual
// ---------------------------------------------------------------------------

/// VGCH-2 Part A — C diamond. 1.45 eV residual, opposite-sign split
/// across one-electron and Hartree (Δ one-e = +1.76, Δ E_H = −0.59
/// per Phase 1b observation).
#[test]
#[ignore = "VGCH-2 Part A diagnostic: C diamond per-component trace; use `-- --ignored`"]
fn vgch2_c_diamond_per_component() {
    let crystal = fcc_crystal(
        3.567,
        vec![
            Atom::new(6, [0.00, 0.00, 0.00]),
            Atom::new(6, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_c = load_pp("C");

    let params = ScfParams {
        smearing_sigma: 0.01 * RY_TO_EV,
        ..default_heavy_params(8, MixingMode::Broyden { kerker: true }, 1)
    };

    // QE: qe_validation/c_diamond_scf.out — total=-23.84343910 Ry,
    // one-e=8.50287341, hartree=1.82838120, xc=-8.60280600, ewald=-25.57188769.
    let qe = QeReference::from_ry(
        -23.843_439_10,
        8.502_873_41,
        1.828_381_20,
        -8.602_806_00,
        -25.571_887_69,
    );
    run_and_assert_sum(
        "C_diamond",
        &crystal,
        &[&pp_c],
        &params,
        30.0 * RY_TO_EV,
        4,
        &qe,
    );
}

// ---------------------------------------------------------------------------
// Al FCC — VGCH light-atom 75 meV residual
// ---------------------------------------------------------------------------

/// VGCH-2 Part A — Al FCC. 75 meV residual at basis-converged
/// ecut=24 Ry (8×8×8).
#[test]
#[ignore = "VGCH-2 Part A diagnostic: Al per-component trace; use `-- --ignored`"]
fn vgch2_al_per_component() {
    let crystal = fcc_crystal(4.05, vec![Atom::new(13, [0.0, 0.0, 0.0])]);
    let pp_al = load_pp("Al");

    let params = default_heavy_params(6, MixingMode::Kerker { q_tf: None }, 1);

    // QE: qe_validation/al_fcc_scf.out — total=-4.72724484 Ry,
    // one-e=2.88539588, hartree=0.00731666, xc=-2.22048340, ewald=-5.39205235.
    let qe = QeReference::from_ry(
        -4.727_244_84,
        2.885_395_88,
        0.007_316_66,
        -2.220_483_40,
        -5.392_052_35,
    );
    run_and_assert_sum(
        "Al",
        &crystal,
        &[&pp_al],
        &params,
        24.0 * RY_TO_EV,
        8,
        &qe,
    );
}

// ---------------------------------------------------------------------------
// Heavy-atom systems
// ---------------------------------------------------------------------------

/// VGCH Phase 1a — Cu FCC per-component diagnostic.
///
/// Matches `qe_validation/cu_fcc_scf.in` cell (a=3.610 Å FCC, nspin=1)
/// but runs at a reduced 4×4×4 k-grid + ecut=25 Ry to keep the test
/// tractable. Prints the per-component decomposition alongside the QE
/// reference. Used by VGCH Phase 1b + VGCH-2 to localize the ~16 eV
/// Cu residual.
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

    let params = default_heavy_params(14, MixingMode::Kerker { q_tf: None }, 1);

    // QE reference at 8×8×8 (from qe_validation/cu_fcc_scf.out):
    //   total = -356.73602869 Ry; one_e = -149.48760888 Ry;
    //   hartree = 76.47616676 Ry; xc = -41.09518682 Ry; ewald = -242.62085475 Ry.
    // Note: this test's 4×4×4 sampling will differ from QE's 8×8×8 by
    // k-sampling noise (~O(100 meV)), which is small against the ~16 eV
    // residual this diagnostic is meant to localize.
    let qe = QeReference::from_ry(
        -356.736_028_69,
        -149.487_608_88,
        76.476_166_76,
        -41.095_186_82,
        -242.620_854_75,
    );
    run_and_assert_sum(
        "Cu_FCC",
        &crystal,
        &[&pp_cu],
        &params,
        25.0 * RY_TO_EV,
        4,
        &qe,
    );
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

    let params = default_heavy_params(12, MixingMode::Kerker { q_tf: None }, 1);

    let qe = QeReference::from_ry(
        -224.917_449_34,
        -50.856_512_13,
        26.641_158_55,
        -28.904_021_34,
        -171.779_065_80,
    );
    run_and_assert_sum(
        "Fe_BCC_FM",
        &crystal,
        &[&pp_fe],
        &params,
        15.0 * RY_TO_EV,
        8,
        &qe,
    );
}

/// VGCH-2 Part A — GaAs zincblende. 33.6 eV residual (Ga Z=31, As Z=33).
#[test]
#[ignore = "VGCH-2 Part A diagnostic: GaAs per-component trace; use `-- --ignored`"]
fn vgch2_gaas_per_component() {
    let crystal = fcc_crystal(
        5.653,
        vec![
            Atom::new(31, [0.00, 0.00, 0.00]),
            Atom::new(33, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_ga = load_pp("Ga");
    let pp_as = load_pp("As");

    let params = ScfParams {
        smearing_sigma: 0.01 * RY_TO_EV,
        ..default_heavy_params(18, MixingMode::Plain, 1)
    };

    let qe = QeReference::from_ry(
        -307.928_895_02,
        -107.675_419_91,
        74.102_634_65,
        -75.978_778_17,
        -198.372_229_22,
    );
    run_and_assert_sum(
        "GaAs",
        &crystal,
        &[&pp_ga, &pp_as],
        &params,
        20.0 * RY_TO_EV,
        4,
        &qe,
    );
}

/// VGCH-2 Part A — NaCl rocksalt. 7.7 eV residual (Cl Z=17).
#[test]
#[ignore = "VGCH-2 Part A diagnostic: NaCl per-component trace; use `-- --ignored`"]
fn vgch2_nacl_per_component() {
    let crystal = fcc_crystal(
        5.614,
        vec![
            Atom::new(11, [0.00, 0.00, 0.00]),
            Atom::new(17, [0.50, 0.50, 0.50]),
        ],
    );
    let pp_na = load_pp("Na");
    let pp_cl = load_pp("Cl");

    let params = ScfParams {
        smearing_sigma: 0.01 * RY_TO_EV,
        ..default_heavy_params(12, MixingMode::Plain, 1)
    };

    let qe = QeReference::from_ry(
        -119.779_703_03,
        -64.294_393_42,
        34.920_824_87,
        -21.274_136_00,
        -69.131_998_46,
    );
    run_and_assert_sum(
        "NaCl",
        &crystal,
        &[&pp_na, &pp_cl],
        &params,
        25.0 * RY_TO_EV,
        4,
        &qe,
    );
}

/// VGCH-2 Part A — MgO rocksalt. 10.1 eV residual (Mg 2s/2p semicore PP).
#[test]
#[ignore = "VGCH-2 Part A diagnostic: MgO per-component trace; use `-- --ignored`"]
fn vgch2_mgo_per_component() {
    let crystal = fcc_crystal(
        4.212,
        vec![
            Atom::new(12, [0.00, 0.00, 0.00]),
            Atom::new(8, [0.50, 0.50, 0.50]),
        ],
    );
    let pp_mg = load_pp("Mg");
    let pp_o = load_pp("O");

    let params = ScfParams {
        smearing_sigma: 0.01 * RY_TO_EV,
        ..default_heavy_params(10, MixingMode::Plain, 1)
    };

    let qe = QeReference::from_ry(
        -147.235_477_68,
        -65.020_630_59,
        35.329_918_47,
        -22.765_634_85,
        -94.779_130_59,
    );
    run_and_assert_sum(
        "MgO",
        &crystal,
        &[&pp_mg, &pp_o],
        &params,
        30.0 * RY_TO_EV,
        4,
        &qe,
    );
}
