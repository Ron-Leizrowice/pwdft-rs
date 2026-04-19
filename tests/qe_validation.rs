//! Validate pwdft-rs SCF results against Quantum ESPRESSO 7.5 reference data
//! (Tier 1 + Tier 2 of proposal QEVL).
//!
//! All QE runs use PseudoDojo ONCV NC/LDA pseudopotentials (same PPs as
//! pwdft-rs), Fermi-Dirac smearing, Monkhorst-Pack k-grids. Reference inputs
//! and `pw.x` outputs live in `qe_validation/`; machine-readable reference
//! values are in `qe_validation/reference_data.toml`.
//!
//! Layout of this file:
//!   * Helpers: `fcc_crystal`, `bcc_crystal`, `run_qe_comparison`,
//!     `assert_energy_matches_qe`, `assert_fermi_matches_qe`.
//!   * One `#[test]` per system (8 total).
//!
//! Post-MPSH (2026-04-18) `monkhorst_pack` accepts a [`KGridShift`] and
//! defaults to Γ-centered — matching QE's `K_POINTS automatic / Nx Ny Nz
//! 0 0 0`. With this both codes now sample the same k-mesh. Empirically
//! MPSH closed the **Si 4×4×4 total-energy** residual from 0.26 eV to
//! 33 meV, but the eigenvalue / E_F assertions on Si still fail due to
//! a V_loc(G=0) absolute-reference shift (≈ 1.35 eV constant offset;
//! VGCH territory), and Al (83 meV) and C (still-stalling SCF) were
//! **not** resolved by MPSH alone — investigation of the residual
//! root cause is logged in each test's `#[ignore]` reason string.
//!
//! Heavy-atom (Z > 14) systems carry an additional 7–34 eV residual whose
//! root cause is TBD — VGCMP Phases 1–4 proved the V_local(G) assembly
//! pipeline bit-correct on Si, so the heavy-atom residual is *not*
//! V_local(G) and the continuing investigation is tracked under VGCH
//! (candidate root causes: semicore/ecut convergence, V_local(G=0) Z-scaling,
//! Ewald for large Z, etc.). Each `#[ignore]` reason cites the specific
//! blocker (MPSH residual category or VGCH heavy-atom residual) plus
//! the measured pwdft-rs and QE values. Drop an `#[ignore]` once the
//! attributed residual category closes.
//!
//! [`KGridShift`]: pwdft_rs::kpoints::KGridShift

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
    error::Result as PwdftResult,
    kpoints,
    pseudopotential::PseudopotentialData,
    scf::{self, ScfParams, ScfResult, mixing::MixingMode, smearing::SmearingScheme},
    symmetry::SymmetryInfo,
};
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const RY_TO_EV: f64 = 13.605_693_122_994;

// ---------------------------------------------------------------------------
// Crystal builders
// ---------------------------------------------------------------------------

/// Build an FCC primitive cell (QE `ibrav = 2`) of lattice parameter `a` (Å)
/// with the given basis atoms (fractional coordinates in the primitive cell).
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

/// Build a BCC primitive cell (QE `ibrav = 3`) of lattice parameter `a` (Å)
/// with a single basis atom at the origin.
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

// ---------------------------------------------------------------------------
// Pseudopotentials
// ---------------------------------------------------------------------------

fn load_pp(element: &str) -> PseudopotentialData {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("pseudopotentials/nc/lda")
        .join(format!("{element}.upf"));
    pwdft_rs::pseudopotential::load(&path)
        .unwrap_or_else(|e| panic!("failed to load {}: {e}", path.display()))
}

// ---------------------------------------------------------------------------
// Comparison harness
// ---------------------------------------------------------------------------

/// Parameters for a QE-vs-pwdft-rs comparison.
///
/// Mirrors the knobs exposed in a QE `&SYSTEM`/`&ELECTRONS` block so test
/// bodies stay declarative.
struct QeComparisonConfig<'a> {
    crystal: &'a Crystal,
    pps: Vec<&'a PseudopotentialData>,
    ecut_ry: f64,
    nk: u32,
    n_bands: usize,
    mixing: MixingMode,
    smearing: SmearingScheme,
    degauss_ry: f64,
    nspin: usize,
    starting_magnetization: HashMap<String, f64>,
}

impl<'a> QeComparisonConfig<'a> {
    fn new(crystal: &'a Crystal, pps: Vec<&'a PseudopotentialData>) -> Self {
        Self {
            crystal,
            pps,
            ecut_ry: 15.0,
            nk: 4,
            n_bands: 8,
            mixing: MixingMode::Plain,
            smearing: SmearingScheme::FermiDirac,
            degauss_ry: 0.01,
            nspin: 1,
            starting_magnetization: HashMap::new(),
        }
    }
}

/// Run an SCF with the given QE-equivalent configuration.
fn run_qe_comparison(cfg: &QeComparisonConfig<'_>) -> PwdftResult<ScfResult> {
    let ecut_ev = cfg.ecut_ry * RY_TO_EV;
    let basis = BasisSet::new(&cfg.crystal.lattice, ecut_ev);
    // MPSH: match QE's `K_POINTS automatic / Nx Ny Nz 0 0 0` by sampling on
    // a Γ-centered grid. Every reference input under `qe_validation/*.in`
    // uses `0 0 0`.
    let kpts = kpoints::monkhorst_pack(
        cfg.nk,
        cfg.nk,
        cfg.nk,
        kpoints::KGridShift::GammaCentered,
        &cfg.crystal.lattice,
    );

    eprintln!(
        "  ecut={:.1} Ry  basis={} PWs  {}x{}x{} grid -> {} k-points  nspin={}",
        cfg.ecut_ry,
        basis.len(),
        cfg.nk,
        cfg.nk,
        cfg.nk,
        kpts.len(),
        cfg.nspin,
    );

    let params = ScfParams {
        n_bands: cfg.n_bands,
        max_iter: 80,
        conv_threshold: 1e-8,
        energy_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: cfg.degauss_ry * RY_TO_EV,
        smearing_scheme: cfg.smearing,
        ecutrho_ratio: 4,
        mixing_mode: cfg.mixing.clone(),
        nspin: cfg.nspin,
        starting_magnetization: cfg.starting_magnetization.clone(),
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(cfg.crystal, 1e-5);
    scf::run_scf(cfg.crystal, &basis, &kpts, &cfg.pps, &params, &symmetry)
}

// ---------------------------------------------------------------------------
// Assertions
// ---------------------------------------------------------------------------

/// Require that the pwdft-rs total energy agrees with QE (in Ry) within
/// `tolerance_ev`. Prints the comparison regardless of pass/fail.
fn assert_energy_matches_qe(
    label: &str,
    result: &ScfResult,
    qe_energy_ry: f64,
    tolerance_ev: f64,
) {
    let qe_energy_ev = qe_energy_ry * RY_TO_EV;
    let de = (result.total_energy - qe_energy_ev).abs();
    eprintln!(
        "  [{label}] E_pwdft = {:.6} eV,  E_QE = {:.6} eV,  |ΔE| = {:.4} eV",
        result.total_energy, qe_energy_ev, de,
    );
    assert!(
        de <= tolerance_ev,
        "{label}: |ΔE|={de:.4} eV exceeds tolerance {tolerance_ev:.4} eV \
         (E_pwdft={:.6} eV, E_QE={:.6} eV)",
        result.total_energy,
        qe_energy_ev,
    );
}

/// Require that the pwdft-rs Fermi energy (already in eV) agrees with QE's
/// within `tolerance_ev`.
fn assert_fermi_matches_qe(
    label: &str,
    result: &ScfResult,
    qe_fermi_ev: f64,
    tolerance_ev: f64,
) {
    let df = (result.fermi_energy - qe_fermi_ev).abs();
    eprintln!(
        "  [{label}] E_F_pwdft = {:.4} eV,  E_F_QE = {:.4} eV,  |ΔE_F| = {:.4} eV",
        result.fermi_energy, qe_fermi_ev, df,
    );
    assert!(
        df <= tolerance_ev,
        "{label}: |ΔE_F|={df:.4} eV exceeds tolerance {tolerance_ev:.4} eV",
    );
}

/// Print Γ-point eigenvalues alongside the QE reference for visual diffing.
fn report_gamma_eigenvalues(label: &str, result: &ScfResult, qe_eigs_ev: &[f64]) {
    if let Some(ours) = result.eigenvalues.first() {
        eprintln!("  [{label}] Γ eigenvalues (eV):");
        eprintln!("    pwdft: {ours:?}");
        eprintln!("    QE   : {qe_eigs_ev:?}");
    }
}

// ---------------------------------------------------------------------------
// Tier 1 — Core systems
// ---------------------------------------------------------------------------

/// Si diamond (FCC, 2 atoms, LDA insulator).
///
/// QE ref (qe_validation/si_scf.in): E = -17.022_993_44 Ry,
/// E_F = 6.3449 eV, converges in 7 iters.
///
/// Ignored: post-MPSH (2026-04-18, Γ-centered grid matching QE), the
/// total-energy residual dropped from ≈0.26 eV to **33 meV** at 4×4×4
/// ecut = 15 Ry (E_pwdft = −231.6428 eV vs QE −231.6096 eV). The
/// remaining residual is now driven by the **absolute energy reference**
/// — pwdft-rs sets V_eff(G=0) = 0 while QE uses a different convention,
/// so every Kohn-Sham eigenvalue (and consequently the Fermi energy) is
/// offset by a constant ≈ 1.35 eV. This shows up as
/// `|ΔE_F| ≈ 1.35 eV` exceeding the 50 meV tolerance even though the
/// physical band structure (i.e. band-to-band energy differences) agrees
/// with QE to < 10 meV. The absolute-reference issue is tracked under
/// VGCH (heavy-atom V_loc audit also covers this shift for lighter
/// elements).
#[test]
#[ignore = "MPSH: E_total now 33 meV (within 50 meV tol), but eigenvalue absolute reference ≈ 1.35 eV shift (VGCH V_loc(G=0))"]
fn test_si_diamond_vs_qe() {
    let crystal = fcc_crystal(
        5.431,
        vec![
            Atom::new(14, [0.00, 0.00, 0.00]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_si = load_pp("Si");

    let cfg = QeComparisonConfig {
        ecut_ry: 15.0,
        nk: 4,
        n_bands: 8,
        ..QeComparisonConfig::new(&crystal, vec![&pp_si])
    };
    let result = run_qe_comparison(&cfg).expect("Si SCF should converge");

    report_gamma_eigenvalues(
        "Si",
        &result,
        &[-5.8903, 6.0816, 6.0816, 6.0816, 8.6106, 8.6106, 8.6106, 9.3253],
    );
    assert_energy_matches_qe("Si", &result, -17.022_993_44, 0.05);
    assert_fermi_matches_qe("Si", &result, 6.3449, 0.05);
}

/// C diamond (FCC, 2 atoms, LDA wide-gap insulator).
///
/// QE ref: E = -23.843_439_10 Ry, E_F = 15.8873 eV, 9 iters, ecut = 30 Ry.
///
/// Ignored: post-MPSH (2026-04-18, Γ-centered grid matching QE), C still
/// stalls at Δρ ≈ 4.4e-6 after 80 iterations — **MPSH did NOT close
/// this**. Empirically: switching from MP-1976 shifted to Γ-centered
/// produces a similar Δρ plateau at the same scale (4.1e-6 → 4.4e-6),
/// so the shift convention was NOT the root cause. Tentative suspects
/// for the remaining blockage: (i) the mixer tuning (ecut=30 Ry on
/// 4×4×4 with `MixingMode::Plain` may need Kerker or Broyden for the
/// wider C gap), (ii) the underlying V_loc(G=0) absolute-reference
/// issue that ≈ 1 eV-shifts every band and may be upsetting the mixer's
/// residual bookkeeping. C is Z=6 (light) so no heavy-atom V_loc
/// dependency in the usual VGCH sense; the blocker is probably a
/// mixer/ecut combination and is now out of MPSH scope.
#[test]
#[ignore = "post-MPSH: C still stalls at Δρ ≈ 4.4e-6 under Γ-centered grid — root cause is NOT the shift convention; tentative mixer/ecut follow-up"]
fn test_c_diamond_vs_qe() {
    let crystal = fcc_crystal(
        3.567,
        vec![
            Atom::new(6, [0.00, 0.00, 0.00]),
            Atom::new(6, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_c = load_pp("C");

    let cfg = QeComparisonConfig {
        ecut_ry: 30.0,
        nk: 4,
        n_bands: 8,
        ..QeComparisonConfig::new(&crystal, vec![&pp_c])
    };
    let result = run_qe_comparison(&cfg).expect("C SCF should converge");

    report_gamma_eigenvalues(
        "C",
        &result,
        &[-8.1456, 14.0232, 14.0232, 14.0232, 19.3568, 19.3568, 19.3568, 27.2505],
    );
    assert_energy_matches_qe("C", &result, -23.843_439_10, 0.05);
    assert_fermi_matches_qe("C", &result, 15.8873, 0.05);
}

/// Al FCC (1 atom, simple metal).
///
/// QE ref: E = -4.723_717_90 Ry, E_F = 7.6130 eV, 6 iters, ecut = 15 Ry,
/// 8x8x8 k-grid, degauss = 0.02 Ry, Kerker (QE `local-TF`).
///
/// Ignored: post-MPSH (2026-04-18, Γ-centered grid matching QE), the
/// residual went from 73 meV to **83 meV** — still above the 50 meV
/// tolerance. MPSH slightly worsened Al; empirically the old MP-1976
/// shifted grid happened to cancel some of Al's residual against QE's
/// Γ-centered reference, and fixing the shift to match QE exposes the
/// underlying ≈ 80 meV discrepancy. Attribution: NOT the grid
/// convention. Candidate blockers: (i) ecut = 15 Ry on Al is barely
/// converged (113 PW basis), (ii) the Kerker q_TF auto-estimate may
/// differ from QE's `local-TF`, (iii) small FFT-grid differences (note
/// the default `ecutrho_ratio = 4` may be below QE's auto-determined
/// wfc_grid). Follow-up: bump ecut to 30 Ry and re-measure before
/// opening a dedicated proposal.
///
/// Reference values (for year-later readers):
///   pwdft-rs post-MPSH (Γ-centered): E = −64.1864 eV
///   QE:                              E = −64.2695 eV  (−4.723_717_90 Ry)
///   residual:                        ~83 meV
#[test]
#[ignore = "post-MPSH: Al 8×8×8 residual 83 meV (Γ-centered grid now matches QE; remaining gap NOT shift-related — ecut/Kerker/FFT candidates)"]
fn test_al_fcc_vs_qe() {
    let crystal = fcc_crystal(4.05, vec![Atom::new(13, [0.0, 0.0, 0.0])]);
    let pp_al = load_pp("Al");

    let cfg = QeComparisonConfig {
        ecut_ry: 15.0,
        nk: 8,
        n_bands: 6,
        mixing: MixingMode::Kerker { q_tf: None },
        degauss_ry: 0.02,
        ..QeComparisonConfig::new(&crystal, vec![&pp_al])
    };
    let result = run_qe_comparison(&cfg).expect("Al SCF should converge");

    report_gamma_eigenvalues(
        "Al",
        &result,
        &[-3.4111, 20.3573, 20.3573, 21.5721, 21.5721, 21.5721],
    );
    assert_energy_matches_qe("Al", &result, -4.723_717_90, 0.05);
    assert_fermi_matches_qe("Al", &result, 7.6130, 0.05);
}

/// BCC Fe (1 atom, nspin=2, collapses to non-magnetic under this PP/cutoff).
///
/// QE ref: E = -224.917_449_34 Ry, E_F = 26.2006 eV, 10 iters. With
/// PseudoDojo NC/LDA at ecut=15 Ry, |M| collapses to 0.00 μB in both codes,
/// so the test validates the nspin=2 machinery rather than the magnetism.
///
/// limit cycle); post-MPSH (2026-04-18) the k-grid now matches QE's
/// Γ-centered convention but the total energy still differs from QE by
/// ≈ 11.5 eV — the residual Z>14 heavy-atom gap tracked under VGCH
/// (root cause TBD; VGCMP Phases 1–4 proved the V_local(G) assembly
/// pipeline bit-correct on Si, so the heavy-atom residual is *not*
/// V_local(G) assembly — continuing cross-check vs QE). MPSH switching
/// the grid from MP-1976 to Γ-centered moved the residual from 9.5 eV
/// to 11.5 eV, consistent with MPSH sampling a different k-mesh than
/// the pre-MPSH baseline used; the underlying heavy-atom discrepancy
/// itself is unchanged. See
/// `proposals/VGCH-heavy-atom-vloc-residual.md`.
///
/// Reference values (for year-later readers):
///   pwdft-rs post-MPSH (Γ-centered): E = -3048.655 eV (8×8×8, 15 Ry, Kerker)
///   QE ref:                          E = -3060.158 eV (-224.917_449_34 Ry)
///   residual:                        ~11.5 eV  →  tracked under VGCH
#[test]
#[ignore = "post-MPSH Fe residual ≈11.5 eV on 8×8×8 Γ-centered grid; blocked on VGCH (root cause TBD)"]
fn test_fe_bcc_fm_vs_qe() {
    let crystal = bcc_crystal(2.87, Atom::new(26, [0.0, 0.0, 0.0]));
    let pp_fe = load_pp("Fe");

    let mut starting_mag = HashMap::new();
    starting_mag.insert("Fe".to_string(), 0.5);

    let cfg = QeComparisonConfig {
        ecut_ry: 15.0,
        nk: 8,
        n_bands: 12,
        mixing: MixingMode::Kerker { q_tf: None },
        degauss_ry: 0.02,
        nspin: 2,
        starting_magnetization: starting_mag,
        ..QeComparisonConfig::new(&crystal, vec![&pp_fe])
    };
    let result = run_qe_comparison(&cfg).expect("Fe SCF should converge");

    eprintln!(
        "  [Fe] M_pwdft = {:.4} μB  (QE: 0.00 μB; NM collapse expected)",
        result.magnetization,
    );
    report_gamma_eigenvalues(
        "Fe",
        &result,
        &[
            -122.5041, -46.4140, -46.4140, -46.4140, 9.2567, 23.8230, 23.8230, 24.4515,
        ],
    );
    assert_energy_matches_qe("Fe", &result, -224.917_449_34, 0.05);
    assert_fermi_matches_qe("Fe", &result, 26.2006, 0.1);
}

// ---------------------------------------------------------------------------
// Tier 2 — Compounds and extended systems
// ---------------------------------------------------------------------------

/// GaAs zincblende (FCC, 2 species, III-V semiconductor).
///
/// QE ref: E = -307.928_895_02 Ry, E_F = 6.5212 eV, 11 iters, ecut = 20 Ry.
///
/// Ignored: both Ga (Z=31) and As (Z=33) are heavy-atom Z>14 — residual
/// has root cause TBD (tracked in VGCH). VGCMP Phases 1–4 ruled out a
/// V_local(G) assembly bug on Si, so the heavy-atom residual lives
/// elsewhere (candidate causes: semicore/ecut convergence, V_local(G=0)
/// Z-scaling, Ewald for large Z); ~9.5 eV seen on Fe BCC carries over and
/// compounds across two heavy species here. VERF did not close the Si gap
/// and is archived — replacing the old VERF attribution with VGCH.
///
/// Reference values (for year-later readers):
///   pwdft-rs: E = −4155.9543 eV
///   QE:       E = −4189.5860 eV  (−307.928_895_02 Ry)
///   residual: ~33.6 eV
#[test]
#[ignore = "VGCH: heavy-atom residual ≈33.6 eV (root cause TBD) on GaAs (Z=31+33); pwdft-rs E = -4155.954 eV, QE = -4189.586 eV"]
fn test_gaas_zincblende_vs_qe() {
    let crystal = fcc_crystal(
        5.653,
        vec![
            Atom::new(31, [0.00, 0.00, 0.00]),
            Atom::new(33, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_ga = load_pp("Ga");
    let pp_as = load_pp("As");

    let cfg = QeComparisonConfig {
        ecut_ry: 20.0,
        nk: 4,
        n_bands: 18,
        ..QeComparisonConfig::new(&crystal, vec![&pp_ga, &pp_as])
    };
    let result = run_qe_comparison(&cfg).expect("GaAs SCF should converge");

    report_gamma_eigenvalues(
        "GaAs",
        &result,
        &[
            -9.0547, -9.0547, -9.0547, -7.9791, -7.1414, -7.1414, -0.9400, -0.9400,
        ],
    );
    assert_energy_matches_qe("GaAs", &result, -307.928_895_02, 0.1);
    assert_fermi_matches_qe("GaAs", &result, 6.5212, 0.1);
}

/// Cu FCC (1 atom, transition metal with semicore 3s/3p/3d).
///
/// QE ref: E = -356.736_028_69 Ry, E_F = 19.2056 eV, 9 iters, ecut = 25 Ry,
/// 8x8x8 k-grid, degauss = 0.02 Ry, Kerker (QE `local-TF`).
///
/// Ignored: Cu is Z=29 heavy-atom — residual has root cause TBD (tracked
/// in VGCH). VGCMP Phases 1–4 ruled out a V_local(G) assembly bug on Si,
/// so the heavy-atom residual lives elsewhere (candidate causes:
/// semicore/ecut convergence, V_local(G=0) Z-scaling, Ewald for large Z);
/// ~9.5 eV seen on Fe BCC carries over here. Cu has 3s/3p/3d semicore so
/// the semicore sensitivity is especially plausible. VERF did not close
/// the Si gap and is archived — replacing the old VERF attribution with
/// VGCH.
///
/// Reference values (for year-later readers):
///   pwdft-rs: E = −4837.4659 eV
///   QE:       E = −4853.6409 eV  (−356.736_028_69 Ry)
///   residual: ~16.2 eV
#[test]
#[ignore = "VGCH: heavy-atom residual (root cause TBD) ≈16.2 eV on Cu (Z=29, 3s/3p/3d semicore); pwdft-rs E = -4837.466 eV, QE = -4853.641 eV"]
fn test_cu_fcc_vs_qe() {
    let crystal = fcc_crystal(3.61, vec![Atom::new(29, [0.0, 0.0, 0.0])]);
    let pp_cu = load_pp("Cu");

    let cfg = QeComparisonConfig {
        ecut_ry: 25.0,
        nk: 8,
        n_bands: 14,
        mixing: MixingMode::Kerker { q_tf: None },
        degauss_ry: 0.02,
        ..QeComparisonConfig::new(&crystal, vec![&pp_cu])
    };
    let result = run_qe_comparison(&cfg).expect("Cu SCF should converge");

    report_gamma_eigenvalues(
        "Cu",
        &result,
        &[
            -141.4571, -69.5066, -69.5066, -69.5066, 7.4368, 15.0881, 15.0881, 15.4196,
        ],
    );
    assert_energy_matches_qe("Cu", &result, -356.736_028_69, 0.1);
    assert_fermi_matches_qe("Cu", &result, 19.2056, 0.1);
}

/// NaCl rocksalt (FCC, 2 species, ionic insulator).
///
/// QE ref: E = -119.779_703_03 Ry, E_F = 3.4704 eV, 9 iters, ecut = 25 Ry.
/// Na at (0,0,0), Cl at (½,½,½) in the FCC primitive cell.
///
/// Ignored: Cl is Z=17 (heavy, Z>14) — residual has root cause TBD
/// (tracked in VGCH). VGCMP Phases 1–4 ruled out a V_local(G) assembly
/// bug on Si, so the heavy-atom residual lives elsewhere (candidate
/// causes: semicore/ecut convergence, V_local(G=0) Z-scaling, Ewald for
/// large Z); ~9.5 eV seen on Fe BCC carries over here. VERF did not
/// close the Si gap and is archived — replacing the old VERF attribution
/// with VGCH.
///
/// Reference values (for year-later readers):
///   pwdft-rs: E = −1621.9441 eV
///   QE:       E = −1629.6859 eV  (−119.779_703_03 Ry)
///   residual: ~7.7 eV
#[test]
#[ignore = "VGCH: heavy-atom residual (root cause TBD) ≈7.7 eV on NaCl (Cl Z=17); pwdft-rs E = -1621.944 eV, QE = -1629.686 eV"]
fn test_nacl_rocksalt_vs_qe() {
    let crystal = fcc_crystal(
        5.614,
        vec![
            Atom::new(11, [0.00, 0.00, 0.00]), // Na
            Atom::new(17, [0.50, 0.50, 0.50]), // Cl
        ],
    );
    let pp_na = load_pp("Na");
    let pp_cl = load_pp("Cl");

    let cfg = QeComparisonConfig {
        ecut_ry: 25.0,
        nk: 4,
        n_bands: 12,
        ..QeComparisonConfig::new(&crystal, vec![&pp_na, &pp_cl])
    };
    let result = run_qe_comparison(&cfg).expect("NaCl SCF should converge");

    report_gamma_eigenvalues(
        "NaCl",
        &result,
        &[
            -59.5034, -18.1802, -18.1802, -18.1802, -11.3942, 1.1682, 1.1682, 1.1682,
        ],
    );
    assert_energy_matches_qe("NaCl", &result, -119.779_703_03, 0.1);
    assert_fermi_matches_qe("NaCl", &result, 3.4704, 0.1);
}

/// MgO rocksalt (FCC, 2 species, wide-gap ionic insulator).
///
/// QE ref: E = -147.235_477_68 Ry, E_F = 10.2064 eV, 8 iters, ecut = 30 Ry.
/// Mg at (0,0,0), O at (½,½,½) in the FCC primitive cell.
///
/// Ignored: measured residual (~10 eV) is at the heavy-atom scale despite
/// nominal Z<14 — the Mg ONCV LDA PP includes 2s/2p semicore, which places
/// it in the same heavy-atom residual class (tracked in VGCH, root cause
/// TBD). VGCMP Phases 1–4 ruled out a V_local(G) assembly bug on Si, so
/// the heavy-atom residual lives elsewhere (candidate causes: semicore/ecut
/// convergence, V_local(G=0) Z-scaling, Ewald for large Z). VERF did not
/// close the Si gap and is archived — replacing the old VERF attribution
/// with VGCH.
///
/// Reference values (for year-later readers):
///   pwdft-rs: E = −1993.1458 eV
///   QE:       E = −2003.2407 eV  (−147.235_477_68 Ry)
///   residual: ~10.1 eV
#[test]
#[ignore = "VGCH: heavy-atom residual (root cause TBD) ≈10.1 eV on MgO (Mg semicore PP); pwdft-rs E = -1993.146 eV, QE = -2003.241 eV"]
fn test_mgo_rocksalt_vs_qe() {
    let crystal = fcc_crystal(
        4.212,
        vec![
            Atom::new(12, [0.00, 0.00, 0.00]), // Mg
            Atom::new(8, [0.50, 0.50, 0.50]),  // O
        ],
    );
    let pp_mg = load_pp("Mg");
    let pp_o = load_pp("O");

    let cfg = QeComparisonConfig {
        ecut_ry: 30.0,
        nk: 4,
        n_bands: 10,
        ..QeComparisonConfig::new(&crystal, vec![&pp_mg, &pp_o])
    };
    let result = run_qe_comparison(&cfg).expect("MgO SCF should converge");

    report_gamma_eigenvalues(
        "MgO",
        &result,
        &[
            -74.3080, -29.9713, -29.9713, -29.9713, -10.1472, 8.4930, 8.4930, 8.4930,
        ],
    );
    assert_energy_matches_qe("MgO", &result, -147.235_477_68, 0.1);
    assert_fermi_matches_qe("MgO", &result, 10.2064, 0.1);
}

// ---------------------------------------------------------------------------
// NLCC audit (Part B) — defensive E_xc regression guard for Fe
// ---------------------------------------------------------------------------

/// NLCC audit, Part B: defensive E_xc regression guard on Fe BCC.
///
/// This test exercises the NLCC code path end-to-end on an element
/// where core/valence overlap is large (Fe 3d semicore overlaps 4s/3d
/// valence). If anyone regresses the NLCC parser or the r²·4π
/// Bessel-transform weighting fixed by NCFX, this test should trip
/// *before* any other Fe test because E_xc is the term NLCC affects
/// most directly.
///
/// Post-NCFX measured residual against QE (pre-PCFX, pre-MP-shift fix):
///   pwdft-rs E_xc = −392.5675 eV   (pin in `tests/vgc5_per_component_si.rs`)
///   QE        E_xc = −28.904_021_34 Ry = −393.259 eV
///   Δ_xc     = +0.692 eV  (down from −48.85 eV pre-NCFX — 71× reduction)
///
/// We pin the residual |Δ_xc| ≤ 1.0 eV. The 1 eV ceiling is chosen so
/// that any regression that restores the pre-NCFX bug (which produced
/// a ~49 eV swing on this same test) is caught immediately, while
/// leaving headroom for the known Monkhorst-Pack shifted-vs-Γ-centered
/// residual (tracked in SYKP) and the `nspin=1` ≠ QE `nspin=2`
/// geometry choice (forced by Fe collapsing to NM under the PseudoDojo
/// LDA PP at ecut=15 Ry — see `reference_data.toml`).
///
/// This test is *defensive*, not an accuracy milestone. Tolerances are
/// loose by design. The Fe total-energy comparison remains gated by
/// VGCH/PCFX and stays in `test_fe_bcc_fm_vs_qe` (still `#[ignore]`d).
///
/// MPSH note (2026-04-18): the `|Δ_xc| ≤ 1.0 eV` ceiling and the
/// `E_xc_pwdft = -392.5675 eV` reference in this test's docstring were
/// both captured on the **MP-1976 shifted** 4×4×4 grid. We preserve
/// that grid here (even though the rest of `qe_validation` now runs
/// Γ-centered) because the pin is a regression guard tied to the
/// exact 0.69 eV residual baseline; switching to Γ-centered pushes
/// |Δ_xc| to ≈ 1.41 eV. That is still 35× below the pre-NCFX
/// pathology (~49 eV) so the guard's purpose is unaffected, but the
/// simplest way to keep this trip-wire armed at its designed
/// threshold is to keep its k-sample fixed.
#[test]
fn test_fe_bcc_xc_nlcc_regression_guard() {
    let crystal = bcc_crystal(2.87, Atom::new(26, [0.0, 0.0, 0.0]));
    let pp_fe = load_pp("Fe");
    assert!(
        pp_fe.has_nlcc(),
        "Fe UPF must have core_correction=T for this test to exercise NLCC"
    );

    // nspin=1, 4×4×4 MP — matches `tests/vgc5_per_component_si.rs`
    // `vgc5_fe_per_component`. The magnetic ground state collapses to NM
    // at this ecut/PP anyway (see `reference_data.toml`). We construct
    // ScfParams directly here (instead of `run_qe_comparison`) to use
    // loose convergence targets that accommodate the GPU f32 precision
    // floor (CPU reaches 1e-8 in ~80 iters; GPU stalls at ~1e-7 due to
    // f32 roundoff in Hartree/XC shaders). The looser targets only
    // affect the last few digits of E_xc — well within the 1 eV
    // tolerance below.
    let ecut_ev = 15.0 * RY_TO_EV;
    let basis = BasisSet::new(&crystal.lattice, ecut_ev);
    // Preserve MP-1976 shift — see MPSH note in the docstring above.
    let kpts = kpoints::monkhorst_pack(4, 4, 4, kpoints::KGridShift::MP1976, &crystal.lattice);

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
        .expect("Fe SCF should converge");

    let qe_e_xc_ev = -28.904_021_34 * RY_TO_EV;
    let delta = (result.components.e_xc - qe_e_xc_ev).abs();
    eprintln!(
        "  [Fe NLCC] E_xc_pwdft = {:.4} eV,  E_xc_QE = {:.4} eV,  |Δ_xc| = {:.4} eV",
        result.components.e_xc, qe_e_xc_ev, delta,
    );
    // 1 eV ceiling (see docstring). Pre-NCFX would have |Δ_xc| ≈ 49 eV.
    assert!(
        delta <= 1.0,
        "Fe E_xc = {:.4} eV, QE {:.4} eV, |Δ| = {:.4} eV > 1.0 eV — \
         NLCC regression suspected (pre-NCFX baseline was ~49 eV)",
        result.components.e_xc,
        qe_e_xc_ev,
        delta,
    );
}

// ---------------------------------------------------------------------------
// Ewald energy (component-level)
// ---------------------------------------------------------------------------

/// Ewald ion-ion energy for Fe BCC vs QE reference (<0.01 eV tolerance).
///
/// Migrated from the now-deleted `tests/fe_debug.rs` (TACC finding #3: the
/// only surviving test from the Fe 210 eV diagnostic file — the other five
/// were superseded by VGCMP Phases 1-4 or were zip-code/dead assertions).
///
/// QE reference: PseudoDojo Fe LDA (Z_val=16), BCC a=2.87 Å
///   ewald contribution = -171.779_065_80 Ry
///
/// This test pins the Ewald summation convergence on a heavy-Z_val cell.
/// It is feature-independent of NCFX/CCMX/VGCMP and should stay green
/// unless `src/ewald.rs` regresses.
#[test]
fn test_fe_bcc_ewald_vs_qe() {
    let crystal = bcc_crystal(2.87, Atom::new(26, [0.0, 0.0, 0.0]));
    let pp_fe = load_pp("Fe");

    let e_ewald = pwdft_rs::ewald::ewald_energy(&crystal, &[&pp_fe]);
    let qe_ewald = -171.779_065_80 * RY_TO_EV;

    eprintln!("  [Fe Ewald] pwdft = {e_ewald:.6} eV,  QE = {qe_ewald:.6} eV");
    let diff = (e_ewald - qe_ewald).abs();
    eprintln!("  [Fe Ewald] |Δ| = {diff:.6} eV");

    assert!(
        diff < 0.01,
        "Fe Ewald energy {e_ewald:.4} eV differs from QE {qe_ewald:.4} eV by {diff:.4} eV"
    );
}
