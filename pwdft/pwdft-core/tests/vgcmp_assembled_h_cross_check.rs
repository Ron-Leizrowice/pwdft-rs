//! VGCMP Phase 4 — cross-check pwdft-core' assembled Hamiltonian diagonal
//! against an independent Python reference for Si (ONCVPSP LDA) at k=Γ.
//!
//! Context
//! -------
//! Phases 1–3 have closed the per-term pseudopotential form-factor checks:
//! V_local(G), β_l(q), and D_ij all match QE to machine precision. The
//! residual 13.4 eV Si total-energy gap must therefore come from either the
//! *assembly* of those form factors into H_{G,G'} (structure factor,
//! (2l+1)/(4π) angular factor, 1/Ω prefactor, D_ij summation pattern) or
//! something outside the PP pipeline entirely (kinetic convention, V_eff,
//! Ewald, symmetry, SCF mixing).
//!
//! Phase 4 tests the assembly directly. We build a Si FCC Γ-point
//! Hamiltonian with V_eff set to **zero**, leaving only kinetic + V_NL on
//! the diagonal, and compare to
//! `data/csv/vgcmp_phase4_reference.csv` which was produced by an
//! independent Python implementation that re-uses the Phase 1/2/3 validated
//! form factors and re-derives the assembly formula from scratch.
//!
//! Setting V_eff = 0 removes the Hartree/XC/V_local(≠0) contribution to
//! off-diagonal elements and the *constant* V_xc(G=0) shift to the diagonal.
//! This isolates the G-dependent part of the diagonal — kinetic + V_NL —
//! which is where any assembly-stage bug would manifest. V_xc(G=0) alone
//! cannot break Γ-point degeneracies or produce per-shell discrepancies.
//!
//! Equivalence: `hamiltonian::build_kinetic` is identical to
//! `scf::potentials::build_hamiltonian_with_v_eff` evaluated at V_eff = 0
//! on the FFT grid (both compute `HBAR2_OVER_2M * |k+G|²` on the diagonal
//! and zero elsewhere). Since `build_hamiltonian_with_v_eff` is crate-private
//! we use the equivalent public entry point in this test.
//!
//! Units
//! -----
//! CSV reference values are in Ry (QE native). pwdft-core Hamiltonian entries
//! are in eV (internal). We convert the reference Ry → eV and require
//! per-shell agreement on both the kinetic and V_NL terms to
//! **< 1e-4 Ry (≈ 1.36e-3 eV) absolute**, matching the Phase 1/2 tolerance.
//!
//! If the CSV is missing, the test prints a skip message and passes; the
//! CSV is committed as a generated artifact so CI does not depend on Python.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use std::{fs, path::PathBuf};

use nalgebra::Vector3;
use pwdft_core::{
    basis::BasisSet,
    consts::{BOHR_TO_ANG, HBAR2_OVER_2M, RY_TO_EV},
    crystal::{Atom, Crystal, Lattice},
    hamiltonian::build_kinetic,
    potential::nonlocal::NonlocalPotential,
    pseudopotential::UpfPseudoPotential,
};

const CSV_REL_PATH: &str = "data/csv/vgcmp_phase4_reference.csv";
const UPF_REL_PATH: &str = "pseudopotentials/nc/lda/Si.upf";

/// Tolerance: < 1e-4 Ry absolute (matches Phases 1 and 2 pass criterion).
const TOL_RY: f64 = 1.0e-4;
const TOL_EV: f64 = TOL_RY * RY_TO_EV;

/// One row of the Phase 4 reference CSV.
struct RefRow {
    shell_index: usize,
    miller: [i32; 3],
    g_bohr_inv: f64,
    g2_int_tpba2: i32,
    kinetic_ry: f64,
    v_nl_diag_ry: f64,
    h_diag_ry: f64,
}

fn load_reference_csv(path: &PathBuf) -> Option<Vec<RefRow>> {
    let content = fs::read_to_string(path).ok()?;
    let mut rows = Vec::new();
    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line_no == 0 {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        assert!(
            fields.len() == 10,
            "malformed CSV row at line {}: {line:?} (got {} fields, expected 10)",
            line_no + 1,
            fields.len()
        );
        let shell_index: usize = fields[0].trim().parse().unwrap();
        let n1: i32 = fields[1].trim().parse().unwrap();
        let n2: i32 = fields[2].trim().parse().unwrap();
        let n3: i32 = fields[3].trim().parse().unwrap();
        let g_bohr_inv: f64 = fields[4].trim().parse().unwrap();
        let g2_int_tpba2: i32 = fields[5].trim().parse().unwrap();
        // fields[6] is q_bohr_inv, equal to g_bohr_inv at k=Γ — not used here.
        let kinetic_ry: f64 = fields[7].trim().parse().unwrap();
        let v_nl_diag_ry: f64 = fields[8].trim().parse().unwrap();
        let h_diag_ry: f64 = fields[9].trim().parse().unwrap();
        rows.push(RefRow {
            shell_index,
            miller: [n1, n2, n3],
            g_bohr_inv,
            g2_int_tpba2,
            kinetic_ry,
            v_nl_diag_ry,
            h_diag_ry,
        });
    }
    Some(rows)
}

fn load_si_pp() -> UpfPseudoPotential {
    let path = PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join(UPF_REL_PATH);
    UpfPseudoPotential::load(&path).expect("failed to parse Si.upf")
}

/// Primitive Si FCC cell at a = 5.431 Å, matching Phase 1/2/3 reference.
/// Two atoms at fractional (0, 0, 0) and (1/4, 1/4, 1/4).
fn build_si_fcc() -> Crystal {
    let a: f64 = 5.431;
    let lattice = Lattice::new(
        (a / 2.0) * Vector3::new(0.0, 1.0, 1.0),
        (a / 2.0) * Vector3::new(1.0, 0.0, 1.0),
        (a / 2.0) * Vector3::new(1.0, 1.0, 0.0),
    );
    let atoms = vec![Atom::new(14, [0.0, 0.0, 0.0]), Atom::new(14, [0.25, 0.25, 0.25])];
    Crystal { atoms, lattice }
}

/// Find the basis index corresponding to a given Miller triple.
/// Panics if the triple is not in the basis (which would mean ecut is too
/// small).
fn basis_index_of(basis: &BasisSet, miller: [i32; 3]) -> usize {
    basis.index_of(miller[0], miller[1], miller[2]).unwrap_or_else(|| {
        panic!(
            "Miller index ({}, {}, {}) not found in basis of size {}. \
                 Increase ecut in the test.",
            miller[0],
            miller[1],
            miller[2],
            basis.len()
        )
    })
}

#[test]
fn vgcmp_phase4_hbar2_over_2m_consistency() {
    // Pre-flight check: the SI-derived HBAR2_OVER_2M in pwdft-core (eV·Å²) must
    // match the QE-native Hartree-atomic-unit value derived from RY_TO_EV and
    // BOHR_TO_ANG. If these disagree, the kinetic term is wrong independent
    // of anything else, and every diagonal H entry is shifted in a
    // G-dependent way.
    //
    //   ħ²/(2m_e) = 0.5 Ha · a₀² = RY_TO_EV · BOHR_TO_ANG²  [eV·Å²]
    //
    // (because Ha = 2 Ry, so 0.5 Ha = 1 Ry = RY_TO_EV eV;
    // a₀² = BOHR_TO_ANG² Å².)
    let qe_derived_ev_a2 = RY_TO_EV * BOHR_TO_ANG * BOHR_TO_ANG;
    let rel = (HBAR2_OVER_2M - qe_derived_ev_a2).abs() / qe_derived_ev_a2;
    eprintln!("\nVGCMP Phase 4 — ħ²/(2m) consistency check:");
    eprintln!("  pwdft-core HBAR2_OVER_2M (SI-derived)     = {HBAR2_OVER_2M:.10} eV·Å²");
    eprintln!("  RY_TO_EV * BOHR_TO_ANG² (QE convention) = {qe_derived_ev_a2:.10} eV·Å²");
    eprintln!("  relative difference: {rel:.3e}");
    assert!(
        rel < 1.0e-6,
        "HBAR2_OVER_2M disagrees with QE convention by {rel:.3e} (rel); > 1e-6 is suspicious."
    );
}

#[test]
fn vgcmp_phase4_assembled_h_matches_python_reference() {
    struct RowOut {
        shell: usize,
        miller: [i32; 3],
        g_bohr_inv: f64,
        g2_int: i32,
        kin_py_ry: f64,
        kin_rust_ry: f64,
        kin_diff_ry: f64,
        vnl_py_ry: f64,
        vnl_rust_ry: f64,
        vnl_diff_ry: f64,
        h_py_ry: f64,
        h_rust_ry: f64,
        h_diff_ry: f64,
    }

    let csv_path = PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join(CSV_REL_PATH);
    let Some(rows) = load_reference_csv(&csv_path) else {
        eprintln!(
            "VGCMP Phase 4: reference CSV not found at {} — skipping cross-check.\n\
             (Run `uv run pwdft-validate reference hamiltonian` to regenerate.)",
            csv_path.display()
        );
        return;
    };
    assert!(!rows.is_empty(), "reference CSV contains no data rows");
    assert_eq!(rows.len(), 5, "expected 5 shells in reference CSV");

    // --- Set up Si FCC Γ-point calculation ----------------------------------
    let crystal = build_si_fcc();
    let pp = load_si_pp();
    let pp_refs: Vec<&UpfPseudoPotential> = vec![&pp];

    // ecut = 200 eV: comfortably covers all 5 reference shells.
    // Shell 4 has |G|² = 11·(2π/a)² ≈ 14.7 Å⁻² → kinetic T ≈ 56 eV.
    let ecut_ev: f64 = 200.0;
    let basis = BasisSet::new(&crystal.lattice, ecut_ev);
    eprintln!(
        "\nVGCMP Phase 4 — Si Γ-point assembled H cross-check\n\
         ecut = {ecut_ev} eV, basis size = {}",
        basis.len()
    );

    let k_gamma: Vector3<f64> = Vector3::zeros();

    // Kinetic-only diagonal: `build_kinetic` is equivalent to
    // `build_hamiltonian_with_v_eff` with V_eff_fft = 0.
    let h_kinetic = build_kinetic(&basis, &k_gamma);

    // Full Hamiltonian = kinetic + V_NL (V_eff still zero).
    let mut h_full = build_kinetic(&basis, &k_gamma);
    let vnl =
        NonlocalPotential::new(&crystal, &basis, &k_gamma, &pp_refs).expect("NonlocalPotential construction failed");
    vnl.add_to_hamiltonian(&mut h_full, &crystal, &basis, &k_gamma);

    // Build per-row output.
    let mut table: Vec<RowOut> = Vec::with_capacity(rows.len());
    let mut max_abs_h_diff_ry: f64 = 0.0;
    let mut max_abs_kin_diff_ry: f64 = 0.0;
    let mut max_abs_vnl_diff_ry: f64 = 0.0;

    for row in &rows {
        let idx = basis_index_of(&basis, row.miller);

        // Convert Rust H (eV) → Ry for comparison with Python reference.
        let h_rust_ry = h_full[(idx, idx)].re / RY_TO_EV;
        let kin_rust_ry = h_kinetic[(idx, idx)].re / RY_TO_EV;
        let vnl_rust_ry = h_rust_ry - kin_rust_ry;

        let kin_diff_ry = kin_rust_ry - row.kinetic_ry;
        let vnl_diff_ry = vnl_rust_ry - row.v_nl_diag_ry;
        let h_diff_ry = h_rust_ry - row.h_diag_ry;

        max_abs_h_diff_ry = max_abs_h_diff_ry.max(h_diff_ry.abs());
        max_abs_kin_diff_ry = max_abs_kin_diff_ry.max(kin_diff_ry.abs());
        max_abs_vnl_diff_ry = max_abs_vnl_diff_ry.max(vnl_diff_ry.abs());

        // Assert non-zero imaginary part is negligible (H is Hermitian on-diagonal
        // ⇒ diagonal must be real).
        let im = h_full[(idx, idx)].im.abs();
        assert!(
            im < 1e-10,
            "H[{idx},{idx}] has non-zero imaginary part {im:.3e} eV at shell {} (Miller {:?})",
            row.shell_index,
            row.miller
        );

        table.push(RowOut {
            shell: row.shell_index,
            miller: row.miller,
            g_bohr_inv: row.g_bohr_inv,
            g2_int: row.g2_int_tpba2,
            kin_py_ry: row.kinetic_ry,
            kin_rust_ry,
            kin_diff_ry,
            vnl_py_ry: row.v_nl_diag_ry,
            vnl_rust_ry,
            vnl_diff_ry,
            h_py_ry: row.h_diag_ry,
            h_rust_ry,
            h_diff_ry,
        });
    }

    // Print the full table unconditionally so passing runs also document the
    // numerics. The per-term split (kinetic vs V_NL) pinpoints which sub-term
    // carries any discrepancy.
    eprintln!("\n(tol = {TOL_RY:.2e} Ry abs on each of kinetic, V_NL, H_diag)");
    eprintln!(
        "\n{:>5}  {:>14}  {:>4}  {:>10}  {:>12}  {:>12}  {:>11}  {:>12}  {:>12}  {:>11}  {:>12}  {:>12}  {:>11}",
        "shell",
        "n1,n2,n3",
        "|G|²",
        "|G|(Bohr⁻¹)",
        "T_py(Ry)",
        "T_rust(Ry)",
        "ΔT (Ry)",
        "V_py(Ry)",
        "V_rust(Ry)",
        "ΔV (Ry)",
        "H_py(Ry)",
        "H_rust(Ry)",
        "ΔH (Ry)"
    );
    eprintln!("{}", "-".repeat(168));
    for r in &table {
        eprintln!(
            "{:>5}  ({:>2},{:>2},{:>2})  {:>4}  {:>10.6}  \
             {:>12.5e}  {:>12.5e}  {:>11.3e}  \
             {:>12.5e}  {:>12.5e}  {:>11.3e}  \
             {:>12.5e}  {:>12.5e}  {:>11.3e}",
            r.shell,
            r.miller[0],
            r.miller[1],
            r.miller[2],
            r.g2_int,
            r.g_bohr_inv,
            r.kin_py_ry,
            r.kin_rust_ry,
            r.kin_diff_ry,
            r.vnl_py_ry,
            r.vnl_rust_ry,
            r.vnl_diff_ry,
            r.h_py_ry,
            r.h_rust_ry,
            r.h_diff_ry,
        );
    }
    eprintln!("\nSummary:");
    eprintln!(
        "  max |Δ kinetic| = {max_abs_kin_diff_ry:.3e} Ry ({:.3e} eV)",
        max_abs_kin_diff_ry * RY_TO_EV
    );
    eprintln!(
        "  max |Δ V_NL|    = {max_abs_vnl_diff_ry:.3e} Ry ({:.3e} eV)",
        max_abs_vnl_diff_ry * RY_TO_EV
    );
    eprintln!(
        "  max |Δ H_diag|  = {max_abs_h_diff_ry:.3e} Ry ({:.3e} eV)",
        max_abs_h_diff_ry * RY_TO_EV
    );

    if max_abs_kin_diff_ry >= TOL_RY {
        panic!(
            "Kinetic diagonal disagrees with Python reference beyond tolerance.\n\
             max |Δ kinetic| = {max_abs_kin_diff_ry:.3e} Ry ({:.3e} eV), tol = {TOL_RY:.3e} Ry ({TOL_EV:.3e} eV)",
            max_abs_kin_diff_ry * RY_TO_EV
        );
    }
    if max_abs_vnl_diff_ry >= TOL_RY {
        panic!(
            "V_NL diagonal disagrees with Python reference beyond tolerance.\n\
             max |Δ V_NL| = {max_abs_vnl_diff_ry:.3e} Ry ({:.3e} eV), tol = {TOL_RY:.3e} Ry ({TOL_EV:.3e} eV)",
            max_abs_vnl_diff_ry * RY_TO_EV
        );
    }
    if max_abs_h_diff_ry >= TOL_RY {
        panic!(
            "H_diag disagrees with Python reference beyond tolerance.\n\
             max |Δ H_diag| = {max_abs_h_diff_ry:.3e} Ry ({:.3e} eV), tol = {TOL_RY:.3e} Ry ({TOL_EV:.3e} eV)",
            max_abs_h_diff_ry * RY_TO_EV
        );
    }
}
