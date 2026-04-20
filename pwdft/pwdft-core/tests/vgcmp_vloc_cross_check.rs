//! VGCMP Phase 1 — cross-check pwdft-core V_local(G) against an independent
//! Python reference for Si (ONCVPSP LDA).
//!
//! This test consumes `pwdft/pwdft-validation/pwdft_validation/scripts/vloc_g_si_reference.csv`, generated
//! by `pwdft/pwdft-validation/pwdft_validation/scripts/vloc_g_reference.py`, which computes V_local(G) for
//! the first 20 distinct |G| shells of Si FCC (a = 5.431 Å) directly from
//! `pseudopotentials/nc/lda/Si.upf` using QE's erf-subtracted formula.
//!
//! The reference is in QE native units (Ry, Bohr). pwdft-core's
//! `PseudopotentialData::v_local_of_g` returns values in internal units
//! (eV, Å). We convert the reference to eV for comparison and require
//! agreement to < 1e-4 Ry (≈ 1.36e-3 eV) absolute per shell.
//!
//! If the CSV is missing (e.g. running in an environment without Python),
//! the test prints a skip message and passes. The CSV is committed as a
//! generated artifact so CI runs do not depend on Python at test time.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use std::fs;
use std::path::PathBuf;

use nalgebra::Vector3;
use pwdft_core::{
    consts::{BOHR_TO_ANG, RY_TO_EV},
    crystal::Lattice,
    pseudopotential::{PseudopotentialData, load},
};

const CSV_REL_PATH: &str = "data/csv/vloc_g_si_reference.csv";
const UPF_REL_PATH: &str = "pseudopotentials/nc/lda/Si.upf";

/// Tolerance: < 1e-4 Ry absolute (matches QE UPF interpolation precision;
/// see VGCMP proposal Phase 1 pass criterion).
const TOL_EV: f64 = 1.0e-4 * RY_TO_EV; // ≈ 1.36e-3 eV

/// One row of the reference CSV.
struct RefRow {
    shell_index: usize,
    g_bohr_inv: f64,
    g2_int: i32,
    v_loc_ry: f64,
}

fn load_reference_csv(path: &PathBuf) -> Option<Vec<RefRow>> {
    let content = fs::read_to_string(path).ok()?;
    let mut rows = Vec::new();
    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line_no == 0 {
            // Skip header + blank lines.
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() != 4 {
            panic!("malformed CSV row at line {}: {line:?}", line_no + 1);
        }
        let shell_index: usize = fields[0].trim().parse().unwrap();
        let g_bohr_inv: f64 = fields[1].trim().parse().unwrap();
        let g2_int: i32 = fields[2].trim().parse().unwrap();
        let v_loc_ry: f64 = fields[3].trim().parse().unwrap();
        rows.push(RefRow { shell_index, g_bohr_inv, g2_int, v_loc_ry });
    }
    Some(rows)
}

fn load_si_pp() -> PseudopotentialData {
    let path = PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join(UPF_REL_PATH);
    load(&path).expect("failed to parse Si.upf")
}

/// Primitive Si FCC cell at a = 5.431 Å.
fn si_cell_volume_ang3() -> f64 {
    let a: f64 = 5.431;
    let lat = Lattice::new(
        (a / 2.0) * Vector3::new(0.0, 1.0, 1.0),
        (a / 2.0) * Vector3::new(1.0, 0.0, 1.0),
        (a / 2.0) * Vector3::new(1.0, 1.0, 0.0),
    );
    lat.volume()
}

#[test]
fn vgcmp_phase1_v_local_g_matches_python_reference() {
    struct RowOut {
        shell: usize,
        g2_int: i32,
        g_bohr_inv: f64,
        py_ry: f64,
        rust_ry: f64,
        diff_ry: f64,
    }

    let csv_path = PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join(CSV_REL_PATH);
    let Some(rows) = load_reference_csv(&csv_path) else {
        eprintln!(
            "VGCMP: reference CSV not found at {} — skipping cross-check.\n\
             (Run `uv run pwdft/pwdft-validation/pwdft_validation/scripts/vloc_g_reference.py` to regenerate.)",
            csv_path.display()
        );
        return;
    };
    assert!(!rows.is_empty(), "reference CSV contains no data rows");
    assert_eq!(rows.len(), 20, "expected 20 shells in reference CSV");

    let pp = load_si_pp();
    let omega_ang3 = si_cell_volume_ang3();

    // Per-shell comparison. Build a table so test failure prints the full
    // structure (small-G vs large-G, growing vs isolated discrepancies).
    let mut max_abs_diff_ev: f64 = 0.0;
    let mut max_abs_diff_ry: f64 = 0.0;
    let mut worst: Option<(usize, i32, f64, f64, f64, f64)> = None;

    let mut table: Vec<RowOut> = Vec::with_capacity(rows.len());

    for row in &rows {
        // Convert |G| from Bohr⁻¹ → Å⁻¹: divide by BOHR_TO_ANG (since
        // 1 Bohr⁻¹ = (1/BOHR_TO_ANG) Å⁻¹ when BOHR_TO_ANG < 1 is Å/Bohr).
        //   g [Å⁻¹] = g [Bohr⁻¹] / BOHR_TO_ANG
        let g_ang_inv = row.g_bohr_inv / BOHR_TO_ANG;
        let v_ev = pp.v_local_of_g(g_ang_inv, omega_ang3);
        let v_rust_ry = v_ev / RY_TO_EV;
        let diff_ry = v_rust_ry - row.v_loc_ry;
        let diff_ev = v_ev - row.v_loc_ry * RY_TO_EV;

        if diff_ry.abs() > max_abs_diff_ry {
            max_abs_diff_ry = diff_ry.abs();
            max_abs_diff_ev = diff_ev.abs();
            worst = Some((
                row.shell_index,
                row.g2_int,
                row.g_bohr_inv,
                row.v_loc_ry,
                v_rust_ry,
                diff_ry,
            ));
        }

        table.push(RowOut {
            shell: row.shell_index,
            g2_int: row.g2_int,
            g_bohr_inv: row.g_bohr_inv,
            py_ry: row.v_loc_ry,
            rust_ry: v_rust_ry,
            diff_ry,
        });
    }

    // Always print the table so successful runs still document the numbers.
    eprintln!(
        "\nVGCMP Phase 1 — Si V_local(G) cross-check (tol = {:.2e} Ry / {:.2e} eV abs):",
        1.0e-4, TOL_EV
    );
    eprintln!(
        "{:>5}  {:>8}  {:>14}  {:>18}  {:>18}  {:>14}",
        "shell", "|G|²", "|G| (Bohr⁻¹)", "Python (Ry)", "Rust (Ry)", "Δ (Ry)"
    );
    eprintln!("{}", "-".repeat(86));
    for r in &table {
        eprintln!(
            "{:>5}  {:>8}  {:>14.6}  {:>18.10e}  {:>18.10e}  {:>14.3e}",
            r.shell, r.g2_int, r.g_bohr_inv, r.py_ry, r.rust_ry, r.diff_ry
        );
    }
    eprintln!(
        "\nmax |Δ| = {max_abs_diff_ry:.3e} Ry  ({max_abs_diff_ev:.3e} eV)"
    );

    if max_abs_diff_ev >= TOL_EV {
        let (shell, g2, g_bohr, py, rust, diff) = worst.unwrap();
        panic!(
            "V_local(G) disagrees with Python reference beyond tolerance.\n\
             Worst shell: i={shell}, |G|²={g2}, |G|={g_bohr:.6} Bohr⁻¹\n\
             Python = {py:+.8e} Ry, Rust = {rust:+.8e} Ry, Δ = {diff:+.3e} Ry\n\
             max |Δ| = {max_abs_diff_ry:.3e} Ry ({max_abs_diff_ev:.3e} eV), tol = {TOL_EV:.3e} eV"
        );
    }
}
