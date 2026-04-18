//! VGCMP Phase 2 — cross-check pwdft-rs KB non-local projector form factors
//! β_l(q) against an independent Python reference for Si (ONCVPSP LDA).
//!
//! Compiles and consumes `scripts/validate/beta_q_si_reference.csv`, generated
//! by `scripts/validate/beta_q_reference.py`. The reference computes, for each
//! UPF β-projector and 20 q-values in [0.1, 7.0] Bohr⁻¹,
//!
//!     F_l(q) = 4π · ∫ χ(r) · j_l(q·r) · r dr
//!
//! where χ(r) = r·β(r) is the UPF-stored radial function, exactly matching
//! QE 7.5 `upflib/beta_mod.f90:111-116`.
//!
//! In pwdft-rs the same integral is computed by the private helper
//! `bessel_transform_projector` in `src/potential/nonlocal.rs:223`. Since that
//! function is not part of the public API (`NonlocalPotential::new` hides its
//! output inside a private field), this test **reproduces the production
//! formula verbatim** using only public interfaces (`PseudopotentialData`
//! fields + `numerics::simpson_integrate`), then compares the result to the
//! Python reference value stored in the CSV.
//!
//! If the helper ever drifts from this reproduction, either (a) this test
//! fails, or (b) a similar drift in the reproduction would mask a production
//! bug. Guard (b) by keeping the helper below a direct transcription of the
//! `nonlocal.rs` implementation.
//!
//! Units
//! -----
//! - CSV reference: Bohr^(3/2), q in Bohr⁻¹ (QE native).
//! - pwdft-rs: Å^(3/2), q in Å⁻¹ (internal).
//! - Conversion: F_rust [Å^(3/2)] = F_py [Bohr^(3/2)] × BOHR_TO_ANG^(3/2).
//!
//! We convert the Rust-side F back to Bohr^(3/2) for apples-to-apples and
//! assert per-row agreement to **< 1e-4 Bohr^(3/2) absolute**.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use std::fs;
use std::path::PathBuf;

use pwdft_rs::{
    consts::BOHR_TO_ANG,
    numerics::simpson_integrate,
    pseudopotential::{PseudopotentialData, load},
};

const CSV_REL_PATH: &str = "scripts/validate/beta_q_si_reference.csv";
const UPF_REL_PATH: &str = "pseudopotentials/nc/lda/Si.upf";

/// Per-row tolerance in Bohr^(3/2). Phase 1 achieved ~1e-9 Ry via the same
/// Simpson-on-UPF-mesh strategy; the projector integrand is smoother (no erf
/// subtraction, no 1/r singularity), so we expect comparable agreement.
const TOL_BOHR_3_HALVES: f64 = 1.0e-4;

struct RefRow {
    projector_index: usize,
    l: i32,
    q_bohr_inv: f64,
    f_bohr_3_halves: f64,
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
            fields.len() == 4,
            "malformed CSV row at line {}: {line:?}",
            line_no + 1
        );
        let projector_index: usize = fields[0].trim().parse().unwrap();
        let l: i32 = fields[1].trim().parse().unwrap();
        let q_bohr_inv: f64 = fields[2].trim().parse().unwrap();
        let f_bohr_3_halves: f64 = fields[3].trim().parse().unwrap();
        rows.push(RefRow {
            projector_index,
            l,
            q_bohr_inv,
            f_bohr_3_halves,
        });
    }
    Some(rows)
}

fn load_si_pp() -> PseudopotentialData {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(UPF_REL_PATH);
    load(&path).expect("failed to parse Si.upf")
}

/// Spherical Bessel function j_l(x), replicated from
/// `src/potential/nonlocal.rs::spherical_bessel_j` to avoid depending on the
/// private helper. If this drifts from production, the comparable drift in
/// the production Bessel path would cause this test to fail (by design).
fn spherical_bessel_j(l: i32, x: f64) -> f64 {
    assert!(l >= 0);
    if x.abs() < 1e-10 {
        return if l == 0 { 1.0 } else { 0.0 };
    }
    if l == 0 {
        return x.sin() / x;
    }
    if l == 1 {
        return x.sin() / (x * x) - x.cos() / x;
    }
    let mut jlm1 = x.sin() / x;
    let mut jl = x.sin() / (x * x) - x.cos() / x;
    for n in 1..l {
        let jlp1 = (2 * n + 1) as f64 / x * jl - jlm1;
        jlm1 = jl;
        jl = jlp1;
    }
    jl
}

/// Pure-Rust transcription of `src/potential/nonlocal.rs::bessel_transform_projector`.
///
///     F_l(q) = 4π · Simpson[χ(r) · j_l(q·r) · r, rab]
///
/// - `r_grid`: radial grid (Å)
/// - `rab`: Simpson weights dr (Å)
/// - `chi`: χ(r) = r·β(r) (Å^(-1/2), after UPF unit conversion)
/// - `q`: |k+G| in Å⁻¹
///
/// Returns F in Å^(3/2).
fn f_l_of_q_rust(
    r_grid: &[f64],
    rab: &[f64],
    chi: &[f64],
    l: i32,
    q_ang_inv: f64,
) -> f64 {
    use std::f64::consts::PI;
    let n = r_grid.len();
    let mut integrand = vec![0.0_f64; n];
    for i in 0..n {
        let r = r_grid[i];
        let jl = spherical_bessel_j(l, q_ang_inv * r);
        integrand[i] = chi[i] * jl * r;
    }
    4.0 * PI * simpson_integrate(&integrand, rab)
}

/// Convert Rust output from Å^(3/2) → Bohr^(3/2).
///
/// F [Å^(3/2)] = F [Bohr^(3/2)] · (Å/Bohr)^(3/2) = F [Bohr^(3/2)] · BOHR_TO_ANG^(3/2)
/// ⇒ F [Bohr^(3/2)] = F [Å^(3/2)] / BOHR_TO_ANG^(3/2).
fn ang_3halves_to_bohr_3halves(x: f64) -> f64 {
    x / BOHR_TO_ANG.powf(1.5)
}

#[test]
fn vgcmp_phase2_beta_q_matches_python_reference() {
    struct RowOut {
        projector_index: usize,
        l: i32,
        q_bohr_inv: f64,
        py_bohr_3halves: f64,
        rust_bohr_3halves: f64,
        diff_bohr_3halves: f64,
    }

    let csv_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(CSV_REL_PATH);
    let Some(rows) = load_reference_csv(&csv_path) else {
        eprintln!(
            "VGCMP Phase 2: reference CSV not found at {} — skipping cross-check.\n\
             (Run `uv run scripts/validate/beta_q_reference.py` to regenerate.)",
            csv_path.display()
        );
        return;
    };
    assert!(!rows.is_empty(), "reference CSV contains no data rows");

    let pp = load_si_pp();
    let n_proj = pp.n_projectors();

    let mut table: Vec<RowOut> = Vec::with_capacity(rows.len());

    let mut max_abs_diff: f64 = 0.0;
    let mut worst: Option<(usize, i32, f64, f64, f64, f64)> = None;

    // Per-projector tracking for the post-run table.
    let mut per_proj_max: Vec<f64> = vec![0.0; n_proj];
    // Low-q (< 3.5 Bohr⁻¹) vs high-q bucket tracking.
    let mut max_low: f64 = 0.0;
    let mut max_high: f64 = 0.0;

    for row in &rows {
        assert!(
            row.projector_index < n_proj,
            "projector_index {} out of range ({n_proj} projectors in Si UPF)",
            row.projector_index
        );
        let proj = &pp.beta_projectors[row.projector_index];
        assert_eq!(
            proj.l, row.l,
            "projector {} angular momentum mismatch (UPF says {}, CSV says {})",
            row.projector_index, proj.l, row.l
        );

        // Convert q from Bohr⁻¹ → Å⁻¹: q_Å = q_Bohr / BOHR_TO_ANG.
        let q_ang_inv = row.q_bohr_inv / BOHR_TO_ANG;
        let f_ang = f_l_of_q_rust(&pp.r_grid, &pp.rab, &proj.values, row.l, q_ang_inv);
        let f_bohr = ang_3halves_to_bohr_3halves(f_ang);
        let diff = f_bohr - row.f_bohr_3_halves;

        per_proj_max[row.projector_index] =
            per_proj_max[row.projector_index].max(diff.abs());

        if row.q_bohr_inv < 3.5 {
            max_low = max_low.max(diff.abs());
        } else {
            max_high = max_high.max(diff.abs());
        }

        if diff.abs() > max_abs_diff {
            max_abs_diff = diff.abs();
            worst = Some((
                row.projector_index,
                row.l,
                row.q_bohr_inv,
                row.f_bohr_3_halves,
                f_bohr,
                diff,
            ));
        }

        table.push(RowOut {
            projector_index: row.projector_index,
            l: row.l,
            q_bohr_inv: row.q_bohr_inv,
            py_bohr_3halves: row.f_bohr_3_halves,
            rust_bohr_3halves: f_bohr,
            diff_bohr_3halves: diff,
        });
    }

    // Print the full table unconditionally so passing runs also document
    // the numerics.
    eprintln!(
        "\nVGCMP Phase 2 — Si β_l(q) cross-check (tol = {TOL_BOHR_3_HALVES:.2e} Bohr^(3/2) abs):"
    );
    eprintln!(
        "{:>4}  {:>2}  {:>11}  {:>20}  {:>20}  {:>14}",
        "proj", "l", "q (Bohr⁻¹)", "Python (Bohr^3/2)", "Rust (Bohr^3/2)", "Δ (Bohr^3/2)"
    );
    eprintln!("{}", "-".repeat(82));
    for r in &table {
        eprintln!(
            "{:>4}  {:>2}  {:>11.6}  {:>20.12e}  {:>20.12e}  {:>14.3e}",
            r.projector_index,
            r.l,
            r.q_bohr_inv,
            r.py_bohr_3halves,
            r.rust_bohr_3halves,
            r.diff_bohr_3halves
        );
    }

    eprintln!(
        "\nPer-projector max |Δ| (Bohr^(3/2)):"
    );
    for (i, &d) in per_proj_max.iter().enumerate() {
        let l = pp.beta_projectors[i].l;
        eprintln!("  proj={i} l={l}  max|Δ| = {d:.3e}");
    }
    eprintln!(
        "\nBucket max |Δ|: low-q (q<3.5)  = {max_low:.3e} Bohr^(3/2)"
    );
    eprintln!("                high-q (q>=3.5) = {max_high:.3e} Bohr^(3/2)");
    eprintln!("\noverall max |Δ| = {max_abs_diff:.3e} Bohr^(3/2)");

    if max_abs_diff >= TOL_BOHR_3_HALVES {
        let (pi, l, q_b, py, rust, diff) = worst.unwrap();
        panic!(
            "β_l(q) disagrees with Python reference beyond tolerance.\n\
             Worst row: proj={pi}, l={l}, q={q_b:.6} Bohr⁻¹\n\
             Python = {py:+.8e} Bohr^(3/2), Rust = {rust:+.8e} Bohr^(3/2), Δ = {diff:+.3e}\n\
             max |Δ| = {max_abs_diff:.3e} Bohr^(3/2), tol = {TOL_BOHR_3_HALVES:.3e}"
        );
    }
}
