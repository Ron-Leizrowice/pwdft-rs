//! VGCH Phase 1c Hypothesis 2 — Superposition-of-Atomic-Densities (SAD)
//! initial-density cross-check against an independent Python reference
//! (`pwdft-validate reference sad`) for every VGCH-class system.
//!
//! The hypothesis under test: does pwdft-core' SAD (its PP_RHOATOM Bessel
//! transform + structure-factor sum + IFFT + clamp/renormalize pipeline)
//! disagree with the QE-convention reference on *any* of the seven
//! VGCH-scope cells (C, Al, Fe, Cu, GaAs, NaCl, MgO)?
//!
//! **If C diamond matches bit-perfect and still carries its 1.45 eV
//! E_total residual, Hypothesis 2 is cleared** — the residual lives
//! outside SAD (mixer basin / energy assembly). Si and Al are informative
//! controls; Si matches QE to <40 meV and Al matches to ~100 meV, so if
//! SAD is the bug it ought to show up there too, modestly.
//!
//! The Python reference reproduces QE's `atomic_rho.f90` recipe:
//!   ρ(G) = Σ_species strf(G,nt) · (1/Ω) · Simpson( ρ_at(r) · j₀(Gr) ; rab )
//!   ρ(r) = Σ_G ρ(G) exp(+iG·r)                (pwdft-core unnormalized IFFT)
//! plus QE's G=0 renormalization (`charge = Ω · ρ(G=0); ρ ← ρ · N_el /
//! charge`). The Rust side calls `build_sad_density_for_diagnostic`, which
//! returns the SCF driver's exact starting ρ(r) — post-clamp and post-renorm.
//! Both sides shell-average around each atom on identical grid dims.
//!
//! Tolerance is per-bin in e/Å³; see `TOL_RHO_E_PER_ANG3` and per-system
//! overrides below. The dominant expected discrepancy at a light-atom
//! regime is the negative-clamp step in pwdft-core (QE skips it; see
//! `qe-7.5/PW/src/atomic_rho.f90:186-188`) that trims a few percent of
//! mass off the core and renormalizes the rest. If the diagnostic
//! surfaces clamp-sized Δρ on a heavy atom, the clamp is the bug.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use std::{
    collections::{BTreeMap, HashMap},
    fs,
    path::PathBuf,
};

use nalgebra::Vector3;
use pwdft_core::{
    crystal::{Atom, Crystal, Lattice},
    pseudopotential::UpfPseudoPotential,
    scf::initial_density::{InitialDensityConfig, build_sad_density_for_diagnostic_verbose},
};

const CSV_REL_PATH: &str = "data/csv/vgch_sad_heavy.csv";
const SAMPLES_CSV_REL_PATH: &str = "data/csv/vgch_sad_heavy_samples.csv";
const BOHR_TO_ANG: f64 = 0.529_177_210_903;
const RY_TO_EV: f64 = 13.605_693_122_994;

/// Default per-bin absolute tolerance in e/Å³ for the *final* (post-
/// clamp + post-renorm) shell-averaged density compared against the
/// Python reference's post-renorm density. The pre-clamp tolerance
/// is pinned separately (see `TOL_PRE_CLAMP_RHO`) — it is where the
/// physics lives, since clamp/renorm is a production-only deviation
/// from QE's `atomic_rho.f90` (line 186-188 of QE explicitly skips
/// the clamp).
///
/// 5e-2 comfortably covers:
///   - Bit-perfect (≤ 1e-9) agreement on systems where no density sample goes negative (Al, Fe, Cu,
///     NaCl, MgO, C — post-wrap-fix).
///   - The tiny (≤ 1.5e-2 e/Å³ local) deviation on GaAs near the As core, where pwdft-core clamps
///     O(1e-5 e) of negative Gibbs ringing that QE keeps. That 1.5e-2 is diagnostic, not a physics
///     bug; the clamp is confined to ~1 bin near the core atom.
const TOL_RHO_E_PER_ANG3: f64 = 5.0e-2;

/// Tolerance for the *pre-clamp* ρ(r) shell averages: this is the
/// tight physics-correctness check that the Bessel transform +
/// structure factor + IFFT pipeline is bit-perfect relative to QE's
/// `atomic_rho.f90` convention. Observed residuals on all seven VGCH
/// systems are ≤ 8e-6 e/Å³; we budget 1e-4 for long-term robustness
/// against trivially-equivalent refactors.
const TOL_PRE_CLAMP_RHO: f64 = 1.0e-4;

/// Tolerance for the raw-sample pointwise ρ(r) diff. Same physics as
/// the pre-clamp tolerance above — asserts that the full ρ(r) array
/// (post-clamp + post-renorm in both Rust and Python) agrees to
/// machine precision except where the clamp intervenes. GaAs shows
/// the largest residual here: 1.2e-5 at one of 200 sampled points,
/// sourced to the same As-core Gibbs-ringing clamp bin. Budget 1e-4.
const TOL_RAW_SAMPLE_RHO: f64 = 1.0e-4;

/// Per-bin tolerance override table. Key = system name; value =
/// tolerance in e/Å³. Any system not present here uses
/// `TOL_RHO_E_PER_ANG3`. The reason a system lands here is *always*
/// documented — a bin-by-bin mismatch is diagnostic, not noise.
///
/// No overrides yet — 5e-2 covers all seven systems including GaAs.
/// The clamp effect is confined to ≤1 bin per atom; widening further
/// would weaken the regression guard on the rest of the density.
fn tolerance_override(_system: &str) -> f64 {
    TOL_RHO_E_PER_ANG3
}

/// One shell-averaged row from the Python reference CSV.
struct RefRow {
    system: String,
    atom_label: String,
    r_bin_center: f64,
    rho_avg_ref: f64,
    #[allow(dead_code, reason = "exposed for future debug prints")]
    n_bin_points_ref: u64,
}

struct SampleRow {
    system: String,
    i1: usize,
    i2: usize,
    i3: usize,
    rho_ref: f64,
}

fn load_samples_csv(path: &PathBuf) -> Option<Vec<SampleRow>> {
    let content = fs::read_to_string(path).ok()?;
    let mut rows = Vec::new();
    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line_no == 0 {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        assert!(fields.len() == 5, "malformed samples row: {line:?}");
        rows.push(SampleRow {
            system: fields[0].trim().to_string(),
            i1: fields[1].trim().parse().unwrap(),
            i2: fields[2].trim().parse().unwrap(),
            i3: fields[3].trim().parse().unwrap(),
            rho_ref: fields[4].trim().parse().unwrap(),
        });
    }
    Some(rows)
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
            fields.len() == 5,
            "malformed CSV row at line {}: {line:?} ({} fields)",
            line_no + 1,
            fields.len()
        );
        rows.push(RefRow {
            system: fields[0].trim().to_string(),
            atom_label: fields[1].trim().to_string(),
            r_bin_center: fields[2].trim().parse().unwrap(),
            rho_avg_ref: fields[3].trim().parse().unwrap(),
            n_bin_points_ref: fields[4].trim().parse().unwrap(),
        });
    }
    Some(rows)
}

// ---------------------------------------------------------------------------
// System registry — MUST match the Python script's SYSTEMS list exactly.
// ---------------------------------------------------------------------------

struct SystemSpec {
    name: &'static str,
    lattice_type: LatticeType,
    celldm1_bohr: f64,
    atoms: Vec<(&'static str, [f64; 3], u32)>, // (element, frac, z)
    grid_dims: [usize; 3],
    ecutwfc_ry: f64,
    ecutrho_ratio: u32,
}

#[derive(Clone, Copy)]
enum LatticeType {
    Fcc,
    Bcc,
    /// Rocksalt uses FCC primitive lattice (Bravais), so same vectors as Fcc.
    /// Distinguished only for documentation; geometry identical.
    RocksaltFcc,
}

fn systems() -> Vec<SystemSpec> {
    vec![
        SystemSpec {
            name: "c_diamond",
            lattice_type: LatticeType::Fcc,
            celldm1_bohr: 6.7409,
            atoms: vec![("C", [0.00, 0.00, 0.00], 6), ("C", [0.25, 0.25, 0.25], 6)],
            // 32³ aligns the two C atoms with integer grid points
            // (0,0,0) and (8,8,8); see the Python twin file for the
            // discretization rationale.
            grid_dims: [32, 32, 32],
            ecutwfc_ry: 30.0,
            ecutrho_ratio: 4,
        },
        SystemSpec {
            name: "al_fcc",
            lattice_type: LatticeType::Fcc,
            celldm1_bohr: 7.6527,
            atoms: vec![("Al", [0.00, 0.00, 0.00], 13)],
            grid_dims: [24, 24, 24],
            ecutwfc_ry: 24.0,
            ecutrho_ratio: 4,
        },
        SystemSpec {
            name: "fe_bcc",
            lattice_type: LatticeType::Bcc,
            celldm1_bohr: 5.4235,
            atoms: vec![("Fe", [0.00, 0.00, 0.00], 26)],
            grid_dims: [24, 24, 24],
            ecutwfc_ry: 15.0,
            ecutrho_ratio: 4,
        },
        SystemSpec {
            name: "cu_fcc",
            lattice_type: LatticeType::Fcc,
            celldm1_bohr: 6.8219,
            atoms: vec![("Cu", [0.00, 0.00, 0.00], 29)],
            grid_dims: [24, 24, 24],
            ecutwfc_ry: 25.0,
            ecutrho_ratio: 4,
        },
        SystemSpec {
            name: "gaas",
            lattice_type: LatticeType::Fcc,
            celldm1_bohr: 10.6829,
            atoms: vec![("Ga", [0.00, 0.00, 0.00], 31), ("As", [0.25, 0.25, 0.25], 33)],
            grid_dims: [32, 32, 32],
            ecutwfc_ry: 20.0,
            ecutrho_ratio: 4,
        },
        SystemSpec {
            name: "nacl",
            lattice_type: LatticeType::RocksaltFcc,
            celldm1_bohr: 10.6078,
            atoms: vec![("Na", [0.00, 0.00, 0.00], 11), ("Cl", [0.50, 0.50, 0.50], 17)],
            grid_dims: [32, 32, 32],
            ecutwfc_ry: 25.0,
            ecutrho_ratio: 4,
        },
        SystemSpec {
            name: "mgo",
            lattice_type: LatticeType::RocksaltFcc,
            celldm1_bohr: 7.9586,
            atoms: vec![("Mg", [0.00, 0.00, 0.00], 12), ("O", [0.50, 0.50, 0.50], 8)],
            grid_dims: [24, 24, 24],
            ecutwfc_ry: 30.0,
            ecutrho_ratio: 4,
        },
    ]
}

fn build_crystal(sys: &SystemSpec) -> Crystal {
    let a_ang = sys.celldm1_bohr * BOHR_TO_ANG;
    let lattice = match sys.lattice_type {
        LatticeType::Fcc | LatticeType::RocksaltFcc => Lattice::new(
            (a_ang / 2.0) * Vector3::new(0.0, 1.0, 1.0),
            (a_ang / 2.0) * Vector3::new(1.0, 0.0, 1.0),
            (a_ang / 2.0) * Vector3::new(1.0, 1.0, 0.0),
        ),
        LatticeType::Bcc => Lattice::new(
            (a_ang / 2.0) * Vector3::new(-1.0, 1.0, 1.0),
            (a_ang / 2.0) * Vector3::new(1.0, -1.0, 1.0),
            (a_ang / 2.0) * Vector3::new(1.0, 1.0, -1.0),
        ),
    };
    let atoms = sys.atoms.iter().map(|(_, frac, z)| Atom::new(*z, *frac)).collect();
    Crystal { lattice, atoms }
}

fn load_pp(element: &str) -> UpfPseudoPotential {
    let path = PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join(format!("pseudopotentials/nc/lda/{element}.upf"));
    UpfPseudoPotential::load(&path).unwrap_or_else(|e| panic!("failed to load {}: {e}", path.display()))
}

/// Shell-average ρ(r) around one atom. MUST mirror Python's
/// `shell_average`: fractional-coord min-image wrap, then bin by
/// cartesian distance.
fn shell_average(
    rho_r: &[f64],
    dims: [usize; 3],
    lattice: &Lattice,
    tau_cart: Vector3<f64>,
    r_edges: &[f64],
) -> (Vec<f64>, Vec<u64>) {
    let [nx, ny, nz] = dims;
    let n_bins = r_edges.len() - 1;

    // Build A_rows matrix and its inverse (reciprocal / 2π).
    // cart = frac · A_rows, so frac = cart · A⁻¹. The reciprocal lattice
    // rows satisfy a_i · b_j = 2π δ_ij → b_j = 2π · (A⁻¹)^T row j.
    // Equivalently frac_k = (b_k · cart) / (2π).
    let recip = lattice.reciprocal();

    let tau_frac_x = recip.a.dot(&tau_cart) / (2.0 * std::f64::consts::PI);
    let tau_frac_y = recip.b.dot(&tau_cart) / (2.0 * std::f64::consts::PI);
    let tau_frac_z = recip.c.dot(&tau_cart) / (2.0 * std::f64::consts::PI);

    let mut rho_sum = vec![0.0_f64; n_bins];
    let mut counts = vec![0_u64; n_bins];

    // Min-image wrap into [-0.5, +0.5) via (d - floor(d + 0.5)).
    // `.round()` in Rust rounds half-away-from-zero; numpy's np.round
    // rounds half-to-even. Using either here would make the Rust/
    // Python shell-average buckets disagree at `d = ±0.5`, producing
    // a spurious 20% asymmetry on 2-atom FCC cells where an atom sits
    // at (0.25, 0.25, 0.25) (C diamond, GaAs). The floor-based wrap
    // is platform-independent and places `d = ±0.5` consistently at
    // the `-0.5` edge.
    let wrap = |d: f64| d - (d + 0.5).floor();
    for i1 in 0..nx {
        let fx = (i1 as f64) / (nx as f64);
        let dfx = wrap(fx - tau_frac_x);
        for i2 in 0..ny {
            let fy = (i2 as f64) / (ny as f64);
            let dfy = wrap(fy - tau_frac_y);
            for i3 in 0..nz {
                let fz = (i3 as f64) / (nz as f64);
                let dfz = wrap(fz - tau_frac_z);

                // Cartesian displacement: dc = dfx·a + dfy·b + dfz·c
                let dc = dfx * lattice.a + dfy * lattice.b + dfz * lattice.c;
                let dist = dc.norm();

                let idx = i1 * ny * nz + i2 * nz + i3;
                let rho = rho_r[idx];

                // Binary search into edges: bin k has [r_edges[k], r_edges[k+1]).
                // For small n_bins, linear search is faster; keep simple.
                let mut bin = None;
                for k in 0..n_bins {
                    if dist >= r_edges[k] && dist < r_edges[k + 1] {
                        bin = Some(k);
                        break;
                    }
                }
                if let Some(k) = bin {
                    rho_sum[k] += rho;
                    counts[k] += 1;
                }
            }
        }
    }

    let mut rho_avg = vec![0.0_f64; n_bins];
    for k in 0..n_bins {
        if counts[k] > 0 {
            rho_avg[k] = rho_sum[k] / counts[k] as f64;
        }
    }
    (rho_avg, counts)
}

fn r_edges_ang() -> Vec<f64> {
    // 50 bins from 0 → 2.0 Å, 0.04 Å wide — matches Python.
    let n_bins: usize = 50;
    let r_max = 2.0_f64;
    (0..=n_bins).map(|k| r_max * (k as f64) / (n_bins as f64)).collect()
}

/// Result bundle for one system: shell-averaged ρ(r) around each atom
/// using both the final (post-clamp, post-renorm) and the pre-clamp
/// intermediate density, plus the SAD diagnostic statistics.
struct SystemDiagnosis {
    #[allow(dead_code, reason = "held for future per-bin debug prints")]
    crystal: Crystal,
    /// atom_label → shell-averaged ρ(r) from the SCF driver's final
    /// initial density (post-clamp + post-renorm).
    rho_final_per_atom: HashMap<String, Vec<f64>>,
    /// atom_label → shell-averaged ρ(r) from the pre-clamp intermediate.
    rho_pre_clamp_per_atom: HashMap<String, Vec<f64>>,
    /// atom_label → bin count (same for both pre and post since the
    /// grid and atom positions are identical).
    counts_per_atom: HashMap<String, Vec<u64>>,
    stats: pwdft_core::scf::initial_density::SadDiagnosticStats,
    n_electrons: f64,
}

/// Compute pwdft-core SAD on the fixed grid and return shell-averaged
/// ρ(r) for each atom — both pre-clamp and final.
fn diagnose_system(sys: &SystemSpec, pp_cache: &HashMap<&str, UpfPseudoPotential>) -> SystemDiagnosis {
    let crystal = build_crystal(sys);
    let mut unique_elements: Vec<&str> = sys.atoms.iter().map(|(e, _, _)| *e).collect();
    unique_elements.sort();
    unique_elements.dedup();
    let pps: Vec<&UpfPseudoPotential> = unique_elements.iter().map(|e| pp_cache.get(e).unwrap()).collect();

    // ecutwfc in eV (pwdft-core internal)
    let n_electrons: f64 = sys
        .atoms
        .iter()
        .map(|(e, _, _)| pp_cache.get(e).unwrap().z_valence)
        .sum();

    let ecutwfc_ev = sys.ecutwfc_ry * RY_TO_EV;

    let config = InitialDensityConfig::non_magnetic(crystal.atoms.len());
    let (dims, rho_pre, rho_final, stats) = build_sad_density_for_diagnostic_verbose(
        &crystal,
        &pps,
        n_electrons,
        ecutwfc_ev,
        sys.ecutrho_ratio,
        Some(sys.grid_dims),
        &config,
    );
    assert_eq!(dims, sys.grid_dims, "grid dim mismatch for {}", sys.name);

    let edges = r_edges_ang();
    let mut rho_final_per_atom = HashMap::new();
    let mut rho_pre_per_atom = HashMap::new();
    let mut counts_per_atom = HashMap::new();
    for (atom_idx, atom) in crystal.atoms.iter().enumerate() {
        let (elem, _, _) = sys.atoms[atom_idx];
        let tau = atom.cart_position(&crystal.lattice);
        let label = format!("{elem}{atom_idx}");
        let (rho_avg, counts) = shell_average(&rho_final, dims, &crystal.lattice, tau, &edges);
        rho_final_per_atom.insert(label.clone(), rho_avg);
        let (rho_avg_pre, _) = shell_average(&rho_pre, dims, &crystal.lattice, tau, &edges);
        rho_pre_per_atom.insert(label.clone(), rho_avg_pre);
        counts_per_atom.insert(label, counts);
    }
    SystemDiagnosis {
        crystal,
        rho_final_per_atom,
        rho_pre_clamp_per_atom: rho_pre_per_atom,
        counts_per_atom,
        stats,
        n_electrons,
    }
}

#[test]
#[ignore = "TSPL Tier-2: VGCH Phase 1c SAD diagnostic — needs CSV from pwdft-validate reference sad. Run with `-- --ignored`."]
fn test_vgch_sad_matches_qe_convention_reference() {
    let csv_path = PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join(CSV_REL_PATH);
    let Some(ref_rows) = load_reference_csv(&csv_path) else {
        panic!(
            "reference CSV not found at {}: regenerate via `uv run pwdft-validate reference sad`",
            csv_path.display()
        );
    };

    // Bucket reference rows: system → atom → Vec<(r_center, rho_ref)>
    let mut ref_buckets: BTreeMap<String, BTreeMap<String, Vec<(f64, f64)>>> = BTreeMap::new();
    for row in &ref_rows {
        ref_buckets
            .entry(row.system.clone())
            .or_default()
            .entry(row.atom_label.clone())
            .or_default()
            .push((row.r_bin_center, row.rho_avg_ref));
    }

    // Load all unique PPs once
    let mut pp_cache: HashMap<&str, UpfPseudoPotential> = HashMap::new();
    for sys in &systems() {
        for (elem, _, _) in &sys.atoms {
            pp_cache.entry(elem).or_insert_with(|| load_pp(elem));
        }
    }

    // ---- Raw-sample point-wise diff: bypasses shell-averaging to
    // directly compare ρ(r) at 200 deterministic grid points per
    // system. This is the authoritative test that the Rust ρ(r) array
    // matches the Python reference at machine precision for systems
    // where it should. ----
    let samples_path = PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join(SAMPLES_CSV_REL_PATH);
    let samples =
        load_samples_csv(&samples_path).unwrap_or_else(|| panic!("samples CSV missing at {}", samples_path.display()));

    // Group samples by system
    let mut sample_by_system: HashMap<String, Vec<&SampleRow>> = HashMap::new();
    for s in &samples {
        sample_by_system.entry(s.system.clone()).or_default().push(s);
    }

    println!("Raw-sample point-wise ρ(r) diff (post-clamp+renorm, 200 pts/system):");
    println!(
        "{:>12} {:>6} {:>14} {:>14} {:>7}",
        "system", "n", "max |Δρ|", "mean |Δρ|", "verdict"
    );
    println!("{}", "-".repeat(58));
    let mut worst_raw: Vec<(String, f64)> = Vec::new();
    let mut any_raw_fail = false;
    for sys in &systems() {
        let diag = diagnose_system(sys, &pp_cache);
        let samples_sys = sample_by_system
            .get(sys.name)
            .unwrap_or_else(|| panic!("no samples rows for system {}", sys.name));
        let [_nx, ny, nz] = sys.grid_dims;
        let mut max_d = 0.0_f64;
        let mut sum_d = 0.0_f64;
        let mut n_checked = 0_usize;
        // We need raw rho_final arrays here, not shell-averaged. Re-run.
        let (_dims, _pre, rho_final, _stats) =
            pwdft_core::scf::initial_density::build_sad_density_for_diagnostic_verbose(
                &diag.crystal,
                &{
                    let mut u: Vec<&str> = sys.atoms.iter().map(|(e, _, _)| *e).collect();
                    u.sort();
                    u.dedup();
                    u.iter().map(|e| pp_cache.get(e).unwrap()).collect::<Vec<_>>()
                },
                diag.n_electrons,
                sys.ecutwfc_ry * RY_TO_EV,
                sys.ecutrho_ratio,
                Some(sys.grid_dims),
                &InitialDensityConfig::non_magnetic(diag.crystal.atoms.len()),
            );
        for s in samples_sys {
            let idx = s.i1 * ny * nz + s.i2 * nz + s.i3;
            let d = (rho_final[idx] - s.rho_ref).abs();
            if d > max_d {
                max_d = d;
            }
            sum_d += d;
            n_checked += 1;
        }
        let mean_d = if n_checked > 0 { sum_d / n_checked as f64 } else { 0.0 };
        let verdict = if max_d <= TOL_RAW_SAMPLE_RHO { "ok" } else { "FAIL" };
        if max_d > TOL_RAW_SAMPLE_RHO {
            any_raw_fail = true;
        }
        println!(
            "{:>12} {:>6} {:>14.4e} {:>14.4e} {:>7}",
            sys.name, n_checked, max_d, mean_d, verdict
        );
        worst_raw.push((sys.name.to_string(), max_d));
    }
    assert!(
        !any_raw_fail,
        "VGCH Phase 1c: raw-sample ρ(r) diff exceeded {TOL_RAW_SAMPLE_RHO:.0e} on one or more \
         systems — the Bessel + structure-factor + IFFT pipeline no longer matches the \
         QE-convention Python reference. See table above and `pwdft-validate reference sad`."
    );

    println!();
    println!("SAD pipeline diagnostics:");
    println!(
        "{:>12} {:>8} {:>16} {:>16} {:>12}",
        "system", "N_el", "∫ρ pre-clamp", "neg mass clamp", "renorm fac"
    );
    println!("{}", "-".repeat(72));
    let mut system_stats: HashMap<String, pwdft_core::scf::initial_density::SadDiagnosticStats> = HashMap::new();
    for sys in &systems() {
        let diag = diagnose_system(sys, &pp_cache);
        system_stats.insert(sys.name.to_string(), diag.stats);
        println!(
            "{:>12} {:>8.3} {:>16.6} {:>16.3e} {:>12.6}",
            sys.name,
            diag.n_electrons,
            diag.stats.integrated_pre_clamp,
            diag.stats.negative_mass_clamped,
            diag.stats.renorm_scale,
        );
    }

    println!();
    println!("Shell-average Δρ(r) vs QE-convention reference (final = post-clamp + renorm):");
    println!(
        "{:>12} {:>4} {:>6} {:>14} {:>14} {:>7}",
        "system", "atom", "bins", "max |Δρ|", "mean |Δρ|", "verdict"
    );
    println!("{}", "-".repeat(80));

    let mut any_fail = false;
    let mut any_pre_clamp_fail = false;
    let mut per_atom_worst: Vec<(String, String, f64, f64)> = Vec::new();
    for sys in &systems() {
        let diag = diagnose_system(sys, &pp_cache);
        let rust_rho_avgs = &diag.rho_final_per_atom;
        let rust_rho_pre = &diag.rho_pre_clamp_per_atom;
        let rust_counts = &diag.counts_per_atom;

        let Some(atom_buckets) = ref_buckets.get(sys.name) else {
            println!("{:>12}  MISSING  reference rows", sys.name);
            any_fail = true;
            continue;
        };

        let tol = tolerance_override(sys.name);

        for (atom_label, ref_bins) in atom_buckets {
            let rust_rho = rust_rho_avgs.get(atom_label).unwrap_or_else(|| {
                panic!(
                    "system {} atom {} present in CSV but Rust did not compute it",
                    sys.name, atom_label
                )
            });
            let rust_rho_pre_atom = rust_rho_pre.get(atom_label).unwrap();
            let counts = rust_counts.get(atom_label).unwrap();
            let edges = r_edges_ang();

            let mut max_delta = 0.0_f64;
            let mut max_delta_r = 0.0_f64;
            let mut max_delta_rust = 0.0_f64;
            let mut max_delta_ref = 0.0_f64;
            let mut sum_delta = 0.0_f64;
            let mut n_checked = 0_usize;
            let mut local_fail = false;

            // Also track the pre-clamp-vs-ref delta — isolates Bessel+FFT
            // bit-perfection from clamp/renorm effects.
            let mut max_delta_pre = 0.0_f64;

            for (r_center, rho_ref) in ref_bins {
                let bin_idx = edges.windows(2).position(|w| *r_center >= w[0] && *r_center < w[1]);
                let Some(bin) = bin_idx else { continue };
                if counts[bin] == 0 {
                    continue;
                }
                let delta = (rust_rho[bin] - rho_ref).abs();
                let delta_pre = (rust_rho_pre_atom[bin] - rho_ref).abs();
                if delta > max_delta {
                    max_delta = delta;
                    max_delta_r = *r_center;
                    max_delta_rust = rust_rho[bin];
                    max_delta_ref = *rho_ref;
                }
                if delta_pre > max_delta_pre {
                    max_delta_pre = delta_pre;
                }
                sum_delta += delta;
                n_checked += 1;
                if delta > tol {
                    local_fail = true;
                }
            }

            let mean_delta = if n_checked > 0 {
                sum_delta / n_checked as f64
            } else {
                0.0
            };
            let verdict = if local_fail { "FAIL" } else { "ok" };
            println!(
                "{:>12} {:>4} {:>6} {:>14.4e} {:>14.4e} {:>7} | pre-clamp max|Δρ|={:.4e} @ r={:.2}Å rust={:+.4e} ref={:+.4e}",
                sys.name,
                atom_label,
                n_checked,
                max_delta,
                mean_delta,
                verdict,
                max_delta_pre,
                max_delta_r,
                max_delta_rust,
                max_delta_ref,
            );
            per_atom_worst.push((sys.name.to_string(), atom_label.clone(), max_delta, max_delta_pre));
            if local_fail {
                any_fail = true;
            }
            if max_delta_pre > TOL_PRE_CLAMP_RHO {
                any_pre_clamp_fail = true;
            }
        }
    }

    println!("\nPer-atom max|Δρ| summary (final | pre-clamp):");
    for (sys, atom, d, d_pre) in &per_atom_worst {
        println!("  {sys:>12} / {atom:>4} : final={d:.4e}   pre-clamp={d_pre:.4e} e/Å³");
    }

    assert!(
        !any_pre_clamp_fail,
        "VGCH Phase 1c: pre-clamp shell-average Δρ exceeded {TOL_PRE_CLAMP_RHO:.0e} on one or \
         more atoms — the Bessel + structure-factor + IFFT pipeline no longer matches the \
         QE-convention reference. This is the physics-correctness assertion; widening the \
         tolerance without a physics reason hides a bug."
    );
    assert!(
        !any_fail,
        "VGCH Phase 1c: SAD diagnostic — one or more (system, atom) shell-averaged \
         bins differed from the QE-convention reference by more than {TOL_RHO_E_PER_ANG3:.2e} e/Å³. \
         See per-atom table above. If this surfaces new physics, refine the \
         tolerance in tolerance_override() with a physics-first comment \
         rather than simply raising the default."
    );
}
