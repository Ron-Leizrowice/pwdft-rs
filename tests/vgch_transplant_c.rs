//! VGCH-2E — C diamond transplant iter-1 diagnostic.
//!
//! Per VGCH-MECH's Class C taxonomy, C diamond is the 1.45 eV light-atom
//! outlier with the same opposite-sign one-electron / Hartree signature
//! as the Class A heavy-atom cells, but at light-atom magnitude.
//! Hypotheses: (a) NLCC Bessel-transform — C LDA PP has
//! `core_correction="T"`, identical code path to Fe/Cu; (b) Kerker q_TF
//! metallic default on a wide-gap insulator; (c) a Class-A-equivalent
//! mechanism at smaller magnitude.
//!
//! Mirrors the Cu transplant (`vgch_transplant_cu.rs`) byte-for-byte in
//! structure; only the crystal, PP, reference values, and FFT-grid pin
//! differ. See that file's header comments for the regeneration recipe.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use num_complex::Complex64;
use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    kpoints,
    pseudopotential::PseudopotentialData,
    scf::{
        mixing::MixingMode,
        smearing::SmearingScheme,
        transplant::{TransplantIter1Result, run_scf_iter1_from_rho_g_fft},
        ScfParams,
    },
    symmetry::SymmetryInfo,
};
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

const RY_TO_EV: f64 = 13.605_693_122_994;
const BOHR_TO_ANG: f64 = 0.529_177_210_903;
const BOHR3_TO_ANG3: f64 = BOHR_TO_ANG * BOHR_TO_ANG * BOHR_TO_ANG;

/// QE density bundle emitted by `scripts/validate/vgch2_parse_qe_density.py`.
/// Layout matches `tests/vgch_transplant_cu.rs::QeDensity`; see that file
/// for the canonical specification.
struct QeDensity {
    nspin: i32,
    ngm: usize,
    mill: Vec<i32>,
    rho_g: Vec<Complex64>,
}

fn read_qe_density(path: &Path) -> QeDensity {
    let bytes = std::fs::read(path)
        .unwrap_or_else(|e| panic!("failed to read {path:?}: {e}"));
    assert!(
        bytes.len() >= 8 + 4 + 1 + 4 + 4 + 72,
        "VGCH2BIN file too short ({} bytes)",
        bytes.len()
    );
    assert_eq!(
        &bytes[..8],
        b"VGCH2BIN",
        "magic mismatch: file is not a VGCH2BIN bundle"
    );
    let mut off = 8_usize;
    let version = u32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
    off += 4;
    assert_eq!(version, 1, "unsupported VGCH2BIN version {version}");
    let gamma_only = bytes[off] != 0;
    off += 1;
    assert!(
        !gamma_only,
        "VGCH-2E: gamma_only=True densities are not supported"
    );
    let nspin = i32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
    off += 4;
    let ngm_i32 = i32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
    assert!(ngm_i32 > 0, "VGCH-2E: ngm must be positive, got {ngm_i32}");
    #[expect(
        clippy::cast_sign_loss,
        reason = "assert!(ngm_i32 > 0) guards against negative cast"
    )]
    let ngm = ngm_i32 as usize;
    off += 4;
    off += 9 * 8; // skip b1..b3

    let mill_bytes = 3 * ngm * 4;
    assert!(off + mill_bytes <= bytes.len(), "truncated mill block");
    let mill: Vec<i32> = bytes[off..off + mill_bytes]
        .chunks_exact(4)
        .map(|c| i32::from_le_bytes(c.try_into().unwrap()))
        .collect();
    off += mill_bytes;

    assert!(nspin > 0, "VGCH-2E: nspin must be positive, got {nspin}");
    #[expect(
        clippy::cast_sign_loss,
        reason = "assert!(nspin > 0) guards against negative cast"
    )]
    let rho_bytes = (nspin as usize) * ngm * 16;
    assert!(off + rho_bytes <= bytes.len(), "truncated rho_g block");
    let rho_g: Vec<Complex64> = bytes[off..off + rho_bytes]
        .chunks_exact(16)
        .map(|c| {
            let re = f64::from_le_bytes(c[..8].try_into().unwrap());
            let im = f64::from_le_bytes(c[8..].try_into().unwrap());
            Complex64::new(re, im)
        })
        .collect();

    QeDensity {
        nspin,
        ngm,
        mill,
        rho_g,
    }
}

fn fcc_crystal(a_ang: f64, atoms: Vec<Atom>) -> Crystal {
    use nalgebra::Vector3;
    Crystal {
        lattice: Lattice::new(
            a_ang / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a_ang / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a_ang / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        ),
        atoms,
    }
}

fn load_pp(element: &str) -> PseudopotentialData {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("pseudopotentials/nc/lda")
        .join(format!("{element}.upf"));
    pwdft_rs::pseudopotential::load(&path)
        .unwrap_or_else(|e| panic!("failed to load {}: {e}", path.display()))
}

fn scatter_rho_onto_fft_grid(
    mill: &[i32],
    rho_g_bohr3: &[Complex64],
    dims: [usize; 3],
    ngm: usize,
) -> Vec<Complex64> {
    let total = dims[0] * dims[1] * dims[2];
    let mut rho_g_fft = vec![Complex64::new(0.0, 0.0); total];
    let unit_scale = 1.0 / BOHR3_TO_ANG3;
    for ig in 0..ngm {
        let n1 = mill[3 * ig];
        let n2 = mill[3 * ig + 1];
        let n3 = mill[3 * ig + 2];
        let idx = miller_to_idx(dims, n1, n2, n3);
        rho_g_fft[idx] = rho_g_bohr3[ig] * unit_scale;
    }
    rho_g_fft
}

#[expect(
    clippy::cast_possible_wrap,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "identical pattern to src/scf/grid.rs::miller_to_idx"
)]
fn miller_to_idx(dims: [usize; 3], n1: i32, n2: i32, n3: i32) -> usize {
    let d0 = dims[0] as i32;
    let d1 = dims[1] as i32;
    let d2 = dims[2] as i32;
    let i1 = (((n1 % d0) + d0) % d0) as usize;
    let i2 = (((n2 % d1) + d1) % d1) as usize;
    let i3 = (((n3 % d2) + d2) % d2) as usize;
    i1 * dims[1] * dims[2] + i2 * dims[2] + i3
}

fn integrated_charge(rho_g_fft: &[Complex64], omega: f64) -> f64 {
    rho_g_fft[0].re * omega
}

// ---------------------------------------------------------------------------
// QE reference (qe_validation/c_diamond_scf.out, 4×4×4 Γ-centered, ecut=30 Ry).
// ---------------------------------------------------------------------------
const QE_C_TOTAL_RY: f64 = -23.843_439_10;
const QE_C_ONE_E_RY: f64 = 8.502_873_41;
const QE_C_HARTREE_RY: f64 = 1.828_381_20;
const QE_C_XC_RY: f64 = -8.602_806_00;
const QE_C_EWALD_RY: f64 = -25.571_887_69;
const QE_C_FERMI_EV: f64 = 15.8873;

fn c_bin_path() -> Option<PathBuf> {
    let candidates = [
        PathBuf::from("/tmp/vgch2e_c/c_rho_qe.bin"),
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("qe_validation/c_rho_qe.bin"),
    ];
    candidates.into_iter().find(|p| p.is_file())
}

#[test]
#[ignore = "TSPL Tier-2: VGCH-2E transplant diagnostic — seed C diamond SCF from QE density, compare iter-1 per-term energies; regenerate c_rho_qe.bin via scripts/validate/vgch2_parse_qe_density.py"]
fn test_c_diamond_transplant_iter1() {
    let bin_path = c_bin_path().unwrap_or_else(|| {
        panic!(
            "VGCH-2E: expected QE density at /tmp/vgch2e_c/c_rho_qe.bin or \
             qe_validation/c_rho_qe.bin — regenerate via \
             `uv run scripts/validate/vgch2_parse_qe_density.py \
                --rho <prefix>.save/charge-density.dat \
                --out /tmp/vgch2e_c/c_rho_qe.bin`"
        )
    });
    let qe = read_qe_density(&bin_path);
    assert_eq!(qe.nspin, 1, "VGCH-2E: C transplant expects nspin=1 density");
    let ngm = qe.ngm;
    let mill = &qe.mill;
    let rho_g_qe = &qe.rho_g;
    assert_eq!(mill.len(), 3 * ngm);
    assert_eq!(rho_g_qe.len(), ngm);

    // Mirror `tests/qe_validation.rs::test_c_diamond_vs_qe` geometry.
    let crystal = fcc_crystal(
        3.567,
        vec![
            Atom::new(6, [0.00, 0.00, 0.00]),
            Atom::new(6, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_c = load_pp("C");
    let ecut_ev = 30.0 * RY_TO_EV;
    let basis = BasisSet::new(&crystal.lattice, ecut_ev);
    let nk = 4_u32;
    let kpts = kpoints::monkhorst_pack(
        nk,
        nk,
        nk,
        kpoints::KGridShift::GammaCentered,
        &crystal.lattice,
    );
    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);

    let params = ScfParams {
        n_bands: 8,
        max_iter: 1,
        conv_threshold: 1e-6,
        energy_threshold: 1e-5,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.01 * RY_TO_EV,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Broyden { kerker: true },
        nspin: 1,
        starting_magnetization: HashMap::new(),
        // QE uses 20×20×20 (see qe_validation/c_diamond_scf.out
        // "Dense grid: 1687 G-vectors FFT dimensions: (20, 20, 20)").
        fft_grid: Some([20, 20, 20]),
        ..Default::default()
    };

    let dims = params.fft_grid.expect("explicit dims pinned above");
    let omega = crystal.lattice.volume();
    eprintln!(
        "  [C transplant] FFT dims = {dims:?}, ngm_qe = {ngm}, Ω = {omega:.3} Å³"
    );

    let rho_g_fft = scatter_rho_onto_fft_grid(mill, rho_g_qe, dims, ngm);

    let n_el_est = integrated_charge(&rho_g_fft, omega);
    eprintln!(
        "  [C transplant] ρ(G=0) · Ω = {n_el_est:.4} (expected 8.0 for C diamond)"
    );
    assert!(
        (n_el_est - 8.0).abs() < 0.01,
        "transplanted ρ(G=0) integrates to {n_el_est}, expected 8.0 — \
         units or Miller-index mapping broken"
    );

    let result: TransplantIter1Result = run_scf_iter1_from_rho_g_fft(
        &crystal, &basis, &kpts, &[&pp_c], &params, &symmetry, &rho_g_fft,
    )
    .unwrap_or_else(|e| panic!("transplant iter-1 failed: {e}"));

    let c = &result.components;
    let ours_one_electron = c.e_kinetic + c.e_local + c.e_local_g0_shift + c.e_nonlocal;
    let qe_total_ev = QE_C_TOTAL_RY * RY_TO_EV;
    let qe_one_e_ev = QE_C_ONE_E_RY * RY_TO_EV;
    let qe_hartree_ev = QE_C_HARTREE_RY * RY_TO_EV;
    let qe_xc_ev = QE_C_XC_RY * RY_TO_EV;
    let qe_ewald_ev = QE_C_EWALD_RY * RY_TO_EV;

    eprintln!("\n==== C diamond transplant iter-1 per-term vs QE (eV) ====");
    eprintln!(
        "  {:<20}  {:>14}  {:>14}  {:>14}",
        "term", "pwdft (iter1)", "QE (converged)", "Δ (ours−QE)"
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "one-electron",
        ours_one_electron,
        qe_one_e_ev,
        ours_one_electron - qe_one_e_ev
    );
    eprintln!(
        "    [breakdown] E_kin={:.4}  E_loc(G≠0)={:.4}  E_loc(G=0)·N_el={:.4}  E_nl={:.4}",
        c.e_kinetic, c.e_local, c.e_local_g0_shift, c.e_nonlocal
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "Hartree",
        c.e_hartree,
        qe_hartree_ev,
        c.e_hartree - qe_hartree_ev
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "XC (bare)",
        c.e_xc,
        qe_xc_ev,
        c.e_xc - qe_xc_ev
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "Ewald",
        c.e_ewald,
        qe_ewald_ev,
        c.e_ewald - qe_ewald_ev
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "Total (E_KS)",
        result.total_energy,
        qe_total_ev,
        result.total_energy - qe_total_ev
    );
    eprintln!(
        "  {:<20}  {:>14.6}",
        "Fermi (iter1)", result.fermi_energy
    );
    eprintln!("  {:<20}  {:>14.6e}", "Δρ (iter in-out)", result.delta_rho);
    eprintln!(
        "  E_HF = {:.6} eV, |E_HF − E_KS| = {:.6} eV",
        result.harris_foulkes_energy,
        (result.harris_foulkes_energy - result.total_energy).abs()
    );

    assert!(result.total_energy.is_finite(), "E_total not finite");

    if let Some(gamma_eigs) = result.eigenvalues.first() {
        let n_print = 8.min(gamma_eigs.len());
        eprintln!("  Γ eigenvalues (iter-1, eV, first {n_print}):");
        for (nb, &e) in gamma_eigs.iter().take(n_print).enumerate() {
            eprintln!("    band {nb}: {e:>10.4} eV");
        }
        // QE reference (qe_validation/c_diamond_scf.out at Γ, eV):
        //   -8.1456  14.0232  14.0232  14.0232  19.3568  19.3568  19.3568  27.2505
        eprintln!("  QE Γ eigenvalues (converged, eV):");
        let qe_eigs = [-8.1456, 14.0232, 14.0232, 14.0232, 19.3568, 19.3568, 19.3568, 27.2505];
        for (nb, &e) in qe_eigs.iter().take(n_print).enumerate() {
            eprintln!("    band {nb}: {e:>10.4} eV");
        }
    }

    let dvol_total = omega / (dims[0] * dims[1] * dims[2]) as f64;
    let n_in: f64 = result.rho_r_in.iter().copied().sum::<f64>() * dvol_total;
    let n_out: f64 = result.rho_r_out.iter().copied().sum::<f64>() * dvol_total;
    eprintln!(
        "  [density sanity] ∫ρ_in = {n_in:.6} e,  ∫ρ_out = {n_out:.6} e  (expect {})",
        8.0_f64,
    );

    let df = (result.fermi_energy - QE_C_FERMI_EV).abs();
    eprintln!(
        "  [Fermi] pwdft iter-1 E_F = {:.4} eV, QE = {:.4} eV, |Δ| = {:.4} eV",
        result.fermi_energy, QE_C_FERMI_EV, df
    );
}
