//! VGCH-2 Part B — Cu FCC transplant iter-1 diagnostic.
//!
//! Seeds pwdft-core' SCF with QE's converged Cu FCC density and runs
//! exactly one iteration. Compares per-term energies against QE's
//! converged output. If the per-term values match QE at iter-1 within
//! ~1 meV/term, the H3 mixer-basin hypothesis is confirmed: pwdft-core and
//! QE converge to different densities from their respective starting
//! points but are algebraically equivalent at QE's fixed point.
//!
//! If the per-term values do NOT match, H3 is cleared and the bug
//! lives in a finer-grained Hamiltonian-assembly / ψ-reconstruction
//! difference.
//!
//! ## How to regenerate the input density
//!
//! 1. Run QE Cu FCC SCF with `disk_io='medium'` (the QE input at
//!    `validation/reference/qe/cu_fcc_scf.in` has `disk_io='low'` and does NOT
//!    write `charge-density.dat`; copy it to `/tmp/vgch2b_cu/cu.in`,
//!    flip `disk_io` to `'medium'`, and run with the machine lock).
//! 2. Parse the resulting `<outdir>/cu.save/charge-density.dat` with
//!    `validation/src/pwdft_validation/scripts/vgch2_parse_qe_density.py --verbose` → produces
//!    `cu_rho_qe.bin` with `{mill, rho_g (e/Bohr³), b1, b2, b3}` in a
//!    flat little-endian binary bundle (VGCH2BIN magic — chosen over
//!    `.npz` to avoid pulling a ZIP crate into the test harness).
//! 3. This test loads that `.bin` via a tiny inline reader and maps
//!    the Miller-indexed density onto pwdft-core' FFT grid, converting
//!    units (e/Bohr³ → e/Å³) at the boundary.
//!
//! The test is Tier-2 gated (TSPL) — it pulls in the full SCF context
//! setup + per-term diagnostics on a heavy transition metal at the QE
//! 8×8×8 Γ-centered k-grid.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use num_complex::Complex64;
use pwdft_core::{
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

/// QE density bundle emitted by `validation/src/pwdft_validation/scripts/vgch2_parse_qe_density.py`.
///
/// Layout (all little-endian, see Python script for the canonical
/// specification):
///
///     magic:      8 bytes = b"VGCH2BIN"
///     version:    u32 = 1
///     gamma_only: u8
///     nspin:      i32
///     ngm:        i32
///     b1..b3:     3 × (3 × f64)   — 1/Bohr
///     mill:       ngm × 3 × i32   — C-order
///     rho_g:      ngm × nspin × c16   — e/Bohr³
struct QeDensity {
    nspin: i32,
    ngm: usize,
    /// Flat Miller indices: `mill[3*ig + c]` for coord c ∈ {0,1,2} of
    /// G-vector ig.
    mill: Vec<i32>,
    /// `rho_g[ig + ispin * ngm]` in e/Bohr³.
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
        "VGCH-2B: gamma_only=True densities are not supported; \
         re-run QE with a k-grid to write the full G-sphere"
    );
    let nspin = i32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
    off += 4;
    let ngm_i32 = i32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
    assert!(ngm_i32 > 0, "VGCH-2B: ngm must be positive, got {ngm_i32}");
    #[expect(
        clippy::cast_sign_loss,
        reason = "assert!(ngm_i32 > 0) guards against negative cast"
    )]
    let ngm = ngm_i32 as usize;
    off += 4;
    // Skip b1, b2, b3 (9 doubles) — not needed by the consumer once we
    // map via Miller indices. Keep the offset cursor sync'd.
    off += 9 * 8;

    let mill_bytes = 3 * ngm * 4;
    assert!(
        off + mill_bytes <= bytes.len(),
        "truncated mill block"
    );
    let mill: Vec<i32> = bytes[off..off + mill_bytes]
        .chunks_exact(4)
        .map(|c| i32::from_le_bytes(c.try_into().unwrap()))
        .collect();
    off += mill_bytes;

    assert!(nspin > 0, "VGCH-2B: nspin must be positive, got {nspin}");
    #[expect(
        clippy::cast_sign_loss,
        reason = "assert!(nspin > 0) guards against negative cast"
    )]
    let rho_bytes = (nspin as usize) * ngm * 16;
    assert!(
        off + rho_bytes <= bytes.len(),
        "truncated rho_g block"
    );
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

// ---------------------------------------------------------------------------
// Crystal + PP helpers
// ---------------------------------------------------------------------------

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
    let path = PathBuf::from(env!("CARGO_WORKSPACE_DIR"))
        .join("pseudopotentials/nc/lda")
        .join(format!("{element}.upf"));
    pwdft_core::pseudopotential::load(&path)
        .unwrap_or_else(|e| panic!("failed to load {}: {e}", path.display()))
}

/// Map (mill_index, rho_g_qe) pairs onto pwdft-core' FFT grid.
///
/// Uses the same Miller-index → FFT-index mapping as `scf::grid::miller_to_idx`
/// (inlined because that helper is `pub(crate)`). Converts units from
/// e/Bohr³ to e/Å³ at the boundary by dividing by `BOHR_TO_ANG³`.
/// Returns a dense `Vec<Complex64>` of length `dims[0]*dims[1]*dims[2]`
/// with zeros at every G not in the QE record — consistent with the
/// spherical cutoff QE writes.
fn scatter_rho_onto_fft_grid(
    mill: &[i32],
    rho_g_bohr3: &[Complex64],
    dims: [usize; 3],
    ngm: usize,
) -> Vec<Complex64> {
    let total = dims[0] * dims[1] * dims[2];
    let mut rho_g_fft = vec![Complex64::new(0.0, 0.0); total];
    // Unit conversion: rho(G) in e/Bohr³ → e/Å³. Multiply by 1/BOHR_TO_ANG³.
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
    reason = "dims bounded by FftGrid::MAX_FFT_DIM (1024); Miller indices are \
              bounded by ecut; `((n % d) + d) % d` is mathematically non-negative \
              so `as usize` loses no sign — identical pattern to \
              `src/scf/grid.rs::miller_to_idx`, duplicated here only because \
              that helper is `pub(crate)`."
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

/// Sanity: integrate ρ on the real-space grid after IFFT and check ≈ N_el.
/// ρ(G=0) under pwdft-core' `1/N`-forward convention is `⟨ρ⟩ = N_el / Ω`, so
/// `N_el = ρ(G=0) · Ω` — the cheapest consistency check on the transplant.
fn integrated_charge(rho_g_fft: &[Complex64], omega: f64, _dims: [usize; 3]) -> f64 {
    rho_g_fft[0].re * omega
}

// ---------------------------------------------------------------------------
// QE reference (validation/reference/qe/cu_fcc_scf.out, 8×8×8 Γ-centered, ecut=25 Ry).
// ---------------------------------------------------------------------------
const QE_CU_TOTAL_RY: f64 = -356.736_028_69;
const QE_CU_ONE_E_RY: f64 = -149.487_608_88;
const QE_CU_HARTREE_RY: f64 = 76.476_166_76;
const QE_CU_XC_RY: f64 = -41.095_186_82;
const QE_CU_EWALD_RY: f64 = -242.620_854_75;
const QE_CU_FERMI_EV: f64 = 19.2056;

// ---------------------------------------------------------------------------
// Test
// ---------------------------------------------------------------------------

fn cu_bin_path() -> Option<PathBuf> {
    // Look in a few well-known places. Default: /tmp/vgch2b_cu/cu_rho_qe.bin
    // (matches the QE-run + parse pipeline in the module header).
    let candidates = [
        PathBuf::from("/tmp/vgch2b_cu/cu_rho_qe.bin"),
        PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("validation/reference/qe/cu_rho_qe.bin"),
    ];
    candidates.into_iter().find(|p| p.is_file())
}

#[test]
#[ignore = "TSPL Tier-2: VGCH-2B transplant diagnostic — seed Cu FCC SCF from QE density, compare iter-1 per-term energies; regenerate `cu_rho_qe.npz` via validation/src/pwdft_validation/scripts/vgch2_parse_qe_density.py"]
fn test_cu_fcc_transplant_iter1() {
    let bin_path = cu_bin_path().unwrap_or_else(|| {
        panic!(
            "VGCH-2B: expected QE density at /tmp/vgch2b_cu/cu_rho_qe.bin or \
             validation/reference/qe/cu_rho_qe.bin — regenerate via \
             `uv run validation/src/pwdft_validation/scripts/vgch2_parse_qe_density.py \
                --rho <prefix>.save/charge-density.dat \
                --out /tmp/vgch2b_cu/cu_rho_qe.bin`"
        )
    });
    let qe = read_qe_density(&bin_path);
    assert_eq!(qe.nspin, 1, "VGCH-2B: Cu transplant expects nspin=1 density");
    let ngm = qe.ngm;
    let mill = &qe.mill;
    let rho_g_qe = &qe.rho_g;
    assert_eq!(mill.len(), 3 * ngm);
    assert_eq!(rho_g_qe.len(), ngm);

    // Mirror `tests/qe_validation.rs::test_cu_fcc_vs_qe` geometry.
    let crystal = fcc_crystal(3.61, vec![Atom::new(29, [0.0, 0.0, 0.0])]);
    let pp_cu = load_pp("Cu");
    let ecut_ev = 25.0 * RY_TO_EV;
    let basis = BasisSet::new(&crystal.lattice, ecut_ev);
    let nk = 8_u32;
    let kpts = kpoints::monkhorst_pack(
        nk,
        nk,
        nk,
        kpoints::KGridShift::GammaCentered,
        &crystal.lattice,
    );
    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);

    let params = ScfParams {
        n_bands: 14,
        max_iter: 1,
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
        // Force the FFT grid to match QE's 15×15×15 (see validation/reference/qe/cu_fcc_scf.out
        // "Dense grid: 1363 G-vectors FFT dimensions: (15, 15, 15)"). Without
        // this override, pwdft-core picks a 2/3/5-smooth grid that may differ,
        // and the transplant would have to interpolate.
        fft_grid: Some([15, 15, 15]),
        ..Default::default()
    };

    let dims = params.fft_grid.expect("explicit dims pinned above");
    let omega = crystal.lattice.volume();
    eprintln!(
        "  [Cu transplant] FFT dims = {dims:?}, ngm_qe = {ngm}, Ω = {omega:.3} Å³"
    );

    // Map QE's rho_g (e/Bohr³) onto pwdft-core' FFT grid (e/Å³).
    let rho_g_fft = scatter_rho_onto_fft_grid(mill, rho_g_qe, dims, ngm);

    // Sanity: ρ(G=0) · Ω ≈ N_el (= 19 for Cu).
    let n_el_est = integrated_charge(&rho_g_fft, omega, dims);
    eprintln!(
        "  [Cu transplant] ρ(G=0) · Ω = {n_el_est:.4} (expected 19.0 for Cu)"
    );
    assert!(
        (n_el_est - 19.0).abs() < 0.01,
        "transplanted ρ(G=0) integrates to {n_el_est}, expected 19.0 — \
         units or Miller-index mapping broken"
    );

    // Run iter-1 from the transplant.
    let result: TransplantIter1Result = run_scf_iter1_from_rho_g_fft(
        &crystal, &basis, &kpts, &[&pp_cu], &params, &symmetry, &rho_g_fft,
    )
    .unwrap_or_else(|e| panic!("transplant iter-1 failed: {e}"));

    // --- Side-by-side print ---
    let c = &result.components;
    let ours_one_electron = c.e_kinetic + c.e_local + c.e_local_g0_shift + c.e_nonlocal;
    let qe_total_ev = QE_CU_TOTAL_RY * RY_TO_EV;
    let qe_one_e_ev = QE_CU_ONE_E_RY * RY_TO_EV;
    let qe_hartree_ev = QE_CU_HARTREE_RY * RY_TO_EV;
    let qe_xc_ev = QE_CU_XC_RY * RY_TO_EV;
    let qe_ewald_ev = QE_CU_EWALD_RY * RY_TO_EV;

    eprintln!("\n==== Cu FCC transplant iter-1 per-term vs QE (eV) ====");
    eprintln!(
        "  {:<20}  {:>14}  {:>14}  {:>14}",
        "term", "pwdft (iter1)", "QE (converged)", "Δ (ours−QE)"
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "one-electron", ours_one_electron, qe_one_e_ev,
        ours_one_electron - qe_one_e_ev
    );
    // Breakdown of one-electron = kinetic + local(G≠0) + G=0 shift + nonlocal.
    // The G=0 shift = V_loc(G=0) · N_el captures the uniform background
    // re-added to compensate the G=0 zeroing of V_local in the Hamiltonian;
    // without it, per-band eigenvalues would differ from QE by V_loc(G=0).
    eprintln!(
        "    [breakdown] E_kin={:.4}  E_loc(G≠0)={:.4}  E_loc(G=0)·N_el={:.4}  E_nl={:.4}",
        c.e_kinetic, c.e_local, c.e_local_g0_shift, c.e_nonlocal
    );
    eprintln!(
        "    [v_loc_g0] = E_loc(G=0)·N_el / N_el = {:.4} eV",
        c.e_local_g0_shift / 19.0
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "Hartree", c.e_hartree, qe_hartree_ev, c.e_hartree - qe_hartree_ev
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "XC (bare)", c.e_xc, qe_xc_ev, c.e_xc - qe_xc_ev
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "Ewald", c.e_ewald, qe_ewald_ev, c.e_ewald - qe_ewald_ev
    );
    eprintln!(
        "  {:<20}  {:>14.6}  {:>14.6}  {:>+14.6}",
        "Total (E_KS)", result.total_energy, qe_total_ev,
        result.total_energy - qe_total_ev
    );
    eprintln!(
        "  {:<20}  {:>14.6}",
        "Fermi (iter1)", result.fermi_energy
    );
    eprintln!(
        "  {:<20}  {:>14.6e}",
        "Δρ (iter in-out)", result.delta_rho
    );
    eprintln!(
        "  E_HF = {:.6} eV, |E_HF − E_KS| = {:.6} eV",
        result.harris_foulkes_energy,
        (result.harris_foulkes_energy - result.total_energy).abs()
    );

    // The direct-sum identity `e_band = T + V_loc + V_nl + 2·E_H + E_vxc`
    // (and therefore Σ(components) == E_total) only holds at
    // self-consistency; at iter-1 after transplant, `e_band` pairs ψ from
    // H[ρ_in] with V_H[ρ_in]·ψ implicitly, while `e_H` and `e_vxc` are
    // built on ρ_out — so the gap is O(‖Δρ‖²) and can be large. We print
    // the mismatch as a diagnostic but do NOT assert — the PCFX
    // converged-only invariant is not applicable here.
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
        "  [diagnostic: Σ − E_total at iter-1 = {sum_err:.3e} eV; expected non-zero \
         because the iter-1 estimator mixes ρ_in eigenvalues with ρ_out \
         double-counting — identity only closes at self-consistency]"
    );

    // Diagnostic-only: do NOT fail on total-energy match. The whole point
    // of this test is to report the iter-1 delta so VGCH-2 Part B can
    // distinguish H3 (mixer basin) from a finer assembly bug. Guard only
    // against numerical nonsense (NaN/Inf).
    assert!(result.total_energy.is_finite(), "E_total not finite");

    // Top 5 Γ eigenvalues of iter-1 — these are the eigenvalues of
    // H[ρ_QE], not QE's own eigenvalues (QE might have been running
    // iterative Davidson with a different warm-start). Print for
    // visual diagnosis.
    if let Some(gamma_eigs) = result.eigenvalues.first() {
        let n_print = 8.min(gamma_eigs.len());
        eprintln!("  Γ eigenvalues (iter-1, eV, first {n_print}):");
        for (nb, &e) in gamma_eigs.iter().take(n_print).enumerate() {
            eprintln!("    band {nb}: {e:>10.4} eV");
        }
    }

    // Density integral sanity: ρ_in and ρ_out should both integrate to
    // N_el = 19 electrons. ρ_in is the real-space IFFT of the QE-seeded
    // ρ(G); ρ_out is the PCFX-symmetrized band sum.
    let dvol_total = omega / (dims[0] * dims[1] * dims[2]) as f64;
    let n_in: f64 = result.rho_r_in.iter().copied().sum::<f64>() * dvol_total;
    let n_out: f64 = result.rho_r_out.iter().copied().sum::<f64>() * dvol_total;
    eprintln!(
        "  [density sanity] ∫ρ_in = {n_in:.6} e,  ∫ρ_out = {n_out:.6} e  (expect {})",
        19.0_f64,
    );

    // QE Fermi level sanity: transplanted iter-1 Fermi should be close
    // to QE's converged Fermi for this density (within smearing noise).
    let df = (result.fermi_energy - QE_CU_FERMI_EV).abs();
    eprintln!(
        "  [Fermi] pwdft iter-1 E_F = {:.4} eV, QE = {:.4} eV, |Δ| = {:.4} eV",
        result.fermi_energy, QE_CU_FERMI_EV, df
    );
}
