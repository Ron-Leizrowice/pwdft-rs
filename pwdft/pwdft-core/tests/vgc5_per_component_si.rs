//! VGC5 (VGCMP Phase 5) — Per-component energy accounting.
//!
//! Runs Si diamond and Fe BCC SCF with the VGC5 parameters (matching
//! `data/qe/si_scf.in` and `fe_bcc_fm_scf.in`) and:
//!
//! 1. Prints the per-component decomposition (`EnergyComponents`) alongside
//!    the QE reference values parsed from `data/qe/*.out`.
//! 2. Pins the pwdft-core per-component values as regression guards.
//!
//! This is Phase 5 of the VGCMP audit: the pseudopotential -> Hamiltonian
//! assembly pipeline is bit-correct (Phases 1-4). This test localizes the
//! residual Si 13.4 eV gap to specific term(s) of the per-component sum.
//!
//! IMPORTANT: Tolerances are loose (0.05 eV per term) because these are
//! regression pins, not QE-match assertions. When the eventual fix
//! proposal (VGFX / EWFX / similar) closes the gap, bump the pinned
//! numbers and keep the tolerance loose; QE-match assertions go into
//! `tests/qe_validation.rs`.
//!
//! Pin values are refreshed post-NCFX (NLCC core-density fix, PR #?).
//! Pre-NCFX baselines are retained as `// PRE-NCFX:` comments on each
//! pin line for context — they are useful to see the impact of NCFX at
//! a glance and help anyone re-running the audit in the future.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use nalgebra::Vector3;
use pwdft_core::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    kpoints,
    pseudopotential::PseudopotentialData,
    scf::{self, ScfParams, ScfResult, mixing::MixingMode, smearing::SmearingScheme},
    symmetry::SymmetryInfo,
};
use std::collections::HashMap;

const RY_TO_EV: f64 = 13.605_693_122_994;

// ----------------------------------------------------------------------------
// Fixtures
// ----------------------------------------------------------------------------

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
    let path = std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR"))
        .join("pseudopotentials/nc/lda")
        .join(format!("{element}.upf"));
    pwdft_core::pseudopotential::load(&path)
        .unwrap_or_else(|e| panic!("failed to load {}: {e}", path.display()))
}

// ----------------------------------------------------------------------------
// Print helpers
// ----------------------------------------------------------------------------

/// QE per-term reference (eV) for side-by-side reporting.
struct QeReference {
    total:        f64,
    one_electron: f64,
    hartree:      f64,
    xc:           f64,
    ewald:        f64,
    fermi:        f64,
}

impl QeReference {
    fn si() -> Self {
        // From data/qe/si_scf.out (see pwdft/pwdft-validation/pwdft_validation/scripts/vgc5_qe_si_components.csv).
        // Parsed: one-electron=4.86744632 Ry, hartree=1.11010670 Ry,
        //         xc=-6.20301507 Ry, ewald=-16.79667313 Ry,
        //         total=-17.02299344 Ry (includes -TS).
        Self {
            total:        -17.022_993_44 * RY_TO_EV,
            one_electron:   4.867_446_32 * RY_TO_EV,
            hartree:        1.110_106_70 * RY_TO_EV,
            xc:            -6.203_015_07 * RY_TO_EV,
            ewald:        -16.796_673_13 * RY_TO_EV,
            fermi:          6.3449,
        }
    }

    fn fe() -> Self {
        // From data/qe/fe_bcc_fm_scf.out.
        // one-electron=-50.85651213 Ry, hartree=26.64115855 Ry,
        // xc=-28.90402134 Ry, ewald=-171.77906580 Ry,
        // total=-224.91744934 Ry.
        Self {
            total:       -224.917_449_34 * RY_TO_EV,
            one_electron: -50.856_512_13 * RY_TO_EV,
            hartree:       26.641_158_55 * RY_TO_EV,
            xc:           -28.904_021_34 * RY_TO_EV,
            ewald:       -171.779_065_80 * RY_TO_EV,
            fermi:         26.2006,
        }
    }
}

/// Print pwdft-core components vs QE side-by-side, log the deltas.
fn print_side_by_side(label: &str, result: &ScfResult, qe: &QeReference) {
    let c = &result.components;
    // pwdft-core "one-electron" equivalent = E_kin + E_loc + E_NL + V_loc(G=0)*N_el
    let ours_one_electron = c.e_kinetic + c.e_local + c.e_local_g0_shift + c.e_nonlocal;

    eprintln!("\n===== {label} : per-component decomposition (eV) =====");
    eprintln!(
        "  {:<24}  {:>14}  {:>14}  {:>12}",
        "term", "pwdft-core", "QE", "Δ (ours−QE)"
    );
    eprintln!(
        "  {:<24}  {:>14}  {:>14}  {:>12}",
        "----", "--------", "--", "-----------"
    );
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14}  {:>12}",
        "E_band", c.e_band, "(n/a: QE)", ""
    );
    eprintln!("  ---- decomposition ----");
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
    eprintln!(
        "  {:<24}  {:>14.6}  {:>14.6}  {:>+12.6}",
        "E_F",
        result.fermi_energy,
        qe.fermi,
        result.fermi_energy - qe.fermi
    );

    // Consistency: components should sum to total_energy.
    //
    // Post-PCFX (`proposals/completed/PCFX-symmetrize-rho-g-space.md`),
    // density symmetrization is exact in G-space for any fractional
    // translation τ, so the identity closes to machine precision
    // regardless of grid–τ commensurability. Pre-PCFX the real-space
    // symmetrizer rounded τ=(¼,¼,¼) to nearest-integer on an 18³ grid
    // (18·¼ = 4.5 ∉ ℤ), smearing density into wrong grid points and
    // producing a plateau residual of 1.204 eV on Si that was invariant
    // under conv_threshold tightening — see the PCRS investigation.
    //
    // Post-TSEN (2026-04-19), `total_energy` additionally carries
    // `−TS` (= `c.e_smearing`); that field is included in the sum so
    // the identity still closes to machine precision on metals.
    let e_sum = c.e_kinetic + c.e_local + c.e_local_g0_shift + c.e_nonlocal
        + c.e_hartree + c.e_xc + c.e_ewald + c.e_smearing;
    let sum_err = e_sum - result.total_energy;
    eprintln!(
        "  [self-check] Σ(components) = {e_sum:.6} eV, E_total = {:.6} eV, Δ = {sum_err:.2e} eV",
        result.total_energy
    );
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

/// VGC5 Si diamond per-component audit.
///
/// Matches `data/qe/si_scf.in` parameters:
///   a = 5.431 Å, ecutwfc = 15 Ry = 204.085 eV, 4×4×4 MP, FD smearing
///   degauss = 0.01 Ry, conv_thr = 1e-8.
///
/// Pins per-component values; prints side-by-side vs QE for localization
/// of the 13.4 eV gap.
#[test]
#[ignore = "TSPL Tier-2: Si diamond 4×4×4 SCF at ecut=15 Ry with conv=1e-8 (per-component regression pins); run with cargo test -- --ignored when touching scf/, potential/, pseudopotential/, or symmetry/ paths"]
fn vgc5_si_per_component() {
    let crystal = fcc_crystal(
        5.431,
        vec![
            Atom::new(14, [0.00, 0.00, 0.00]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_si = load_pp("Si");

    let basis = BasisSet::new(&crystal.lattice, 15.0 * RY_TO_EV);
    // VGC5 pins were captured on the MP-1976 shifted grid; preserve that
    // convention to keep the per-component regression pins valid. Moving
    // to Γ-centered (the MPSH default elsewhere) would shift every pin by
    // k-sampling noise, which is out of scope for this regression guard.
    let kpts = kpoints::monkhorst_pack(4, 4, 4, kpoints::KGridShift::MP1976, &crystal.lattice);

    let params = ScfParams {
        n_bands: 8,
        max_iter: 80,
        conv_threshold: 1e-8,
        energy_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.01 * RY_TO_EV,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Plain,
        nspin: 1,
        starting_magnetization: HashMap::new(),
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(
        &crystal, &basis, &kpts, &[&pp_si], &params, &symmetry,
    )
    .expect("Si SCF should converge");

    eprintln!(
        "  Si SCF: {} iters, |HF-KS|={:.2e} eV, E_total={:.6} eV",
        result.n_iterations,
        (result.harris_foulkes_energy - result.total_energy).abs(),
        result.total_energy,
    );

    print_side_by_side("Si diamond", &result, &QeReference::si());

    // -------- Regression pins (pwdft-core, NOT QE-match) --------
    // Post-PCFX values (see `proposals/completed/PCFX-symmetrize-rho-g-space.md`).
    // PCFX moved density symmetrization from real-space (rounding-sensitive,
    // wrong for non-symmorphic τ on incompatible grids) to G-space (exact
    // via phase factors). That shifted the Si total by ~23 meV
    // (-231.865 → -231.843 eV) — this is the bug being fixed; the
    // pre-PCFX pin was numerically stable but off by the symmetrization
    // artefact documented in the PCRS residual scan.
    //
    // Per-component tolerances are 0.05 eV (CI / machine noise budget).
    let tol = 0.05;
    let pin = |name: &str, got: f64, expected: f64| {
        let d = (got - expected).abs();
        assert!(
            d <= tol,
            "VGC5 Si pin: {name} = {got:.4} eV, pinned {expected:.4} eV, |Δ|={d:.4} eV > {tol} eV",
        );
    };

    let c = &result.components;
    // PCFX correction: τ now applied exactly in G-space; pre-PCFX was E_total=-231.8653 eV.
    //
    // Post-VGCH-SiEF-B1 (2026-04-19): the `V_loc(G=0)` DC offset
    // (`10.7447 eV` = 2 Si atoms × N_el/atom scaling) now lives on
    // the Hamiltonian diagonal, so both `e_band` and `e_local` have
    // moved up by `V_loc(G=0)·N_el` while `e_local_g0_shift` is now
    // `0.0`. The per-component sum and `E_total` are algebraically
    // unchanged.
    pin("E_band",             c.e_band,              7.7748); // pre-B1:  -2.9699
    pin("E_kinetic",          c.e_kinetic,          83.4120); // pre-PCFX: 83.7203
    pin("E_local",            c.e_local,           -52.9772); // pre-B1: -63.7219 (G≠0 only)
    pin("E_local(G=0)*N_el",  c.e_local_g0_shift,    0.0);    // pre-B1:  10.7447
    pin("E_nonlocal",         c.e_nonlocal,         35.7641); // pre-PCFX: 35.9570
    pin("E_hartree",          c.e_hartree,          14.8249); // pre-PCFX: 14.3111
    pin("E_xc",               c.e_xc,              -84.3474); // pre-PCFX: -84.7026
    pin("E_ewald",            c.e_ewald,          -228.5192); // unchanged (lattice-only)
    pin("E_total",            result.total_energy, -231.8429); // pre-PCFX: -231.8653; post-TSEN: −TS on Si is ~10 meV

    // PCFX regression guard: the per-component identity closes to machine
    // precision post-fix (was 1.204 eV plateau pre-PCFX). Target was
    // ≤ 1e-5 eV; observed ~3.5e-11 eV. Post-TSEN the identity includes
    // `c.e_smearing` (= −TS).
    let e_sum = c.e_kinetic + c.e_local + c.e_local_g0_shift + c.e_nonlocal
        + c.e_hartree + c.e_xc + c.e_ewald + c.e_smearing;
    let sum_residual = (e_sum - result.total_energy).abs();
    assert!(
        sum_residual < 1e-5,
        "PCFX per-component self-check regressed: Σ - E_total = {sum_residual:.2e} eV, \
         pre-PCFX plateau was 1.204 eV"
    );
}

/// VGC5 Fe BCC per-component audit.
///
/// Matches `data/qe/fe_bcc_fm_scf.in` parameters:
///   a = 2.87 Å, ecutwfc = 15 Ry, 8×8×8 MP, FD smearing degauss = 0.02 Ry,
///   nspin = 2 with starting_magnetization(Fe) = 0.5.
///
/// Fe BCC matches QE to ~0.02 eV at the present state; this test pins the
/// per-component values to confirm the 13.4 eV Si gap is geometry-specific.
#[test]
#[ignore = "TSPL Tier-2: Fe BCC 4×4×4 SCF at ecut=15 Ry with max_iter=150 (per-component regression pins); run with cargo test -- --ignored when touching scf/, potential/, pseudopotential/, or NLCC paths"]
fn vgc5_fe_per_component() {
    let crystal = bcc_crystal(2.87, Atom::new(26, [0.0, 0.0, 0.0]));
    let pp_fe = load_pp("Fe");

    let basis = BasisSet::new(&crystal.lattice, 15.0 * RY_TO_EV);
    // 4x4x4 MP (not 8x8x8 as in QE ref): SCF at 8x8x8 fails to converge
    // in 80 iters with current mixer defaults. The 4x4x4 case captures
    // enough information to localize the per-component discrepancy
    // relative to QE without inheriting k-sampling noise >~10 meV.
    //
    // Pinned on the MP-1976 shifted grid; preserve it to keep the Fe
    // per-component regression pins valid (MPSH note: swapping to
    // Γ-centered would shift every pin by k-sampling noise).
    let kpts = kpoints::monkhorst_pack(4, 4, 4, kpoints::KGridShift::MP1976, &crystal.lattice);

    // NB: nspin=1 used here (not nspin=2 as in QE ref). The PseudoDojo Fe PP
    // at ecut=15 Ry collapses to non-magnetic anyway (see reference_data.toml
    // note on fe_bcc_fm), and the nspin=2 run fails to converge at these
    // parameters. Non-spin reproduces the Fe energy landscape we want to
    // audit for per-component decomposition.
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
    let result = scf::run_scf(
        &crystal, &basis, &kpts, &[&pp_fe], &params, &symmetry,
    )
    .expect("Fe SCF should converge");

    print_side_by_side("Fe BCC (FM)", &result, &QeReference::fe());

    // -------- Regression pins (pwdft-core) --------
    // Pinned from the first run under VGC5. Update as needed.
    let tol = 0.10; // Fe runs are slightly noisier; 0.1 eV safe.
    let pin = |name: &str, got: f64, expected: f64| {
        let d = (got - expected).abs();
        assert!(
            d <= tol,
            "VGC5 Fe pin: {name} = {got:.4} eV, pinned {expected:.4} eV, |Δ|={d:.4} eV > {tol} eV",
        );
    };

    let c = &result.components;
    assert!(result.total_energy.is_finite(), "Fe total energy is NaN");
    assert!(c.e_band.is_finite() && c.e_kinetic.is_finite(), "Fe components are NaN");

    // Post-NCFX pins: nspin=1, 4×4×4 MP, ecut=15 Ry, Kerker, 150 iters.
    // Not byte-matched to QE (nspin=2, 8×8×8); these are regression guards
    // only. See `proposals/completed/NCFX-nlcc-core-density-fix.md` for the
    // impact of NCFX on these numbers.
    // Post-VGCH-SiEF-B1 (2026-04-19): `e_band` and `e_local` each
    // gained `V_loc(G=0)·N_el = +82.7774 eV` (Fe, 1 atom, N_el = 16),
    // `e_local_g0_shift` is now 0. `E_total` is algebraically
    // unchanged.
    pin("E_band",             c.e_band,           -336.8252); // pre-B1: -419.6026
    pin("E_kinetic",          c.e_kinetic,         942.0082); // PRE-NCFX:  942.2052
    pin("E_local",            c.e_local,         -1666.0622); // pre-B1: -1748.8396 (G≠0 only)
    pin("E_local(G=0)*N_el",  c.e_local_g0_shift,    0.0);    // pre-B1:   82.7774
    pin("E_nonlocal",         c.e_nonlocal,         38.7678); // PRE-NCFX:   38.9673
    pin("E_hartree",          c.e_hartree,         363.1071); // PRE-NCFX:  363.4925
    pin("E_xc",               c.e_xc,             -392.5675); // PRE-NCFX: -442.1090 (Δ_QE: −48.85 → +0.69)
    pin("E_ewald",            c.e_ewald,         -2337.1672); // PRE-NCFX: -2337.1672 (unchanged)
    // Post-TSEN (2026-04-19) `total_energy` includes `−TS`. Fe at σ =
    // 0.02 Ry (0.272 eV) and 4×4×4 MP (nspin=1) carries a metallic
    // `−TS` of ≈ −430 meV, shifting the pin down accordingly. If this
    // pin needs to move by more than the 0.1 eV tolerance, revisit
    // the smearing entropy formula (`src/scf/smearing.rs::entropy_ts`),
    // not the per-component decomposition.
    pin("E_total",            result.total_energy, -3052.3209); // PRE-TSEN: -3051.8909; PRE-NCFX: -3101.2389

    // PCFX regression guard. Fe BCC Im-3m is symmorphic (τ=0), so the
    // pre-PCFX residual was already small (~0.045 eV at conv=1e-6) and
    // PCFX doesn't change this case materially — but keep the guard so
    // any future regression in the G-space symmetrizer surfaces here.
    // Post-TSEN the sum includes `c.e_smearing`.
    let e_sum = c.e_kinetic + c.e_local + c.e_local_g0_shift + c.e_nonlocal
        + c.e_hartree + c.e_xc + c.e_ewald + c.e_smearing;
    let sum_residual = (e_sum - result.total_energy).abs();
    assert!(
        sum_residual < 0.1,
        "PCFX per-component self-check: Σ - E_total = {sum_residual:.2e} eV on Fe"
    );
}

// ----------------------------------------------------------------------------
// MADOC band-sum identity (TRV2 Finding #1)
// ----------------------------------------------------------------------------
//
// MADOC Phase A (PR #87) documented in `src/scf/energy.rs` the Kohn-Sham
// double-counting identity that every converged SCF must satisfy:
//
//     E_band = E_kinetic + E_local + E_nonlocal + 2·E_hartree + E_vxc
//
// Derivation: each band eigenvalue is the expectation value of the full
// Kohn-Sham Hamiltonian,
//     ε_{n,k} = ⟨ψ_{n,k}| T + V_ext + V_H + V_xc |ψ_{n,k}⟩,
// so summing with Fermi-Dirac occupations and k-point weights,
//     E_band = E_kin + E_loc + E_nl + E_H⟨ψ|ψ⟩-coupling + ∫ρ·V_xc dr.
// V_H is self-linear in ρ (picks up a factor of 2 upon integration), V_xc
// is not — hence the 2·E_H + E_vxc pattern.
//
// This identity is independent of the direct-sum identity
// E_total = Σ_components (which is trivially true by construction at the
// site where `EnergyComponents` is built). A factor-2 Hartree bug or a
// sign flip in V_xc would break band-sum while leaving direct-sum intact,
// so VGC5's existing self-check cannot catch this bug class.
//
// The residual is O(Δρ) off the fixed point (because E_H and E_vxc in
// `EnergyComponents` are evaluated on the OUTPUT density while ε_{n,k}
// came from diagonalizing H built on the INPUT density). At tight SCF
// convergence the identity must close to well within that noise floor.

/// MADOC band-sum identity for Si (nspin=1).
///
/// Si diamond at the same ecut/k-mesh/conv as `vgc5_si_per_component`
/// above; reuses that fixture rather than standing up a new one.
/// Asserts the double-counting identity
/// `E_band = E_kin + E_loc + E_nl + 2·E_H + E_vxc` closes to well
/// below 1 μeV on a conv=1e-8 SCF.
#[test]
#[ignore = "TSPL Tier-2: Si diamond 4×4×4 SCF at ecut=15 Ry with conv=1e-8 (MADOC band-sum identity pin); run with cargo test -- --ignored when touching scf/energy, potential/, or density paths"]
fn test_madoc_band_sum_identity_si() {
    let crystal = fcc_crystal(
        5.431,
        vec![
            Atom::new(14, [0.00, 0.00, 0.00]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    );
    let pp_si = load_pp("Si");

    let basis = BasisSet::new(&crystal.lattice, 15.0 * RY_TO_EV);
    // Γ-centered (MPSH default) — band-sum identity is exact regardless
    // of shift, so use the new default for consistency with the rest of
    // qe_validation.
    let kpts = kpoints::monkhorst_pack(4, 4, 4, kpoints::KGridShift::GammaCentered, &crystal.lattice);

    let params = ScfParams {
        n_bands: 8,
        max_iter: 80,
        conv_threshold: 1e-8,
        energy_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.01 * RY_TO_EV,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Plain,
        nspin: 1,
        starting_magnetization: HashMap::new(),
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(
        &crystal, &basis, &kpts, &[&pp_si], &params, &symmetry,
    )
    .expect("Si SCF should converge at conv=1e-8");

    let c = &result.components;
    let e_band_from_identity =
        c.e_kinetic + c.e_local + c.e_nonlocal + 2.0 * c.e_hartree + c.e_vxc;
    let residual = (e_band_from_identity - c.e_band).abs();

    eprintln!(
        "MADOC band-sum identity (Si, nspin=1): E_band={:.10} eV  \
         E_kin+E_loc+E_nl+2·E_H+E_vxc={e_band_from_identity:.10} eV  \
         residual={residual:.3e} eV  ({} iters)",
        c.e_band, result.n_iterations,
    );

    // Tolerance rationale.
    //
    // CPU (default features, all f64): at conv_threshold = 1e-8 on Si the
    // residual is dominated by the remaining O(Δρ) mismatch between E_band
    // (built from eigenvalues of H[ρ_in]) and E_H/E_vxc (evaluated on
    // ρ_out). Observed residual on this fixture: ~1e-9 eV (release build,
    // LTO) to ~5e-11 eV (dev build with faer optimized). 1e-7 eV leaves
    // ~100× empirical headroom over the release-build observation.
    //
    // GPU (`--features gpu`): the Hartree / XC / V_eff grid operations
    // run in f32 inside WGSL kernels, with conversion at the host
    // boundary. The f32 precision ceiling (~1e-7 relative) accumulated
    // over ~18³ grid points and ~O(10³) eV of per-grid potential adds up
    // to an identity residual of ~5.58e-6 eV on this fixture (observed
    // on Apple M2 Metal, 2026-04-18 under ALOC-F5 audit). 5e-5 eV leaves
    // ~10× headroom over that observation.
    //
    // Both bounds still catch the bug classes this test exists to pin:
    //   - factor-2 in E_hartree: |ΔE_H| ≈ 14 eV → 2·10⁸× (CPU) /
    //     2·10⁵× (GPU) the tolerance.
    //   - sign-flipped e_vxc: |ΔE_vxc| ≈ 10 eV → 10⁸× (CPU) /
    //     10⁵× (GPU) the tolerance.
    //   - spin-channel swap (on the Fe counterpart): > 0.5 eV, still
    //     10⁴× the GPU tolerance.
    //
    // We use a feature-gated tolerance rather than a single loose bound
    // so the CPU path keeps its nano-eV pin — tightening that silently
    // on GPU would hide a future real regression in the Hartree/XC host
    // path, which runs in f64 even under `--features gpu`.
    let tol = if cfg!(feature = "gpu") { 5e-5 } else { 1e-7 };
    assert!(
        residual < tol,
        "MADOC band-sum identity violated on Si (nspin=1):\n  \
         E_band                                = {:.10} eV\n  \
         E_kin + E_loc + E_nl + 2·E_H + E_vxc  = {e_band_from_identity:.10} eV\n  \
         residual                              = {residual:.3e} eV  (tol {tol:.1e})\n\
         This identity is required by the Kohn-Sham eigenvalue equation; a\n\
         nonzero residual indicates a factor-2 bug in E_hartree, a sign\n\
         flip in e_vxc, or a missing term in the Hamiltonian assembly.",
        c.e_band,
    );
}

/// MADOC band-sum identity for Fe BCC (nspin=2, LSDA).
///
/// Exercises the spin-polarized driver and the coupled-channel mixer
/// (CCMX). Setup mirrors `test_ccmx_fe_free_magnetization_converges`
/// in `tests/spin_polarization.rs` (nc/lda/Fe.upf, 4×4×4, Kerker) but
/// with a tighter `conv_threshold` so the identity residual is well
/// below the physics-discrimination threshold. The identity for LSDA is
/// identical in form to the non-spin case —
/// `E_band = E_kin + E_loc + E_nl + 2·E_H + E_vxc` — with the caveat
/// that `E_kin`, `E_nl`, and `E_vxc` are each sums over both spin
/// channels; `E_H` is still built from the total density only (Hartree
/// does not couple spin).
#[test]
#[ignore = "TSPL Tier-2: Fe BCC nspin=2 4×4×4 SCF at ecut=15 Ry (MADOC band-sum identity pin, LSDA + CCMX); run with cargo test -- --ignored when touching scf/ spin driver, mixing, energy, or XC paths"]
fn test_madoc_band_sum_identity_fe_bcc() {
    let crystal = bcc_crystal(2.87, Atom::new(26, [0.0, 0.0, 0.0]));
    let pp_fe = load_pp("Fe");

    let ecut = 15.0 * RY_TO_EV;
    let basis = BasisSet::new(&crystal.lattice, ecut);
    // Γ-centered (MPSH default) — band-sum identity is exact regardless
    // of shift.
    let kpts = kpoints::monkhorst_pack(4, 4, 4, kpoints::KGridShift::GammaCentered, &crystal.lattice);

    let mut starting_mag = HashMap::new();
    starting_mag.insert("Fe".to_string(), 0.5);

    let params = ScfParams {
        n_bands: 12,
        max_iter: 150,
        // Fe BCC nspin=2 at ecut=15 Ry 4×4×4 converges to `delta ~ 5e-7`
        // reliably under CCMX + Kerker in ≲ 80 iters. Tighter than 1e-6
        // does not converge inside a reasonable iter budget on this PP;
        // looser than 1e-6 leaves the O(Δρ) systematic drift between
        // `V_H[ρ_in]` (inside H, determines ε_{n,k}) and `E_H[ρ_out]`
        // (stored in `components.e_hartree`) large enough to dominate
        // the band-sum identity residual. 1e-6 is the minimum-workable
        // setting that both converges and yields an informative
        // identity residual of O(1e-2 eV) on this system.
        conv_threshold: 1e-6,
        energy_threshold: 1e-5,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.02 * RY_TO_EV,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Kerker { q_tf: None },
        nspin: 2,
        starting_magnetization: starting_mag,
        ..Default::default()
    };

    let symmetry = SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(
        &crystal, &basis, &kpts, &[&pp_fe], &params, &symmetry,
    )
    .expect("Fe BCC nspin=2 SCF should converge under CCMX+Kerker at conv=1e-6");

    let c = &result.components;
    let e_band_from_identity =
        c.e_kinetic + c.e_local + c.e_nonlocal + 2.0 * c.e_hartree + c.e_vxc;
    let residual = (e_band_from_identity - c.e_band).abs();

    eprintln!(
        "MADOC band-sum identity (Fe BCC, nspin=2): E_band={:.10} eV  \
         E_kin+E_loc+E_nl+2·E_H+E_vxc={e_band_from_identity:.10} eV  \
         residual={residual:.3e} eV  ({} iters, M={:.4} μB)",
        c.e_band, result.n_iterations, result.magnetization,
    );

    // Tolerance rationale: at conv_threshold = 1e-6 the identity residual
    // on Fe BCC is dominated by the O(Δρ) systematic drift between
    // H[ρ_in] (used by the eigenvalue equation, hence by E_band) and
    // V_H[ρ_out] + V_xc[ρ_out] (the quantities stored in `components`).
    // Because the Hartree coupling on Fe is ~360 eV (|E_H|) and the XC
    // coupling is ~400 eV (|E_vxc|), even a Δρ_RMS of 1e-6 leaves an
    // identity residual of order Σ_G ρ(G)·4π·Δρ(G)/|G|² ~ O(1e-2) eV.
    // We cannot close the identity tighter without converging the SCF
    // itself tighter — which is not feasible on this PP/ecut inside a
    // CI time budget.
    //
    // 5e-2 eV keeps ~2× empirical headroom over the observed ~2.3e-2 eV
    // and still catches the major bug classes this test exists to pin:
    //   - factor-2 in E_H: error = |E_H| ≈ 360 eV → 7·10³× the tolerance.
    //   - sign-flipped spin-channel vxc: error = 2·|E_vxc_σ| ≈ 400 eV
    //     → 8·10³× the tolerance.
    //   - kernel swap (V_xc_up/down confused): error depends on |ρ↑-ρ↓|
    //     × V_xc spread; conservatively ≳ 0.5 eV on Fe.
    //
    // For an order-of-magnitude tighter regression guard on this class,
    // convert Fe to a converged external-dataset fixture (cache a
    // reference `EnergyComponents`) rather than running SCF inside the
    // test. Tracked as a follow-up once a regression in this space
    // actually surfaces.
    let tol = 5e-2;
    assert!(
        residual < tol,
        "MADOC band-sum identity violated on Fe BCC (nspin=2):\n  \
         E_band                                = {:.10} eV\n  \
         E_kin + E_loc + E_nl + 2·E_H + E_vxc  = {e_band_from_identity:.10} eV\n  \
         residual                              = {residual:.3e} eV  (tol {tol:.1e})\n\
         This identity is required by the LSDA Kohn-Sham eigenvalue\n\
         equation summed over both spin channels; a nonzero residual\n\
         indicates a factor-2 bug in E_hartree, a sign/scale error in\n\
         one of the spin vxc channels, or a missing term in the\n\
         spin-driver Hamiltonian assembly.",
        c.e_band,
    );
}
