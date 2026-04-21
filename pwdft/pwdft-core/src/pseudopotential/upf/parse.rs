//! Unit conversion from UPF native (Ry, Bohr) to internal (eV, Å) units.
//!
//! The UPF format stores quantities in Rydberg atomic units with various
//! on-disk conventions (`r·β(r)`, `4πr²·ρ_at(r)`, bare ρ_core(r)). This
//! module reads those blocks via [`super::xml`] and converts each one to
//! the engine's internal representation.

use std::str::FromStr;

use elements_rs::Element;

use super::xml::{extract_attr, extract_beta_angular_momentum, extract_data_block};
use crate::{
    consts::{BOHR_TO_ANG, BOHR3_TO_ANG3, RY_TO_EV},
    error::{PwdftError, Result},
    pseudopotential::{BetaProjector, UpfPseudoPotential},
};

/// Parse the full UPF body and return a unit-converted UpfPseudoPotential
pub(super) fn parse_upf_body(content: &str) -> Result<UpfPseudoPotential> {
    let element_str = extract_attr(content, "element")
        .ok_or_else(|| PwdftError::Parse("missing element in PP_HEADER".into()))?
        .trim();
    let element = Element::from_str(element_str).map_err(|_| PwdftError::UnknownElement {
        symbol: element_str.to_owned(),
    })?;

    let z_valence: f64 = extract_attr(content, "z_valence")
        .ok_or_else(|| PwdftError::Parse("missing z_valence".into()))?
        .trim()
        .parse()
        .map_err(|e| PwdftError::Parse(format!("z_valence parse error: {e}")))?;

    let l_max: i32 = extract_attr(content, "l_max")
        .ok_or_else(|| PwdftError::Parse("missing l_max".into()))?
        .trim()
        .parse()
        .map_err(|e| PwdftError::Parse(format!("l_max parse error: {e}")))?;

    let n_proj: usize = extract_attr(content, "number_of_proj")
        .ok_or_else(|| PwdftError::Parse("missing number_of_proj".into()))?
        .trim()
        .parse()
        .map_err(|e| PwdftError::Parse(format!("number_of_proj parse error: {e}")))?;

    let mesh_size: usize = extract_attr(content, "mesh_size")
        .ok_or_else(|| PwdftError::Parse("missing mesh_size".into()))?
        .trim()
        .parse()
        .map_err(|e| PwdftError::Parse(format!("mesh_size parse error: {e}")))?;

    // Parse radial grid (Bohr → Å)
    let r_grid_bohr = extract_data_block(content, "PP_R", mesh_size)?;
    let r_grid: Vec<f64> = r_grid_bohr.iter().map(|&r| r * BOHR_TO_ANG).collect();

    // Parse integration weights (Bohr → Å)
    let rab_bohr = extract_data_block(content, "PP_RAB", mesh_size)?;
    let rab: Vec<f64> = rab_bohr.iter().map(|&dr| dr * BOHR_TO_ANG).collect();

    // Parse local potential (Ry → eV)
    let v_local_ry = extract_data_block(content, "PP_LOCAL", mesh_size)?;
    let v_local: Vec<f64> = v_local_ry.iter().map(|&v| v * RY_TO_EV).collect();

    // Parse non-local projectors
    // UPF stores β(r) · r, in Ry^{1/2} units, on the radial grid.
    // We need β(r) in our units.
    let mut beta_projectors = Vec::with_capacity(n_proj);
    for i in 1..=n_proj {
        let tag = format!("PP_BETA.{i}");
        let l = extract_beta_angular_momentum(content, &tag)?;

        let beta_r_ry = extract_data_block(content, &tag, mesh_size)?;
        // UPF stores χ(r) = r · β(r) in Bohr^{-1/2} (no energy dimension).
        // β(r) is a wavefunction-like quantity in Bohr^{-3/2}, so χ = r·β is in
        // Bohr^{-1/2}. Energy enters only through D_ij (in Ry, converted to
        // eV).
        //
        // Convert Bohr^{-1/2} → Å^{-1/2}: divide by √(BOHR_TO_ANG).
        let values: Vec<f64> = beta_r_ry.iter().map(|&v| v / BOHR_TO_ANG.sqrt()).collect();

        beta_projectors.push(BetaProjector { l, values });
    }

    // Parse D_ij matrix (Ry → eV)
    let dij_ry = extract_data_block(content, "PP_DIJ", n_proj * n_proj)?;
    let dij: Vec<f64> = dij_ry.iter().map(|&d| d * RY_TO_EV).collect();

    // Parse atomic charge density (optional — may be all zeros for HGH)
    // UPF stores 4πr²ρ(r) in e/Bohr on the radial grid.
    // The Bessel transform ∫ [4πr²ρ(r)] j₀(Gr) dr needs consistent units:
    // r_grid is in Å, rab is in Å, G is in 1/Å, so 4πr²ρ(r) must be in e/Å.
    // Convert e/Bohr → e/Å by dividing by BOHR_TO_ANG.
    let rho_atom = if let Ok(rho_raw) = extract_data_block(content, "PP_RHOATOM", mesh_size) {
        rho_raw.iter().map(|&rho| rho / BOHR_TO_ANG).collect()
    } else {
        vec![0.0; mesh_size]
    };

    // Parse nonlinear core correction (NLCC) charge density.
    //
    // PP_NLCC stores the **bare** core density ρ_core(r) in e/Bohr³ (NOT
    // 4πr²·ρ). This differs from PP_RHOATOM which is stored as 4πr²·ρ_at
    // (see the rho_atom block above). QE confirms this convention at
    // `qe-7.5/upflib/rhoc_mod.f90:107`:
    //
    //     aux(ir) = upf%rho_atc(ir) * rgrid%r2(ir) * sin(q r)/(q r)
    //
    // — QE multiplies by r² (and later by 4π/Ω) in its Bessel transform,
    // proving rho_atc is the bare volumetric density.
    //
    // Convert e/Bohr³ → e/Å³ (volumetric, not linear) by dividing by
    // BOHR_TO_ANG³. The downstream Bessel transform in
    // `src/scf/potentials.rs::compute_core_density` supplies the r²·4π
    // factors in Å units.
    let has_nlcc = content.contains("core_correction=\"T\"")
        || content.contains("core_correction=\"t\"")
        || content.contains("nlcc=.true.");
    let core_charge = if has_nlcc {
        if let Ok(nlcc_raw) = extract_data_block(content, "PP_NLCC", mesh_size) {
            nlcc_raw.iter().map(|&rho| rho / BOHR3_TO_ANG3).collect()
        } else {
            vec![0.0; mesh_size]
        }
    } else {
        vec![]
    };

    Ok(UpfPseudoPotential {
        element,
        z_valence,
        l_max,
        r_grid,
        rab,
        v_local,
        beta_projectors,
        dij,
        rho_atom,
        core_charge,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn si_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap()
    }

    #[test]
    fn test_parse_header() {
        let pp = parse_upf_body(&si_content()).unwrap();
        assert_eq!(pp.element, Element::Si);
        assert!((pp.z_valence - 4.0).abs() < 1e-10);
        assert!(pp.l_max >= 1, "Si should have l_max >= 1");
        assert!(pp.n_projectors() > 0, "Si should have projectors");
    }

    #[test]
    fn test_radial_grid_monotonic() {
        let pp = parse_upf_body(&si_content()).unwrap();
        for i in 1..pp.r_grid.len() {
            assert!(
                pp.r_grid[i] > pp.r_grid[i - 1],
                "radial grid not monotonic at index {i}"
            );
        }
    }

    #[test]
    fn test_radial_grid_units() {
        let pp = parse_upf_body(&si_content()).unwrap();
        // First grid point should be small (< 0.01 Å)
        assert!(pp.r_grid[0] < 0.01, "first r = {} Å too large", pp.r_grid[0]);
        // Last grid point should be reasonable (> 1 Å)
        assert!(
            pp.r_grid.last().unwrap() > &1.0,
            "last r = {} Å too small",
            pp.r_grid.last().unwrap()
        );
    }

    #[test]
    fn test_v_local_coulomb_tail() {
        let pp = parse_upf_body(&si_content()).unwrap();
        // At large r, V_local should approach -Z_val e²/r
        // e² = 14.3997 eV·Å, Z_val = 4
        let e2 = crate::consts::E2_COULOMB;
        let n = pp.r_grid.len();
        let r_far = pp.r_grid[n - 10];
        let v_far = pp.v_local[n - 10];
        let v_coulomb = -pp.z_valence * e2 / r_far;
        let relative_err = ((v_far - v_coulomb) / v_coulomb).abs();
        assert!(
            relative_err < 0.01,
            "V_local at r={r_far:.2} Å: {v_far:.6} eV vs Coulomb {v_coulomb:.6} eV (err={relative_err:.4})"
        );
    }

    #[test]
    fn test_dij_matrix_size() {
        let pp = parse_upf_body(&si_content()).unwrap();
        assert_eq!(pp.dij.len(), pp.n_projectors() * pp.n_projectors());
    }

    /// NCFX: verify the partial core charge integrates to a physically
    /// reasonable value after the corrected unit conversion.
    ///
    /// `pp.core_charge` now stores the bare ρ_core(r) in e/Å³. The
    /// integrated partial core charge is
    ///     Q_core = ∫ 4π r² ρ_core(r) dr
    /// with r and rab in Å. For Si ONCVPSP LDA this is ≈ 0.74 e — a
    /// reasonable partial-core charge for a Si atom (confirmed against
    /// independent Python parse of Si.upf PP_NLCC).
    ///
    /// The previous (pre-NCFX) erroneous `/BOHR_TO_ANG` conversion would
    /// give ≈ 0.74·BOHR_TO_ANG² ≈ 0.207 e — implausibly small.
    #[test]
    fn test_si_core_charge_integrates_to_partial_core() {
        let pp = parse_upf_body(&si_content()).unwrap();
        assert!(pp.core_charge.len() == pp.r_grid.len());
        assert!(pp.has_nlcc(), "Si ONCVPSP should have core_correction=T");

        let four_pi = 4.0 * std::f64::consts::PI;
        let q_core: f64 = pp
            .core_charge
            .iter()
            .zip(pp.r_grid.iter())
            .zip(pp.rab.iter())
            .map(|((&rho, &r), &dr)| four_pi * r * r * rho * dr)
            .sum();
        // Expected ~0.74 e for Si ONCVPSP LDA. Allow 0.01 e tolerance for
        // trapezoidal-rule round-off on the log mesh.
        assert!(
            (q_core - 0.74).abs() < 0.05,
            "Si partial core charge = {q_core:.6} e (expected ≈ 0.74 e)"
        );
    }

    // -----------------------------------------------------------------
    // NLCC audit (Part A) — pin ρ_core(G) at G=0 and the first non-zero
    // G shell for Si and Fe. Reference values were computed by an
    // independent Python implementation
    // (`pwdft-validate reference nlcc`, which parses the UPF
    // with regex and integrates with `scipy.integrate.simpson` in QE
    // native units, then converts e/Bohr³ → e/Å³). The Rust integrator
    // below trapezoidal-sums the same formula on the *Å-unit*
    // `pp.core_charge`, `pp.r_grid`, `pp.rab` — so the pins also
    // exercise the post-NCFX unit conversion (`e/Bohr³ → e/Å³` via
    // `/BOHR3_TO_ANG3` in the parser) and the r²·4π weighting used by
    // `src/scf/potentials.rs::compute_core_density`.
    //
    // The Bessel formula:
    //     ρ_core(G) = (4π / Ω) · ∫ ρ_core(r) · r² · j₀(|G|r) dr
    // with j₀(x) = sin(x)/x (series form at x→0). The Rust helper uses
    // the trapezoidal rule (instead of Simpson) to keep it dependency-
    // free; the residual vs. the Python/Simpson reference is O(1e-5
    // e/Å³) on the ONCVPSP log mesh, which sets the pin tolerance.

    fn bessel_j0(gr: f64) -> f64 {
        if gr < 1e-10 { 1.0 - gr * gr / 6.0 } else { gr.sin() / gr }
    }

    /// Compute ρ_core(G) in e/Å³ from the parsed pseudopotential data.
    ///
    /// Inputs in internal units: `pp.r_grid` in Å, `pp.rab` in Å,
    /// `pp.core_charge` in e/Å³. Caller supplies the cell volume `omega`
    /// in Å³ and the wavenumber magnitude `g_norm` in Å⁻¹.
    fn rho_core_of_g_ang(pp: &UpfPseudoPotential, g_norm: f64, omega: f64) -> f64 {
        let four_pi = 4.0 * std::f64::consts::PI;
        let integral: f64 = pp
            .core_charge
            .iter()
            .zip(pp.r_grid.iter())
            .zip(pp.rab.iter())
            .map(|((&rho, &r), &dr)| {
                let j0 = bessel_j0(g_norm * r);
                rho * r * r * j0 * dr
            })
            .sum();
        four_pi * integral / omega
    }

    /// NLCC audit, Part A.1: pin Si ρ_core(G=0) to the known value.
    ///
    /// ρ_core(G=0) = Q_core / Ω (since j₀(0) = 1). For Si ONCVPSP LDA
    /// (a = 5.431 Å, FCC: Ω = a³/4 = 40.032 Å³) with Q_core ≈ 0.7399 e,
    /// we expect ρ_core(0) ≈ 1.8476·10⁻² e/Å³.
    ///
    /// Reference: `data/csv/rho_core_g_reference.csv` row
    /// `(si, shell 0)` = 1.8476428665e-02 e/Å³.
    #[test]
    fn test_si_rho_core_of_g_zero() {
        let pp = parse_upf_body(&si_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 5.431_f64; // Å
        let omega = a * a * a / 4.0; // FCC primitive cell volume
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        let expected = 1.847_642_866_5e-2;
        // Tolerance 1e-5 covers trapezoidal-vs-Simpson difference on the
        // ONCVPSP mesh (dr ≈ 5e-3 Å in the relevant region).
        assert!(
            (rho_g0 - expected).abs() < 1.0e-5,
            "Si ρ_core(G=0) = {rho_g0:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    /// NLCC audit, Part A.2: pin Si ρ_core(G≠0) at |G|² = 3·(2π/a)²
    /// — the first non-zero FCC shell (the {111} family). In Å⁻¹:
    ///     |G| = 2π/a · √3 ≈ 2.003873 Å⁻¹.
    ///
    /// Chosen G is the smallest non-zero |G| for Si FCC; it exercises
    /// the full Bessel-transform integrand (not just j₀ = 1). Reference:
    /// `data/csv/rho_core_g_reference.csv` row `(si, shell 1)`
    /// = 1.5427684529e-02 e/Å³.
    #[test]
    fn test_si_rho_core_of_g_first_shell() {
        let pp = parse_upf_body(&si_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 5.431_f64; // Å
        let omega = a * a * a / 4.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (3.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        let expected = 1.542_768_452_9e-2;
        assert!(
            (rho_g - expected).abs() < 1.0e-5,
            "Si ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    fn fe_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Fe.upf"),
        )
        .unwrap()
    }

    /// NLCC audit, Part A.3: pin Fe ρ_core(G=0). Fe has a substantially
    /// larger core charge than Si (Q_core ≈ 2.9171 e vs. 0.74 e) because
    /// the 3d semicore overlaps the valence 4s/3d — this is the case
    /// that motivates NLCC most strongly.
    ///
    /// BCC primitive Ω = a³/2 = 11.820 Å³ at a = 2.87 Å. Reference:
    /// `data/csv/rho_core_g_reference.csv` row `(fe, shell 0)`
    /// = 2.4679728958e-01 e/Å³.
    #[test]
    fn test_fe_rho_core_of_g_zero() {
        let pp = parse_upf_body(&fe_content()).unwrap();
        assert!(pp.has_nlcc(), "Fe ONCVPSP should have core_correction=T");

        let a = 2.87_f64; // Å
        let omega = a * a * a / 2.0; // BCC primitive cell volume
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        let expected = 2.467_972_895_8e-1;
        // Fe NLCC is ~13× larger than Si's at G=0; tolerance 1e-4 e/Å³
        // ≈ 4·10⁻⁴ relative, consistent with the trapezoidal-vs-Simpson
        // residual scaled by magnitude.
        assert!(
            (rho_g0 - expected).abs() < 1.0e-4,
            "Fe ρ_core(G=0) = {rho_g0:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    /// NLCC audit, Part A.4: pin Fe ρ_core(G≠0) at |G|² = 2·(2π/a)²
    /// — the first non-zero BCC shell (the {110} family). In Å⁻¹:
    ///     |G| = 2π/a · √2 ≈ 3.096355 Å⁻¹.
    ///
    /// Reference: `data/csv/rho_core_g_reference.csv` row
    /// `(fe, shell 1)` = 2.2502662616e-01 e/Å³.
    #[test]
    fn test_fe_rho_core_of_g_first_shell() {
        let pp = parse_upf_body(&fe_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 2.87_f64;
        let omega = a * a * a / 2.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (2.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        let expected = 2.250_266_261_6e-1;
        assert!(
            (rho_g - expected).abs() < 1.0e-4,
            "Fe ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    // -----------------------------------------------------------------
    // TRV2 Finding #3 — Extend NLCC ρ_core(G) regression coverage beyond
    // Si/Fe.  NCFX (PR #40) fixed a universal unit/radial-weight bug,
    // but the original test coverage only pinned Si and Fe; the other
    // NLCC-active pseudopotentials in `pseudopotentials/nc/lda/` were
    // silent to future regressions of the same bug class on a single
    // element's mesh.  Cu exercises the 3s/3p/3d semicore edge case
    // (Z_val=19, larger ρ_core than Fe), and Mn extends to the magnetic
    // reference (Z_val=15, largest Q_core of the four pinned elements).
    // The remaining NLCC-active pseudos are flagged in FLUP.

    fn cu_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Cu.upf"),
        )
        .unwrap()
    }

    /// NLCC audit, Part A.5: pin Cu ρ_core(G=0). Cu ONCVPSP has the full
    /// 3s/3p/3d semicore in the valence (Z_val = 19), Q_core ≈ 3.01 e.
    ///
    /// Cell: FCC a = 6.8219 Bohr = 3.610 Å (matches
    /// `data/qe/cu_fcc_scf.in`).
    /// Reference: `data/csv/rho_core_g_reference.csv` row
    /// `(cu, shell 0)` = 2.5568474262e-01 e/Å³.
    #[test]
    fn test_cu_rho_core_of_g_zero() {
        let pp = parse_upf_body(&cu_content()).unwrap();
        assert!(pp.has_nlcc(), "Cu ONCVPSP should have core_correction=T");

        let a = 6.8219_f64 * crate::consts::BOHR_TO_ANG; // 3.6100 Å
        let omega = a * a * a / 4.0; // FCC primitive cell volume
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        let expected = 2.556_847_426_2e-1;
        // Cu NLCC magnitude ≈ Fe's, so use the same 1e-4 e/Å³ tolerance
        // (≈ 4·10⁻⁴ relative).
        assert!(
            (rho_g0 - expected).abs() < 1.0e-4,
            "Cu ρ_core(G=0) = {rho_g0:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    /// NLCC audit, Part A.6: pin Cu ρ_core(G≠0) at |G|² = 3·(2π/a)²
    /// — the first non-zero FCC shell (the {111} family). In Å⁻¹:
    ///     |G| = 2π/a · √3 ≈ 3.014181 Å⁻¹.
    ///
    /// Reference: `data/csv/rho_core_g_reference.csv` row
    /// `(cu, shell 1)` = 2.3967132540e-01 e/Å³.
    #[test]
    fn test_cu_rho_core_of_g_first_shell() {
        let pp = parse_upf_body(&cu_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 6.8219_f64 * crate::consts::BOHR_TO_ANG;
        let omega = a * a * a / 4.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (3.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        let expected = 2.396_713_254_0e-1;
        assert!(
            (rho_g - expected).abs() < 1.0e-4,
            "Cu ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    fn mn_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Mn.upf"),
        )
        .unwrap()
    }

    /// NLCC audit, Part A.7: pin Mn ρ_core(G=0). Mn ONCVPSP has
    /// Z_val = 15 (3s²3p⁶3d⁵4s²) with a substantial partial core charge
    /// (Q_core ≈ 4.19 e — the largest of the four pinned elements).
    ///
    /// Cell choice: α-Mn has a complex 58-atom cubic ground state, but
    /// for NLCC regression only the cell volume matters. We use a
    /// simple BCC container at a = 2.89 Å (close to Fe) so Mn's
    /// ρ_core(G) appears in the same |G|-shell range as Fe's.
    /// Reference: `data/csv/rho_core_g_reference.csv` row
    /// `(mn, shell 0)` = 3.4733264804e-01 e/Å³.
    #[test]
    fn test_mn_rho_core_of_g_zero() {
        let pp = parse_upf_body(&mn_content()).unwrap();
        assert!(pp.has_nlcc(), "Mn ONCVPSP should have core_correction=T");

        let a = 2.89_f64; // Å
        let omega = a * a * a / 2.0; // BCC primitive cell volume
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        let expected = 3.473_326_480_4e-1;
        // Mn's NLCC is ≈ 1.4× Fe's at G=0; keep the same 1e-4 tolerance
        // — still ≈ 3·10⁻⁴ relative on the ONCVPSP log mesh.
        assert!(
            (rho_g0 - expected).abs() < 1.0e-4,
            "Mn ρ_core(G=0) = {rho_g0:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    /// NLCC audit, Part A.8: pin Mn ρ_core(G≠0) at |G|² = 2·(2π/a)²
    /// — the first non-zero BCC shell (the {110} family). In Å⁻¹:
    ///     |G| = 2π/a · √2 ≈ 3.074921 Å⁻¹.
    ///
    /// Reference: `data/csv/rho_core_g_reference.csv` row
    /// `(mn, shell 1)` = 3.1210083482e-01 e/Å³.
    #[test]
    fn test_mn_rho_core_of_g_first_shell() {
        let pp = parse_upf_body(&mn_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 2.89_f64;
        let omega = a * a * a / 2.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (2.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        let expected = 3.121_008_348_2e-1;
        assert!(
            (rho_g - expected).abs() < 1.0e-4,
            "Mn ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    // -----------------------------------------------------------------
    // VGCH-2F Part C session-2 — H-C4 probe: close the remaining NLCC
    // ρ_core(G) coverage gap.  Cu was pinned by TRV2; VGCH-2 Part B
    // identified a shared-density ΔE_xc = +8.87 eV on the Cu transplant
    // whose fingerprint is consistent with a silent ρ_core unit- or
    // mesh-conversion bug on a heavy element.  The Class A NLCC
    // contributors without a pin pre-VGCH-2F were Ga and As (GaAs
    // zinc-blende), O (MgO rocksalt), and Cl (NaCl rocksalt).  After
    // these four pins every NLCC-active PP that enters a Class A QE
    // reference cell is pinned at G=0 and the first non-zero G-shell.
    //
    // Mg and Na have `core_correction="F"` → no PP_NLCC block → no pin
    // needed.  The Class A NLCC surface is closed.
    //
    // Reference values: `scripts/validate/rho_core_g_reference.csv`,
    // regenerated from `scripts/validate/rho_core_g_reference.py` with
    // cells matching `qe_validation/{gaas,mgo,nacl}_scf.in` (ibrav=2,
    // FCC primitive Ω = a³/4).

    fn ga_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Ga.upf"),
        )
        .unwrap()
    }

    /// VGCH-2F A.9 — pin Ga ρ_core(G=0).  Ga ONCVPSP has the full
    /// 3d semicore in the valence (Z_val = 13, 3d¹⁰4s²4p¹) with
    /// Q_core ≈ 7.94 e — the largest of the Class A pinned elements.
    ///
    /// Cell: GaAs zinc-blende, a = 10.6829 Bohr = 5.6530 Å (matches
    /// `qe_validation/gaas_scf.in`).  Ω = a³/4 Å³.
    /// Reference: `rho_core_g_reference.csv` row `(ga, shell 0)`
    /// = 1.7582117492e-01 e/Å³.
    #[test]
    fn test_ga_rho_core_of_g_zero() {
        let pp = parse_upf_body(&ga_content()).unwrap();
        assert!(pp.has_nlcc(), "Ga ONCVPSP should have core_correction=T");

        let a = 10.6829_f64 * crate::consts::BOHR_TO_ANG; // 5.6530 Å
        let omega = a * a * a / 4.0;
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        let expected = 1.758_211_749_2e-1;
        // Ga's NLCC is ≈ Fe's at G=0 in absolute magnitude; 1e-4
        // tolerance ≈ 6·10⁻⁴ relative on the ONCVPSP log mesh.
        assert!(
            (rho_g0 - expected).abs() < 1.0e-4,
            "Ga ρ_core(G=0) = {rho_g0:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    /// VGCH-2F A.10 — pin Ga ρ_core(G≠0) at |G|² = 3·(2π/a)² — the
    /// first non-zero FCC shell (the {111} family).  In Å⁻¹:
    ///     |G| = 2π/a · √3 ≈ 1.9254 Å⁻¹.
    ///
    /// Reference: `rho_core_g_reference.csv` row `(ga, shell 1)`
    /// = 1.6454252388e-01 e/Å³.
    #[test]
    fn test_ga_rho_core_of_g_first_shell() {
        let pp = parse_upf_body(&ga_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 10.6829_f64 * crate::consts::BOHR_TO_ANG;
        let omega = a * a * a / 4.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (3.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        let expected = 1.645_425_238_8e-1;
        assert!(
            (rho_g - expected).abs() < 1.0e-4,
            "Ga ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    fn as_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/As.upf"),
        )
        .unwrap()
    }

    /// VGCH-2F A.11 — pin As ρ_core(G=0).  As ONCVPSP (Z_val = 5,
    /// 4s²4p³) + 3d semicore → Q_core ≈ 7.99 e.
    ///
    /// Cell: GaAs zinc-blende, a = 5.6530 Å (same as Ga).
    /// Reference: `rho_core_g_reference.csv` row `(as, shell 0)`
    /// = 1.7690334999e-01 e/Å³.
    #[test]
    fn test_as_rho_core_of_g_zero() {
        let pp = parse_upf_body(&as_content()).unwrap();
        assert!(pp.has_nlcc(), "As ONCVPSP should have core_correction=T");

        let a = 10.6829_f64 * crate::consts::BOHR_TO_ANG;
        let omega = a * a * a / 4.0;
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        let expected = 1.769_033_499_9e-1;
        assert!(
            (rho_g0 - expected).abs() < 1.0e-4,
            "As ρ_core(G=0) = {rho_g0:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    /// VGCH-2F A.12 — pin As ρ_core(G≠0) at |G|² = 3·(2π/a)².
    ///
    /// Reference: `rho_core_g_reference.csv` row `(as, shell 1)`
    /// = 1.6708351731e-01 e/Å³.
    #[test]
    fn test_as_rho_core_of_g_first_shell() {
        let pp = parse_upf_body(&as_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 10.6829_f64 * crate::consts::BOHR_TO_ANG;
        let omega = a * a * a / 4.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (3.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        let expected = 1.670_835_173_1e-1;
        assert!(
            (rho_g - expected).abs() < 1.0e-4,
            "As ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    fn o_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/O.upf"),
        )
        .unwrap()
    }

    /// VGCH-2F A.13 — pin O ρ_core(G=0).  O ONCVPSP (Z_val = 6,
    /// 2s²2p⁴) has a small core (no semicore) → Q_core ≈ 0.574 e.
    /// Smallest ρ_core magnitude of the pinned Class A elements —
    /// sensitive to trapezoidal-vs-Simpson mesh differences, so the
    /// tolerance is tightened to 1e-5 to match Si's regime.
    ///
    /// Cell: MgO rocksalt, a = 7.9586 Bohr = 4.2115 Å (matches
    /// `qe_validation/mgo_scf.in`).  Reference: `rho_core_g_reference.csv`
    /// row `(o, shell 0)` = 3.0760542725e-02 e/Å³.
    #[test]
    fn test_o_rho_core_of_g_zero() {
        let pp = parse_upf_body(&o_content()).unwrap();
        assert!(pp.has_nlcc(), "O ONCVPSP should have core_correction=T");

        let a = 7.9586_f64 * crate::consts::BOHR_TO_ANG; // 4.2115 Å
        let omega = a * a * a / 4.0;
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        let expected = 3.076_054_272_5e-2;
        // O NLCC is ≈ 2× Si's at G=0; keep the same 1e-5 tolerance
        // as Si (O and Si are both on the no-semicore regime).
        assert!(
            (rho_g0 - expected).abs() < 1.0e-5,
            "O ρ_core(G=0) = {rho_g0:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    /// VGCH-2F A.14 — pin O ρ_core(G≠0) at |G|² = 3·(2π/a)².
    ///
    /// Reference: `rho_core_g_reference.csv` row `(o, shell 1)`
    /// = 2.9521130484e-02 e/Å³.
    #[test]
    fn test_o_rho_core_of_g_first_shell() {
        let pp = parse_upf_body(&o_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 7.9586_f64 * crate::consts::BOHR_TO_ANG;
        let omega = a * a * a / 4.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (3.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        let expected = 2.952_113_048_4e-2;
        assert!(
            (rho_g - expected).abs() < 1.0e-5,
            "O ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    fn cl_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Cl.upf"),
        )
        .unwrap()
    }

    /// VGCH-2F A.15 — pin Cl ρ_core(G=0).  Cl ONCVPSP (Z_val = 7,
    /// 3s²3p⁵) → Q_core ≈ 1.872 e.  Intermediate magnitude between O
    /// and Si on one side and the heavier Ga/As/Cu on the other.
    ///
    /// Cell: NaCl rocksalt, a = 10.6078 Bohr = 5.6134 Å (matches
    /// `qe_validation/nacl_scf.in`).  Reference:
    /// `rho_core_g_reference.csv` row `(cl, shell 0)` = 4.2341213379e-02 e/Å³.
    #[test]
    fn test_cl_rho_core_of_g_zero() {
        let pp = parse_upf_body(&cl_content()).unwrap();
        assert!(pp.has_nlcc(), "Cl ONCVPSP should have core_correction=T");

        let a = 10.6078_f64 * crate::consts::BOHR_TO_ANG; // 5.6134 Å
        let omega = a * a * a / 4.0;
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        let expected = 4.234_121_337_9e-2;
        assert!(
            (rho_g0 - expected).abs() < 1.0e-5,
            "Cl ρ_core(G=0) = {rho_g0:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    /// VGCH-2F A.16 — pin Cl ρ_core(G≠0) at |G|² = 3·(2π/a)².
    ///
    /// Reference: `rho_core_g_reference.csv` row `(cl, shell 1)`
    /// = 3.9423292394e-02 e/Å³.
    #[test]
    fn test_cl_rho_core_of_g_first_shell() {
        let pp = parse_upf_body(&cl_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 10.6078_f64 * crate::consts::BOHR_TO_ANG;
        let omega = a * a * a / 4.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (3.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        let expected = 3.942_329_239_4e-2;
        assert!(
            (rho_g - expected).abs() < 1.0e-5,
            "Cl ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    // UPFV: parse-time rejection of negative `angular_momentum`.
    //
    // `extract_beta_angular_momentum` returned `Option<i32>` before UPFV,
    // which happily accepted `-1`, `-2`, … — values that are physically
    // meaningless and would later propagate as an out-of-bounds index
    // into the spherical-harmonic tables inside `NonlocalPotential::new`.
    // The runtime `assert!(proj.l >= 0)` in that constructor (landed with
    // CAST, PR #56) stays as defense-in-depth; these regression tests
    // guard the parser-level early exit with a clearer
    // `PwdftError::InvalidPseudopotential` message (ERR2 P1.d migrated
    // this path from the catch-all `InvalidInput`).
    //
    // Positive integers (0, 1, 2, 3, …) are unaffected — `test_parse_header`
    // and `test_dij_matrix_size` above already exercise that path.

    use super::super::xml::extract_beta_angular_momentum;
    use crate::{error::PwdftError, pseudopotential::upf::parse_upf_body};

    #[test]
    fn test_extract_beta_l_rejects_negative() {
        // Minimal PP_BETA-shaped tag with a negative angular_momentum.
        let content = r#"<PP_BETA.1 type="real" index="1" angular_momentum="-1" >
        1.0 2.0 3.0
        </PP_BETA.1>"#;
        let err = extract_beta_angular_momentum(content, "PP_BETA.1")
            .expect_err("negative angular_momentum must be rejected");
        match err {
            PwdftError::InvalidPseudopotential { file, reason } => {
                assert_eq!(file, "PP_BETA.1");
                assert!(
                    reason.contains("angular_momentum") && reason.contains("-1"),
                    "expected reason to mention angular_momentum and -1, got: {reason}"
                );
            },
            other => panic!("expected PwdftError::InvalidPseudopotential, got {other:?}"),
        }
    }

    #[test]
    fn test_extract_beta_l_accepts_zero_and_positive() {
        for (raw, expected) in [("0", 0_i32), ("1", 1), ("2", 2), ("3", 3)] {
            let content = format!(
                r#"<PP_BETA.1 index="1" angular_momentum="{raw}" >
                </PP_BETA.1>"#
            );
            let l =
                extract_beta_angular_momentum(&content, "PP_BETA.1").expect("non-negative angular_momentum must parse");
            assert_eq!(l, expected);
        }
    }

    /// End-to-end regression: a malformed UPF (a real Si ONCV file with
    /// the PP_BETA.1 angular_momentum flipped from 0 to −1) must fail
    /// `parse_body` with `PwdftError::InvalidPseudopotential`, not a
    /// later internal error from unit conversion or the D_ij block.
    #[test]
    fn test_parse_rejects_negative_angular_momentum_in_upf() {
        let raw = si_content();
        assert!(
            raw.contains(r#"angular_momentum="0""#),
            "test assumes Si.upf has at least one angular_momentum=\"0\" projector"
        );
        let corrupted = raw.replacen(r#"angular_momentum="0""#, r#"angular_momentum="-1""#, 1);
        let err = parse_upf_body(&corrupted).expect_err("corrupted UPF must fail to parse");
        match err {
            PwdftError::InvalidPseudopotential { file, reason } => {
                assert!(
                    file.starts_with("PP_BETA"),
                    "expected file context to name a PP_BETA tag, got: {file}"
                );
                assert!(
                    reason.contains("angular_momentum"),
                    "expected reason mentioning angular_momentum, got: {reason}"
                );
            },
            other => panic!("expected PwdftError::InvalidPseudopotential, got {other:?}"),
        }
    }
}
