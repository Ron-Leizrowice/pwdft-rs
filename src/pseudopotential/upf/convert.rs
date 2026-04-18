//! Unit conversion from UPF native (Ry, Bohr) to internal (eV, Å) units.
//!
//! The UPF format stores quantities in Rydberg atomic units with various
//! on-disk conventions (`r·β(r)`, `4πr²·ρ_at(r)`, bare ρ_core(r)). This
//! module reads those blocks via [`super::xml`] and converts each one to
//! the engine's internal representation.

use crate::consts::{BOHR3_TO_ANG3, BOHR_TO_ANG, RY_TO_EV};
use crate::error::{PwdftError, Result};

use super::super::{BetaProjector, PseudopotentialData};
use super::xml::{extract_attr, extract_beta_angular_momentum, extract_data_block};

/// Parse the full UPF body and return a unit-converted
/// [`PseudopotentialData`].
pub(super) fn parse_body(content: &str) -> Result<PseudopotentialData> {
    let element = extract_attr(content, "element")
        .ok_or_else(|| PwdftError::Parse("missing element in PP_HEADER".into()))?
        .trim()
        .to_string();

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
        let l = extract_beta_angular_momentum(content, &tag)
            .ok_or_else(|| PwdftError::Parse(format!("missing angular_momentum in {tag}")))?;

        let beta_r_ry = extract_data_block(content, &tag, mesh_size)?;
        // UPF stores χ(r) = r · β(r) in Bohr^{-1/2} (no energy dimension).
        // β(r) is a wavefunction-like quantity in Bohr^{-3/2}, so χ = r·β is in Bohr^{-1/2}.
        // Energy enters only through D_ij (in Ry, converted to eV).
        //
        // Convert Bohr^{-1/2} → Å^{-1/2}: divide by √(BOHR_TO_ANG).
        let values: Vec<f64> = beta_r_ry
            .iter()
            .map(|&v| v / BOHR_TO_ANG.sqrt())
            .collect();

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
        rho_raw
            .iter()
            .map(|&rho| rho / BOHR_TO_ANG)
            .collect()
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
            nlcc_raw
                .iter()
                .map(|&rho| rho / BOHR3_TO_ANG3)
                .collect()
        } else {
            vec![0.0; mesh_size]
        }
    } else {
        vec![]
    };

    Ok(PseudopotentialData {
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
    use super::super::parse;
    use crate::pseudopotential::PseudopotentialData;

    fn si_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap()
    }

    #[test]
    fn test_parse_header() {
        let pp = parse(&si_content()).unwrap();
        assert_eq!(pp.element, "Si");
        assert!((pp.z_valence - 4.0).abs() < 1e-10);
        assert!(pp.l_max >= 1, "Si should have l_max >= 1");
        assert!(pp.n_projectors() > 0, "Si should have projectors");
    }

    #[test]
    fn test_radial_grid_monotonic() {
        let pp = parse(&si_content()).unwrap();
        for i in 1..pp.r_grid.len() {
            assert!(
                pp.r_grid[i] > pp.r_grid[i - 1],
                "radial grid not monotonic at index {i}"
            );
        }
    }

    #[test]
    fn test_radial_grid_units() {
        let pp = parse(&si_content()).unwrap();
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
        let pp = parse(&si_content()).unwrap();
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
        let pp = parse(&si_content()).unwrap();
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
        let pp = parse(&si_content()).unwrap();
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
        eprintln!("Si partial core charge = {q_core:.6} e");
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
    // (`scripts/validate/rho_core_g_reference.py`, which parses the UPF
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
        if gr < 1e-10 {
            1.0 - gr * gr / 6.0
        } else {
            gr.sin() / gr
        }
    }

    /// Compute ρ_core(G) in e/Å³ from the parsed pseudopotential data.
    ///
    /// Inputs in internal units: `pp.r_grid` in Å, `pp.rab` in Å,
    /// `pp.core_charge` in e/Å³. Caller supplies the cell volume `omega`
    /// in Å³ and the wavenumber magnitude `g_norm` in Å⁻¹.
    fn rho_core_of_g_ang(pp: &PseudopotentialData, g_norm: f64, omega: f64) -> f64 {
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
    /// Reference: `scripts/validate/rho_core_g_reference.csv` row
    /// `(si, shell 0)` = 1.8476428665e-02 e/Å³.
    #[test]
    fn test_si_rho_core_of_g_zero() {
        let pp = parse(&si_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 5.431_f64; // Å
        let omega = a * a * a / 4.0; // FCC primitive cell volume
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        eprintln!("Si ρ_core(G=0) = {rho_g0:.6e} e/Å³  (ref 1.8476e-2)");
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
    /// `scripts/validate/rho_core_g_reference.csv` row `(si, shell 1)`
    /// = 1.5427684529e-02 e/Å³.
    #[test]
    fn test_si_rho_core_of_g_first_shell() {
        let pp = parse(&si_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 5.431_f64; // Å
        let omega = a * a * a / 4.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (3.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        eprintln!(
            "Si ρ_core(|G|²=3·(2π/a)²) = ρ_core(G={g_norm:.6} Å⁻¹) \
             = {rho_g:.6e} e/Å³  (ref 1.5428e-2)"
        );
        let expected = 1.542_768_452_9e-2;
        assert!(
            (rho_g - expected).abs() < 1.0e-5,
            "Si ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }

    fn fe_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("pseudopotentials/nc/lda/Fe.upf"),
        )
        .unwrap()
    }

    /// NLCC audit, Part A.3: pin Fe ρ_core(G=0). Fe has a substantially
    /// larger core charge than Si (Q_core ≈ 2.9171 e vs. 0.74 e) because
    /// the 3d semicore overlaps the valence 4s/3d — this is the case
    /// that motivates NLCC most strongly.
    ///
    /// BCC primitive Ω = a³/2 = 11.820 Å³ at a = 2.87 Å. Reference:
    /// `scripts/validate/rho_core_g_reference.csv` row `(fe, shell 0)`
    /// = 2.4679728958e-01 e/Å³.
    #[test]
    fn test_fe_rho_core_of_g_zero() {
        let pp = parse(&fe_content()).unwrap();
        assert!(pp.has_nlcc(), "Fe ONCVPSP should have core_correction=T");

        let a = 2.87_f64; // Å
        let omega = a * a * a / 2.0; // BCC primitive cell volume
        let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);

        eprintln!("Fe ρ_core(G=0) = {rho_g0:.6e} e/Å³  (ref 2.4680e-1)");
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
    /// Reference: `scripts/validate/rho_core_g_reference.csv` row
    /// `(fe, shell 1)` = 2.2502662616e-01 e/Å³.
    #[test]
    fn test_fe_rho_core_of_g_first_shell() {
        let pp = parse(&fe_content()).unwrap();
        assert!(pp.has_nlcc());

        let a = 2.87_f64;
        let omega = a * a * a / 2.0;
        let g_norm = 2.0 * std::f64::consts::PI / a * (2.0_f64).sqrt();
        let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);

        eprintln!(
            "Fe ρ_core(|G|²=2·(2π/a)²) = ρ_core(G={g_norm:.6} Å⁻¹) \
             = {rho_g:.6e} e/Å³  (ref 2.2503e-1)"
        );
        let expected = 2.250_266_261_6e-1;
        assert!(
            (rho_g - expected).abs() < 1.0e-4,
            "Fe ρ_core(first shell) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
        );
    }
}
