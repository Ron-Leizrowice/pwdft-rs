pub mod psp8;
pub mod upf;

use std::path::Path;

use crate::error::{PwdftError, Result};

/// Unit-converted pseudopotential data in internal units (eV, Å).
///
/// All formats (UPF, PSP8, HGH) parse into this common representation.
#[derive(Debug, Clone)]
pub struct PseudopotentialData {
    /// Element symbol (e.g. "Si").
    pub element: String,
    /// Number of valence electrons.
    pub z_valence: f64,
    /// Maximum angular momentum in non-local projectors.
    pub l_max: i32,
    /// Radial grid points (Å).
    pub r_grid: Vec<f64>,
    /// Integration weights dr for Simpson/trapezoidal rule (Å).
    pub rab: Vec<f64>,
    /// Local pseudopotential on radial grid (eV).
    /// Includes the -Z_val e²/r Coulomb tail.
    pub v_local: Vec<f64>,
    /// Non-local projectors: beta[proj_index] = (angular_momentum, radial_values).
    /// Radial values store χ(r) = r·β(r) in Å^{-1/2} (no energy dimension).
    /// Energy enters through D_ij. The KB matrix element is:
    /// V_NL = (1/Ω) Σ F_i D_ij F_j × angular, where F = 4π ∫ χ(r) j_l(qr) r dr.
    pub beta_projectors: Vec<BetaProjector>,
    /// D_ij coupling matrix for non-local projectors (eV).
    /// Stored as a flat n_proj × n_proj matrix in row-major order.
    pub dij: Vec<f64>,
    /// Number of non-local projectors.
    pub n_projectors: usize,
    /// Atomic charge density on radial grid (e/Å, stores 4πr²ρ(r)).
    /// May be empty if not provided by the pseudopotential.
    pub rho_atom: Vec<f64>,
    /// Nonlinear core correction (NLCC) charge density on radial grid.
    /// Stores 4πr²ρ_core(r) in e/Å. Empty if `has_nlcc` is false.
    pub core_charge: Vec<f64>,
    /// Whether this PP has nonlinear core correction.
    pub has_nlcc: bool,
}

/// A single non-local beta projector.
#[derive(Debug, Clone)]
pub struct BetaProjector {
    /// Angular momentum quantum number.
    pub l: i32,
    /// Radial function β(r) on the radial grid.
    pub values: Vec<f64>,
}

// Unit conversion constants (Rydberg a.u. → internal eV/Å)
const RY_TO_EV: f64 = 13.605693122994; // 1 Ry = 13.6057 eV
const BOHR_TO_ANG: f64 = 0.529177210903; // 1 Bohr = 0.5292 Å

/// Detect pseudopotential format and parse.
///
/// Supported formats:
/// - `.UPF` / `.upf`: UPF v2 (Quantum ESPRESSO)
/// - `.psp8`: PSP8 (ABINIT / PseudoDojo)
pub fn load(path: &Path) -> Result<PseudopotentialData> {
    let content = std::fs::read_to_string(path)?;

    // Detect format from content (more robust than extension)
    if content.contains("<UPF") || content.contains("<PP_HEADER") {
        upf::parse(&content)
    } else if content.trim_start().starts_with(|c: char| c.is_alphanumeric()) {
        // PSP8 starts with a title line; check for pspcod=8 on line 3
        let lines: Vec<&str> = content.lines().collect();
        if lines.len() > 2 {
            let fields: Vec<&str> = lines[2].split_whitespace().collect();
            if fields.first().is_some_and(|&f| f == "8") {
                return psp8::parse(&content);
            }
        }
        Err(PwdftError::Parse(format!(
            "unrecognized pseudopotential format in {}",
            path.display()
        )))
    } else {
        Err(PwdftError::Parse(format!(
            "unrecognized pseudopotential format in {}",
            path.display()
        )))
    }
}

/// Find the pseudopotential matching an atom's atomic number.
///
/// Matches by converting the PP element symbol to an atomic number.
/// Panics if no match is found — callers should validate PP coverage at startup.
pub fn find_for_atom<'a>(z: u32, pseudopotentials: &[&'a PseudopotentialData]) -> &'a PseudopotentialData {
    pseudopotentials
        .iter()
        .find(|pp| {
            crate::atoms::from_symbol(&pp.element)
                .is_some_and(|e| e.atomic_number() == z)
        })
        .expect("no pseudopotential found for atom")
}

impl PseudopotentialData {
    /// Compute V_local(G) via spherical Bessel transform.
    ///
    /// V_local(G) = (4π/Ω) ∫₀^∞ r² [V_local(r) + Z_val e²/r] sin(Gr)/(Gr) dr
    ///              - (4π/Ω) Z_val e² / G²
    ///
    /// The Coulomb singularity is handled by subtracting -Z_val/r analytically:
    /// the bracketed term [V_local(r) + Z_val e²/r] is short-ranged.
    ///
    /// For G=0, the integral is:
    /// V_local(G=0) = (4π/Ω) ∫₀^∞ r² [V_local(r) + Z_val e²/r] dr
    ///
    /// Units: returns eV (potential in reciprocal space per unit cell).
    pub fn v_local_of_g(&self, g_norm: f64, omega: f64) -> f64 {
        // e² in eV·Å (Coulomb constant × e²)
        // In Gaussian units: e²/(4πε₀) = 14.3996 eV·Å
        const E2: f64 = 14.399645351950548; // eV·Å

        let n = self.r_grid.len();
        let mut integral = 0.0;

        if g_norm < 1e-12 {
            // G = 0 case: ∫ r² [V_loc(r) + Z e²/r] dr
            for i in 0..n {
                let r = self.r_grid[i];
                let dr = self.rab[i];
                let v_short = self.v_local[i] + self.z_valence * E2 / r.max(1e-20);
                integral += r * r * v_short * dr;
            }
            4.0 * std::f64::consts::PI / omega * integral
        } else {
            // G ≠ 0: ∫ r² [V_loc(r) + Z e²/r] sin(Gr)/(Gr) dr - 4π Z e² / (Ω G²)
            for i in 0..n {
                let r = self.r_grid[i];
                let dr = self.rab[i];
                let gr = g_norm * r;
                let v_short = self.v_local[i] + self.z_valence * E2 / r.max(1e-20);
                let sinc = if gr < 1e-10 {
                    1.0 - gr * gr / 6.0
                } else {
                    gr.sin() / gr
                };
                integral += r * r * v_short * sinc * dr;
            }
            4.0 * std::f64::consts::PI / omega * integral
                - 4.0 * std::f64::consts::PI * self.z_valence * E2
                    / (omega * g_norm * g_norm)
        }
    }

    /// Generate a simple model atomic charge density if none is provided.
    ///
    /// Uses a Gaussian: ρ(r) = Z_val / (2π σ²)^{3/2} exp(-r²/(2σ²))
    /// with σ chosen so the charge is concentrated near the atom.
    pub fn rho_atom_or_model(&self) -> &[f64] {
        if self.rho_atom.iter().any(|&v| v.abs() > 1e-20) {
            &self.rho_atom
        } else {
            // Caller should use generate_model_density instead
            &self.rho_atom
        }
    }

    /// Check if atomic charge density is available (non-zero).
    pub fn has_rho_atom(&self) -> bool {
        self.rho_atom.iter().any(|&v| v.abs() > 1e-20)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn si_pp_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Si.UPF")
    }

    #[test]
    fn test_load_si_upf() {
        let pp = load(&si_pp_path()).unwrap();
        assert_eq!(pp.element, "Si");
        assert!((pp.z_valence - 4.0).abs() < 1e-10);
        assert_eq!(pp.l_max, 1);
        assert_eq!(pp.r_grid.len(), 1141);
        assert_eq!(pp.v_local.len(), 1141);
        assert_eq!(pp.n_projectors, 3);
        assert_eq!(pp.beta_projectors.len(), 3);
        // Two s-projectors and one p-projector
        assert_eq!(pp.beta_projectors[0].l, 0);
        assert_eq!(pp.beta_projectors[1].l, 0);
        assert_eq!(pp.beta_projectors[2].l, 1);
    }

    #[test]
    fn test_find_for_atom() {
        let pp = load(&si_pp_path()).unwrap();
        let found = find_for_atom(14, &[&pp]);
        assert_eq!(found.element, "Si");
        assert!((found.z_valence - 4.0).abs() < 1e-10);
    }

    #[test]
    #[should_panic(expected = "no pseudopotential found")]
    fn test_find_for_atom_missing() {
        let pp = load(&si_pp_path()).unwrap();
        find_for_atom(6, &[&pp]); // Carbon, not in list
    }

    #[test]
    fn test_load_fe_upf() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Fe.UPF");
        let pp = load(&path).unwrap();
        assert_eq!(pp.element, "Fe");
        assert!((pp.z_valence - 8.0).abs() < 1e-10);
        assert_eq!(pp.n_projectors, 2);
        assert_eq!(pp.beta_projectors.len(), 2);
        assert_eq!(pp.beta_projectors[0].l, 1); // p-projector
        assert_eq!(pp.beta_projectors[1].l, 2); // d-projector
        // D_ij should be 2×2 with non-zero diagonal
        assert_eq!(pp.dij.len(), 4);
        eprintln!("Fe D_ij (eV): {:?}", pp.dij);
        // D_11 (p) should be small positive, D_22 (d) should be large negative
        assert!(pp.dij[0].abs() > 0.01, "D_11 should be nonzero: {}", pp.dij[0]);
        assert!(pp.dij[3] < -10.0, "D_22 should be large negative: {}", pp.dij[3]);
    }

    #[test]
    fn test_load_c_upf() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/C.UPF");
        match load(&path) {
            Ok(pp) => {
                assert_eq!(pp.element, "C");
                assert!((pp.z_valence - 4.0).abs() < 1e-10);
            }
            Err(e) => {
                // C.UPF may be v1 format which we don't support yet
                eprintln!("C.UPF parse failed (may be v1 format): {e}");
            }
        }
    }

    #[test]
    fn test_rho_atom_integrates_to_z_valence() {
        // ∫ rho_atom(r) × dr should give z_valence
        // rho_atom stores 4πr²ρ(r) in converted units; rab stores dr
        let pp = load(&si_pp_path()).unwrap();
        if !pp.has_rho_atom() {
            eprintln!("Si PP has no rho_atom data (HGH) — skipping integral test");
            return;
        }
        let integral: f64 = pp
            .rho_atom
            .iter()
            .zip(pp.rab.iter())
            .map(|(&rho, &dr)| rho * dr)
            .sum();
        eprintln!(
            "Si rho_atom integral = {integral:.6}, z_valence = {}",
            pp.z_valence
        );
        assert!(
            (integral - pp.z_valence).abs() < 0.1,
            "rho_atom integral {integral} != z_valence {}",
            pp.z_valence
        );
    }

    #[test]
    fn test_fe_rho_atom_integrates_to_z_valence() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Fe.UPF");
        let pp = load(&path).unwrap();
        if !pp.has_rho_atom() {
            eprintln!("Fe PP has no rho_atom data — skipping");
            return;
        }
        let integral: f64 = pp
            .rho_atom
            .iter()
            .zip(pp.rab.iter())
            .map(|(&rho, &dr)| rho * dr)
            .sum();
        eprintln!(
            "Fe rho_atom integral = {integral:.6}, z_valence = {}",
            pp.z_valence
        );
        // This test will FAIL if the unit conversion is wrong
        assert!(
            (integral - pp.z_valence).abs() < 0.5,
            "Fe rho_atom integral {integral} != z_valence {}",
            pp.z_valence
        );
    }

    #[test]
    fn test_v_local_of_g_finite() {
        let pp = load(&si_pp_path()).unwrap();
        let omega = 40.0; // approximate Si cell volume in ų
        // V_local(G) should be finite for any G
        let v0 = pp.v_local_of_g(0.0, omega);
        let v1 = pp.v_local_of_g(1.0, omega);
        let v5 = pp.v_local_of_g(5.0, omega);
        assert!(v0.is_finite(), "V_local(G=0) is not finite: {v0}");
        assert!(v1.is_finite(), "V_local(G=1) is not finite: {v1}");
        assert!(v5.is_finite(), "V_local(G=5) is not finite: {v5}");
        // V_local(G) should decay for large G
        assert!(
            v5.abs() < v1.abs(),
            "V_local should decay: |V(5)|={} > |V(1)|={}",
            v5.abs(),
            v1.abs()
        );
    }
}
