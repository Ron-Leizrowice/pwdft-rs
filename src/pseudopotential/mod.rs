pub mod upf;

use std::path::Path;

use crate::error::{PwdftError, Result};

/// Unit-converted pseudopotential data in internal units (eV, Å).
///
/// Parsed from UPF v2 format (Quantum ESPRESSO).
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
    /// Atomic charge density on radial grid (e/Å, stores 4πr²ρ(r)).
    /// May be empty if not provided by the pseudopotential.
    pub rho_atom: Vec<f64>,
    /// Nonlinear core correction (NLCC) charge density on radial grid.
    /// Stores 4πr²ρ_core(r) in e/Å. Empty if NLCC is not present.
    pub core_charge: Vec<f64>,
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
use crate::consts::{RY_TO_EV, BOHR_TO_ANG};

/// Load and parse a UPF v2 pseudopotential file.
pub fn load(path: &Path) -> Result<PseudopotentialData> {
    let content = std::fs::read_to_string(path)?;

    if content.contains("<UPF") || content.contains("<PP_HEADER") {
        upf::parse(&content)
    } else {
        Err(PwdftError::Parse(format!(
            "unrecognized pseudopotential format in {} (only UPF v2 is supported)",
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
    /// Number of non-local beta projectors.
    pub fn n_projectors(&self) -> usize {
        self.beta_projectors.len()
    }

    /// Whether this PP has nonlinear core correction.
    pub fn has_nlcc(&self) -> bool {
        !self.core_charge.is_empty()
    }

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
        use crate::consts::E2_COULOMB as E2;

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
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf")
    }

    #[test]
    fn test_load_si_upf() {
        let pp = load(&si_pp_path()).unwrap();
        assert_eq!(pp.element, "Si");
        assert!((pp.z_valence - 4.0).abs() < 1e-10);
        assert!(pp.l_max >= 1, "Si should have l_max >= 1");
        assert!(pp.r_grid.len() > 100, "Radial grid too small");
        assert_eq!(pp.v_local.len(), pp.r_grid.len());
        assert!(pp.n_projectors() > 0, "Should have projectors");
    }

    #[test]
    fn test_find_for_atom() {
        let pp = load(&si_pp_path()).unwrap();
        let found = find_for_atom(14, &[&pp]);
        assert_eq!(found.element, "Si");
        assert!((found.z_valence - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_load_fe_upf() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Fe.upf");
        let pp = load(&path).unwrap();
        assert_eq!(pp.element, "Fe");
        assert!(pp.z_valence >= 8.0, "Fe should have >= 8 valence electrons");
        assert!(pp.n_projectors() > 0, "Fe should have projectors");
        // D_ij should be non-trivial (nonzero diagonal)
        let dij_max: f64 = pp.dij.iter().map(|d| d.abs()).fold(0.0, f64::max);
        assert!(dij_max > 0.01, "D_ij should have nonzero entries, max={dij_max}");
    }

    #[test]
    fn test_load_c_upf() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/C.upf");
        let pp = load(&path).unwrap();
        assert_eq!(pp.element, "C");
        assert!((pp.z_valence - 4.0).abs() < 1e-10);
        assert!(pp.n_projectors() > 0);
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
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Fe.upf");
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
