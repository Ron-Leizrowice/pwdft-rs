//! Pseudopotential data and loading.
//!
//! Holds the unit-converted norm-conserving pseudopotential in internal
//! units (eV, Å) as [`PseudopotentialData`]: radial grid, local potential
//! `V_local(r)`, Kleinman–Bylander `β` projectors with the `D_ij` coupling
//! matrix, atomic density `ρ_atom`, and the optional NLCC core density.
//! [`load`] parses a UPF v2 file via [`upf::parse`];
//! [`PseudopotentialData::v_local_of_g`] returns the spherical-Bessel
//! transform `V_local(G)` used by the SCF driver to assemble the local
//! part of the effective potential.

pub mod upf;
use std::{path::PathBuf, sync::LazyLock};

pub use upf::{BetaProjector, UpfPseudoPotential};

static PP_LIBRARY_ROOT: LazyLock<PathBuf> =
    LazyLock::new(|| PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials"));

pub enum Pseudopotential {
    NormConserving(UpfPseudoPotential),
}

#[cfg(test)]
mod tests {
    use elements_rs::Element;

    use super::*;

    #[test]
    fn test_load_si_upf() {
        let pp = UpfPseudoPotential::load(Element::Si).unwrap();
        assert_eq!(pp.element, Element::Si);
        assert!((pp.z_valence - 4.0).abs() < 1e-10);
        assert!(pp.l_max >= 1, "Si should have l_max >= 1");
        assert!(pp.r_grid.len() > 100, "Radial grid too small");
        assert_eq!(pp.v_local.len(), pp.r_grid.len());
        assert!(pp.n_projectors() > 0, "Should have projectors");
    }

    #[test]
    fn test_load_fe_upf() {
        let pp = UpfPseudoPotential::load(Element::Fe).unwrap();
        assert_eq!(pp.element, Element::Fe);
        assert!(pp.z_valence >= 8.0, "Fe should have >= 8 valence electrons");
        assert!(pp.n_projectors() > 0, "Fe should have projectors");
        // D_ij should be non-trivial (nonzero diagonal)
        let dij_max: f64 = pp.dij.iter().map(|d| d.abs()).fold(0.0, f64::max);
        assert!(dij_max > 0.01, "D_ij should have nonzero entries, max={dij_max}");
    }

    #[test]
    fn test_load_c_upf() {
        let pp = UpfPseudoPotential::load(Element::C).unwrap();
        assert_eq!(pp.element, Element::C);
        assert!((pp.z_valence - 4.0).abs() < 1e-10);
        assert!(pp.n_projectors() > 0);
    }

    #[test]
    fn test_rho_atom_integrates_to_z_valence() {
        // ∫ rho_atom(r) × dr should give z_valence
        // rho_atom stores 4πr²ρ(r) in converted units; rab stores dr
        let pp = UpfPseudoPotential::load(Element::Si).unwrap();
        if !pp.has_rho_atom() {
            log::debug!("Si PP has no rho_atom data (HGH) — skipping integral test");
            return;
        }
        let integral: f64 = pp.rho_atom.iter().zip(pp.rab.iter()).map(|(&rho, &dr)| rho * dr).sum();
        assert!(
            (integral - pp.z_valence).abs() < 0.1,
            "rho_atom integral {integral} != z_valence {}",
            pp.z_valence
        );
    }

    #[test]
    fn test_fe_rho_atom_integrates_to_z_valence() {
        let pp = UpfPseudoPotential::load(Element::Fe).unwrap();
        if !pp.has_rho_atom() {
            log::debug!("Fe PP has no rho_atom data — skipping");
            return;
        }
        let integral: f64 = pp.rho_atom.iter().zip(pp.rab.iter()).map(|(&rho, &dr)| rho * dr).sum();
        // This test will FAIL if the unit conversion is wrong
        assert!(
            (integral - pp.z_valence).abs() < 0.5,
            "Fe rho_atom integral {integral} != z_valence {}",
            pp.z_valence
        );
    }

    #[test]
    fn test_v_local_of_g_finite() {
        let pp = UpfPseudoPotential::load(Element::Si).unwrap();
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
