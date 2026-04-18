//! Local pseudopotential in reciprocal space.
//!
//! V_local(G) = Σ_atoms S_atom(G) · v_local(|G|)
//!
//! where S_atom(G) = exp(-i G · τ_atom) is the structure factor
//! and v_local(|G|) is the spherical Bessel transform of the radial local potential.

use num_complex::Complex64;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    error::{PwdftError, Result},
    pseudopotential::PseudopotentialData,
};

/// Precomputed local pseudopotential for efficient Hamiltonian construction.
pub struct LocalPotential {
    /// V_local(G) for each G-vector in the basis (eV).
    v_g: Vec<Complex64>,
}

impl LocalPotential {
    /// Build the local pseudopotential in reciprocal space.
    ///
    /// For each G-vector, computes:
    ///   V_ps(G) = Σ_atoms S_atom(G) · v_local(|G|)
    ///
    /// where S_atom(G) = exp(-i G · τ_atom) and v_local(|G|) comes from
    /// the pseudopotential's spherical Bessel transform.
    ///
    /// # Errors
    /// Returns `PwdftError::MissingPseudopotential` if any atom lacks a loaded PP.
    pub fn new(
        crystal: &Crystal,
        basis: &BasisSet,
        pseudopotentials: &[&PseudopotentialData],
    ) -> Result<Self> {
        let omega = crystal.lattice.volume();
        let n = basis.len();
        let mut v_g = vec![Complex64::new(0.0, 0.0); n];

        for (ig, g) in basis.g_vectors().iter().enumerate() {
            let g_norm = g.norm();

            for atom in &crystal.atoms {
                let pp = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials)
                    .ok_or_else(|| PwdftError::MissingPseudopotential(
                        format!("Z={} not found in loaded pseudopotentials", atom.z)
                    ))?;

                // Structure factor: S(G) = exp(-i G · τ)
                let tau = atom.cart_position(&crystal.lattice);
                let phase = -g.dot(&tau);
                let structure_factor = Complex64::cis(phase);

                // Form factor: v_local(|G|)
                let v_form = pp.v_local_of_g(g_norm, omega);

                v_g[ig] += structure_factor * v_form;
            }
        }

        Ok(Self { v_g })
    }

    /// Get V_local(G) for a given G-vector index.
    #[must_use]
    pub fn v_of_g(&self, ig: usize) -> Complex64 {
        self.v_g[ig]
    }

    /// Get the full V_local(G) array.
    #[must_use]
    pub fn as_slice(&self) -> &[Complex64] {
        &self.v_g
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "ERR2 § Phase 0: in-src test modules are allowed to panic")]
mod tests {
    use super::*;
    use crate::crystal::{Atom, Lattice};
    use nalgebra::Vector3;

    fn si_crystal() -> Crystal {
        let a = 5.431;
        Crystal {
            lattice: Lattice::new(
                a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
            ),
            atoms: vec![
                Atom::new(14, [0.0, 0.0, 0.0]),
                Atom::new(14, [0.25, 0.25, 0.25]),
            ],
        }
    }

    #[test]
    fn test_local_potential_size() {
        let crystal = si_crystal();
        let basis = BasisSet::new(&crystal.lattice, 100.0);
        let pp = crate::pseudopotential::load(
            &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap();

        let v_loc = LocalPotential::new(&crystal, &basis, &[&pp]).unwrap();
        assert_eq!(v_loc.as_slice().len(), basis.len());
    }

    #[test]
    fn test_local_potential_g0_real() {
        // V_local(G=0) should be real for a crystal with inversion symmetry
        let crystal = si_crystal();
        let basis = BasisSet::new(&crystal.lattice, 100.0);
        let pp = crate::pseudopotential::load(
            &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap();

        let v_loc = LocalPotential::new(&crystal, &basis, &[&pp]).unwrap();
        let g0_idx = basis.index_of(0, 0, 0).unwrap();
        let v_g0 = v_loc.v_of_g(g0_idx);
        // With 2 atoms at (0,0,0) and (1/4,1/4,1/4), the structure factor
        // for G=0 is 2 (both exp(0) = 1), so V(G=0) should be real.
        assert!(
            v_g0.im.abs() < 1e-6,
            "V_local(G=0) should be real, got im={}",
            v_g0.im
        );
    }
}
