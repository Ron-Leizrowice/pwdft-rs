//! Local pseudopotential assembled in reciprocal space.
//!
//! Each atomic species contributes a spherically symmetric radial potential
//! `v_local^(α)(r)` (eV) whose spherical Bessel transform
//!
//! ```text
//!      v_local^(α)(|G|) = (4π/Ω) · ∫₀^∞ r² · v_local^(α)(r) · j₀(|G|r) dr
//! ```
//!
//! (in eV, per-cell units; the erf-regularized form used in production is
//! in [`crate::pseudopotential::UpfPseudoPotential::v_local_of_g`]) is
//! multiplied by the per-site structure factor `S_α(G) = exp(−iG·τ_α)`
//! and summed over all atoms in the cell:
//!
//! ```text
//!      V_local(G) = Σ_α S_α(G) · v_local^(α)(|G|).
//! ```
//!
//! Units: `r` in Å, `G` in 1/Å, `τ_α` in Å, `Ω` in Å³, `V_local(G)` in eV.

use std::collections::HashMap;

use elements_rs::Element;
use num_complex::Complex64;
use rayon::prelude::*;

use crate::{basis::BasisSet, crystal::Crystal, pseudopotential::UpfPseudoPotential};

/// Precomputed local pseudopotential in reciprocal space for a fixed
/// crystal geometry and basis cutoff.
pub struct LocalPotential {
    /// `V_local(G)` at each basis G-vector, in eV. Indexed 1:1 with
    /// [`BasisSet::g_vectors`].
    v_g: Vec<Complex64>,
}

impl LocalPotential {
    /// Build `V_local(G)` at every G-vector of `basis` for the supplied
    /// crystal geometry.
    ///
    /// # Panics
    /// Panics if any atom in `crystal` does not have a matching pseudopotential
    /// in the provided map. This invariant should be checked at the start of
    /// the simulation.
    pub fn new(crystal: &Crystal, basis: &BasisSet, pseudopotentials: &HashMap<Element, UpfPseudoPotential>) -> Self {
        let omega = crystal.lattice.volume();

        // Pre-collect atom positions and their corresponding PPs to avoid
        // repeated HashMap lookups in the parallel loop.
        let atom_data: Vec<_> = crystal
            .atoms
            .iter()
            .map(|atom| {
                let pp = pseudopotentials
                    .get(&atom.symbol)
                    .expect("BUG: Missing pseudopotential during LocalPotential assembly.");
                (atom.cart_position(&crystal.lattice), pp)
            })
            .collect();

        // Map every G-vector in the basis to its summed potential contribution.
        // We use into_par_iter because N_basis can be very large.
        let v_g = basis
            .g_vectors()
            .into_par_iter()
            .map(|g| {
                let g_norm = g.norm();

                atom_data.iter().fold(Complex64::default(), |acc, (tau, pp)| {
                    // Structure factor: S(G) = exp(-i G · τ)
                    let phase = -g.dot(tau);
                    let structure_factor = Complex64::cis(phase);

                    // Form factor: v_local(|G|) (includes 4π/Ω factor)
                    let v_form = pp.v_local_of_g(g_norm, omega);

                    acc + structure_factor * v_form
                })
            })
            .collect();

        Self { v_g }
    }

    /// Local pseudopotential in reciprocal space at the `ig`-th basis
    /// G-vector.
    pub fn v_of_g(&self, ig: usize) -> Complex64 {
        self.v_g[ig]
    }

    /// Borrow the full `V_local(G)` array as a slice.
    pub fn as_slice(&self) -> &[Complex64] {
        &self.v_g
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    use super::*;
    use crate::crystal::{Atom, Lattice};

    fn setup_si_test() -> (Crystal, HashMap<Element, UpfPseudoPotential>) {
        let a = 5.431;
        let crystal = Crystal {
            lattice: Lattice::new(
                a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
            ),
            atoms: vec![
                Atom::new(Element::Si, [0.0, 0.0, 0.0]),
                Atom::new(Element::Si, [0.25, 0.25, 0.25]),
            ],
        };

        let mut pps = HashMap::new();
        pps.insert(Element::Si, UpfPseudoPotential::load("Si").unwrap());

        (crystal, pps)
    }

    #[test]
    fn test_local_potential_size() {
        let (crystal, pps) = setup_si_test();
        let basis = BasisSet::new(&crystal.lattice, 100.0);

        let v_loc = LocalPotential::new(&crystal, &basis, &pps);
        assert_eq!(v_loc.as_slice().len(), basis.len());
    }

    #[test]
    fn test_local_potential_g0_real() {
        let (crystal, pps) = setup_si_test();
        let basis = BasisSet::new(&crystal.lattice, 100.0);

        let v_loc = LocalPotential::new(&crystal, &basis, &pps);
        let g0_idx = basis.index_of(0, 0, 0).expect("G=0 must exist in basis");
        let v_g0 = v_loc.v_of_g(g0_idx);

        // For a crystal with inversion symmetry around (0,0,0), V(G=0)
        // must be real. Here, structure factor for G=0 is exactly 2.0.
        assert_relative_eq!(v_g0.im, 0.0, epsilon = 1e-10);
    }
}
