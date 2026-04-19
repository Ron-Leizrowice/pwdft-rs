//! Local pseudopotential assembled in reciprocal space.
//!
//! Each atomic species contributes a spherically symmetric radial potential
//! `v_local^(α)(r)` (eV) whose spherical Bessel transform
//!
//! ```text
//!     v_local^(α)(|G|) = (4π/Ω) · ∫₀^∞ r² · v_local^(α)(r) · j₀(|G|r) dr
//! ```
//!
//! (in eV, per-cell units; the erf-regularized form used in production is
//! in [`crate::pseudopotential::PseudopotentialData::v_local_of_g`]) is
//! multiplied by the per-site structure factor `S_α(G) = exp(−iG·τ_α)`
//! and summed over all atoms in the cell:
//!
//! ```text
//!     V_local(G) = Σ_α S_α(G) · v_local^(α)(|G|).
//! ```
//!
//! Units: `r` in Å, `G` in 1/Å, `τ_α` in Å, `Ω` in Å³, `V_local(G)` in eV.
//! The `G = 0` component is finite (the Coulomb tail is subtracted inside
//! `v_local_of_g`); it is zeroed out inside the Hamiltonian to keep the
//! one-body operator diagonal-finite and re-added to the total energy by
//! the `scf::energy::with_g0_shift` helper.
//!
//! Reference: Martin, *Electronic Structure*, §11.4, Eq. (11.15)
//! (reciprocal-space form of a superposition of spherical atomic
//! potentials); Pickett, *Comput. Phys. Rep.* **9**, 115 (1989),
//! §III.B for the pseudopotential conventions.

use num_complex::Complex64;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    error::{PwdftError, Result},
    pseudopotential::PseudopotentialData,
};

/// Precomputed local pseudopotential in reciprocal space for a fixed
/// crystal geometry and basis cutoff.
///
/// Storage:
///
/// ```text
///     v_g[ig] = V_local(G_ig) = Σ_α S_α(G_ig) · v_local^(α)(|G_ig|)
/// ```
///
/// in eV, indexed 1:1 with [`BasisSet::g_vectors`]. The structure-factor
/// sum is finished at construction time; the SCF hot path only reads
/// this slice, so the per-iteration cost is a single contiguous scan.
///
/// Conventions:
/// - `V_local(G)` is complex in general; it is real whenever the crystal
///   has inversion symmetry about the origin and the atoms live at
///   inversion-paired positions (see `test_local_potential_g0_real`).
/// - The `G = 0` entry carries the finite residual after the Coulomb
///   `−4πZ_α·e²/|G|²` tail has been subtracted inside
///   [`crate::pseudopotential::PseudopotentialData::v_local_of_g`]. That
///   residual enters the total energy as `V_local(G=0) · N_el` via the
///   `scf::energy::with_g0_shift` helper, while the one-body
///   Hamiltonian keeps `V_local(G=0) = 0` so the diagonal stays
///   finite under the neutral-background convention.
///
/// Reference: Kleinman & Bylander, *Phys. Rev. Lett.* **48**, 1425
/// (1982) for the local+separable split that keeps this struct
/// independent of the angular-momentum-dependent
/// [`crate::potential::nonlocal::NonlocalPotential`].
pub struct LocalPotential {
    /// `V_local(G)` at each basis G-vector, in eV. Indexed 1:1 with
    /// [`BasisSet::g_vectors`]; length equals
    /// [`BasisSet::len`].
    v_g: Vec<Complex64>,
}

impl LocalPotential {
    /// Build `V_local(G)` at every G-vector of `basis` for the supplied
    /// crystal geometry.
    ///
    /// For each G the constructor evaluates
    ///
    /// ```text
    ///     V_local(G) = Σ_α exp(−iG·τ_α) · v_local^(α)(|G|)
    /// ```
    ///
    /// where
    /// - `G` runs over `basis.g_vectors()` (each in 1/Å);
    /// - `τ_α` is the Cartesian position of atom α in Å
    ///   (`Atom::cart_position`);
    /// - `v_local^(α)(|G|)` is in eV and comes from
    ///   [`crate::pseudopotential::PseudopotentialData::v_local_of_g`],
    ///   which applies the erf-regularized spherical Bessel transform
    ///   (see its docstring for the Coulomb-tail subtraction conventions);
    /// - `Ω = crystal.lattice.volume()` (Å³) is folded into
    ///   `v_local^(α)` via the `4π/Ω` prefactor of the Bessel transform.
    ///
    /// Invariants on input:
    /// - `pseudopotentials` must contain at least one
    ///   [`PseudopotentialData`] whose `z` matches each `atom.z` in the
    ///   crystal; otherwise the constructor errors (see below).
    /// - `basis` must be built from the same lattice as `crystal`; the
    ///   G-vectors come directly from `basis.g_vectors()` so no check is
    ///   possible here — the caller owns that invariant.
    ///
    /// Returns a [`LocalPotential`] whose `v_g` has length `basis.len()`
    /// and units of eV.
    ///
    /// Reference: Martin, *Electronic Structure*, §11.4 (separable
    /// local+non-local split of the pseudopotential).
    ///
    /// # Errors
    /// Returns [`PwdftError::MissingPseudopotential`] if any atom's
    /// atomic number `z` has no corresponding
    /// [`PseudopotentialData`] in `pseudopotentials`.
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

    /// Local pseudopotential in reciprocal space at the `ig`-th basis
    /// G-vector.
    ///
    /// Returns
    ///
    /// ```text
    ///     V_local(G_ig) = Σ_α exp(−iG_ig·τ_α) · v_local^(α)(|G_ig|)
    /// ```
    ///
    /// in eV, where
    /// - `S_α(G) = exp(−iG·τ_α)` is the structure factor for atom α at
    ///   Cartesian position `τ_α` (Å);
    /// - `v_local^(α)(|G|)` is the erf-regularized spherical Bessel
    ///   transform
    ///   ```text
    ///     v_local^(α)(|G|) = (4π/Ω) · ∫ r² · [v_loc(r) + Z_α e² · erf(r)/r]
    ///                                     · sin(|G|r)/(|G|r) dr
    ///                       − 4π · Z_α e² · exp(−|G|²/4) / (Ω · |G|²)
    ///   ```
    ///   in eV, computed by
    ///   [`crate::pseudopotential::PseudopotentialData::v_local_of_g`]
    ///   (see that docstring for the decomposition; the Gaussian width
    ///   is fixed at 1 Å and cancels identically in the sum).
    /// - The `G = 0` component carries a finite residual; the divergent
    ///   `−4π·Z_α·e²/|G|²` piece has been subtracted by the two
    ///   branches of `v_local_of_g`. That residual is paid back into
    ///   the total energy as `V_local(G=0) · N_el` via the
    ///   `scf::energy::with_g0_shift` helper; the Hamiltonian's
    ///   diagonal keeps its own G=0 entry zero under the neutral-
    ///   background convention (stashed by `ScfContext::new` in
    ///   `scf::context`).
    ///
    /// Reference: Martin, *Electronic Structure*, §11.4 for the
    /// reciprocal-space structure-factor form; Kleinman & Bylander,
    /// *Phys. Rev. Lett.* **48**, 1425 (1982) for the local/non-local
    /// split the rest of the pseudopotential machinery relies on.
    ///
    /// # Panics
    /// Panics if `ig >= self.as_slice().len()` (plain slice-index
    /// out-of-bounds).
    #[must_use]
    pub fn v_of_g(&self, ig: usize) -> Complex64 {
        self.v_g[ig]
    }

    /// Borrow the full `V_local(G)` array.
    ///
    /// Returns a slice of length [`BasisSet::len`] (the value passed to
    /// [`LocalPotential::new`]), indexed 1:1 with
    /// [`BasisSet::g_vectors`]; each entry is in eV and is the same
    /// quantity as the single-index result of [`Self::v_of_g`].
    ///
    /// Used by Hamiltonian assembly (e.g. `scf::potentials`) to add the
    /// local-PP column to `V_eff(G)` in bulk, and by the total-energy
    /// kernel for the `⟨ρ|V_local⟩` convolution.
    #[must_use]
    pub fn as_slice(&self) -> &[Complex64] {
        &self.v_g
    }
}

#[cfg(test)]
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
