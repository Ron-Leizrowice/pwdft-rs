//! Crystal geometry: lattice vectors, atomic positions, reciprocal basis.
//!
//! [`Crystal`] aggregates a [`Lattice`] (the three primitive vectors a, b, c
//! in Å) with a list of [`Atom`]s (atomic number + fractional coordinates).
//! The lattice knows how to compute its reciprocal-space partners, unit
//! cell volume, and nearest-image distances, all in the engine's native
//! Å / Å⁻¹ units.
//!
//! This module is the geometric input to essentially everything else:
//! [`crate::basis`] enumerates G-vectors against the reciprocal lattice,
//! [`crate::kpoints`] builds Monkhorst-Pack grids in its Brillouin zone,
//! and [`crate::ewald`] sums ion-ion Coulomb interactions over its atoms.

use std::f64::consts::PI;

use elements_rs::Element;
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

use crate::settings::{AtomInput, CrystalInput};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Crystal {
    pub atoms: Vec<Atom>,
    pub lattice: Lattice,
}

impl From<CrystalInput> for Crystal {
    fn from(input: CrystalInput) -> Self {
        Self {
            lattice: Lattice::from(input.lattice),
            atoms: input.atoms.into_iter().map(Atom::from).collect(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub symbol: Element,
    pub z: u8,
    pub position: [f64; 3], // fractional coordinates
}

impl From<AtomInput> for Atom {
    fn from(a: AtomInput) -> Self {
        Self {
            z: a.symbol.into(),
            symbol: a.symbol,
            position: a.position,
        }
    }
}

impl Atom {
    /// Creates a new Atom. Accepts Element, u8 (atomic number), or &str
    /// (symbol).
    ///
    /// # Panics
    /// Panics if the provided element identifier is invalid.
    pub fn new<E>(element_like: E, frac: [f64; 3]) -> Self
    where
        E: TryInto<Element> + std::fmt::Debug + Copy,
        <E as TryInto<Element>>::Error: std::fmt::Debug,
    {
        let symbol: Element = element_like.try_into().unwrap_or_else(|_| {
            panic!("Failed to create Atom: '{element_like:?}' is not a valid element symbol or atomic number.")
        });

        Self {
            symbol,
            z: symbol.into(),
            position: frac,
        }
    }

    /// Convert fractional coordinates to Cartesian (Å).
    pub fn cart_position(&self, lattice: &Lattice) -> Vector3<f64> {
        let [f1, f2, f3] = self.position;
        f1 * lattice.a + f2 * lattice.b + f3 * lattice.c
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Lattice {
    pub a: Vector3<f64>,
    pub b: Vector3<f64>,
    pub c: Vector3<f64>,
}

impl From<[[f64; 3]; 3]> for Lattice {
    fn from([a, b, c]: [[f64; 3]; 3]) -> Self {
        Self::new(a, b, c)
    }
}

impl Lattice {
    pub fn new<V: Into<Vector3<f64>>>(a: V, b: V, c: V) -> Self {
        Self {
            a: a.into(),
            b: b.into(),
            c: c.into(),
        }
    }

    /// Cell volume Ω = |a · (b × c)|.
    ///
    /// Always positive regardless of lattice vector handedness.
    pub fn volume(&self) -> f64 {
        self.a.cross(&self.b).dot(&self.c).abs()
    }

    pub fn reciprocal(&self) -> Self {
        // Use signed triple product to get correct reciprocal vector directions
        let triple = self.a.cross(&self.b).dot(&self.c);
        let factor = 2.0 * PI / triple;
        Self {
            a: self.b.cross(&self.c) * factor,
            b: self.c.cross(&self.a) * factor,
            c: self.a.cross(&self.b) * factor,
        }
    }

    /// 3×3 matrix whose columns are the lattice vectors.
    pub fn matrix(&self) -> Matrix3<f64> {
        Matrix3::from_columns(&[self.a, self.b, self.c])
    }
}

#[cfg(test)]
mod tests {
    use approx::relative_eq;

    use super::*;

    fn si_lattice() -> Lattice {
        let si_a = 5.431;
        Lattice::new(
            si_a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            si_a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            si_a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        )
    }

    #[test]
    fn test_reciprocal() {
        let reciprocal = si_lattice().reciprocal();
        assert!(relative_eq!(
            reciprocal.a,
            Vector3::new(-1.157, 1.157, 1.157),
            epsilon = 1e-3
        ));
        assert!(relative_eq!(
            reciprocal.b,
            Vector3::new(1.157, -1.157, 1.157),
            epsilon = 1e-3
        ));
        assert!(relative_eq!(
            reciprocal.c,
            Vector3::new(1.157, 1.157, -1.157),
            epsilon = 1e-3
        ));
    }

    #[test]
    fn test_cart_position() {
        let lat = si_lattice();
        // Atom at (0.25, 0.25, 0.25) in fractional coords
        let atom = Atom::new(Element::Si, [0.25, 0.25, 0.25]);
        let cart = atom.cart_position(&lat);
        // Should be at (a/4)(0+1+1, 1+0+1, 1+1+0) = (a/4)(2,2,2) = a/2 * (1,1,1) * 0.5
        let expected = 0.25 * (lat.a + lat.b + lat.c);
        assert!(relative_eq!(cart, expected, epsilon = 1e-10));
    }

    #[test]
    fn test_serde_roundtrip() {
        let lat = si_lattice();
        let crystal = Crystal {
            lattice: lat,
            atoms: vec![
                Atom::new(Element::Si, [0.0, 0.0, 0.0]),
                Atom::new(Element::Si, [0.25, 0.25, 0.25]),
            ],
        };
        let serialized = serde_yaml_ng::to_string(&crystal).unwrap();
        let deserialized: Crystal = serde_yaml_ng::from_str(&serialized).unwrap();
        assert_eq!(deserialized.atoms.len(), 2);
        assert_eq!(deserialized.atoms[0].z, 14);
    }
}
