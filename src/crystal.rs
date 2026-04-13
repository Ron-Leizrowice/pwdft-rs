use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

use crate::consts::PI;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Crystal {
    pub atoms: Vec<Atom>,
    pub lattice: Lattice,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub z: u32,
    pub position: [f64; 3], // fractional coordinates
}

impl Atom {
    pub fn new(z: u32, frac: [f64; 3]) -> Self {
        Self { z, position: frac }
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

impl Lattice {
    pub fn new(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> Self {
        Self { a, b, c }
    }

    pub fn volume(&self) -> f64 {
        self.a.cross(&self.b).dot(&self.c)
    }

    pub fn reciprocal(&self) -> Self {
        let factor = 2.0 * PI / self.volume();
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
    use super::*;
    use approx::relative_eq;

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
        let atom = Atom::new(14, [0.25, 0.25, 0.25]);
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
                Atom::new(14, [0.0, 0.0, 0.0]),
                Atom::new(14, [0.25, 0.25, 0.25]),
            ],
        };
        let serialized = toml::to_string(&crystal).unwrap();
        let deserialized: Crystal = toml::from_str(&serialized).unwrap();
        assert_eq!(deserialized.atoms.len(), 2);
        assert_eq!(deserialized.atoms[0].z, 14);
    }
}
