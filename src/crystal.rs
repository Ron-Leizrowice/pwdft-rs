use nalgebra::Vector3;

use crate::consts::PI;

pub struct Crystal {
    pub atoms: Vec<Atom>,
    pub lattice: Lattice,
}

pub struct Atom {
    pub z: u32,
    pub position: Vector3<f64>, // fractional coordinates
}

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
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

    #[test]
    fn test_reciprocal() {
        let si_a = 5.431;
        let si = Lattice::new(
            si_a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            si_a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            si_a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        );
        let reciprocal = si.reciprocal();
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
}
