use crate::{consts::HBAR2_OVER_2M, crystal::Lattice};
use nalgebra::Vector3;

pub struct BasisSet {
    pw: Vec<Vector3<f64>>,
    ecut: f64,
}
impl BasisSet {
    pub fn new(lattice: &Lattice, ecut: f64) -> Self {
        let recip = lattice.reciprocal();
        let g_max_sq = ecut / HBAR2_OVER_2M;

        let g_max = g_max_sq.sqrt();
        let (b1, b2, b3) = (recip.a, recip.b, recip.c);

        let n1_max = (g_max / b1.norm()).ceil() as i32;
        let n2_max = (g_max / b2.norm()).ceil() as i32;
        let n3_max = (g_max / b3.norm()).ceil() as i32;

        let mut pw = Vec::new();
        for n1 in -n1_max..=n1_max {
            for n2 in -n2_max..=n2_max {
                for n3 in -n3_max..=n3_max {
                    let g = n1 as f64 * b1 + n2 as f64 * b2 + n3 as f64 * b3;
                    if g.norm_squared() <= g_max_sq {
                        pw.push(g);
                    }
                }
            }
        }

        Self { pw, ecut }
    }

    pub fn len(&self) -> usize {
        self.pw.len()
    }

    pub fn kinetic_energy(&self) -> Vec<f64> {
        self.pw.iter().map(|g| HBAR2_OVER_2M * g.norm_squared()).collect()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn si_basis(ecut: f64) -> BasisSet {
        let a = 5.431;
        let si = Lattice::new(
            a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        );
        BasisSet::new(&si, ecut)
    }

    #[test]
    fn test_si_basis_count() {
        let basis = si_basis(200.0);
        println!("G-vectors at 200 eV: {}", basis.len());
        assert!(basis.len() == 259, "expected 259 G-vectors at 200 eV, got {}", basis.len());
    }

    #[test]
    fn test_kinetic_energy_contains_zero() {
        let basis = si_basis(200.0);
        let ke = basis.kinetic_energy();
        let min_ke = ke.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(min_ke.abs() < 1e-12, "G=0 should have zero kinetic energy, got {min_ke}");
    }

    #[test]
    fn test_kinetic_energy_within_cutoff() {
        let ecut = 200.0;
        let basis = si_basis(ecut);
        let ke = basis.kinetic_energy();
        let max_ke = ke.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_ke <= ecut + 1e-10, "max KE {max_ke} exceeds cutoff {ecut}");
    }
}
