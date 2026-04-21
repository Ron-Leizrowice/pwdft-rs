//! Kinetic part of the Kohn-Sham Hamiltonian in the plane-wave basis.
//!
//! Exposes [`build_kinetic`], which returns the diagonal matrix
//! `H⁰_{G,G'}(k) = δ_{GG'} · (ℏ²/2m) |k + G|²` at a given k-point.
//! This is the kinetic-only (free-electron) Hamiltonian used by the
//! band-structure diagnostic in [`crate::bandstructure`] and by tests
//! that isolate the kinetic contribution.
//!
//! The full SCF Hamiltonian (kinetic + local + Hartree + XC + non-local)
//! is assembled inline inside the SCF driver; this module only owns the
//! kinetic block.

use nalgebra::Vector3;
use num_complex::Complex64;

use crate::{basis::BasisSet, consts::HBAR2_OVER_2M};

/// Build the free-electron (kinetic-only) Hamiltonian at k-point k.
///
/// H_{G,G'}(k) = δ_{GG'} · (ℏ²/2m) |k + G|²
///
/// This is a diagonal matrix in the plane-wave basis.
pub fn build_kinetic(basis: &BasisSet, k: &Vector3<f64>) -> faer::Mat<Complex64> {
    let n = basis.len();
    let mut h = faer::Mat::<Complex64>::zeros(n, n);

    for (i, g) in basis.g_vectors().iter().enumerate() {
        let kpg = k + g;
        let ke = HBAR2_OVER_2M * kpg.norm_squared();
        h[(i, i)] = Complex64::new(ke, 0.0);
    }

    h
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::Lattice;
    use approx::relative_eq;

    fn si_basis() -> BasisSet {
        let a = 5.431;
        let si = Lattice::new(
            a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        );
        BasisSet::new(&si, 200.0)
    }

    #[test]
    fn test_kinetic_diagonal() {
        let basis = si_basis();
        let k = Vector3::zeros();
        let h = build_kinetic(&basis, &k);
        let n = basis.len();

        for i in 0..n {
            for j in 0..n {
                if i != j {
                    assert!(
                        h[(i, j)].norm() < 1e-15,
                        "off-diagonal H[{i},{j}] = {} should be zero",
                        h[(i, j)]
                    );
                }
            }
        }
    }

    #[test]
    fn test_kinetic_at_gamma() {
        let basis = si_basis();
        let k = Vector3::zeros();
        let h = build_kinetic(&basis, &k);

        let ke = basis.kinetic_energy();
        for (i, &expected) in ke.iter().enumerate() {
            assert!(
                relative_eq!(h[(i, i)].re, expected, epsilon = 1e-10),
                "H[{i},{i}] = {}, expected {expected}",
                h[(i, i)].re
            );
            assert!(h[(i, i)].im.abs() < 1e-15);
        }
    }

    #[test]
    fn test_kinetic_hermitian() {
        let basis = si_basis();
        let k = Vector3::new(0.1, 0.2, 0.3);
        let h = build_kinetic(&basis, &k);
        let n = basis.len();

        for i in 0..n {
            for j in 0..n {
                let diff = (h[(i, j)] - h[(j, i)].conj()).norm();
                assert!(diff < 1e-14, "H not Hermitian at ({i},{j}): diff={diff}");
            }
        }
    }
}
