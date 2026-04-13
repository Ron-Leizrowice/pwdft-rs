//! Hartree potential: the classical electrostatic potential from the electron density.
//!
//! In reciprocal space: V_H(G) = 4π e² ρ(G) / |G|²  for G ≠ 0
//!                      V_H(G=0) = 0  (neutralizing background)
//!
//! The Hartree energy is: E_H = (Ω/2) Σ_G |ρ(G)|² · 4π e² / |G|²

use num_complex::Complex64;

/// Coulomb constant e² in eV·Å (= e²/(4πε₀) in Gaussian units).
pub const E2: f64 = 14.399645351950548;

/// Compute Hartree potential V_H(G) from charge density ρ(G).
///
/// V_H(G) = 4π e² ρ(G) / |G|²  for G ≠ 0
/// V_H(G=0) = 0
///
/// `rho_g`: charge density in reciprocal space (e/ų in G-space convention).
/// `g_vectors`: Cartesian G-vectors (1/Å).
///
/// Returns V_H(G) in eV.
pub fn hartree_potential(
    rho_g: &[Complex64],
    g_vectors: &[nalgebra::Vector3<f64>],
) -> Vec<Complex64> {
    let fourpi_e2 = 4.0 * std::f64::consts::PI * E2;

    rho_g
        .iter()
        .zip(g_vectors.iter())
        .map(|(&rho, g)| {
            let g2 = g.norm_squared();
            if g2 < 1e-20 {
                Complex64::new(0.0, 0.0) // G=0: neutralizing background
            } else {
                rho * fourpi_e2 / g2
            }
        })
        .collect()
}

/// Compute Hartree energy from charge density in reciprocal space.
///
/// E_H = (Ω/2) Σ_{G≠0} |ρ(G)|² · 4π e² / |G|²
///
/// Returns energy in eV.
pub fn hartree_energy(
    rho_g: &[Complex64],
    g_vectors: &[nalgebra::Vector3<f64>],
    omega: f64,
) -> f64 {
    let fourpi_e2 = 4.0 * std::f64::consts::PI * E2;

    let sum: f64 = rho_g
        .iter()
        .zip(g_vectors.iter())
        .map(|(&rho, g)| {
            let g2 = g.norm_squared();
            if g2 < 1e-20 {
                0.0
            } else {
                rho.norm_sqr() * fourpi_e2 / g2
            }
        })
        .sum();

    0.5 * omega * sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_hartree_g0_zero() {
        let rho_g = vec![Complex64::new(1.0, 0.0)];
        let g_vectors = vec![Vector3::zeros()];
        let v_h = hartree_potential(&rho_g, &g_vectors);
        assert!((v_h[0].norm()) < 1e-15, "V_H(G=0) should be zero");
    }

    #[test]
    fn test_hartree_positive_definite() {
        // E_H should always be non-negative
        let g_vectors = vec![
            Vector3::zeros(),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let rho_g = vec![
            Complex64::new(5.0, 0.0),
            Complex64::new(0.3, 0.1),
            Complex64::new(-0.1, 0.2),
        ];
        let e_h = hartree_energy(&rho_g, &g_vectors, 40.0);
        assert!(e_h >= 0.0, "Hartree energy should be non-negative: {e_h}");
    }
}
