//! Ewald summation for ion-ion electrostatic energy.
//!
//! E_ewald = E_real + E_recip + E_self + E_background
//!
//! Uses the standard Ewald decomposition of the Coulomb sum between point charges
//! in a periodic crystal.

use std::f64::consts::PI;

use nalgebra::Vector3;
use num_complex::Complex64;

use crate::{
    crystal::Crystal,
    potential::hartree::E2,
    pseudopotential::PseudopotentialData,
};

/// Compute the Ewald ion-ion energy for a crystal. Returns energy in eV.
///
/// Decomposes the Coulomb sum of periodic point charges into four terms:
///
///   E_real  = (e²/2) Σ'_{i,j,T} Z_i Z_j erfc(η|r_ij+T|) / |r_ij+T|
///   E_recip = (2πe²/Ω) Σ_{G≠0} |S(G)|² exp(-|G|²/(4η²)) / |G|²
///   E_self  = -(η/√π) e² Σ_i Z_i²
///   E_bg    = -πe² (Σ Z_i)² / (2Ωη²)
///
/// where S(G) = Σ_i Z_i exp(iG·r_i) is the charge-weighted structure factor,
/// η = (N_atoms π/Ω)^{1/3} balances real/reciprocal cost, and the primed sum
/// excludes i=j when T=0 (self-interaction).
///
/// Cutoffs: g_max = 10η (reciprocal), r_max = 10/η (real).
pub fn ewald_energy(crystal: &Crystal, pseudopotentials: &[&PseudopotentialData]) -> f64 {
    let omega = crystal.lattice.volume();

    // Get charges
    let charges: Vec<f64> = crystal
        .atoms
        .iter()
        .map(|a| crate::pseudopotential::find_for_atom(a.z, pseudopotentials).z_valence)
        .collect();

    let positions: Vec<Vector3<f64>> = crystal
        .atoms
        .iter()
        .map(|a| a.cart_position(&crystal.lattice))
        .collect();

    let n_atoms = positions.len();

    // Choose Ewald parameter η (controls real/reciprocal space partition)
    // η ∝ (N_atoms / Ω)^{1/3}
    let eta = (n_atoms as f64 * PI / omega).powf(1.0 / 3.0);
    let eta2 = eta * eta;

    // Reciprocal space sum
    let recip = crystal.lattice.reciprocal();
    let g_max = 10.0 * eta; // cutoff for G-vectors
    let n1_max = (g_max / recip.a.norm()).ceil() as i32;
    let n2_max = (g_max / recip.b.norm()).ceil() as i32;
    let n3_max = (g_max / recip.c.norm()).ceil() as i32;

    let mut e_recip = 0.0;
    for n1 in -n1_max..=n1_max {
        for n2 in -n2_max..=n2_max {
            for n3 in -n3_max..=n3_max {
                if n1 == 0 && n2 == 0 && n3 == 0 {
                    continue;
                }
                let g = n1 as f64 * recip.a + n2 as f64 * recip.b + n3 as f64 * recip.c;
                let g2 = g.norm_squared();
                if g2 < 1e-12 {
                    continue; // guard against near-zero G² from floating-point noise
                }

                // Structure factor S(G) = Σ_i Z_i exp(iG·r_i)
                let s: Complex64 = positions.iter().enumerate()
                    .map(|(i, pos)| charges[i] * Complex64::cis(g.dot(pos)))
                    .sum();

                e_recip += s.norm_sqr() * (-g2 / (4.0 * eta2)).exp() / g2;
            }
        }
    }
    e_recip *= 2.0 * PI * E2 / omega;

    // Real space sum
    let r_max = 10.0 / eta;
    let l1_max = (r_max / crystal.lattice.a.norm()).ceil() as i32;
    let l2_max = (r_max / crystal.lattice.b.norm()).ceil() as i32;
    let l3_max = (r_max / crystal.lattice.c.norm()).ceil() as i32;

    let mut e_real = 0.0;
    for l1 in -l1_max..=l1_max {
        for l2 in -l2_max..=l2_max {
            for l3 in -l3_max..=l3_max {
                let t = l1 as f64 * crystal.lattice.a
                    + l2 as f64 * crystal.lattice.b
                    + l3 as f64 * crystal.lattice.c;

                for i in 0..n_atoms {
                    for j in 0..n_atoms {
                        let r = positions[i] - positions[j] + t;
                        let r_norm = r.norm();
                        if r_norm < 1e-10 {
                            continue; // skip self-interaction in same cell
                        }
                        e_real += charges[i] * charges[j] * erfc(eta * r_norm) / r_norm;
                    }
                }
            }
        }
    }
    e_real *= 0.5 * E2;

    // Self energy correction
    let z2_sum: f64 = charges.iter().map(|&z| z * z).sum();
    let e_self = -eta / PI.sqrt() * z2_sum * E2;

    // Background charge correction (for charged cells — typically zero)
    let z_sum: f64 = charges.iter().sum();
    let e_bg = -PI * z_sum * z_sum * E2 / (2.0 * omega * eta2);

    e_real + e_recip + e_self + e_bg
}

/// Complementary error function erfc(x) = 1 - erf(x).
///
/// Delegates to `puruspe::erfc` — validated, full f64 precision.
/// puruspe also provides erf, gamma, beta for future GGA/PAW use.
fn erfc(x: f64) -> f64 {
    puruspe::erfc(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::{Atom, Lattice};
    use approx::relative_eq;

    #[test]
    fn test_erfc_values() {
        // erfc(0) = 1 exactly
        assert!(relative_eq!(erfc(0.0), 1.0, epsilon = 1e-12));
        // erfc(1) = 0.15729920705... (NIST DLMF)
        assert!(relative_eq!(erfc(1.0), 0.157_299_207_050_285_13, epsilon = 1e-7));
        // erfc(2) = 0.00467773498...
        assert!(relative_eq!(erfc(2.0), 0.004_677_734_981_047_266, epsilon = 1e-7));
        // erfc(5) ≈ 1.537e-12
        assert!(erfc(5.0) < 1e-10);
        // Symmetry: erfc(-x) = 2 - erfc(x)
        assert!(relative_eq!(erfc(-1.0), 2.0 - erfc(1.0), epsilon = 1e-14));
        assert!(relative_eq!(erfc(-2.0), 2.0 - erfc(2.0), epsilon = 1e-14));
    }

    #[test]
    fn test_ewald_nacl() {
        // NaCl structure: Madelung constant M = 1.747_565
        // E_ewald = -M × e² / a₀ per ion pair
        // For a cubic cell with a = 5.64 Å (NaCl lattice constant)
        let a = 5.64;
        let crystal = Crystal {
            lattice: Lattice::new(
                Vector3::new(a, 0.0, 0.0),
                Vector3::new(0.0, a, 0.0),
                Vector3::new(0.0, 0.0, a),
            ),
            atoms: vec![
                // Na+ at corners and face centers (FCC)
                Atom::new(11, [0.0, 0.0, 0.0]),
                Atom::new(11, [0.5, 0.5, 0.0]),
                Atom::new(11, [0.5, 0.0, 0.5]),
                Atom::new(11, [0.0, 0.5, 0.5]),
                // Cl- at edge centers and body center (FCC shifted by a/2)
                Atom::new(17, [0.5, 0.0, 0.0]),
                Atom::new(17, [0.0, 0.5, 0.0]),
                Atom::new(17, [0.0, 0.0, 0.5]),
                Atom::new(17, [0.5, 0.5, 0.5]),
            ],
        };

        // Mock pseudopotentials with Z_val = +1 (Na+) and -1 (Cl-)
        let pp_na = mock_pp("Na", 1.0);
        let pp_cl = mock_pp("Cl", -1.0);

        let e = ewald_energy(&crystal, &[&pp_na, &pp_cl]);
        // Expected: E = -M × e² × 4 (ion pairs) / (a/2)
        // M = 1.747_565, e² = 14.3997 eV·Å, nearest-neighbor distance = a/2
        let e_expected = -1.747_565 * E2 / (a / 2.0) * 4.0;
        let relative_err = ((e - e_expected) / e_expected).abs();
        assert!(
            relative_err < 0.01,
            "Ewald energy: {e:.6} eV, expected {e_expected:.6} eV (err={relative_err:.4})"
        );
    }

    fn mock_pp(element: &str, z_valence: f64) -> PseudopotentialData {
        PseudopotentialData {
            element: element.into(),
            z_valence,
            l_max: 0,
            r_grid: vec![],
            rab: vec![],
            v_local: vec![],
            beta_projectors: vec![],
            dij: vec![],
            n_projectors: 0,
            rho_atom: vec![],
            core_charge: vec![],
            has_nlcc: false,
        }
    }

    #[test]
    fn test_ewald_zero_charges() {
        let a = 5.0;
        let crystal = Crystal {
            lattice: Lattice::new(
                Vector3::new(a, 0.0, 0.0),
                Vector3::new(0.0, a, 0.0),
                Vector3::new(0.0, 0.0, a),
            ),
            atoms: vec![Atom::new(14, [0.0, 0.0, 0.0])],
        };
        let pp = mock_pp("Si", 0.0);
        let e = ewald_energy(&crystal, &[&pp]);
        assert!(e.abs() < 1e-12, "Zero charges should give zero energy: {e}");
    }

    #[test]
    fn test_ewald_single_atom() {
        // Single atom: only self-energy and background, no real-space pairs
        let a = 5.0;
        let crystal = Crystal {
            lattice: Lattice::new(
                Vector3::new(a, 0.0, 0.0),
                Vector3::new(0.0, a, 0.0),
                Vector3::new(0.0, 0.0, a),
            ),
            atoms: vec![Atom::new(14, [0.0, 0.0, 0.0])],
        };
        let pp = mock_pp("Si", 4.0);
        let e = ewald_energy(&crystal, &[&pp]);
        assert!(e.is_finite(), "Single atom Ewald should be finite: {e}");
        assert!(e < 0.0, "Single atom Ewald should be negative: {e}");
    }

    #[test]
    fn test_ewald_anisotropic_cell() {
        // Slab-like geometry: very short c-axis
        let crystal = Crystal {
            lattice: Lattice::new(
                Vector3::new(10.0, 0.0, 0.0),
                Vector3::new(0.0, 10.0, 0.0),
                Vector3::new(0.0, 0.0, 2.0),
            ),
            atoms: vec![
                Atom::new(14, [0.0, 0.0, 0.0]),
                Atom::new(14, [0.5, 0.5, 0.5]),
            ],
        };
        let pp = mock_pp("Si", 4.0);
        let e = ewald_energy(&crystal, &[&pp]);
        assert!(e.is_finite(), "Anisotropic cell energy should be finite: {e}");
    }
}
