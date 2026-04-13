//! Occupation number smearing for metallic or finite-temperature calculations.
//!
//! Fermi-Dirac smearing: f(ε) = 1 / (1 + exp((ε - E_F) / σ))

use crate::kpoints::KPoint;

/// Fermi-Dirac occupation function.
///
/// Returns occupation (0 to 1) for a single state.
/// For spin-unpolarized: multiply by 2 for the total occupation.
pub fn fermi_dirac(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 {
        // Zero-temperature limit
        return if energy < fermi_energy {
            2.0
        } else if (energy - fermi_energy).abs() < 1e-12 {
            1.0
        } else {
            0.0
        };
    }

    let x = (energy - fermi_energy) / sigma;
    if x > 40.0 {
        0.0
    } else if x < -40.0 {
        2.0 // spin factor of 2
    } else {
        2.0 / (1.0 + x.exp()) // spin factor of 2
    }
}

/// Find the Fermi energy by bisection such that:
///
/// N_electrons = Σ_{n,k} w_k f(ε_{n,k}, E_F, σ)
pub fn find_fermi_energy(
    eigenvalues: &[Vec<f64>],
    kpoints: &[KPoint],
    n_electrons: f64,
    sigma: f64,
) -> f64 {
    // Find bounds
    let mut e_min = f64::INFINITY;
    let mut e_max = f64::NEG_INFINITY;
    for evs in eigenvalues {
        for &e in evs {
            e_min = e_min.min(e);
            e_max = e_max.max(e);
        }
    }
    e_min -= 10.0 * sigma.max(0.1);
    e_max += 10.0 * sigma.max(0.1);

    // Bisection
    for _ in 0..200 {
        let e_mid = 0.5 * (e_min + e_max);
        let n = electron_count(eigenvalues, kpoints, e_mid, sigma);
        if n < n_electrons {
            e_min = e_mid;
        } else {
            e_max = e_mid;
        }
        if (e_max - e_min).abs() < 1e-14 {
            break;
        }
    }

    0.5 * (e_min + e_max)
}

/// Count electrons for a given Fermi energy.
fn electron_count(
    eigenvalues: &[Vec<f64>],
    kpoints: &[KPoint],
    fermi_energy: f64,
    sigma: f64,
) -> f64 {
    eigenvalues
        .iter()
        .zip(kpoints.iter())
        .map(|(evs, kp)| {
            evs.iter()
                .map(|&e| kp.weight * fermi_dirac(e, fermi_energy, sigma))
                .sum::<f64>()
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_fermi_dirac_limits() {
        // Well below Fermi energy: occupation = 2
        assert!((fermi_dirac(-1.0, 0.0, 0.01) - 2.0).abs() < 1e-10);
        // Well above: occupation = 0
        assert!(fermi_dirac(1.0, 0.0, 0.01).abs() < 1e-10);
        // At Fermi energy: occupation = 1
        assert!((fermi_dirac(0.0, 0.0, 0.01) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_find_fermi_energy() {
        // 4 electrons, 4 bands at a single k-point
        let eigenvalues = vec![vec![-2.0, -1.0, 1.0, 2.0]];
        let kpoints = vec![KPoint {
            k: Vector3::zeros(),
            weight: 1.0,
            label: None,
        }];
        let ef = find_fermi_energy(&eigenvalues, &kpoints, 4.0, 0.001);
        // With 4 electrons and spin factor 2: first two bands fully occupied
        // Fermi energy should be between -1.0 and 1.0
        assert!(
            ef > -1.0 && ef < 1.0,
            "Fermi energy should be between -1 and 1: {ef}"
        );
    }
}
