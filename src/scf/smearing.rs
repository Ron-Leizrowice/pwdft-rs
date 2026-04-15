//! Occupation number smearing and electronic entropy.
//!
//! Supports Fermi-Dirac, Gaussian, Methfessel-Paxton (order 1), and
//! Marzari-Vanderbilt cold smearing. Each scheme has an occupation function
//! and an entropy formula for the Mermin free energy F = E - TS.

use std::f64::consts::PI;

use crate::kpoints::KPoint;

/// Available smearing schemes for occupation numbers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SmearingScheme {
    #[default]
    FermiDirac,
    Gaussian,
    MethfesselPaxton,
    Cold,
}

// ---------------------------------------------------------------------------
// Occupation functions (all include spin factor of 2)
// ---------------------------------------------------------------------------

/// Compute occupation for any smearing scheme.
pub fn occupation(scheme: SmearingScheme, energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    match scheme {
        SmearingScheme::FermiDirac => fermi_dirac(energy, fermi_energy, sigma),
        SmearingScheme::Gaussian => gaussian(energy, fermi_energy, sigma),
        SmearingScheme::MethfesselPaxton => methfessel_paxton(energy, fermi_energy, sigma),
        SmearingScheme::Cold => cold(energy, fermi_energy, sigma),
    }
}

/// Fermi-Dirac: f(x) = 2 / (1 + exp(x)), x = (ε - E_F) / σ
pub fn fermi_dirac(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 {
        return zero_temp(energy, fermi_energy);
    }
    let x = (energy - fermi_energy) / sigma;
    if x > 40.0 {
        0.0
    } else if x < -40.0 {
        2.0
    } else {
        2.0 / (1.0 + x.exp())
    }
}

/// Gaussian: f(x) = erfc(x), x = (ε - E_F) / σ
fn gaussian(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 {
        return zero_temp(energy, fermi_energy);
    }
    let x = (energy - fermi_energy) / sigma;
    puruspe::erfc(x) // erfc goes from 2 to 0 — already includes spin factor
}

/// Methfessel-Paxton order 1: f₁(x) = erfc(x)/2 - x·exp(-x²)/√π, then ×2
/// Note: occupations CAN be negative for states well above E_F.
fn methfessel_paxton(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 {
        return zero_temp(energy, fermi_energy);
    }
    let x = (energy - fermi_energy) / sigma;
    let f0 = puruspe::erfc(x) / 2.0;
    let gauss = (-x * x).exp() / PI.sqrt();
    // M-P order-1 correction: A₁·H₁(x)·exp(-x²) where A₁ = -1/(4√π), H₁ = 2x
    let f1 = f0 - 0.5 * x * gauss;
    2.0 * f1
}

/// Marzari-Vanderbilt cold smearing:
/// f(x) = 1 + erf(x + 1/√2) + exp(-(x+1/√2)²)/√(2π), then adjust for convention
fn cold(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 {
        return zero_temp(energy, fermi_energy);
    }
    let x = (energy - fermi_energy) / sigma;
    let sq2_inv = std::f64::consts::FRAC_1_SQRT_2;
    let arg = x + sq2_inv;
    // erfc goes from 2 (x→-∞, occupied) to 0 (x→+∞, empty)
    // Cold smearing: f = erfc(x + 1/√2)/2 + exp(-(x+1/√2)²)/√(2π)
    let f_half = 0.5 * puruspe::erfc(arg)
        + (1.0 / (2.0 * PI).sqrt()) * (-arg * arg).exp();
    2.0 * f_half // spin factor
}

fn zero_temp(energy: f64, fermi_energy: f64) -> f64 {
    if energy < fermi_energy {
        2.0
    } else if (energy - fermi_energy).abs() < 1e-12 {
        1.0
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Fermi energy search
// ---------------------------------------------------------------------------

/// Find the Fermi energy by bisection such that:
/// N_electrons = Σ_{n,k} w_k f(ε_{n,k}, E_F, σ)
pub fn find_fermi_energy(
    eigenvalues: &[Vec<f64>],
    kpoints: &[KPoint],
    n_electrons: f64,
    sigma: f64,
    scheme: SmearingScheme,
) -> f64 {
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

    for _ in 0..200 {
        let e_mid = 0.5 * (e_min + e_max);
        let n: f64 = eigenvalues
            .iter()
            .zip(kpoints.iter())
            .map(|(evs, kp)| {
                evs.iter()
                    .map(|&e| kp.weight * occupation(scheme, e, e_mid, sigma))
                    .sum::<f64>()
            })
            .sum();
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

// ---------------------------------------------------------------------------
// Entropy: T*S for each smearing scheme
// ---------------------------------------------------------------------------

/// Compute the electronic entropy contribution T*S in eV.
///
/// For the Mermin free energy: F = E - T*S.
/// The sigma→0 extrapolated energy: E₀ = (E + F) / 2.
pub fn entropy_ts(
    eigenvalues: &[Vec<f64>],
    kpoints: &[KPoint],
    fermi_energy: f64,
    sigma: f64,
    scheme: SmearingScheme,
) -> f64 {
    if sigma < 1e-15 {
        return 0.0;
    }

    let mut s = 0.0;
    for (evs, kp) in eigenvalues.iter().zip(kpoints.iter()) {
        for &e in evs {
            let x = (e - fermi_energy) / sigma;
            s += kp.weight * entropy_weight(scheme, x);
        }
    }

    s * sigma
}

/// Per-state entropy weight for a given scheme and reduced variable x = (ε-E_F)/σ.
fn entropy_weight(scheme: SmearingScheme, x: f64) -> f64 {
    match scheme {
        SmearingScheme::FermiDirac => {
            // S = -2 [f/2 ln(f/2) + (1-f/2) ln(1-f/2)] where f = 2/(1+exp(x))
            // Guard against overflow: for |x| > 40, entropy is negligible
            if x.abs() > 30.0 {
                return 0.0;
            }
            let f_half = 1.0 / (1.0 + x.exp());
            let f_half = f_half.clamp(1e-30, 1.0 - 1e-30);
            -2.0 * (f_half * f_half.ln() + (1.0 - f_half) * (1.0 - f_half).ln())
        }
        SmearingScheme::Gaussian => {
            // S = (2/√π) exp(-x²)
            2.0 * (-x * x).exp() / PI.sqrt()
        }
        SmearingScheme::MethfesselPaxton => {
            // S = (2/√π) (1/2 - x²) exp(-x²) [includes M-P order-1 correction]
            2.0 * (0.5 - x * x) * (-x * x).exp() / PI.sqrt()
        }
        SmearingScheme::Cold => {
            // S = (2/√π) (x + 1/√2) exp(-(x + 1/√2)²)  [simplified]
            let sq2_inv = std::f64::consts::FRAC_1_SQRT_2;
            let arg = x + sq2_inv;
            2.0 * arg * (-arg * arg).exp() / PI.sqrt()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;
    use nalgebra::Vector3;

    fn single_kpoint() -> Vec<KPoint> {
        vec![KPoint {
            k: Vector3::zeros(),
            weight: 1.0,
            label: None,
        }]
    }

    #[test]
    fn test_fermi_dirac_limits() {
        assert!((fermi_dirac(-1.0, 0.0, 0.01) - 2.0).abs() < 1e-10);
        assert!(fermi_dirac(1.0, 0.0, 0.01).abs() < 1e-10);
        assert!((fermi_dirac(0.0, 0.0, 0.01) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gaussian_limits() {
        // Well below → 2, well above → 0, at E_F → 1
        assert!((gaussian(-1.0, 0.0, 0.01) - 2.0).abs() < 1e-6);
        assert!(gaussian(1.0, 0.0, 0.01).abs() < 1e-6);
        assert!((gaussian(0.0, 0.0, 0.01) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_methfessel_paxton_limits() {
        // Well below → 2, at E_F → 1
        assert!((methfessel_paxton(-2.0, 0.0, 0.1) - 2.0).abs() < 0.01);
        assert!((methfessel_paxton(0.0, 0.0, 0.1) - 1.0).abs() < 1e-10);
        // M-P CAN produce negative occupations well above E_F
        let f_high = methfessel_paxton(0.5, 0.0, 0.1);
        // This is expected to be close to 0 or slightly negative
        assert!(f_high < 0.1, "M-P should be near 0 above E_F, got {f_high}");
    }

    #[test]
    fn test_cold_smearing_limits() {
        assert!((cold(-2.0, 0.0, 0.1) - 2.0).abs() < 0.01);
        assert!(cold(2.0, 0.0, 0.1).abs() < 0.01);
        // Cold smearing occupations are always non-negative
        for x_10 in -50..50 {
            let x = x_10 as f64 * 0.1;
            assert!(
                cold(x, 0.0, 0.1) >= -1e-10,
                "Cold smearing should be non-negative at x={x}"
            );
        }
    }

    #[test]
    fn test_find_fermi_energy_all_schemes() {
        let eigenvalues = vec![vec![-2.0, -1.0, 1.0, 2.0]];
        let kpoints = single_kpoint();

        for scheme in [
            SmearingScheme::FermiDirac,
            SmearingScheme::Gaussian,
            SmearingScheme::MethfesselPaxton,
            SmearingScheme::Cold,
        ] {
            let ef = find_fermi_energy(&eigenvalues, &kpoints, 4.0, 0.05, scheme);
            assert!(
                ef > -1.5 && ef < 1.5,
                "{scheme:?}: Fermi energy {ef} outside expected range"
            );
            // Verify electron count
            let n: f64 = eigenvalues[0]
                .iter()
                .map(|&e| occupation(scheme, e, ef, 0.05))
                .sum();
            assert!(
                relative_eq!(n, 4.0, epsilon = 1e-6),
                "{scheme:?}: electron count {n} != 4.0"
            );
        }
    }

    #[test]
    fn test_entropy_zero_at_zero_sigma() {
        let eigenvalues = vec![vec![-1.0, 0.0, 1.0]];
        let kpoints = single_kpoint();
        for scheme in [
            SmearingScheme::FermiDirac,
            SmearingScheme::Gaussian,
            SmearingScheme::MethfesselPaxton,
            SmearingScheme::Cold,
        ] {
            let ts = entropy_ts(&eigenvalues, &kpoints, 0.0, 0.0, scheme);
            assert!(
                ts.abs() < 1e-15,
                "{scheme:?}: entropy should be 0 at sigma=0, got {ts}"
            );
        }
    }

    #[test]
    fn test_entropy_positive_for_fermi_dirac() {
        // F-D entropy is always positive (- Σ f ln f ≥ 0)
        let eigenvalues = vec![vec![-0.5, -0.1, 0.1, 0.5]];
        let kpoints = single_kpoint();
        let ef = find_fermi_energy(&eigenvalues, &kpoints, 4.0, 0.1, SmearingScheme::FermiDirac);
        let ts = entropy_ts(
            &eigenvalues, &kpoints, ef, 0.1, SmearingScheme::FermiDirac,
        );
        assert!(ts > 0.0, "F-D entropy should be positive, got {ts}");
    }

    #[test]
    fn test_zero_temp_occupations() {
        // All schemes should agree at zero temperature
        for scheme in [
            SmearingScheme::FermiDirac,
            SmearingScheme::Gaussian,
            SmearingScheme::MethfesselPaxton,
            SmearingScheme::Cold,
        ] {
            assert_eq!(occupation(scheme, -1.0, 0.0, 0.0), 2.0);
            assert_eq!(occupation(scheme, 1.0, 0.0, 0.0), 0.0);
            assert_eq!(occupation(scheme, 0.0, 0.0, 0.0), 1.0);
        }
    }
}
