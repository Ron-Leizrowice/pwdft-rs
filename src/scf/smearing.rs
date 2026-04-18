//! Occupation number smearing and electronic entropy.
//!
//! Supports Fermi-Dirac, Gaussian, Methfessel-Paxton (order 1), and
//! Marzari-Vanderbilt cold smearing. Each scheme has an occupation function
//! and an entropy formula for the Mermin free energy F = E - TS.
//!
//! All internal functions return occupation in [0, 1] per state.
//! The public `occupation()` multiplies by `spin_factor` (2/nspin):
//! - nspin=1: spin_factor=2 (each state holds 2 electrons)
//! - nspin=2: spin_factor=1 (each state holds 1 electron)

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

/// Available smearing schemes for occupation numbers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SmearingScheme {
    #[default]
    FermiDirac,
    Gaussian,
    MethfesselPaxton,
    Cold,
    /// Fixed occupations (no smearing, insulator mode).
    /// Behaves like Fermi-Dirac when passed through occupation functions.
    Fixed,
}

// ---------------------------------------------------------------------------
// Public occupation API
// ---------------------------------------------------------------------------

/// Compute occupation for a given smearing scheme.
///
/// `spin_factor`: 2.0/nspin (2.0 for unpolarized, 1.0 for spin-polarized).
/// Returns occupation in [0, spin_factor].
#[must_use]
pub fn occupation(
    scheme: SmearingScheme,
    energy: f64,
    fermi_energy: f64,
    sigma: f64,
    spin_factor: f64,
) -> f64 {
    let f01 = occupation_01(scheme, energy, fermi_energy, sigma);
    f01 * spin_factor
}

/// Find the Fermi energy by bisection such that:
/// N_electrons = Σ_{n,k} w_k f(ε_{n,k}, E_F, σ)
///
/// `eigenvalues_flat`: eigenvalues for all (spin, k) pairs.
/// `kpoint_weights_flat`: k-point weight for each (spin, k) pair.
/// `spin_factor`: 2.0/nspin.
#[must_use]
pub fn find_fermi_energy(
    eigenvalues_flat: &[Vec<f64>],
    kpoint_weights_flat: &[f64],
    n_electrons: f64,
    sigma: f64,
    scheme: SmearingScheme,
    spin_factor: f64,
) -> f64 {
    let mut e_min = f64::INFINITY;
    let mut e_max = f64::NEG_INFINITY;
    for evs in eigenvalues_flat {
        for &e in evs {
            e_min = e_min.min(e);
            e_max = e_max.max(e);
        }
    }
    e_min -= 10.0 * sigma.max(0.1);
    e_max += 10.0 * sigma.max(0.1);

    for _ in 0..200 {
        let e_mid = 0.5 * (e_min + e_max);
        let n: f64 = eigenvalues_flat
            .iter()
            .zip(kpoint_weights_flat.iter())
            .map(|(evs, &w)| {
                evs.iter()
                    .map(|&e| w * occupation(scheme, e, e_mid, sigma, spin_factor))
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
// Internal: occupation in [0, 1] per state
// ---------------------------------------------------------------------------

/// Occupation in [0, 1] for any scheme (before spin factor).
fn occupation_01(scheme: SmearingScheme, energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    match scheme {
        SmearingScheme::FermiDirac | SmearingScheme::Fixed => {
            fermi_dirac_01(energy, fermi_energy, sigma)
        }
        SmearingScheme::Gaussian => gaussian_01(energy, fermi_energy, sigma),
        SmearingScheme::MethfesselPaxton => methfessel_paxton_01(energy, fermi_energy, sigma),
        SmearingScheme::Cold => cold_01(energy, fermi_energy, sigma),
    }
}

/// Fermi-Dirac occupation (before spin factor).
///
/// f(ε) = 1 / (1 + exp(x))  where x = (ε - E_F) / σ
///
/// At T=0 (σ→0): step function θ(E_F - ε), with f(E_F) = 1/2.
/// Overflow-protected for |x| > 40.
fn fermi_dirac_01(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 {
        return if energy < fermi_energy { 1.0 }
               else if (energy - fermi_energy).abs() < 1e-12 { 0.5 }
               else { 0.0 };
    }
    let x = (energy - fermi_energy) / sigma;
    if x > 40.0 { 0.0 } else if x < -40.0 { 1.0 } else { 1.0 / (1.0 + x.exp()) }
}

/// Gaussian smearing occupation (before spin factor).
///
/// f(ε) = erfc(x) / 2  where x = (ε - E_F) / σ
fn gaussian_01(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 {
        return if energy < fermi_energy { 1.0 }
               else if (energy - fermi_energy).abs() < 1e-12 { 0.5 }
               else { 0.0 };
    }
    let x = (energy - fermi_energy) / sigma;
    puruspe::erfc(x) / 2.0
}

/// Methfessel-Paxton order-1 occupation (before spin factor).
///
/// f(ε) = erfc(x)/2 - (x/2) exp(-x²) / √π
///
/// Reference: Methfessel & Paxton, Phys. Rev. B 40, 3616 (1989).
fn methfessel_paxton_01(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 {
        return if energy < fermi_energy { 1.0 }
               else if (energy - fermi_energy).abs() < 1e-12 { 0.5 }
               else { 0.0 };
    }
    let x = (energy - fermi_energy) / sigma;
    let f0 = puruspe::erfc(x) / 2.0;
    let gauss = (-x * x).exp() / PI.sqrt();
    (0.5 * x).mul_add(-gauss, f0)
}

/// Marzari-Vanderbilt "cold" smearing occupation (before spin factor).
///
/// f(ε) = (1/2) erfc(x + 1/√2) + exp(-(x + 1/√2)²) / √(2π)
///
/// Designed to give positive-definite entropy. Argument shifted by 1/√2
/// so that f(E_F) = 1/2 exactly.
///
/// Reference: Marzari, Vanderbilt, De Vita, Payne, Phys. Rev. Lett. 82, 3296 (1999).
fn cold_01(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 {
        return if energy < fermi_energy { 1.0 }
               else if (energy - fermi_energy).abs() < 1e-12 { 0.5 }
               else { 0.0 };
    }
    let x = (energy - fermi_energy) / sigma;
    let sq2_inv = std::f64::consts::FRAC_1_SQRT_2;
    let arg = x + sq2_inv;
    0.5f64.mul_add(puruspe::erfc(arg), (1.0 / (2.0 * PI).sqrt()) * (-arg * arg).exp())
}

// ---------------------------------------------------------------------------
// Entropy: T*S for each smearing scheme
// ---------------------------------------------------------------------------

/// Compute the electronic entropy contribution T*S in eV.
///
/// `spin_factor`: 2.0/nspin.
#[must_use]
pub fn entropy_ts(
    eigenvalues_flat: &[Vec<f64>],
    kpoint_weights_flat: &[f64],
    fermi_energy: f64,
    sigma: f64,
    scheme: SmearingScheme,
    spin_factor: f64,
) -> f64 {
    if sigma < 1e-15 {
        return 0.0;
    }

    let mut s = 0.0;
    for (evs, &w) in eigenvalues_flat.iter().zip(kpoint_weights_flat.iter()) {
        for &e in evs {
            let x = (e - fermi_energy) / sigma;
            s += w * entropy_weight(scheme, x);
        }
    }

    s * sigma * spin_factor
}

/// Per-state entropy weight s(x) for reduced variable x = (ε - E_F)/σ.
///
/// The total entropy is TS = σ × spin_factor × Σ_{n,k} w_k s(x_{n,k}).
///
/// Formulas by scheme:
/// - Fermi-Dirac: s = -[f ln f + (1-f) ln(1-f)]
/// - Gaussian:    s = exp(-x²) / √π
/// - Methfessel-Paxton: s = (1/2 - x²) exp(-x²) / √π
/// - Cold:        s = (x + 1/√2) exp(-(x + 1/√2)²) / √π
fn entropy_weight(scheme: SmearingScheme, x: f64) -> f64 {
    match scheme {
        SmearingScheme::FermiDirac | SmearingScheme::Fixed => {
            if x.abs() > 30.0 {
                return 0.0;
            }
            let f = 1.0 / (1.0 + x.exp());
            let f = f.clamp(1e-30, 1.0 - 1e-30);
            -(f * f.ln() + (1.0 - f) * (1.0 - f).ln())
        }
        SmearingScheme::Gaussian => {
            (-x * x).exp() / PI.sqrt()
        }
        SmearingScheme::MethfesselPaxton => {
            x.mul_add(-x, 0.5) * (-x * x).exp() / PI.sqrt()
        }
        SmearingScheme::Cold => {
            let sq2_inv = std::f64::consts::FRAC_1_SQRT_2;
            let arg = x + sq2_inv;
            arg * (-arg * arg).exp() / PI.sqrt()
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "ERR2 § Phase 0: in-src test modules are allowed to panic")]
mod tests {
    use super::*;
    use approx::relative_eq;

    fn single_kpoint_weights() -> Vec<f64> {
        vec![1.0]
    }

    // -----------------------------------------------------------------------
    // Spin factor tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_occupation_spin_factor_2() {
        // nspin=1: spin_factor=2, max occupation = 2.0
        let f = occupation(SmearingScheme::FermiDirac, -1.0, 0.0, 0.01, 2.0);
        assert!((f - 2.0).abs() < 1e-10, "Expected 2.0, got {f}");
        let f = occupation(SmearingScheme::FermiDirac, 1.0, 0.0, 0.01, 2.0);
        assert!(f.abs() < 1e-10, "Expected ~0, got {f}");
        let f = occupation(SmearingScheme::FermiDirac, 0.0, 0.0, 0.01, 2.0);
        assert!((f - 1.0).abs() < 1e-10, "Expected 1.0, got {f}");
    }

    #[test]
    fn test_occupation_spin_factor_1() {
        // nspin=2: spin_factor=1, max occupation = 1.0
        let f = occupation(SmearingScheme::FermiDirac, -1.0, 0.0, 0.01, 1.0);
        assert!((f - 1.0).abs() < 1e-10, "Expected 1.0, got {f}");
        let f = occupation(SmearingScheme::FermiDirac, 1.0, 0.0, 0.01, 1.0);
        assert!(f.abs() < 1e-10, "Expected ~0, got {f}");
        let f = occupation(SmearingScheme::FermiDirac, 0.0, 0.0, 0.01, 1.0);
        assert!((f - 0.5).abs() < 1e-10, "Expected 0.5, got {f}");
    }

    #[test]
    fn test_all_schemes_spin_factor_consistency() {
        // For all schemes: occupation(spin_factor=2) = 2 * occupation(spin_factor=1)
        for scheme in [
            SmearingScheme::FermiDirac,
            SmearingScheme::Gaussian,
            SmearingScheme::MethfesselPaxton,
            SmearingScheme::Cold,
        ] {
            for &e in &[-0.5, -0.1, 0.0, 0.1, 0.5] {
                let f1 = occupation(scheme, e, 0.0, 0.1, 1.0);
                let f2 = occupation(scheme, e, 0.0, 0.1, 2.0);
                assert!(
                    relative_eq!(f2, 2.0 * f1, epsilon = 1e-12),
                    "{scheme:?} at e={e}: f2={f2}, 2*f1={}", 2.0 * f1
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Fermi energy search with both spin factors
    // -----------------------------------------------------------------------

    #[test]
    fn test_find_fermi_energy_nspin1() {
        // 4 electrons, spin_factor=2 (nspin=1): 2 bands fully occupied
        let eigenvalues = vec![vec![-2.0, -1.0, 1.0, 2.0]];
        let weights = single_kpoint_weights();
        let ef = find_fermi_energy(&eigenvalues, &weights, 4.0, 0.05, SmearingScheme::FermiDirac, 2.0);
        assert!(ef > -1.5 && ef < 1.5, "E_F={ef} outside range");
        let n: f64 = eigenvalues[0].iter()
            .map(|&e| occupation(SmearingScheme::FermiDirac, e, ef, 0.05, 2.0))
            .sum();
        assert!(relative_eq!(n, 4.0, epsilon = 1e-6), "N={n} != 4.0");
    }

    #[test]
    fn test_find_fermi_energy_nspin2() {
        // 4 electrons total, spin_factor=1 (nspin=2): need eigenvalues for both spins
        // Spin up: [-2, -1, 1, 2], spin down: [-2, -1, 1, 2]
        // With spin_factor=1, each state holds 1 electron → need 4 occupied states total
        let eigenvalues = vec![
            vec![-2.0, -1.0, 1.0, 2.0],  // spin up at k=0
            vec![-2.0, -1.0, 1.0, 2.0],  // spin down at k=0
        ];
        let weights = vec![1.0, 1.0]; // same k-point for both spins
        let ef = find_fermi_energy(&eigenvalues, &weights, 4.0, 0.05, SmearingScheme::FermiDirac, 1.0);
        assert!(ef > -1.5 && ef < 1.5, "E_F={ef} outside range");
        let n: f64 = eigenvalues.iter().zip(weights.iter())
            .flat_map(|(evs, &w)| evs.iter().map(move |&e| w * occupation(SmearingScheme::FermiDirac, e, ef, 0.05, 1.0)))
            .sum();
        assert!(relative_eq!(n, 4.0, epsilon = 1e-6), "N={n} != 4.0");
    }

    #[test]
    fn test_find_fermi_nspin1_vs_nspin2_unpolarized() {
        // For identical spin channels, nspin=2 should give same E_F as nspin=1
        let evs = vec![-2.0, -1.0, 1.0, 2.0];
        let n_el = 4.0;
        let sigma = 0.05;

        // nspin=1
        let ef1 = find_fermi_energy(
            std::slice::from_ref(&evs), &[1.0], n_el, sigma, SmearingScheme::FermiDirac, 2.0,
        );

        // nspin=2 with identical channels
        let ef2 = find_fermi_energy(
            &[evs.clone(), evs], &[1.0, 1.0], n_el, sigma, SmearingScheme::FermiDirac, 1.0,
        );

        assert!(
            (ef1 - ef2).abs() < 1e-10,
            "nspin=1 E_F={ef1} != nspin=2 E_F={ef2}"
        );
    }

    // -----------------------------------------------------------------------
    // Entropy with spin factor
    // -----------------------------------------------------------------------

    #[test]
    fn test_entropy_zero_at_zero_sigma() {
        let eigenvalues = vec![vec![-1.0, 0.0, 1.0]];
        let weights = single_kpoint_weights();
        for scheme in [
            SmearingScheme::FermiDirac,
            SmearingScheme::Gaussian,
            SmearingScheme::MethfesselPaxton,
            SmearingScheme::Cold,
        ] {
            let ts = entropy_ts(&eigenvalues, &weights, 0.0, 0.0, scheme, 2.0);
            assert!(ts.abs() < 1e-15, "{scheme:?}: TS should be 0 at sigma=0, got {ts}");
        }
    }

    #[test]
    fn test_entropy_positive_for_fermi_dirac() {
        let eigenvalues = vec![vec![-0.5, -0.1, 0.1, 0.5]];
        let weights = single_kpoint_weights();
        let ef = find_fermi_energy(&eigenvalues, &weights, 4.0, 0.1, SmearingScheme::FermiDirac, 2.0);
        let ts = entropy_ts(&eigenvalues, &weights, ef, 0.1, SmearingScheme::FermiDirac, 2.0);
        assert!(ts > 0.0, "F-D entropy should be positive, got {ts}");
    }

    #[test]
    fn test_entropy_nspin1_vs_nspin2_unpolarized() {
        // Entropy for identical spin channels with nspin=2 should equal nspin=1
        let evs = vec![-0.5, -0.1, 0.1, 0.5];
        let ef = 0.0;
        let sigma = 0.1;

        let ts1 = entropy_ts(
            std::slice::from_ref(&evs), &[1.0], ef, sigma, SmearingScheme::FermiDirac, 2.0,
        );
        let ts2 = entropy_ts(
            &[evs.clone(), evs], &[1.0, 1.0], ef, sigma, SmearingScheme::FermiDirac, 1.0,
        );

        assert!(
            (ts1 - ts2).abs() < 1e-10,
            "nspin=1 TS={ts1} != nspin=2 TS={ts2}"
        );
    }

    // -----------------------------------------------------------------------
    // Zero temperature
    // -----------------------------------------------------------------------

    #[test]
    #[allow(clippy::float_cmp)] // Zero-temp occupations return exact values by design
    fn test_zero_temp_all_schemes() {
        for scheme in [
            SmearingScheme::FermiDirac,
            SmearingScheme::Gaussian,
            SmearingScheme::MethfesselPaxton,
            SmearingScheme::Cold,
        ] {
            // nspin=1
            assert_eq!(occupation(scheme, -1.0, 0.0, 0.0, 2.0), 2.0);
            assert_eq!(occupation(scheme, 1.0, 0.0, 0.0, 2.0), 0.0);
            assert_eq!(occupation(scheme, 0.0, 0.0, 0.0, 2.0), 1.0);
            // nspin=2
            assert_eq!(occupation(scheme, -1.0, 0.0, 0.0, 1.0), 1.0);
            assert_eq!(occupation(scheme, 1.0, 0.0, 0.0, 1.0), 0.0);
            assert_eq!(occupation(scheme, 0.0, 0.0, 0.0, 1.0), 0.5);
        }
    }

}
