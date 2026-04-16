pub const PI: f64 = std::f64::consts::PI;

// Atomic unit conversions
/// Hartree to electronvolt.
pub const HA_TO_EV: f64 = 27.211_386_245_988;
/// Rydberg to electronvolt.
pub const RY_TO_EV: f64 = 13.605_693_122_994;
/// Bohr radius in Ångströms.
pub const BOHR_TO_ANG: f64 = 0.529_177_210_903;
/// Bohr³ in ų (volume conversion).
pub const BOHR3_TO_ANG3: f64 = BOHR_TO_ANG * BOHR_TO_ANG * BOHR_TO_ANG;
/// Coulomb constant e² in eV·Å (Gaussian units).
pub const E2_COULOMB: f64 = 14.399_645_351_950_548;

// Numerical thresholds
/// Threshold for treating |G|² as zero (skip G=0 in Coulomb sums).
pub const G2_ZERO_THRESHOLD: f64 = 1e-12;
/// Minimum electron density for XC evaluation (e/ų).
pub const RHO_FLOOR: f64 = 1e-20;

// SI constants
pub const H_SI: f64 = 6.626_070_15e-34; // J·s (exact, SI definition)
pub const HBAR_SI: f64 = H_SI / (2.0 * PI); // J·s
pub const M_E: f64 = 9.109_383_701_5e-31; // kg

// Conversion factors
pub const EV_PER_J: f64 = 6.241_509_074e18; // 1 J = this many eV
pub const ANG_PER_M: f64 = 1e10; // 1 m = 1e10 Å

// Derived: ħ²/2m in eV·Å²
// = (HBAR_SI² / (2 * M_E)) [J·m²] * EV_PER_J [eV/J] * ANG_PER_M² [Å²/m²]
pub const HBAR2_OVER_2M: f64 =
    HBAR_SI * HBAR_SI / (2.0 * M_E) * EV_PER_J * ANG_PER_M * ANG_PER_M;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

    #[test]
    fn test_constants() {
        assert!(relative_eq!(HBAR2_OVER_2M, 3.81, epsilon = 1e-2), "HBAR2_OVER_2M = {HBAR2_OVER_2M}, expected 3.81");
    }
}
