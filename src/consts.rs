pub const PI: f64 = std::f64::consts::PI;

// SI constants
pub const H_SI: f64 = 6.62607015e-34; // J·s (exact, SI definition)
pub const HBAR_SI: f64 = H_SI / (2.0 * PI); // J·s
pub const M_E: f64 = 9.1093837015e-31; // kg

// Conversion factors
pub const EV_PER_J: f64 = 6.241509074e18; // 1 J = this many eV
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
