//! Physical constants and unit conversions.
//!
//! The engine works internally in eV for energies, Å for lengths, and
//! e/Å³ for densities; this module holds the conversion factors from
//! CODATA atomic units (Hartree, Rydberg, Bohr) used to normalize
//! pseudopotential tables at the UPF boundary, plus SI constants for
//! deriving quantities like `ℏ²/2m` in the engine's native units.
//!
//! Also defines the small numerical floors
//! ([`G2_ZERO_THRESHOLD`], [`G_ZERO_THRESHOLD`], [`RHO_FLOOR`]) that gate
//! Coulomb and XC evaluations against divide-by-zero.

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
/// Threshold for treating |G| as zero (branch on G=0 in reciprocal-space
/// transforms that compare against `|G|` directly rather than `|G|²`).
///
/// Distinct from [`G2_ZERO_THRESHOLD`]: the two are numerically different
/// floors because they compare against different quantities — using
/// `1e-12` as a `|G|²` floor is equivalent to using `1e-6` as a `|G|`
/// floor. Callers should pick the one that matches the quantity they
/// already have in hand.
pub const G_ZERO_THRESHOLD: f64 = 1e-12;
/// Minimum electron density for XC evaluation (e/ų).
pub const RHO_FLOOR: f64 = 1e-30;

// SI constants
pub const H_SI: f64 = 6.626_070_15e-34; // J·s (exact, SI definition)
pub const HBAR_SI: f64 = H_SI / (2.0 * std::f64::consts::PI); // J·s
pub const M_E: f64 = 9.109_383_701_5e-31; // kg

// Conversion factors
pub const EV_PER_J: f64 = 6.241_509_074e18; // 1 J = this many eV
pub const ANG_PER_M: f64 = 1e10; // 1 m = 1e10 Å

// Derived: ħ²/2m in eV·Å²
pub const HBAR2_OVER_2M: f64 =
    HBAR_SI * HBAR_SI / (2.0 * M_E) * EV_PER_J * ANG_PER_M * ANG_PER_M;
