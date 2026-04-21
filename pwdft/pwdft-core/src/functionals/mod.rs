pub mod lda;
pub mod pbe;

/// Minimum electron density for XC evaluation (e/ų).
pub const RHO_FLOOR: f64 = 1e-30;
