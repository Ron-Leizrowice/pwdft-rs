//! Crate-wide error type and result alias.
//!
//! [`PwdftError`] is the single fallible boundary for pwdft-core — I/O
//! failures, YAML / UPF parse errors, SCF non-convergence, eigensolver
//! breakdowns, GPU failures, and the [`PwdftError::NotImplemented`]
//! guard that keeps syntactically valid but unimplemented YAML options
//! from silently picking up the wrong physics path.
//!
//! [`Result<T>`] is shorthand for `std::result::Result<T, PwdftError>`.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum PwdftError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("SCF did not converge after {iterations} iterations (delta = {delta:.2e}); raise `scf.max_iter` or tighten mixing/smearing in the input YAML")]
    ConvergenceFailure { iterations: usize, delta: f64 },

    #[error("invalid input: {0}")]
    InvalidInput(String),

    #[error("missing pseudopotential for element {0}")]
    MissingPseudopotential(String),

    #[error("parse error: {0}")]
    Parse(String),

    #[error("eigendecomposition failed for {size}x{size} matrix: {detail}")]
    Eigensolver { size: usize, detail: String },

    #[error("GPU error: {0}")]
    Gpu(String),

    /// A configuration path that is syntactically valid but has no working
    /// implementation yet. Used to guard against silent-wrong-physics when
    /// the YAML parser accepts an option that the SCF pipeline does not
    /// dispatch on.
    #[error("{what} is not yet implemented")]
    NotImplemented { what: String },

    /// A user-supplied parameter failed a range or shape check. `name`
    /// is a compile-time-constant parameter label (a YAML key or
    /// `ScfParams` field name). `reason` interpolates the offending
    /// value when useful. Prefer this variant over the catch-all
    /// `InvalidInput` for per-parameter validation errors.
    #[error("invalid parameter {name}: {reason}")]
    InvalidParam {
        name: &'static str,
        reason: String,
    },

    /// A crystal, lattice, or k-point structural precondition failed.
    /// Used for whole-input shape errors (empty atom list, empty k-point
    /// list, degenerate lattice) rather than per-parameter range checks.
    #[error("invalid crystal input: {reason}")]
    InvalidCrystal { reason: &'static str },

    /// YAML or CLI input referenced an element symbol or atomic number
    /// that is not in the `crate::atoms` periodic table.
    #[error("unknown element symbol: {symbol}")]
    UnknownElement { symbol: String },

    /// A pseudopotential file parsed structurally but failed a
    /// physical-validity check (for example, negative angular momentum
    /// or a missing required block). Distinct from `Parse`, which is
    /// reserved for syntax-level UPF errors: here the XML parsed
    /// fine, but the *value* it carried was out of physical range.
    /// `file` identifies the file path or tag context; `reason`
    /// describes the specific invariant that was violated. The field
    /// is named `file` rather than `source` to avoid thiserror's
    /// reserved-field treatment of `source` (which implies an inner
    /// `std::error::Error` cause).
    #[error("invalid pseudopotential {file}: {reason}")]
    InvalidPseudopotential { file: String, reason: String },
}

pub type Result<T> = std::result::Result<T, PwdftError>;

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: test bodies are allowed to panic"
)]
mod tests {
    use super::*;

    #[test]
    fn invalid_param_fmt() {
        let err = PwdftError::InvalidParam {
            name: "n_bands",
            reason: "must be > 0".into(),
        };
        assert_eq!(err.to_string(), "invalid parameter n_bands: must be > 0");
    }

    #[test]
    fn invalid_param_fmt_with_value() {
        // Exercises the pattern used by PARAM-cluster callers where
        // the reason string interpolates the offending numeric value.
        let err = PwdftError::InvalidParam {
            name: "mixing_beta",
            reason: format!("must be in (0, 1], got {}", 1.5_f64),
        };
        assert_eq!(
            err.to_string(),
            "invalid parameter mixing_beta: must be in (0, 1], got 1.5"
        );
    }

    #[test]
    fn invalid_crystal_fmt() {
        let err = PwdftError::InvalidCrystal {
            reason: "at least one atom is required",
        };
        assert_eq!(
            err.to_string(),
            "invalid crystal input: at least one atom is required"
        );
    }

    #[test]
    fn unknown_element_fmt() {
        let err = PwdftError::UnknownElement {
            symbol: "Xz".into(),
        };
        assert_eq!(err.to_string(), "unknown element symbol: Xz");
    }

    #[test]
    fn invalid_pseudopotential_fmt() {
        let err = PwdftError::InvalidPseudopotential {
            file: "PP_BETA.1".into(),
            reason: "angular_momentum must be non-negative (got -1)".into(),
        };
        assert_eq!(
            err.to_string(),
            "invalid pseudopotential PP_BETA.1: angular_momentum must be non-negative (got -1)"
        );
    }

    #[test]
    fn display_dispatches_through_format_macro() {
        // Confirms the `Display` impl generated by `thiserror` is
        // wired up for the new variants — `format!("{}", err)` is the
        // path most call sites (and the CLI main-binary error path)
        // take.
        let err = PwdftError::InvalidParam {
            name: "conv_threshold",
            reason: "must be positive".into(),
        };
        let rendered = format!("{err}");
        assert_eq!(
            rendered,
            "invalid parameter conv_threshold: must be positive"
        );
    }
}
