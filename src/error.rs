use thiserror::Error;

#[derive(Debug, Error)]
pub enum PwdftError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("SCF did not converge after {iterations} iterations (delta = {delta:.2e})")]
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
}

pub type Result<T> = std::result::Result<T, PwdftError>;
