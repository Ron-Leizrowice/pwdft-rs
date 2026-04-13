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
}

pub type Result<T> = std::result::Result<T, PwdftError>;
