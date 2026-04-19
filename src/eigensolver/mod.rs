//! Eigensolvers for the Kohn-Sham Hamiltonian.
//!
//! Each SCF iteration diagonalizes the Hermitian Hamiltonian at every
//! k-point. Two backends are available:
//!
//! - [`dense`] — full `faer::SelfAdjointEigen` decomposition (all n
//!   eigenpairs, LAPACK-equivalent O(n³)); also home of the WFRX
//!   subspace-rotation warm-start path used between SCF iterations.
//! - [`iterative`] — implicitly-restarted Arnoldi / Krylov-Schur
//!   partial solver that computes only the lowest `n_bands` eigenpairs.
//!
//! Pick between them via [`EigensolverKind`] on
//! [`crate::scf::ScfParams`]. Dense is the default and the only fully
//! validated path at present.

pub mod dense;
pub mod iterative;

pub use dense::EigenResult;

/// Which dense-Hamiltonian eigensolver backend the SCF loop should use.
///
/// `Dense` is the reference `faer::SelfAdjointEigen` full decomposition
/// (LAPACK-equivalent O(n³)). `Iterative` uses faer's implicitly-restarted
/// Arnoldi / Krylov-Schur partial solver ([`iterative::diagonalize_lowest_iterative`])
/// which computes only the lowest `n_bands` eigenpairs.
///
/// The iterative solver is typically 3-10× faster at `n_pw ≥ 200` and
/// 10-50× faster at `n_pw ≥ 700`, but falls back transparently to dense
/// when convergence within the restart budget fails. See proposal ITEV.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EigensolverKind {
    /// Full dense Hermitian eigendecomposition (default for now).
    #[default]
    Dense,
    /// Iterative partial Hermitian eigensolver (ITEV).
    Iterative,
}
