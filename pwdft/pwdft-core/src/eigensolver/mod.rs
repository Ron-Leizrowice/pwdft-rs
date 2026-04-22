//! Eigensolvers for the Kohn-Sham Hamiltonian.
//!
//! Each SCF iteration diagonalizes the Hermitian Hamiltonian at every
//! k-point. Two backends are available:
//!
//! - [`dense`] — full `faer::SelfAdjointEigen` decomposition (all n eigenpairs, LAPACK-equivalent
//!   O(n³)); also home of the WFRX subspace-rotation warm-start path used between SCF iterations.
//! - [`iterative`] — implicitly-restarted Arnoldi / Krylov-Schur partial solver that computes only
//!   the lowest `n_bands` eigenpairs.
//!
//! Pick between them via [`EigensolverKind`] on
//! [`crate::scf::ScfParams`]. Dense is the default and the only fully
//! validated path at present.

pub mod dense;
pub mod iterative;

pub use dense::EigenResult;

/// Which Hermitian-eigensolver backend the SCF loop should use for the
/// Kohn-Sham Hamiltonian.
///
/// `Dense` is the reference `faer::SelfAdjointEigen` full decomposition
/// (LAPACK-equivalent O(n³)) and is the default and the only fully
/// validated path. `Iterative` uses faer's implicitly-restarted Arnoldi
/// / Krylov-Schur partial solver
/// ([`iterative::diagonalize_lowest_iterative`]), which computes only
/// the lowest `n_bands` eigenpairs via a shift-and-flip of the spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EigensolverKind {
    /// Full dense Hermitian eigendecomposition. Default.
    #[default]
    Dense,
    /// Partial Hermitian eigensolver returning only the lowest `n_bands`
    /// eigenpairs.
    ///
    /// Opt-in via `scf.eigensolver: iterative` in YAML. Experimental —
    /// correctness and performance are still under investigation:
    ///
    /// - On a realistic Si Kohn-Sham Hamiltonian at `n_pw = 725`, single-shot iterative is ~0.48×
    ///   the wall-time of `Dense` (i.e. slower). Earlier projections of a 3–10× speedup were built
    ///   from synthetic matrices and do not survive contact with real clustered/degenerate spectra.
    /// - The iterative dispatch path does not yet consume the subspace-rotation warm start used by
    ///   `Dense`, so SCF wall-time parity depends on an end-to-end benchmark that has not been run.
    /// - Size-independent `n_request` padding can drop 3-fold-degenerate valence clusters at larger
    ///   `n_pw`.
    ///
    /// Prefer `Dense` until the end-to-end warm-started SCF benchmark
    /// lands. See `proposals/ITEV-faer-partial-eigen.md` for the full
    /// status and the open follow-ups.
    Iterative,
}
