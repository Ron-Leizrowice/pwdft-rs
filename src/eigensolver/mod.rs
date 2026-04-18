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
