//! Small dense linear solver used by the DIIS/Broyden mixers.
//!
//! Both Anderson (DIIS) and modified Broyden mixing need to solve a small
//! (typically 2–8 × 2–8) linear system every iteration to extract the
//! mixing coefficients. At these sizes the cost is negligible; we delegate
//! to `faer`'s partial-pivoting LU so we don't have to maintain a
//! hand-rolled Gauss-elimination primitive. Bit-equivalent to textbook
//! Gauss-elim with partial pivoting for well-conditioned systems.
//!
//! When the factorization surfaces a tiny pivot (|U[i,i]| < 1e-15), the
//! overlap matrix is effectively singular; we fall back to uniform
//! coefficients so the mixer degrades gracefully rather than emitting
//! NaNs downstream.
//!
//! See MODR's "Flagged for follow-up" and GLUS (FLUP) for the history.

use faer::Mat;
use faer::prelude::Solve;

/// Solve A x = b for small systems via `faer` partial-pivoting LU.
///
/// `a` is a row-major flattened `n × n` matrix. Returns `x` such that
/// `A x = b`, or a uniform coefficient vector `[1/(n+1); n]` if `A` is
/// numerically singular (smallest `|U[i,i]| < 1e-15` after LU).
pub(super) fn solve_linear_system(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }

    let mat_a = Mat::from_fn(n, n, |i, j| a[i * n + j]);
    let lu = mat_a.partial_piv_lu();

    // Mirror the hand-rolled Gauss-elim singular-pivot guard. After
    // partial-pivoting LU, the reduced-matrix pivots appear on the
    // diagonal of U (up to row permutation); a tiny |U[i,i]| marks a
    // rank-deficient system.
    let u = lu.U();
    for i in 0..n {
        if u[(i, i)].abs() < 1e-15 {
            log::warn!(
                "Anderson mixer: singular overlap matrix (pivot={:.2e}), \
                 falling back to uniform coefficients",
                u[(i, i)]
            );
            return vec![1.0 / (n + 1) as f64; n];
        }
    }

    let mut rhs = Mat::from_fn(n, 1, |i, _| b[i]);
    lu.solve_in_place(rhs.as_mut());

    (0..n).map(|i| rhs[(i, 0)]).collect()
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "ERR2 § Phase 0: in-src test modules are allowed to panic")]
mod tests {
    use super::*;

    #[test]
    fn test_solve_linear_system() {
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let b = vec![5.0, 7.0];
        let x = solve_linear_system(&a, &b, 2);
        assert!((x[0] - 1.6).abs() < 1e-10);
        assert!((x[1] - 1.8).abs() < 1e-10);
    }
}
