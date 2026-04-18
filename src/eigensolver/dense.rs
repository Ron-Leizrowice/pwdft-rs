use num_complex::Complex64;

use crate::error::{PwdftError, Result};

/// Result of diagonalizing a Hermitian matrix.
pub struct EigenResult {
    /// Eigenvalues in ascending order.
    pub eigenvalues: Vec<f64>,
    /// Eigenvectors as columns of a matrix (column i corresponds to eigenvalue i).
    pub eigenvectors: faer::Mat<Complex64>,
}

/// Full Hermitian eigendecomposition of H via faer.
///
/// Solves Hψ = εψ for all eigenvalues and eigenvectors.
/// Returns eigenvalues in ascending order (ε₁ ≤ ε₂ ≤ ... ≤ εₙ).
///
/// Uses faer's `self_adjoint_eigen` (dense, O(n³) LAPACK-equivalent).
/// Only the lower triangle of H is read.
///
/// # Errors
/// Returns `PwdftError::Eigensolver` if the matrix is not square or if faer
/// fails to compute the eigendecomposition.
pub fn diagonalize_hermitian(h: &faer::Mat<Complex64>) -> Result<EigenResult> {
    let n = h.nrows();
    if n != h.ncols() {
        return Err(PwdftError::Eigensolver {
            size: n,
            detail: format!("matrix is not square: {}x{}", n, h.ncols()),
        });
    }

    if n == 0 {
        return Ok(EigenResult {
            eigenvalues: vec![],
            eigenvectors: faer::Mat::zeros(0, 0),
        });
    }

    let decomp = h
        .self_adjoint_eigen(faer::Side::Lower)
        .map_err(|_| PwdftError::Eigensolver {
            size: n,
            detail: "faer returned no eigendecomposition".into(),
        })?;

    let s_col = decomp.S().column_vector();
    let eigenvalues: Vec<f64> = (0..n).map(|i| s_col[i].re).collect();
    let eigenvectors = decomp.U().to_owned();

    Ok(EigenResult {
        eigenvalues,
        eigenvectors,
    })
}

/// Diagonalize and return only the lowest `n_bands` eigenvalues/eigenvectors.
///
/// # Errors
/// Returns `PwdftError::Eigensolver` if the eigendecomposition fails.
pub fn diagonalize_lowest(h: &faer::Mat<Complex64>, n_bands: usize) -> Result<EigenResult> {
    let full = diagonalize_hermitian(h)?;
    let n = n_bands.min(full.eigenvalues.len());

    let eigenvalues = full.eigenvalues[..n].to_vec();
    let eigenvectors = full.eigenvectors.subcols(0, n).to_owned();

    Ok(EigenResult {
        eigenvalues,
        eigenvectors,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "ERR2 § Phase 0: in-src test modules are allowed to panic")]
mod tests {
    use super::*;
    use approx::relative_eq;

    fn mat_from_rows(n: usize, data: &[Complex64]) -> faer::Mat<Complex64> {
        faer::Mat::from_fn(n, n, |r, c| data[r * n + c])
    }

    #[test]
    fn test_real_symmetric_2x2() {
        // [3  1]
        // [1  3]  → eigenvalues 2, 4
        let h = mat_from_rows(2, &[
            Complex64::new(3.0, 0.0), Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0), Complex64::new(3.0, 0.0),
        ]);
        let result = diagonalize_hermitian(&h).unwrap();
        assert!(relative_eq!(result.eigenvalues[0], 2.0, epsilon = 1e-10));
        assert!(relative_eq!(result.eigenvalues[1], 4.0, epsilon = 1e-10));
    }

    #[test]
    fn test_complex_hermitian_2x2() {
        let i = Complex64::i();
        // [1   i ]
        // [-i  1 ]  → eigenvalues 0, 2
        let h = mat_from_rows(2, &[
            Complex64::new(1.0, 0.0), i,
            -i, Complex64::new(1.0, 0.0),
        ]);
        let result = diagonalize_hermitian(&h).unwrap();
        assert!(relative_eq!(result.eigenvalues[0], 0.0, epsilon = 1e-10));
        assert!(relative_eq!(result.eigenvalues[1], 2.0, epsilon = 1e-10));
    }

    #[test]
    fn test_eigenvectors_orthonormal() {
        let i = Complex64::i();
        let h = mat_from_rows(3, &[
            Complex64::new(2.0, 0.0), Complex64::new(1.0, 0.0) + i, Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0) - i, Complex64::new(3.0, 0.0), Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0),
        ]);
        let result = diagonalize_hermitian(&h).unwrap();

        // Check Hv = λv for each eigenpair
        for idx in 0..3 {
            let v = result.eigenvectors.col(idx);
            let mut hv = faer::Col::<Complex64>::zeros(3);
            for r in 0..3 {
                let mut sum = Complex64::new(0.0, 0.0);
                for c in 0..3 {
                    sum += h[(r, c)] * v[c];
                }
                hv[r] = sum;
            }
            let lambda = Complex64::new(result.eigenvalues[idx], 0.0);
            let residual: f64 = (0..3).map(|r| (hv[r] - lambda * v[r]).norm_sqr()).sum::<f64>().sqrt();
            assert!(
                residual < 1e-10,
                "eigenvector {idx}: |Hv - λv| = {residual:.2e}"
            );
        }
    }

    #[test]
    fn test_diagonalize_lowest() {
        let h = mat_from_rows(3, &[
            Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(5.0, 0.0), Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(9.0, 0.0),
        ]);
        let result = diagonalize_lowest(&h, 2).unwrap();
        assert_eq!(result.eigenvalues.len(), 2);
        assert!(relative_eq!(result.eigenvalues[0], 1.0, epsilon = 1e-10));
        assert!(relative_eq!(result.eigenvalues[1], 5.0, epsilon = 1e-10));
        assert_eq!(result.eigenvectors.ncols(), 2);
    }
}
