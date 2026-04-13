use nalgebra::DMatrix;
use num_complex::Complex64;

/// Result of diagonalizing a Hermitian matrix.
pub struct EigenResult {
    /// Eigenvalues in ascending order.
    pub eigenvalues: Vec<f64>,
    /// Eigenvectors as columns of a matrix (column i corresponds to eigenvalue i).
    pub eigenvectors: DMatrix<Complex64>,
}

/// Diagonalize a complex Hermitian matrix using LAPACK zheev.
///
/// Returns eigenvalues (ascending) and eigenvectors. Only the upper triangle
/// of `h` is referenced.
///
/// # Panics
/// Panics if the matrix is not square or if LAPACK returns an error.
pub fn diagonalize_hermitian(h: &DMatrix<Complex64>) -> EigenResult {
    let n = h.nrows();
    assert_eq!(n, h.ncols(), "matrix must be square");

    if n == 0 {
        return EigenResult {
            eigenvalues: vec![],
            eigenvectors: DMatrix::zeros(0, 0),
        };
    }

    // LAPACK expects column-major storage — nalgebra stores column-major by default.
    // Clone into a mutable flat buffer.
    let mut a: Vec<Complex64> = Vec::with_capacity(n * n);
    for col in 0..n {
        for row in 0..n {
            a.push(h[(row, col)]);
        }
    }

    let n_i32 = n as i32;
    let mut w = vec![0.0f64; n];
    let mut rwork = vec![0.0f64; (3 * n).max(1) - 2 + 1]; // at least 3n-2

    // Query optimal workspace
    let mut work_query = vec![Complex64::new(0.0, 0.0); 1];
    let mut info = 0i32;
    unsafe {
        lapack::zheev(
            b'V',
            b'U',
            n_i32,
            &mut a,
            n_i32,
            &mut w,
            &mut work_query,
            -1,
            &mut rwork,
            &mut info,
        );
    }
    assert_eq!(info, 0, "zheev workspace query failed with info={info}");

    let lwork = work_query[0].re as i32;
    let mut work = vec![Complex64::new(0.0, 0.0); lwork as usize];

    // Actual diagonalization
    unsafe {
        lapack::zheev(
            b'V',
            b'U',
            n_i32,
            &mut a,
            n_i32,
            &mut w,
            &mut work,
            lwork,
            &mut rwork,
            &mut info,
        );
    }
    assert_eq!(info, 0, "zheev diagonalization failed with info={info}");

    // Reconstruct eigenvector matrix from flat column-major buffer
    let eigenvectors = DMatrix::from_column_slice(n, n, &a);

    EigenResult {
        eigenvalues: w,
        eigenvectors,
    }
}

/// Diagonalize and return only the lowest `n_bands` eigenvalues/eigenvectors.
///
/// This still performs full diagonalization internally (Phase 4 will add
/// iterative solvers), but only returns the requested subset.
pub fn diagonalize_lowest(h: &DMatrix<Complex64>, n_bands: usize) -> EigenResult {
    let full = diagonalize_hermitian(h);
    let n = n_bands.min(full.eigenvalues.len());

    let eigenvalues = full.eigenvalues[..n].to_vec();
    let eigenvectors = full.eigenvectors.columns(0, n).into_owned();

    EigenResult {
        eigenvalues,
        eigenvectors,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

    #[test]
    fn test_real_symmetric_2x2() {
        // [3  1]
        // [1  3]  → eigenvalues 2, 4
        let h = DMatrix::from_row_slice(
            2,
            2,
            &[
                Complex64::new(3.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(3.0, 0.0),
            ],
        );
        let result = diagonalize_hermitian(&h);
        assert!(relative_eq!(result.eigenvalues[0], 2.0, epsilon = 1e-10));
        assert!(relative_eq!(result.eigenvalues[1], 4.0, epsilon = 1e-10));
    }

    #[test]
    fn test_complex_hermitian_2x2() {
        let i = Complex64::i();
        // [1   i ]
        // [-i  1 ]  → eigenvalues 0, 2
        let h = DMatrix::from_row_slice(
            2,
            2,
            &[
                Complex64::new(1.0, 0.0),
                i,
                -i,
                Complex64::new(1.0, 0.0),
            ],
        );
        let result = diagonalize_hermitian(&h);
        assert!(relative_eq!(result.eigenvalues[0], 0.0, epsilon = 1e-10));
        assert!(relative_eq!(result.eigenvalues[1], 2.0, epsilon = 1e-10));
    }

    #[test]
    fn test_eigenvectors_orthonormal() {
        let i = Complex64::i();
        let h = DMatrix::from_row_slice(
            3,
            3,
            &[
                Complex64::new(2.0, 0.0),
                Complex64::new(1.0, 0.0) + i,
                Complex64::new(0.0, 0.0),
                Complex64::new(1.0, 0.0) - i,
                Complex64::new(3.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(2.0, 0.0),
            ],
        );
        let result = diagonalize_hermitian(&h);

        // Check orthonormality: V^H V = I
        let vh = result.eigenvectors.adjoint();
        let prod = &vh * &result.eigenvectors;
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    relative_eq!(prod[(i, j)].norm(), expected, epsilon = 1e-10),
                    "V^H V [{i},{j}] = {}, expected {expected}",
                    prod[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_diagonalize_lowest() {
        let h = DMatrix::from_row_slice(
            3,
            3,
            &[
                Complex64::new(1.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(5.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(9.0, 0.0),
            ],
        );
        let result = diagonalize_lowest(&h, 2);
        assert_eq!(result.eigenvalues.len(), 2);
        assert!(relative_eq!(result.eigenvalues[0], 1.0, epsilon = 1e-10));
        assert!(relative_eq!(result.eigenvalues[1], 5.0, epsilon = 1e-10));
        assert_eq!(result.eigenvectors.ncols(), 2);
    }
}
