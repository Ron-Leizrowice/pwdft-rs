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

/// Diagonalize a complex Hermitian matrix using faer (pure Rust, multi-threaded).
///
/// Drop-in replacement for `diagonalize_hermitian`. Returns eigenvalues in
/// ascending order and eigenvectors as columns of a DMatrix (for compatibility
/// with existing code). The O(n²) DMatrix↔faer::Mat conversion is negligible
/// vs the O(n³) eigendecomposition.
pub fn diagonalize_hermitian_faer(h: &DMatrix<Complex64>) -> EigenResult {
    let n = h.nrows();
    assert_eq!(n, h.ncols(), "matrix must be square");

    if n == 0 {
        return EigenResult {
            eigenvalues: vec![],
            eigenvectors: DMatrix::zeros(0, 0),
        };
    }

    // Convert nalgebra DMatrix → faer Mat (same complex type, different container)
    let h_faer = faer::Mat::<Complex64>::from_fn(n, n, |r, c| h[(r, c)]);

    // faer eigenvalues are in nondecreasing order (same as LAPACK zheev)
    let decomp = h_faer
        .self_adjoint_eigen(faer::Side::Lower)
        .expect("faer eigendecomposition failed");

    // Extract eigenvalues (Hermitian → real, stored as complex with im=0)
    let s = decomp.S();
    let s_col = s.column_vector();
    let eigenvalues: Vec<f64> = (0..n).map(|i| s_col[i].re).collect();

    // Convert eigenvectors back to nalgebra DMatrix
    let u = decomp.U();
    let eigenvectors = DMatrix::from_fn(n, n, |r, c| u[(r, c)]);

    EigenResult {
        eigenvalues,
        eigenvectors,
    }
}

/// Diagonalize with faer and return only the lowest `n_bands` eigenvalues/eigenvectors.
pub fn diagonalize_lowest_faer(h: &DMatrix<Complex64>, n_bands: usize) -> EigenResult {
    let full = diagonalize_hermitian_faer(h);
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
    fn test_faer_matches_lapack_2x2() {
        // Same 2×2 matrices as the LAPACK tests
        let h = DMatrix::from_row_slice(
            2, 2,
            &[
                Complex64::new(3.0, 0.0), Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0), Complex64::new(3.0, 0.0),
            ],
        );
        let lapack = diagonalize_hermitian(&h);
        let faer = diagonalize_hermitian_faer(&h);
        for (i, (&l, &f)) in lapack.eigenvalues.iter().zip(faer.eigenvalues.iter()).enumerate() {
            assert!(
                (l - f).abs() < 1e-12,
                "eigenvalue {i}: lapack={l}, faer={f}"
            );
        }
    }

    #[test]
    fn test_faer_matches_lapack_complex_hermitian() {
        let i = Complex64::i();
        let h = DMatrix::from_row_slice(
            3, 3,
            &[
                Complex64::new(2.0, 0.0), Complex64::new(1.0, 0.0) + i, Complex64::new(0.0, 0.0),
                Complex64::new(1.0, 0.0) - i, Complex64::new(3.0, 0.0), Complex64::new(1.0, 0.0),
                Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0),
            ],
        );
        let lapack = diagonalize_hermitian(&h);
        let faer = diagonalize_hermitian_faer(&h);
        for (idx, (&l, &f)) in lapack.eigenvalues.iter().zip(faer.eigenvalues.iter()).enumerate() {
            assert!(
                (l - f).abs() < 1e-12,
                "eigenvalue {idx}: lapack={l}, faer={f}"
            );
        }

        // Check faer eigenvectors satisfy Hv = λv
        for idx in 0..3 {
            let v = faer.eigenvectors.column(idx);
            let hv = &h * &v;
            let lv = &v * Complex64::new(faer.eigenvalues[idx], 0.0);
            let residual: f64 = hv.iter().zip(lv.iter()).map(|(a, b)| (a - b).norm_sqr()).sum::<f64>().sqrt();
            assert!(
                residual < 1e-10,
                "faer eigenvector {idx}: |Hv - λv| = {residual:.2e}"
            );
        }
    }

    #[test]
    fn test_faer_matches_lapack_realistic_hamiltonian() {
        // Build a real Si Hamiltonian at Γ — this is the matrix size
        // that matters for actual SCF calculations
        use crate::{basis::BasisSet, crystal::Lattice, hamiltonian, potential::nonlocal::NonlocalPotential};
        use nalgebra::Vector3;

        let a = 5.431;
        let lattice = Lattice::new(
            a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        );
        let crystal = crate::crystal::Crystal {
            lattice: lattice.clone(),
            atoms: vec![
                crate::crystal::Atom::new(14, [0.0, 0.0, 0.0]),
                crate::crystal::Atom::new(14, [0.25, 0.25, 0.25]),
            ],
        };
        let basis = BasisSet::new(&lattice, 200.0);
        let pp = crate::pseudopotential::load(
            &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Si.UPF"),
        ).unwrap();
        let k = Vector3::zeros();

        // Build full Hamiltonian with nonlocal potential
        let mut h = hamiltonian::build_kinetic(&basis, &k);
        let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]);
        vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k);

        let n = basis.len();
        eprintln!("Comparing LAPACK vs faer on {n}×{n} Si Hamiltonian at Γ");

        let lapack_result = diagonalize_hermitian(&h);
        let faer_result = diagonalize_hermitian_faer(&h);

        // Eigenvalues must match to high precision
        let mut max_diff = 0.0_f64;
        for (idx, (&l, &f)) in lapack_result.eigenvalues.iter().zip(faer_result.eigenvalues.iter()).enumerate() {
            let diff = (l - f).abs();
            max_diff = max_diff.max(diff);
            assert!(
                diff < 1e-10,
                "eigenvalue {idx}/{n}: lapack={l:.10}, faer={f:.10}, diff={diff:.2e}"
            );
        }
        eprintln!("Max eigenvalue difference: {max_diff:.2e}");

        // Verify faer eigenvectors are unitary: V^H V = I
        let vh = faer_result.eigenvectors.adjoint();
        let prod = &vh * &faer_result.eigenvectors;
        for i_idx in 0..n.min(20) {
            for j_idx in 0..n.min(20) {
                let expected = if i_idx == j_idx { 1.0 } else { 0.0 };
                let got = prod[(i_idx, j_idx)].norm();
                assert!(
                    (got - expected).abs() < 1e-10,
                    "V^H V [{i_idx},{j_idx}] = {got:.6}, expected {expected}"
                );
            }
        }
        eprintln!("Faer eigenvectors: unitary to 1e-10");
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
