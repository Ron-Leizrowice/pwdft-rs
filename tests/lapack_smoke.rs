use approx::relative_eq;
use nalgebra::DMatrix;
use nalgebra_lapack::SymmetricEigen;
use num_complex::Complex64;

/// Verify that nalgebra-lapack links correctly against Apple Accelerate
/// by diagonalizing a known 3×3 real symmetric matrix.
#[test]
fn test_lapack_real_symmetric_eigen() {
    // Symmetric real matrix:
    //  [2  1  0]
    //  [1  3  1]
    //  [0  1  2]
    // Characteristic eq: (2-λ)(λ²-5λ+4) = 0 → eigenvalues: 1, 2, 4
    let m = DMatrix::from_row_slice(3, 3, &[2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0]);

    let eigen = SymmetricEigen::new(m);
    let mut eigenvalues: Vec<f64> = eigen.eigenvalues.iter().cloned().collect();
    eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let expected = [1.0, 2.0, 4.0];

    for (got, want) in eigenvalues.iter().zip(expected.iter()) {
        assert!(
            relative_eq!(got, want, epsilon = 1e-10),
            "eigenvalue mismatch: got {got}, expected {want}"
        );
    }
}

/// Test complex Hermitian matrix diagonalization via direct LAPACK zheev call.
/// This is the codepath needed for plane-wave DFT Hamiltonians.
#[test]
fn test_lapack_complex_hermitian_zheev() {
    // 3×3 Hermitian matrix:
    //  [ 2    1+i   0  ]
    //  [ 1-i   3    1  ]
    //  [ 0     1    2  ]
    let i = Complex64::i();
    let zero = Complex64::new(0.0, 0.0);
    let one = Complex64::new(1.0, 0.0);
    let two = Complex64::new(2.0, 0.0);
    let three = Complex64::new(3.0, 0.0);

    // LAPACK uses column-major order
    // Column 0: [2, 1-i, 0]
    // Column 1: [1+i, 3, 1]
    // Column 2: [0, 1, 2]
    let mut a = vec![
        two, one - i, zero, // col 0
        one + i, three, one, // col 1
        zero, one, two,     // col 2
    ];
    let n = 3i32;
    let mut w = vec![0.0f64; 3]; // eigenvalues output

    // Query optimal workspace size
    let mut work = vec![Complex64::new(0.0, 0.0); 1];
    let mut rwork = vec![0.0f64; 3 * 3 - 2];
    let mut info = 0i32;
    let lwork = -1i32;

    unsafe {
        lapack::zheev(
            b'V', b'U', n, &mut a, n, &mut w, &mut work, lwork, &mut rwork, &mut info,
        );
    }
    assert_eq!(info, 0, "zheev workspace query failed");

    let lwork = work[0].re as i32;
    work.resize(lwork as usize, Complex64::new(0.0, 0.0));

    unsafe {
        lapack::zheev(
            b'V', b'U', n, &mut a, n, &mut w, &mut work, lwork, &mut rwork, &mut info,
        );
    }
    assert_eq!(info, 0, "zheev diagonalization failed");

    // Eigenvalues are returned in ascending order
    // For this matrix, verify they sum to trace = 2 + 3 + 2 = 7
    let trace: f64 = w.iter().sum();
    assert!(
        relative_eq!(trace, 7.0, epsilon = 1e-10),
        "eigenvalue sum {trace} != trace 7.0"
    );

    // All eigenvalues should be real and positive for this matrix
    for &ev in &w {
        assert!(ev > -1e-10, "unexpected negative eigenvalue: {ev}");
    }
}
