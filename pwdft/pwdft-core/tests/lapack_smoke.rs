#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: integration tests are allowed to panic"
)]

use approx::relative_eq;
use num_complex::Complex64;

/// Verify that faer eigendecomposition works correctly on a known
/// 3×3 real symmetric matrix.
#[test]
fn test_faer_real_symmetric_eigen() {
    // Symmetric real matrix:
    //  [2  1  0]
    //  [1  3  1]
    //  [0  1  2]
    // Eigenvalues: 1, 2, 4
    let h = faer::Mat::from_fn(3, 3, |r, c| {
        let data = [2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0];
        data[r * 3 + c]
    });

    let decomp = h.self_adjoint_eigen(faer::Side::Lower).unwrap();
    let s = decomp.S().column_vector();
    let mut eigenvalues: Vec<f64> = (0..3).map(|i| s[i]).collect();
    eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let expected = [1.0, 2.0, 4.0];

    for (got, want) in eigenvalues.iter().zip(expected.iter()) {
        assert!(
            relative_eq!(got, want, epsilon = 1e-10),
            "eigenvalue mismatch: got {got}, expected {want}"
        );
    }
}

/// Test complex Hermitian matrix eigendecomposition via faer.
/// This is the codepath used for plane-wave DFT Hamiltonians.
#[test]
fn test_faer_complex_hermitian_eigen() {
    // 3×3 Hermitian matrix:
    //  [ 2    1+i   0  ]
    //  [ 1-i   3    1  ]
    //  [ 0     1    2  ]
    let i = Complex64::i();
    let one = Complex64::new(1.0, 0.0);
    let two = Complex64::new(2.0, 0.0);
    let three = Complex64::new(3.0, 0.0);
    let zero = Complex64::new(0.0, 0.0);

    let data = [two, one + i, zero, one - i, three, one, zero, one, two];
    let h = faer::Mat::from_fn(3, 3, |r, c| data[r * 3 + c]);

    let decomp = h.self_adjoint_eigen(faer::Side::Lower).unwrap();
    let s = decomp.S().column_vector();
    let eigenvalues: Vec<f64> = (0..3).map(|i| s[i].re).collect();

    // Eigenvalues are returned in ascending order
    // Verify they sum to trace = 2 + 3 + 2 = 7
    let trace: f64 = eigenvalues.iter().sum();
    assert!(
        relative_eq!(trace, 7.0, epsilon = 1e-10),
        "eigenvalue sum {trace} != trace 7.0"
    );

    // All eigenvalues should be real and positive for this matrix
    for &ev in &eigenvalues {
        assert!(ev > -1e-10, "unexpected negative eigenvalue: {ev}");
    }
}
