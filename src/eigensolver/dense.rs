//! Dense Hermitian eigensolver and subspace rotation.
//!
//! [`diagonalize_lowest`] wraps faer's `self_adjoint_eigen` for the full
//! O(n³) decomposition of Hψ = εψ and returns the lowest `n_bands`
//! eigenpairs. [`diagonalize_subspace`] is the WFRX warm-start path: at
//! SCF iterations after the first, it projects the new Hamiltonian onto
//! the previous iteration's occupied+buffer subspace, solves the small
//! `n_sub × n_sub` eigenproblem there, and falls back to the full
//! solver if the residual `‖Hv − εv‖₂` exceeds [`WFRX_RESIDUAL_TOL`].
//!
//! Returns an [`EigenResult`] with eigenvalues sorted ascending and
//! corresponding eigenvectors as matrix columns.

use faer::linalg::matmul::matmul;
use num_complex::Complex64;

use crate::error::{PwdftError, Result};

/// Result of diagonalizing a Hermitian matrix.
pub struct EigenResult {
    /// Eigenvalues in ascending order.
    pub eigenvalues: Vec<f64>,
    /// Eigenvectors as columns of a matrix (column i corresponds to eigenvalue i).
    pub eigenvectors: faer::Mat<Complex64>,
}

/// Subspace-diagonalization residual tolerance (per eigenpair, L2 norm).
///
/// When the max residual `||H v_i − λ_i v_i||₂` across the returned
/// eigenpairs exceeds this threshold, `diagonalize_subspace` considers
/// the warm-start projection insufficient and falls back to the full
/// `diagonalize_lowest` path for that k-point. The tolerance is chosen
/// to be well below the SCF's own density/energy convergence criteria
/// (1e-6 / 1e-5) so that the subspace rotation does not degrade SCF
/// outputs.
pub const WFRX_RESIDUAL_TOL: f64 = 1e-6;

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
pub(crate) fn diagonalize_hermitian(h: &faer::Mat<Complex64>) -> Result<EigenResult> {
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

/// Compute the per-eigenpair residual `||H v_i − λ_i v_i||₂`.
///
/// A small residual (below [`WFRX_RESIDUAL_TOL`]) indicates that the pair
/// (λ_i, v_i) satisfies the eigenvalue equation to the requested accuracy;
/// a large residual means the subspace used to produce them was not
/// invariant under H and the result should be discarded in favour of a
/// full diagonalization.
///
/// The returned `Vec` has one entry per column of `eigenvectors`. It is
/// the caller's responsibility to ensure that `eigenvalues.len() ==
/// eigenvectors.ncols()`; any trailing eigenvalues beyond that column
/// count are silently ignored.
#[must_use]
pub fn eigenvalue_residuals(
    h: &faer::Mat<Complex64>,
    eigenvalues: &[f64],
    eigenvectors: &faer::Mat<Complex64>,
) -> Vec<f64> {
    let n_bands = eigenvectors.ncols().min(eigenvalues.len());
    let n = h.nrows();
    let mut out = Vec::with_capacity(n_bands);
    // Compute H V in one GEMM, then for each column subtract λ_i v_i and take
    // the L2 norm. One matmul is much cheaper than n_bands sequential mat-vec
    // calls on small n_bands × large n shapes.
    let mut hv: faer::Mat<Complex64> = faer::Mat::zeros(n, n_bands);
    let v_sub = eigenvectors.subcols(0, n_bands);
    matmul(
        hv.as_mut(),
        faer::Accum::Replace,
        h.as_ref(),
        v_sub,
        Complex64::new(1.0, 0.0),
        faer::Par::Seq,
    );
    for i in 0..n_bands {
        let lambda = Complex64::new(eigenvalues[i], 0.0);
        let mut acc = 0.0_f64;
        for row in 0..n {
            let diff = hv[(row, i)] - lambda * v_sub[(row, i)];
            acc += diff.norm_sqr();
        }
        out.push(acc.sqrt());
    }
    out
}

/// Subspace (Rayleigh–Ritz) diagonalization using a warm-start basis.
///
/// Instead of solving the full O(n³) eigendecomposition of `H`, project
/// `H` onto the subspace spanned by the columns of `v_prev` (shape `n ×
/// n_bands`), diagonalize the resulting small `n_bands × n_bands` matrix,
/// and rotate the small eigenvectors back to the full basis. For a
/// converging SCF — where `H` changes by only a small V_H / V_xc update
/// between iterations — the previous iteration's eigenvectors span a
/// subspace that is very close to invariant under the new `H`, so the
/// Rayleigh–Ritz values are excellent approximations to the true lowest
/// eigenpairs. See Payne et al., *Rev. Mod. Phys.* **64**, 1045 (1992).
///
/// **Cost:** `O(n² · n_bands)` for the projection `H_sub = V_prev^H · H ·
/// V_prev`, vs `O(n³)` for a full diagonalization. For `n = 725, n_bands
/// = 8`, the projection is ~90× cheaper than the full solve.
///
/// **Accuracy gate.** Because the subspace is only approximately invariant,
/// the returned eigenpairs are approximate. `diagonalize_subspace`
/// therefore computes the per-eigenpair residuals
/// `||H v_new − λ v_new||₂` (see [`eigenvalue_residuals`]) and falls back
/// to the exact `diagonalize_lowest` path if any residual exceeds
/// [`WFRX_RESIDUAL_TOL`]. Early SCF iterations — where the Hamiltonian
/// has changed substantially since the last solve — will typically
/// trigger the fallback; late iterations (near convergence) will not.
///
/// **Fallback paths:**
/// - `v_prev.is_none()` → full `diagonalize_lowest`.
/// - `v_prev` has the wrong shape (rows != `h.nrows()` or fewer columns
///   than `n_bands`) → full `diagonalize_lowest`.
/// - Residual gate trips → full `diagonalize_lowest`.
///
/// In all fallback paths the final result is bit-identical to what
/// `diagonalize_lowest` would have returned on its own, so enabling
/// the subspace warm-start can never degrade numerical accuracy.
///
/// # Errors
/// Returns `PwdftError::Eigensolver` if the (small) projected
/// diagonalization or the full fallback fails.
pub fn diagonalize_subspace(
    h: &faer::Mat<Complex64>,
    n_bands: usize,
    v_prev: Option<&faer::Mat<Complex64>>,
) -> Result<EigenResult> {
    let n = h.nrows();

    // Shape / nullability guards. Any failure here falls back to the
    // reference dense path so WFRX cannot introduce a crash on the first
    // SCF iteration or on a k-point whose band count changed.
    let Some(v_prev) = v_prev else {
        return diagonalize_lowest(h, n_bands);
    };
    if v_prev.nrows() != n || v_prev.ncols() < n_bands {
        return diagonalize_lowest(h, n_bands);
    }
    if n_bands == 0 {
        return Ok(EigenResult {
            eigenvalues: vec![],
            eigenvectors: faer::Mat::zeros(n, 0),
        });
    }

    let v_sub = v_prev.subcols(0, n_bands);

    // Project H into the warm-start subspace:
    //   HV := H · V_prev       (n × n_bands)
    //   H_sub := V_prev^H · HV (n_bands × n_bands, Hermitian in exact arith.)
    let mut hv: faer::Mat<Complex64> = faer::Mat::zeros(n, n_bands);
    matmul(
        hv.as_mut(),
        faer::Accum::Replace,
        h.as_ref(),
        v_sub,
        Complex64::new(1.0, 0.0),
        faer::Par::Seq,
    );
    let mut h_sub: faer::Mat<Complex64> = faer::Mat::zeros(n_bands, n_bands);
    matmul(
        h_sub.as_mut(),
        faer::Accum::Replace,
        v_sub.adjoint(),
        hv.as_ref(),
        Complex64::new(1.0, 0.0),
        faer::Par::Seq,
    );

    // Diagonalize the small matrix. `self_adjoint_eigen` reads only the
    // lower triangle, which sidesteps any numerical Hermitian-ness error
    // introduced by finite-precision arithmetic in the projection above.
    let sub_result = diagonalize_hermitian(&h_sub)?;

    // Rotate the small eigenvectors back to the full basis:
    //   V_new := V_prev · U    (n × n_bands)
    let mut v_new: faer::Mat<Complex64> = faer::Mat::zeros(n, n_bands);
    matmul(
        v_new.as_mut(),
        faer::Accum::Replace,
        v_sub,
        sub_result.eigenvectors.as_ref(),
        Complex64::new(1.0, 0.0),
        faer::Par::Seq,
    );

    // Accuracy gate. If any eigenpair residual exceeds the tolerance,
    // the warm-start subspace wasn't close enough to H-invariant and
    // we fall back to an exact full solve. See WFRX proposal § Step 3.
    let residuals = eigenvalue_residuals(h, &sub_result.eigenvalues, &v_new);
    let max_residual = residuals
        .iter()
        .copied()
        .fold(0.0_f64, f64::max);
    if max_residual > WFRX_RESIDUAL_TOL {
        return diagonalize_lowest(h, n_bands);
    }

    Ok(EigenResult {
        eigenvalues: sub_result.eigenvalues,
        eigenvectors: v_new,
    })
}

#[cfg(test)]
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

    // ---------------------------------------------------------------------
    // WFRX Phase 1 tests
    // ---------------------------------------------------------------------

    #[test]
    fn test_subspace_falls_back_when_v_prev_is_none() {
        // WFRX first-iteration path: v_prev == None must return the same
        // result as diagonalize_lowest with no degradation.
        let h = mat_from_rows(3, &[
            Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(5.0, 0.0), Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(9.0, 0.0),
        ]);
        let result = diagonalize_subspace(&h, 2, None).unwrap();
        assert_eq!(result.eigenvalues.len(), 2);
        assert!(relative_eq!(result.eigenvalues[0], 1.0, epsilon = 1e-10));
        assert!(relative_eq!(result.eigenvalues[1], 5.0, epsilon = 1e-10));
    }

    #[test]
    fn test_subspace_falls_back_on_shape_mismatch() {
        // Wrong-sized v_prev must NOT crash; falls back to full solve.
        let h = mat_from_rows(3, &[
            Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(5.0, 0.0), Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(9.0, 0.0),
        ]);
        // Wrong row count.
        let bad_rows: faer::Mat<Complex64> = faer::Mat::identity(4, 2);
        let r1 = diagonalize_subspace(&h, 2, Some(&bad_rows)).unwrap();
        assert!(relative_eq!(r1.eigenvalues[0], 1.0, epsilon = 1e-10));
        // Too few columns.
        let bad_cols: faer::Mat<Complex64> = faer::Mat::identity(3, 1);
        let r2 = diagonalize_subspace(&h, 2, Some(&bad_cols)).unwrap();
        assert!(relative_eq!(r2.eigenvalues[0], 1.0, epsilon = 1e-10));
        assert!(relative_eq!(r2.eigenvalues[1], 5.0, epsilon = 1e-10));
    }

    #[test]
    fn test_subspace_exact_when_v_prev_is_exact_eigenvectors() {
        // Proposal acceptance criterion #1: when the subspace is invariant
        // (i.e. v_prev already contains the true eigenvectors), subspace
        // diag returns the same eigenvalues as full diag to machine
        // precision.
        let i = Complex64::i();
        let h = mat_from_rows(3, &[
            Complex64::new(2.0, 0.0), Complex64::new(1.0, 0.0) + i, Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0) - i, Complex64::new(3.0, 0.0), Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0),
        ]);
        let full = diagonalize_hermitian(&h).unwrap();
        // Warm-start with the true eigenvectors (all 3 cols) to rotate 2 bands.
        let sub = diagonalize_subspace(&h, 2, Some(&full.eigenvectors)).unwrap();
        for k in 0..2 {
            assert!(
                relative_eq!(sub.eigenvalues[k], full.eigenvalues[k], epsilon = 1e-12),
                "band {k}: full={}, sub={}",
                full.eigenvalues[k],
                sub.eigenvalues[k]
            );
        }
    }

    #[test]
    fn test_subspace_matches_dense_under_small_perturbation() {
        // Simulate the SCF continuity case: H_new = H_old + ε · δ, where
        // δ is a small Hermitian perturbation. Warm-starting with the
        // eigenvectors of H_old should produce eigenvalues of H_new that
        // match a full solve to well under the SCF energy threshold
        // (1e-5 eV), exercising the main-path success.
        let n = 20;
        // H_old: diagonal kinetic + small off-diagonal coupling.
        let mut h_old: faer::Mat<Complex64> = faer::Mat::zeros(n, n);
        for j in 0..n {
            h_old[(j, j)] = Complex64::new(1.0 + j as f64 * 0.5, 0.0);
        }
        for j in 0..n - 1 {
            h_old[(j, j + 1)] = Complex64::new(0.1, 0.05);
            h_old[(j + 1, j)] = Complex64::new(0.1, -0.05);
        }
        let old = diagonalize_lowest(&h_old, 8).unwrap();

        // H_new: same structure, slightly perturbed diagonal (mimics
        // V_xc/V_H update between SCF iterations).
        let mut h_new = h_old.clone();
        for j in 0..n {
            h_new[(j, j)] += Complex64::new(1e-4 * (j as f64 + 1.0).sin(), 0.0);
        }

        let full_new = diagonalize_lowest(&h_new, 8).unwrap();
        let sub_new =
            diagonalize_subspace(&h_new, 8, Some(&old.eigenvectors)).unwrap();

        for k in 0..8 {
            let diff = (full_new.eigenvalues[k] - sub_new.eigenvalues[k]).abs();
            assert!(
                diff < 1e-8,
                "band {k}: full={:.10}, sub={:.10}, diff={:.2e}",
                full_new.eigenvalues[k],
                sub_new.eigenvalues[k],
                diff
            );
        }
    }

    #[test]
    fn test_subspace_falls_back_when_v_prev_is_random() {
        // Stress test: random v_prev that is NOT close to any eigenvector.
        // The residual gate must trip and produce the same result as
        // diagonalize_lowest. This is the "bad warm-start → correct fallback"
        // path that the SCF relies on for correctness in early iterations.
        let n = 10;
        let mut h: faer::Mat<Complex64> = faer::Mat::zeros(n, n);
        for j in 0..n {
            h[(j, j)] = Complex64::new(j as f64, 0.0);
        }
        for j in 0..n - 1 {
            h[(j, j + 1)] = Complex64::new(0.3, 0.0);
            h[(j + 1, j)] = Complex64::new(0.3, 0.0);
        }
        // Deterministic "random" v_prev: columns are simple sinusoidal
        // vectors, manifestly not the eigenvectors of the above tridiag.
        let v_prev: faer::Mat<Complex64> = faer::Mat::from_fn(n, 4, |r, c| {
            Complex64::new(
                ((r as f64 + 1.0) * (c as f64 + 1.0)).sin(),
                ((r as f64 + 1.0) * (c as f64 + 1.0)).cos(),
            )
        });
        let reference = diagonalize_lowest(&h, 4).unwrap();
        let sub = diagonalize_subspace(&h, 4, Some(&v_prev)).unwrap();
        for k in 0..4 {
            assert!(
                relative_eq!(sub.eigenvalues[k], reference.eigenvalues[k], epsilon = 1e-10),
                "band {k}: ref={:.10}, sub={:.10}",
                reference.eigenvalues[k],
                sub.eigenvalues[k]
            );
        }
    }

    #[test]
    fn test_eigenvalue_residuals_zero_for_true_eigenpairs() {
        // Sanity check: residuals of true eigenpairs are at floating-
        // point noise.
        let h = mat_from_rows(3, &[
            Complex64::new(2.0, 0.0), Complex64::new(0.3, 0.0), Complex64::new(0.0, 0.0),
            Complex64::new(0.3, 0.0), Complex64::new(3.0, 0.0), Complex64::new(0.2, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(0.2, 0.0), Complex64::new(4.0, 0.0),
        ]);
        let full = diagonalize_hermitian(&h).unwrap();
        let res = eigenvalue_residuals(&h, &full.eigenvalues, &full.eigenvectors);
        for (k, r) in res.iter().enumerate() {
            assert!(*r < 1e-12, "band {k}: residual = {r:.3e}");
        }
    }

    #[test]
    fn test_eigenvalue_residuals_flags_wrong_pair() {
        // A wrong eigenvalue produces a residual on the order of the shift.
        let h = mat_from_rows(2, &[
            Complex64::new(3.0, 0.0), Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0), Complex64::new(3.0, 0.0),
        ]);
        let full = diagonalize_hermitian(&h).unwrap();
        // Corrupt the first eigenvalue by 0.5.
        let mut evs = full.eigenvalues.clone();
        evs[0] += 0.5;
        let res = eigenvalue_residuals(&h, &evs, &full.eigenvectors);
        assert!(
            res[0] > 0.4,
            "corrupted eigenvalue should give residual ≈ 0.5, got {:.3e}",
            res[0]
        );
    }
}
