//! Iterative partial Hermitian eigensolver via faer's implicitly-restarted
//! Arnoldi / Krylov-Schur method ([`faer::matrix_free::eigen::
//! partial_self_adjoint_eigen`]).
//!
//! ## Status: experimental (opt-in)
//!
//! Faer 0.24's internal `iterate_lanczos` routine contains a
//! reorthogonalization loop that can spin indefinitely when a Krylov
//! vector becomes numerically null during Gram-Schmidt (see
//! `operator/self_adjoint_eigen/mod.rs` line 42-59 in the faer source).
//! Individual calls on well-conditioned Hamiltonians converge correctly;
//! multi-iteration SCF loops on ill-conditioned inputs can sporadically
//! hang. Until the upstream library is fixed, the `EigensolverKind::Dense`
//! backend remains the default; users who enable `Iterative` should be
//! prepared to fall back to `Dense` if they observe non-termination.
//!
//! The SCF loop only needs the `n_bands` lowest eigenpairs of the Kohn-Sham
//! Hamiltonian, but the dense solver (`src/eigensolver/dense.rs`) computes all
//! `n_pw` eigenpairs and discards >99% of them. For production `n_pw` (≥200)
//! the dense solve is the #1 SCF bottleneck (see proposal `ITEV`). This module
//! computes only the lowest `n_bands` eigenpairs using faer's upstream
//! iterative Krylov solver, with a shift-and-flip trick to map
//! "algebraically-lowest of H" to "largest-magnitude of A' = σI − H".
//!
//! ## Shift-and-flip strategy
//!
//! faer's `partial_self_adjoint_eigen` returns the eigenvalues of `A` with
//! the largest magnitude, sorted in descending order. For a Hermitian
//! Kohn-Sham `H`, the quantity of interest is the lowest `n_bands` algebraic
//! eigenvalues. We define the operator
//!
//! ```text
//!   A' = σ · I − H,        σ = max(diag(H)) + ||H − diag(H)||_∞
//! ```
//!
//! which is Hermitian (since `σ I` and `H` are), positive-semidefinite by
//! construction of σ, and whose spectrum is simply `{σ − λᵢ(H)}`. The
//! largest-magnitude eigenvalues of `A'` correspond to the algebraically
//! smallest eigenvalues of `H`. Eigenvectors are preserved (unchanged), and
//! the eigenvalue of `H` is recovered as `λ_H = σ − λ_A'`.

use dyn_stack::{MemBuffer, MemStack};
use faer::{
    ColRef, Mat, MatMut, MatRef, Par,
    matrix_free::{
        LinOp,
        eigen::{PartialEigenParams, partial_self_adjoint_eigen, partial_self_adjoint_eigen_scratch},
    },
};
use faer::linalg::matmul::matmul;
use faer::prelude::{Reborrow, ReborrowMut};
use num_complex::Complex64;

use crate::error::{PwdftError, Result};

use super::dense::EigenResult;

/// Default convergence tolerance for the iterative eigensolver.
///
/// Faer's internal default is `f64::EPSILON * 128 ≈ 2.8e-14`, which is tight
/// enough that SCF energies match the dense solver to machine precision.
pub const DEFAULT_TOL: f64 = f64::EPSILON * 128.0;

/// Default maximum number of Arnoldi restarts before giving up.
///
/// Well within what faer's upstream tests use (up to 1000). If the solver
/// fails to converge within this budget, callers should fall back to the
/// dense path and log a warning.
pub const DEFAULT_MAX_RESTARTS: usize = 500;

/// Shift-and-flip wrapper: applies `(σ · I − H) · v` for a Hermitian `H`.
///
/// The dense matrix `H` is borrowed; the shift `σ` is precomputed once at
/// construction time so that the largest-magnitude eigenvalues of this
/// operator correspond to the algebraically smallest eigenvalues of `H`.
#[derive(Debug)]
struct ShiftedHermitianOp<'a> {
    /// Dense Hermitian matrix. Only the full storage is used (via matmul);
    /// no special symmetric-storage optimization.
    h: MatRef<'a, Complex64>,
    /// Shift σ. Chosen so that `σ I − H` is positive-semidefinite.
    sigma: f64,
}

impl<'a> ShiftedHermitianOp<'a> {
    fn new(h: MatRef<'a, Complex64>, sigma: f64) -> Self {
        Self { h, sigma }
    }

    fn dim(&self) -> usize {
        self.h.nrows()
    }
}

impl LinOp<Complex64> for ShiftedHermitianOp<'_> {
    #[inline]
    fn nrows(&self) -> usize {
        self.dim()
    }

    #[inline]
    fn ncols(&self) -> usize {
        self.dim()
    }

    #[inline]
    fn apply_scratch(
        &self,
        _rhs_ncols: usize,
        _par: faer::Par,
    ) -> dyn_stack::StackReq {
        dyn_stack::StackReq::EMPTY
    }

    fn apply(
        &self,
        out: MatMut<'_, Complex64>,
        rhs: MatRef<'_, Complex64>,
        par: faer::Par,
        _stack: &mut MemStack,
    ) {
        let sigma = Complex64::new(self.sigma, 0.0);
        let mut out = out;
        // out := σ · rhs
        {
            let rhs_ref = rhs.rb();
            let mut oi = out.rb_mut();
            for j in 0..oi.ncols() {
                for i in 0..oi.nrows() {
                    oi[(i, j)] = sigma * rhs_ref[(i, j)];
                }
            }
        }
        // out := out − H · rhs  (accumulate -H·rhs into out)
        matmul(
            out.rb_mut(),
            faer::Accum::Add,
            self.h,
            rhs,
            -Complex64::new(1.0, 0.0),
            par,
        );
    }

    fn conj_apply(
        &self,
        out: MatMut<'_, Complex64>,
        rhs: MatRef<'_, Complex64>,
        par: faer::Par,
        stack: &mut MemStack,
    ) {
        // Self-adjoint: conj_apply = apply. (For a Hermitian operator,
        // (σI − H)^* = σ I − H^* = σ I − H.)
        self.apply(out, rhs, par, stack);
    }
}

/// Estimate a tight upper bound on `max(eigval(H))` via a per-row Gershgorin
/// disk bound: `max_i(diag_i + Σ_{j≠i} |H_{ij}|)`. For a Kohn-Sham
/// Hamiltonian the kinetic diagonal dominates, so this is typically within
/// a small factor of the true max eigenvalue — much tighter than the
/// (max diag) + (max off-diag row sum) approximation, which over-estimates
/// when the maximally-coupled row differs from the max-diagonal row.
///
/// A tight shift is essential: `(σI − H)` with σ far above λ_max has an
/// eigenvalue spectrum almost uniformly clustered near σ, which causes
/// Lanczos reorthogonalization to loop near-infinitely trying to resolve
/// nearly-degenerate Krylov vectors.
fn estimate_shift(h: MatRef<'_, Complex64>) -> f64 {
    let n = h.nrows();
    let mut max_row_bound = f64::NEG_INFINITY;
    for i in 0..n {
        let d = h[(i, i)].re;
        let mut row_sum = 0.0_f64;
        for j in 0..n {
            if i != j {
                row_sum += h[(i, j)].norm();
            }
        }
        let bound = d + row_sum;
        if bound > max_row_bound {
            max_row_bound = bound;
        }
    }
    // Small positive cushion — enough to guarantee A' = σI − H stays PSD
    // under floating-point round-off but not so large that conditioning
    // degrades.
    max_row_bound + max_row_bound.abs().mul_add(1e-6, 1e-6)
}

/// Normalize a complex column in place. Returns the original L2 norm.
fn normalize(v: &mut [Complex64]) -> f64 {
    let norm: f64 = v.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
    if norm > 0.0 {
        let inv = 1.0 / norm;
        for c in v.iter_mut() {
            *c *= inv;
        }
    }
    norm
}

/// Compute the lowest `n_bands` eigenpairs of `h` using faer's
/// implicitly-restarted Arnoldi iterative eigensolver.
///
/// ## Arguments
/// - `h`: Hermitian matrix (the Kohn-Sham Hamiltonian at one k-point). Only
///   read via its dense representation (via faer's generic matmul).
/// - `n_bands`: number of lowest eigenpairs to compute.
/// - `v0`: optional warm-start vector. When the caller has a good guess for
///   the subspace (e.g. prior SCF iteration's ground-state eigenvector),
///   Arnoldi can converge in dramatically fewer restarts. If `None` or
///   empty, a deterministic unit-norm seed is used.
/// - `tol`: convergence tolerance on the residual norm. `DEFAULT_TOL` is
///   safe for SCF.
///
/// ## Returns
/// - `Ok(EigenResult)`: eigenvalues sorted ascending, eigenvectors as
///   columns of a matrix. If fewer than `n_bands` eigenpairs converged,
///   returns an error so the caller can fall back to the dense path.
/// - `Err(PwdftError::Eigensolver)`: convergence failure. The caller (SCF)
///   should log a warning and fall back to dense.
///
/// ## Correctness notes
///
/// Eigenvalues returned by faer's partial solver are those of `σI − H` (the
/// shifted operator), sorted by descending magnitude. We invert the shift
/// via `λ_H = σ − λ_{A'}` and re-sort ascending. Eigenvectors are identical
/// under the shift and are returned as-is.
///
/// # Errors
/// Returns `PwdftError::Eigensolver` if the matrix is not square, empty, or
/// if faer fails to converge at least `n_bands` eigenpairs within the max
/// restart budget.
pub fn diagonalize_lowest_iterative(
    h: &Mat<Complex64>,
    n_bands: usize,
    v0: Option<&[Complex64]>,
    tol: f64,
) -> Result<EigenResult> {
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
            eigenvectors: Mat::zeros(0, 0),
        });
    }
    // Request extra eigenvalues beyond `n_bands` to protect against
    // Krylov-Schur degeneracy collapse. Without this, Lanczos can miss one
    // of a set of (near-)degenerate eigenvalues at high-symmetry
    // Brillouin-zone points (e.g. the 3-fold valence-band degeneracy at
    // Γ in Si diamond), producing a one-eigenvalue gap in the returned
    // spectrum. 50% padding is standard (ARPACK's "ncv" parameter).
    let n_request = (n_bands + n_bands / 2).max(n_bands + 4);

    // Faer's `partial_self_adjoint_eigen_imp` panics when `max_dim >= n`
    // (it has no "full-dense" fallback on the self-adjoint path — line 290
    // in operator/self_adjoint_eigen/mod.rs). The default `max_dim` is
    // `max(2 * MIN_DIM, 2 * n_eigval) = max(64, 2·n_request)`, so we must
    // guarantee `n > max_dim`. In practice we also want a margin before
    // Arnoldi-style iteration is actually cheaper than a direct dense
    // Hessenberg reduction, so we fall back to dense for small `n`.
    let max_dim = core::cmp::max(64, 2 * n_request);
    if n <= max_dim {
        return super::dense::diagonalize_lowest(h, n_bands);
    }

    let sigma = estimate_shift(h.as_ref());
    let op = ShiftedHermitianOp::new(h.as_ref(), sigma);
    let par = Par::Seq; // k-point loop parallelism lives outside this call.

    // Starting vector. faer normalizes internally, but we seed a non-zero
    // vector to avoid hitting faer's "zero v0 → default seed" fallback.
    let mut v0_owned: Vec<Complex64> = match v0 {
        Some(slice) if slice.len() == n => slice.to_vec(),
        _ => (0..n)
            // Deterministic but non-degenerate seed; entries chosen to have
            // support on every basis function so that all eigenvectors
            // enter the Krylov subspace.
            .map(|i| {
                let phase = (i as f64).mul_add(0.618_033_988_749_89, 1.0);
                Complex64::new(phase.sin(), phase.cos())
            })
            .collect(),
    };
    let _ = normalize(&mut v0_owned);
    // SAFETY: v0_owned has length n; form a column reference.
    let v0_mat = MatRef::from_column_major_slice(&v0_owned, n, 1);
    let v0_ref: ColRef<'_, Complex64> = v0_mat.col(0);

    let params = PartialEigenParams {
        max_restarts: DEFAULT_MAX_RESTARTS,
        ..Default::default()
    };

    // Allocate eigenvalue / eigenvector outputs for `n_request` pairs (we
    // over-request to defend against degeneracy collapse; see above).
    // We truncate to `n_bands` after sorting.
    let mut eigvals_shifted: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n_request];
    let mut eigvecs: Mat<Complex64> = Mat::zeros(n, n_request);

    let scratch_req = partial_self_adjoint_eigen_scratch::<Complex64>(
        &op as &dyn LinOp<Complex64>,
        n_request,
        par,
        params,
    );
    let mut mem = MemBuffer::new(scratch_req);
    let stack = MemStack::new(&mut mem);

    let info = partial_self_adjoint_eigen::<Complex64>(
        eigvecs.as_mut(),
        &mut eigvals_shifted,
        &op as &dyn LinOp<Complex64>,
        v0_ref,
        tol,
        par,
        stack,
        params,
    );
    let n_converged = info.n_converged_eigen;
    if n_converged < n_bands {
        return Err(PwdftError::Eigensolver {
            size: n,
            detail: format!(
                "partial_self_adjoint_eigen converged only {n_converged}/{n_bands} eigenpairs \
                 within {DEFAULT_MAX_RESTARTS} restarts"
            ),
        });
    }

    // Eigenvalues of A' = σ I − H are {σ − λᵢ(H)}; faer returns them sorted
    // descending by magnitude. Since A' is PSD-by-construction, magnitude =
    // value, so descending-by-magnitude = descending-by-value. Inverting the
    // shift flips the order to ascending in H, which is what we want.
    // We requested `n_request` pairs; sort the available converged ones
    // ascending in H, then take the lowest `n_bands`.
    let n_available = n_converged.min(n_request);
    let raw_eigenvalues: Vec<f64> = eigvals_shifted[..n_available]
        .iter()
        .map(|c| sigma - c.re)
        .collect();
    let mut order: Vec<usize> = (0..n_available).collect();
    order.sort_by(|&a, &b| {
        raw_eigenvalues[a]
            .partial_cmp(&raw_eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    // Truncate to n_bands after sorting.
    let keep = order.iter().take(n_bands).copied().collect::<Vec<_>>();
    let eigenvalues: Vec<f64> = keep.iter().map(|&i| raw_eigenvalues[i]).collect();
    let mut sorted_eigvecs: Mat<Complex64> = Mat::zeros(n, n_bands);
    for (new_col, &old_col) in keep.iter().enumerate() {
        for row in 0..n {
            sorted_eigvecs[(row, new_col)] = eigvecs[(row, old_col)];
        }
    }

    Ok(EigenResult {
        eigenvalues,
        eigenvectors: sorted_eigvecs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;
    use num_complex::Complex64;

    fn mat_from_rows(n: usize, data: &[Complex64]) -> Mat<Complex64> {
        Mat::from_fn(n, n, |r, c| data[r * n + c])
    }

    #[test]
    fn iterative_falls_back_to_dense_at_small_n() {
        // n=2 is below Arnoldi's minimum viable size (max_dim = 64 > n),
        // so `diagonalize_lowest_iterative` must internally fall back to
        // dense without panicking, producing the correct eigenvalues.
        let h = mat_from_rows(
            2,
            &[
                Complex64::new(3.0, 0.0), Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0), Complex64::new(3.0, 0.0),
            ],
        );
        let r = diagonalize_lowest_iterative(&h, 1, None, DEFAULT_TOL).unwrap();
        assert!(relative_eq!(r.eigenvalues[0], 2.0, epsilon = 1e-10));
    }

    #[test]
    fn iterative_matches_dense_diagonal_100() {
        // A 100×100 diagonal matrix with eigenvalues 1, 2, ..., 100. This
        // is above faer's Arnoldi floor (max_dim ≈ 64), so the Krylov
        // path is actually exercised.
        let n = 100;
        let mut h: Mat<Complex64> = Mat::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = Complex64::new((i + 1) as f64, 0.0);
        }
        let r = diagonalize_lowest_iterative(&h, 4, None, DEFAULT_TOL).unwrap();
        for (got, want) in r.eigenvalues.iter().zip([1.0, 2.0, 3.0, 4.0]) {
            assert!(
                relative_eq!(*got, want, epsilon = 1e-8),
                "expected {want}, got {got}"
            );
        }
    }

    #[test]
    fn iterative_matches_dense_complex_hermitian_100() {
        // Build a non-trivial Hermitian matrix: H = D + B + B^†.
        let n = 100;
        let mut h: Mat<Complex64> = Mat::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = Complex64::new((i as f64).mul_add(0.5, 1.0), 0.0);
        }
        for i in 0..n - 1 {
            let c = Complex64::new(0.1, 0.05);
            h[(i, i + 1)] = c;
            h[(i + 1, i)] = c.conj();
        }
        for i in 0..n - 4 {
            let c = Complex64::new(0.03, -0.02);
            h[(i, i + 4)] = c;
            h[(i + 4, i)] = c.conj();
        }

        let dense = super::super::dense::diagonalize_lowest(&h, 6).unwrap();
        let iter = diagonalize_lowest_iterative(&h, 6, None, DEFAULT_TOL).unwrap();

        for (d, i) in dense.eigenvalues.iter().zip(iter.eigenvalues.iter()) {
            assert!(
                relative_eq!(d, i, epsilon = 1e-8),
                "eigenvalue mismatch: dense={d:.10}, iterative={i:.10}"
            );
        }
    }

    #[test]
    fn iterative_residual_is_small() {
        let n = 100;
        let mut h: Mat<Complex64> = Mat::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = Complex64::new(i as f64, 0.0);
        }
        for i in 0..n - 1 {
            h[(i, i + 1)] = Complex64::new(0.2, 0.0);
            h[(i + 1, i)] = Complex64::new(0.2, 0.0);
        }
        let r = diagonalize_lowest_iterative(&h, 5, None, DEFAULT_TOL).unwrap();
        for k in 0..5 {
            let v = r.eigenvectors.col(k);
            let lambda = Complex64::new(r.eigenvalues[k], 0.0);
            let mut hv = vec![Complex64::new(0.0, 0.0); n];
            for row in 0..n {
                let mut acc = Complex64::new(0.0, 0.0);
                for col in 0..n {
                    acc += h[(row, col)] * v[col];
                }
                hv[row] = acc;
            }
            let residual: f64 = hv
                .iter()
                .enumerate()
                .map(|(i, &x)| (x - lambda * v[i]).norm_sqr())
                .sum::<f64>()
                .sqrt();
            assert!(
                residual < 1e-8,
                "eigenpair {k}: |Hv−λv| = {residual:.3e}"
            );
        }
    }

    #[test]
    fn iterative_handles_degenerate_eigenvalues() {
        // Diamond-structure Si has 4-fold degenerate valence-band triplet
        // at the Γ and L high-symmetry points. Simulate that by building a
        // matrix with explicit 3-fold degenerate eigenvalues and confirm
        // that the iterative solver returns three independent eigenvectors
        // (not copies of one). This is the regression test for the bug
        // where Lanczos collapsed the Γ-point degeneracy during SCF.
        let n = 100;
        let mut h: Mat<Complex64> = Mat::zeros(n, n);
        // First three eigenvalues: 1, 1, 1 (triple degenerate).
        // Next ones: 2, 3, ..., n-2.
        for i in 0..n {
            h[(i, i)] = Complex64::new(if i < 3 { 1.0 } else { (i - 1) as f64 }, 0.0);
        }
        // Light off-diagonal coupling to avoid purely-diagonal-matrix edge case.
        for i in 0..n - 1 {
            h[(i, i + 1)] = Complex64::new(0.01, 0.0);
            h[(i + 1, i)] = Complex64::new(0.01, 0.0);
        }
        let dense = super::super::dense::diagonalize_lowest(&h, 5).unwrap();
        let iter = diagonalize_lowest_iterative(&h, 5, None, DEFAULT_TOL).unwrap();
        for (k, (d, i)) in dense
            .eigenvalues
            .iter()
            .zip(iter.eigenvalues.iter())
            .enumerate()
        {
            assert!(
                (d - i).abs() < 1e-8,
                "band {k}: dense={d:.10}, iterative={i:.10}"
            );
        }
        // Eigenvectors for the degenerate triplet must span a rank-3
        // subspace. Check via SVD-equivalent: the 3 columns should be
        // pairwise orthogonal (within Lanczos tolerance).
        for i in 0..3 {
            for j in (i + 1)..3 {
                let mut ip = Complex64::new(0.0, 0.0);
                for k in 0..n {
                    ip += iter.eigenvectors[(k, i)].conj() * iter.eigenvectors[(k, j)];
                }
                assert!(
                    ip.norm() < 1e-6,
                    "degenerate eigenvectors {i} and {j} not orthogonal: |<v_i,v_j>|={}",
                    ip.norm()
                );
            }
        }
    }

    #[test]
    fn iterative_matches_dense_on_real_si_hamiltonian() {
        // Regression test for the SCF convergence bug: on a real Si
        // diamond Hamiltonian (kinetic + KB non-local, Γ-point, ecut
        // = 100 eV → n_pw = 89, which straddles the Arnoldi breakover),
        // Lanczos collapsed the Γ-point valence degeneracy and produced
        // nonsense eigenvalues. This test reconstructs that exact
        // scenario and requires a match to dense within 1e-6 Ha.
        use crate::basis::BasisSet;
        use crate::crystal::{Atom, Crystal, Lattice};
        use crate::hamiltonian;
        use crate::potential::nonlocal::NonlocalPotential;
        use nalgebra::Vector3;

        let a = 5.431;
        let crystal = Crystal {
            lattice: Lattice::new(
                a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
            ),
            atoms: vec![
                Atom::new(14, [0.0, 0.0, 0.0]),
                Atom::new(14, [0.25, 0.25, 0.25]),
            ],
        };
        let pp_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("pseudopotentials/nc/lda/Si.upf");
        let pp = crate::pseudopotential::load(&pp_path).unwrap();
        let basis = BasisSet::new(&crystal.lattice, 100.0);
        let k_gamma = Vector3::zeros();

        let mut h = hamiltonian::build_kinetic(&basis, &k_gamma);
        let vnl = NonlocalPotential::new(&crystal, &basis, &k_gamma, &[&pp]).unwrap();
        vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k_gamma);
        let n = h.nrows();
        assert!(n > 64, "need n > 64 for this test to exercise Arnoldi (got n={n})");

        let n_bands = 8;
        let dense = super::super::dense::diagonalize_lowest(&h, n_bands).unwrap();
        let iter = diagonalize_lowest_iterative(&h, n_bands, None, DEFAULT_TOL).unwrap();
        for (k, (d, i)) in dense
            .eigenvalues
            .iter()
            .zip(iter.eigenvalues.iter())
            .enumerate()
        {
            assert!(
                (d - i).abs() < 1e-6,
                "band {k}: dense={d:.10} eV, iterative={i:.10} eV, Δ={:.2e}",
                (d - i).abs()
            );
        }
    }

    #[test]
    fn iterative_matches_dense_scf_size_hamiltonian() {
        // Build a "realistic" SCF-size Hermitian via a lattice-of-decaying-
        // amplitude pattern. Purpose: ensure the iterative path works at
        // the n~90-300 regime where real SCF runs live, not just the
        // contrived small tests above. Takes ≈ 50 ms in release mode, well
        // within test-suite budget.
        let n = 120;
        let mut h: Mat<Complex64> = Mat::zeros(n, n);
        for i in 0..n {
            // Kinetic-like diagonal: E ~ |G|² increasing with index.
            h[(i, i)] = Complex64::new(1.0 + i as f64 * 0.2, 0.0);
        }
        // Off-diagonal couplings decaying with distance (mimics V_NL structure).
        for i in 0..n {
            for j in (i + 1)..n.min(i + 8) {
                let amp = 0.3 / (j - i) as f64;
                let phase = Complex64::new(0.0, (i as f64).mul_add(0.1, j as f64 * 0.07)).exp();
                h[(i, j)] = amp * phase;
                h[(j, i)] = h[(i, j)].conj();
            }
        }

        let dense = super::super::dense::diagonalize_lowest(&h, 8).unwrap();
        let iter = diagonalize_lowest_iterative(&h, 8, None, DEFAULT_TOL).unwrap();
        for (k, (d, i)) in dense
            .eigenvalues
            .iter()
            .zip(iter.eigenvalues.iter())
            .enumerate()
        {
            assert!(
                relative_eq!(d, i, epsilon = 1e-8),
                "band {k} mismatch: dense={d:.10}, iterative={i:.10}"
            );
        }
    }

    #[test]
    fn iterative_warm_start_matches_cold() {
        let n = 100;
        let mut h: Mat<Complex64> = Mat::zeros(n, n);
        for i in 0..n {
            h[(i, i)] = Complex64::new((i + 1) as f64, 0.0);
        }
        for i in 0..n - 1 {
            h[(i, i + 1)] = Complex64::new(0.1, 0.0);
            h[(i + 1, i)] = Complex64::new(0.1, 0.0);
        }
        let cold = diagonalize_lowest_iterative(&h, 4, None, DEFAULT_TOL).unwrap();
        let mut v0: Vec<Complex64> = (0..n).map(|i| cold.eigenvectors[(i, 0)]).collect();
        for c in v0.iter_mut().take(5) {
            *c += Complex64::new(0.001, 0.0);
        }
        let warm = diagonalize_lowest_iterative(&h, 4, Some(&v0), DEFAULT_TOL).unwrap();
        for (c, w) in cold.eigenvalues.iter().zip(warm.eigenvalues.iter()) {
            assert!(relative_eq!(c, w, epsilon = 1e-8));
        }
    }
}
