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
//! FIXME(faer-upstream): upstream issue not yet filed. When filed, paste
//! the tracking URL here and in the integration test's `#[ignore]`
//! reason string at `tests/itev_iterative_eigensolver.rs`. Grep for
//! `FIXME(faer-upstream)` to find all three sites.
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

/// Basis-size threshold above which `n_request` widens by a cluster-
/// margin buffer. Below this, the pre-ITEV2 `n_bands + n_bands/2`
/// padding is sufficient; above it, realistic Kohn-Sham Hamiltonians
/// routinely carry clusters whose resolution benefits from the wider
/// request.
const LARGE_BASIS_THRESHOLD: usize = 500;

/// Minimum absolute padding above `n_bands` in the "large-basis"
/// regime.
const MIN_PADDING: usize = 8;

/// Minimum Krylov subspace dimension. Faer's own default is `2 · MIN_DIM
/// = 64`, but empirical defect-1 sweeps show that the cluster-resolution
/// threshold sits at `max_dim = 128` for realistic Kohn-Sham
/// Hamiltonians in the `n_pw ≈ 180–260` regime (Cu FCC ecut=400,
/// Si ecut=200). Lifting the floor from 64 → 128 closes that regime
/// without affecting small-basis performance: if `n_pw < 128 +
/// margin`, the `n ≤ max_dim` check further down falls back to dense,
/// which is faster than Arnoldi at those sizes anyway.
const MIN_KRYLOV_MAX_DIM: usize = 128;

/// Fraction of `n_pw` to use as the baseline Krylov subspace
/// dimension.
///
/// Empirical diagnostic sweep on Si diamond Kohn-Sham Hamiltonians
/// (`itev2_tune_n_request_across_ecuts`, `#[ignore]`'d) at three
/// representative basis sizes:
///
/// | `n_pw` | pass threshold on `max_dim` | `max_dim / n_pw` |
/// |--------|-----------------------------|------------------|
/// | 89     | 64                          | 0.72             |
/// | 259    | 128                         | 0.49             |
/// | 725    | 176                         | 0.24             |
///
/// The ratio falls monotonically with basis size. Setting
/// `max_dim = max(MIN_KRYLOV_MAX_DIM, n_pw / 2)` clears the 0.24–0.49
/// boundary by a comfortable margin at every tabulated `n_pw` and
/// keeps the Arnoldi work bounded by half the dense-eigen budget.
/// The n_pw / 2 factor is encoded as the divisor
/// `KRYLOV_MAX_DIM_DIVISOR = 2`.
const KRYLOV_MAX_DIM_DIVISOR: usize = 2;

/// Choose the number of eigenpairs to request from faer's partial
/// solver (`n_request ≥ n_bands`). `n_request` sets how many Ritz
/// values are *returned*; the Krylov subspace *dimension* (faer's
/// `max_dim`, ARPACK's `NCV`) is chosen separately in
/// [`krylov_max_dim`]. Widening `n_request` without widening `max_dim`
/// can shrink the Arnoldi shift budget per restart and harm
/// convergence on clustered spectra (Lehoucq & Sorensen, *SIAM J.
/// Matrix Anal. Appl.* **17**, 789 (1996), §3.2 on the Implicit QR
/// Shift step).
///
/// Two regimes:
///
/// - **Small basis** (`n_pw < LARGE_BASIS_THRESHOLD`): the base
///   padding `n_bands + n_bands/2` (with a `+4` floor) resolves the
///   cluster. This matches the pre-ITEV2 heuristic exactly; the Si
///   ecut=100 regression test at `n_pw = 89` pins this constant.
/// - **Large basis** (`n_pw ≥ LARGE_BASIS_THRESHOLD`): the finer basis
///   packs more near-degenerate states into the same energy window,
///   so we add `max(n_bands, MIN_PADDING)` more slots to include the
///   full cluster.
///
/// The result is capped by `max_request` so the caller's `max_dim < n`
/// guarantee cannot be invalidated.
fn krylov_n_request(n_bands: usize, n_pw: usize, max_request: usize) -> usize {
    let base = (n_bands + n_bands / 2).max(n_bands + 4);
    let large_margin = if n_pw >= LARGE_BASIS_THRESHOLD {
        core::cmp::max(n_bands, MIN_PADDING)
    } else {
        0
    };
    let requested = base + large_margin;
    core::cmp::min(requested, max_request)
}

/// Choose the Krylov subspace dimension `max_dim` (ARPACK `NCV`,
/// faer's `PartialEigenParams::max_dim`).
///
/// Scales with the basis size: `max_dim = max(MIN_KRYLOV_MAX_DIM,
/// n_pw / KRYLOV_MAX_DIM_DIVISOR)`. The divisor is set from the
/// empirical pass/fail table (see [`KRYLOV_MAX_DIM_DIVISOR`]). The
/// `max_allowed` cap (caller's `n - 2`) guarantees faer's internal
/// assertion `max_dim < n` is never violated.
///
/// Also enforced: `max_dim ≥ 2 · n_request`, which is the ARPACK
/// lower bound and faer's own internal floor (see
/// `partial_self_adjoint_eigen` in `faer/faer/src/operator/eigen/mod.rs`
/// — `max_dim = max(params.max_dim, max(64, 2·n_eigval))`).
fn krylov_max_dim(n_request: usize, n_pw: usize, max_allowed: usize) -> usize {
    let by_basis = n_pw / KRYLOV_MAX_DIM_DIVISOR;
    let by_request = 2 * n_request;
    let wanted = core::cmp::max(by_basis, by_request);
    let floored = core::cmp::max(wanted, MIN_KRYLOV_MAX_DIM);
    core::cmp::min(floored, max_allowed)
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
    // Faer's `partial_self_adjoint_eigen_imp` panics when `max_dim >= n`
    // (no "full-dense" fallback on the self-adjoint path — see
    // `partial_self_adjoint_eigen_imp` in
    // `faer/faer/src/operator/self_adjoint_eigen/mod.rs`). We must size
    // both `n_request` (ARPACK `NEV`) and `max_dim` (ARPACK `NCV`)
    // below the matrix size.
    //
    // Ceiling chain:
    //   max_dim   ≤ n − 2     (faer's assertion + one slot of air)
    //   max_dim   ≥ 2 · n_request  (ARPACK lower bound)
    //   n_request ≥ n_bands
    //
    // When `n` is small enough that those constraints can't hold,
    // fall straight through to the dense solver. The Arnoldi crossover
    // is empirically at n ≈ 64–128 anyway.
    let max_dim_ceiling = n.saturating_sub(2);
    let max_request_ceiling = max_dim_ceiling / 2;
    if max_request_ceiling < n_bands {
        return super::dense::diagonalize_lowest(h, n_bands);
    }
    // Request extra eigenvalues beyond `n_bands` to protect against
    // Krylov-Schur degeneracy collapse.
    let n_request = krylov_n_request(n_bands, n, max_request_ceiling);
    let max_dim = krylov_max_dim(n_request, n, max_dim_ceiling);

    // Dense fallback when the matrix is too small for Arnoldi to pay
    // off: if `max_dim >= n`, faer panics on the self-adjoint path, and
    // even a hair smaller is in the regime where the O(n³) dense solver
    // wins outright.
    if n <= max_dim {
        return super::dense::diagonalize_lowest(h, n_bands);
    }

    run_partial_self_adjoint(h, n_bands, n_request, max_dim, v0, tol)
}

/// Core partial-eigensolver call with explicit Krylov subspace controls.
///
/// Callers of [`diagonalize_lowest_iterative`] get the heuristic-driven
/// (`n_request`, `max_dim`); this routine is exposed for diagnostics and
/// tuning (`#[cfg(test)]`-only callers in the module's unit tests).
#[allow(clippy::too_many_lines, reason = "single monolithic Arnoldi driver; the sort/retruncate tail is part of the same operation")]
fn run_partial_self_adjoint(
    h: &Mat<Complex64>,
    n_bands: usize,
    n_request: usize,
    max_dim: usize,
    v0: Option<&[Complex64]>,
    tol: f64,
) -> Result<EigenResult> {
    let n = h.nrows();
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
        // Pass our chosen Krylov subspace dimension explicitly. Faer
        // defaults `max_dim = max(64, 2 · n_request)`, which silently
        // under-resolves clustered spectra even when we've widened
        // `n_request` to include the whole cluster; setting `max_dim`
        // to `KRYLOV_NCV_RATIO · n_request` gives the IRAM restart
        // enough shifts to purge spurious Ritz pairs without disturbing
        // the cluster.
        max_dim,
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
    fn krylov_n_request_small_basis_regime() {
        // Below the large-basis threshold: base = max(n_bands + n_bands/2,
        // n_bands + 4) (pre-ITEV2 floor preserved exactly).
        // n_bands=8 → base = max(12, 12) = 12; no large margin.
        assert_eq!(krylov_n_request(8, 89, 10_000), 12);
        assert_eq!(krylov_n_request(8, 259, 10_000), 12);
        // n_bands=4 → base = max(6, 8) = 8.
        assert_eq!(krylov_n_request(4, 100, 10_000), 8);
        // n_bands=20 → base = max(30, 24) = 30.
        assert_eq!(krylov_n_request(20, 300, 10_000), 30);
    }

    #[test]
    fn krylov_n_request_large_basis_regime() {
        // At or above LARGE_BASIS_THRESHOLD (500): add max(n_bands, 8).
        // n_bands=8, n_pw=500 → base 12 + margin max(8, 8) = 8 → 20.
        assert_eq!(krylov_n_request(8, 500, 10_000), 20);
        // Canonical Si n_pw=725 case (defect 1 regression):
        // n_bands=8 → 12 + 8 = 20 (was 12 before the fix).
        assert_eq!(krylov_n_request(8, 725, 10_000), 20);
        // Larger n_bands scales the margin:
        // n_bands=12 → base = max(18, 16) = 18, margin = max(12, 8) = 12 → 30.
        assert_eq!(krylov_n_request(12, 1000, 10_000), 30);
    }

    #[test]
    fn krylov_n_request_respects_max_cap() {
        // The `max_request` cap must be honored so callers can guarantee
        // `KRYLOV_NCV_RATIO * n_request < n`. E.g. at a tiny cap of 10,
        // the requested value must not exceed 10.
        assert_eq!(krylov_n_request(8, 50, 10), 10);
        assert_eq!(krylov_n_request(8, 600, 20), 20);
    }

    #[test]
    fn krylov_max_dim_scales_with_basis() {
        // max_dim = max(128, max(2*n_request, n_pw/2)), capped.
        // Si ecut=100 n_pw=89: 89/2=44, floor 128 -> 128 (dense fallback downstream).
        assert_eq!(krylov_max_dim(12, 89, 10_000), 128);
        // Si ecut=200 n_pw=259: 259/2=129, floor 128 -> 129.
        assert_eq!(krylov_max_dim(12, 259, 10_000), 129);
        // Si ecut=400 n_pw=725: 725/2=362.
        assert_eq!(krylov_max_dim(12, 725, 10_000), 362);
    }

    #[test]
    fn krylov_max_dim_respects_arpack_floor() {
        // When 2*n_request > n_pw/2 and also > MIN_KRYLOV_MAX_DIM, the
        // ARPACK floor wins. n_request=80, n_pw=50: 50/2=25,
        // 2·80=160 > MIN_KRYLOV_MAX_DIM=128 -> 160.
        assert_eq!(krylov_max_dim(80, 50, 10_000), 160);
    }

    #[test]
    fn krylov_max_dim_caps_at_allowed() {
        // The `max_allowed` cap (caller's `n - 2`) must override the
        // scaling rule — faer panics otherwise.
        assert_eq!(krylov_max_dim(12, 259, 100), 100);
        assert_eq!(krylov_max_dim(8, 89, 30), 30);
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

    /// Diagnostic sweep: real Si Kohn-Sham Hamiltonian at several
    /// ecut values, finding `n_request` / `max_dim` combos that resolve
    /// the 3-fold-degenerate Γ valence cluster to `|Δ| ≤ 1e-10 eV` vs
    /// dense. Ignored by default — investigation tool, not a gate.
    #[test]
    #[ignore = "diagnostic sweep, not a gate"]
    fn itev2_tune_n_request_across_ecuts() {
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
        let pp_path = std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR"))
            .join("pseudopotentials/nc/lda/Si.upf");
        let pp = crate::pseudopotential::load(&pp_path).unwrap();
        let n_bands = 8;

        for ecut in [100.0_f64, 200.0, 400.0] {
            let basis = BasisSet::new(&crystal.lattice, ecut);
            let k_gamma = Vector3::zeros();
            let mut h = hamiltonian::build_kinetic(&basis, &k_gamma);
            let vnl = NonlocalPotential::new(&crystal, &basis, &k_gamma, &[&pp]).unwrap();
            vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k_gamma);
            let n = h.nrows();
            eprintln!("\n=== Si ecut={ecut} Γ: n_pw = {n} ===");

            let dense = super::super::dense::diagonalize_lowest(&h, n_bands).unwrap();
            for &n_req in &[12usize, 16, 20, 24] {
                for &md in &[64usize, 96, 128, 144, 160, 176, 192, 256] {
                    if md >= n || md < 2 * n_req {
                        continue;
                    }
                    let r = run_partial_self_adjoint(&h, n_bands, n_req, md, None, DEFAULT_TOL);
                    let max_err = match r {
                        Ok(res) => dense
                            .eigenvalues
                            .iter()
                            .zip(res.eigenvalues.iter())
                            .map(|(d, i)| (d - i).abs())
                            .fold(0.0_f64, f64::max),
                        Err(e) => {
                            eprintln!("  n_req={n_req:3} max_dim={md:4}: ERR {e}");
                            continue;
                        }
                    };
                    eprintln!(
                        "  n_req={n_req:3} max_dim={md:4}: max |Δ| = {max_err:.3e} eV{}",
                        if max_err <= 1e-10 { "  <-- PASS" } else { "" }
                    );
                }
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
        let pp_path = std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR"))
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
