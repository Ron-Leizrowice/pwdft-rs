---
id: DVSN
status: active
priority: low
complexity: large
risk: medium
depends_on: []
blocks: [SPRS]
---

# DVSN: Iterative Eigensolver (Davidson / LOBPCG)

> **Note:** Line numbers reference the pre-ScfContext codebase. Proposal 01 (faer) is now completed. Verify locations before implementing.

## Problem

The current eigensolver (`src/eigensolver/dense.rs`) performs full O(n^3) diagonalization via LAPACK `zheev` even though only the lowest `n_bands` eigenvalues are needed. For n_pw = 200 and n_bands = 8, 96% of the computed eigenvalues are discarded (line 95-106: `diagonalize_lowest` truncates the full result).

Production DFT codes use iterative methods that compute only the desired eigenvalues:

| Method | Cost per iteration | Iterations | Total | Memory |
|--------|-------------------|------------|-------|--------|
| Dense (zheev) | O(n^3) | 1 | O(n^3) | O(n^2) |
| Davidson | O(n^2 * m) | 5-15 | O(n^2 * n_bands * 10) | O(n * m_subspace) |
| LOBPCG | O(n^2 * n_bands) | 10-30 | O(n^2 * n_bands * 20) | O(n * 3*n_bands) |
| RMM-DIIS | O(n^2 * n_bands) | 5-10 | O(n^2 * n_bands * 7) | O(n * m*n_bands) |

For n = 200, n_bands = 8: dense = 8M ops, Davidson ~ 3M ops. For n = 1000, n_bands = 20: dense = 1G ops, Davidson ~ 200M ops. The advantage grows with system size.

## References

- Davidson, E.R., J. Comp. Phys. 17, 87 (1975) — original Davidson method
- Knyazev, A.V., SIAM J. Sci. Comput. 23, 517 (2001) — LOBPCG
- Wood, D.M. & Zunger, A., J. Phys. A 18, 1343 (1985) — RMM-DIIS for DFT
- Kresse, G. & Furthmuller, J., Phys. Rev. B 54, 11169 (1996) — VASP methodology
- QE source: `cegterg.f90` — Davidson implementation in QE

## Implementation

### Option A: Davidson (recommended first implementation)

The Davidson method is the most widely used (QE default) and easiest to implement correctly.

**New file: `src/eigensolver/davidson.rs`**

```rust
use nalgebra::DMatrix;
use num_complex::Complex64;

pub struct DavidsonConfig {
    pub max_iter: usize,        // default 20
    pub residual_tol: f64,      // default 1e-8 eV
    pub max_subspace: usize,    // default 4 * n_bands (restart when exceeded)
}

/// Davidson eigensolver for the lowest n_bands eigenvalues of H.
///
/// Algorithm:
/// 1. Start with n_bands trial vectors (random or from previous iteration)
/// 2. Build subspace matrix: H_sub = V^H H V
/// 3. Diagonalize H_sub (small dense solve)
/// 4. Compute Ritz values (approximate eigenvalues) and residuals
/// 5. If all residuals < tol, converged
/// 6. Apply preconditioner to residuals: t_i = P * r_i
/// 7. Orthogonalize t against V, expand subspace: V <- [V, t]
/// 8. If subspace too large, restart with current Ritz vectors
/// 9. Goto 2
pub fn davidson(
    h: &DMatrix<Complex64>,
    n_bands: usize,
    config: &DavidsonConfig,
    initial_guess: Option<&DMatrix<Complex64>>,
) -> super::dense::EigenResult {
    let n = h.nrows();
    let mut subspace_size = n_bands;

    // Initialize trial vectors
    let mut v = match initial_guess {
        Some(guess) => guess.columns(0, n_bands.min(guess.ncols())).into_owned(),
        None => random_orthonormal_vectors(n, n_bands),
    };

    for _iter in 0..config.max_iter {
        // Build and solve subspace problem
        let hv = h * &v;
        let h_sub = v.adjoint() * &hv;
        let sub_result = super::dense::diagonalize_hermitian(&h_sub);

        // Ritz vectors: y_i = V * u_i (where u_i are eigenvectors of H_sub)
        let ritz_vectors = &v * &sub_result.eigenvectors;
        let ritz_values = &sub_result.eigenvalues;

        // Compute residuals: r_i = H|y_i> - e_i|y_i>
        let h_ritz = h * &ritz_vectors;
        let mut all_converged = true;
        let mut new_vectors = Vec::new();

        for i in 0..n_bands {
            let y_i = ritz_vectors.column(i);
            let r_i = h_ritz.column(i) - ritz_values[i] * &y_i;
            let res_norm = r_i.norm();

            if res_norm > config.residual_tol {
                all_converged = false;
                // Preconditioner: t_i = (diag(H) - e_i)^{-1} r_i
                let t = precondition_residual(&r_i, h, ritz_values[i]);
                new_vectors.push(t);
            }
        }

        if all_converged {
            return super::dense::EigenResult {
                eigenvalues: ritz_values[..n_bands].to_vec(),
                eigenvectors: ritz_vectors.columns(0, n_bands).into_owned(),
            };
        }

        // Expand subspace
        for t in new_vectors {
            // Orthogonalize against existing subspace (modified Gram-Schmidt)
            let t_orth = orthogonalize(&t, &v);
            if t_orth.norm() > 1e-14 {
                let t_norm = &t_orth / t_orth.norm();
                v = v.insert_column(v.ncols(), 0.0);
                v.set_column(v.ncols() - 1, &t_norm);
            }
        }

        // Restart if subspace too large
        if v.ncols() > config.max_subspace {
            v = ritz_vectors.columns(0, n_bands).into_owned();
        }
    }

    // Didn't converge in max_iter; return best approximation
    log::warn!("Davidson: did not converge in {} iterations, using best approximation", config.max_iter);
    let hv = h * &v;
    let h_sub = v.adjoint() * &hv;
    let sub_result = super::dense::diagonalize_hermitian(&h_sub);
    super::dense::EigenResult {
        eigenvalues: sub_result.eigenvalues[..n_bands].to_vec(),
        eigenvectors: (&v * &sub_result.eigenvectors).columns(0, n_bands).into_owned(),
    }
}

/// Diagonal preconditioner: t_i = r_i / (H_diag - e_i + shift)
fn precondition_residual(
    residual: &nalgebra::DVectorSlice<Complex64>,
    h: &DMatrix<Complex64>,
    eigenvalue: f64,
) -> nalgebra::DVector<Complex64> {
    let n = residual.len();
    let mut t = nalgebra::DVector::zeros(n);
    for i in 0..n {
        let denom = h[(i, i)].re - eigenvalue;
        let denom = if denom.abs() < 0.1 { denom.signum() * 0.1 } else { denom };
        t[i] = residual[i] / Complex64::new(denom, 0.0);
    }
    t
}
```

### Option B: LOBPCG (better parallel scaling)

LOBPCG is block-based (all n_bands vectors updated simultaneously) and has better parallel scaling than Davidson. It's the default in ABINIT and GPAW. The core iteration is:

```
X_{k+1} = argmin(Rayleigh quotient) over span(X_k, W_k, P_k)
where W_k = preconditioned residuals, P_k = previous search directions
```

This is more complex to implement correctly (requires careful orthogonalization and locking of converged eigenvectors) but has the advantage of a fixed memory footprint (3 * n_bands vectors).

### Integration with SCF

In `src/scf/mod.rs`, replace line 251:

```rust
// Before:
dense::diagonalize_lowest(&h, params.n_bands)

// After:
let prev_v = prev_wavefunctions.as_ref().map(|wfns| &wfns[ik_index]);
davidson::davidson(&h, params.n_bands, &davidson_config, prev_v)
```

### Input configuration

```toml
[scf]
eigensolver = "davidson"        # "dense" | "davidson" | "lobpcg"
davidson_max_iter = 20
davidson_tol = 1e-8             # eV
```

### Relationship to other proposals

- **Proposal 01 (faer):** If faer is adopted, the small subspace diagonalization inside Davidson uses faer instead of LAPACK. The Davidson structure is independent of the dense solver backend.
- **Proposal 15 (wavefunction reuse):** Davidson with warm-start from previous SCF eigenvectors is the standard production approach. Proposal 15's subspace diag is a stepping stone; this proposal is the full solution.
- **Proposal 07 (sparse matrices):** Davidson only needs the matrix-vector product H|v>, not the full matrix. This enables the sparse/matrix-free Hamiltonian from Proposal 07.

## Acceptance Criteria

1. **Correct eigenvalues:** Davidson produces the same lowest n_bands eigenvalues as dense diag (within 1e-10 eV) for the Si test case.
2. **Correct eigenvectors:** Eigenvectors are orthonormal (V^H V = I within 1e-12) and satisfy H|v> = e|v> within `residual_tol`.
3. **Fewer operations than dense:** For n_pw=200, n_bands=8, Davidson uses fewer matrix-vector products than O(n^3) dense diag. Benchmark wall time.
4. **Warm-start acceleration:** With previous eigenvectors, Davidson converges in 2-5 iterations instead of 10-15.
5. **SCF convergence unchanged:** Final SCF energy matches dense-diag SCF to within 1e-8 eV.
6. **Graceful convergence failure:** If Davidson doesn't converge, it returns the best approximation with a warning (not a panic).
7. **Dense fallback:** `eigensolver = "dense"` reproduces current behavior exactly.
8. **Missing eigenvalue detection:** Davidson finds all n_bands lowest eigenvalues, including degenerate ones (test with a system known to have degeneracies at high-symmetry k-points).
