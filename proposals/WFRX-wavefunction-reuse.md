---
id: WFRX
status: active
priority: low
complexity: medium
risk: low
depends_on: []
blocks: [DVSN]
---

# WFRX: Wavefunction Reuse Between SCF Iterations

> **Note:** Line numbers reference the pre-ScfContext codebase (src/scf/mod.rs was ~1127 lines, now ~709). Verify locations before implementing.

## Problem

At every SCF iteration, the Hamiltonian is diagonalized from scratch via `dense::diagonalize_lowest` (`src/scf/mod.rs`, line 251). The eigensolver starts with no knowledge of the previous iteration's eigenvectors, even though the Hamiltonian changes incrementally between iterations (only V_H and V_xc are updated; the kinetic and local terms are fixed).

For a converging SCF, the Hamiltonian at iteration n+1 differs from iteration n by a small perturbation. The eigenvectors change correspondingly little. By using the previous eigenvectors as a starting subspace, the eigensolver can converge in far fewer internal iterations.

Even with dense diagonalization (LAPACK zheev), there are two techniques that exploit continuity:

1. **Subspace diagonalization:** Project the new Hamiltonian into the old eigenvector subspace, diagonalize the small n_bands x n_bands matrix, then optionally refine.
2. **Warm-started iterative solver:** When an iterative eigensolver is available (Proposal 01, faer Phase 4), the previous eigenvectors provide an optimal starting guess.

This proposal covers technique 1 (compatible with the current dense eigensolver) and lays groundwork for technique 2.

## References

- Payne, M.C. et al., Rev. Mod. Phys. 64, 1045 (1992) — subspace rotation in DFT
- Kresse, G. & Furthmuller, J., Comp. Mater. Sci. 6, 15 (1996) — VASP methodology
- Davidson, E.R., J. Comp. Phys. 17, 87 (1975) — Davidson eigensolver with restart

## Implementation

### Step 1: Store previous eigenvectors

In `src/scf/mod.rs`, before the SCF loop:

```rust
// Previous eigenvectors per k-point, for subspace warm-start
let mut prev_wavefunctions: Option<Vec<DMatrix<Complex64>>> = None;
```

After the eigensolve (line 256):

```rust
prev_wavefunctions = Some(all_kpoint_wavefns.clone());
```

### Step 2: Subspace diagonalization

Add to `src/eigensolver/dense.rs`:

```rust
/// Diagonalize H using the previous eigenvectors as an approximate subspace.
///
/// 1. Project: H_sub = V_prev^H * H * V_prev  (n_bands x n_bands)
/// 2. Diagonalize H_sub (small dense solve)
/// 3. Rotate: V_new = V_prev * U  where U are eigenvectors of H_sub
///
/// This gives exact eigenvalues/vectors if the subspace is invariant (at
/// convergence), and excellent approximations when H has changed little.
///
/// For the first iteration (no previous vectors), falls back to full diag.
pub fn diagonalize_subspace(
    h: &DMatrix<Complex64>,
    n_bands: usize,
    v_prev: Option<&DMatrix<Complex64>>,
) -> EigenResult {
    let n = h.nrows();

    let v_prev = match v_prev {
        Some(v) if v.ncols() >= n_bands && v.nrows() == n => v,
        _ => return diagonalize_lowest(h, n_bands),
    };

    let v_sub = v_prev.columns(0, n_bands);

    // Project: H_sub = V^H * H * V
    let hv = h * &v_sub;
    let h_sub = v_sub.adjoint() * &hv;

    // Diagonalize the small matrix
    let sub_result = diagonalize_hermitian(&h_sub);

    // Rotate eigenvectors back to full basis
    let eigenvectors = &v_sub * &sub_result.eigenvectors;

    EigenResult {
        eigenvalues: sub_result.eigenvalues,
        eigenvectors: eigenvectors.into_owned(),
    }
}
```

Cost: O(n^2 * n_bands) for the projection (vs O(n^3) for full diag). For n=200, n_bands=8, this is ~100x faster.

### Step 3: Subspace refinement (optional)

The subspace diag gives approximate eigenvalues. For tighter accuracy, follow with one or two residual-based refinement steps:

```rust
/// Compute eigenvalue residual: r_i = H|v_i> - e_i|v_i>
/// If max(||r_i||) < tol, the subspace result is accurate enough.
pub fn eigenvalue_residuals(
    h: &DMatrix<Complex64>,
    eigenvalues: &[f64],
    eigenvectors: &DMatrix<Complex64>,
) -> Vec<f64> {
    (0..eigenvalues.len())
        .map(|i| {
            let v = eigenvectors.column(i);
            let hv = h * &v;
            let res = &hv - eigenvalues[i] * &v;
            res.norm()
        })
        .collect()
}
```

If residuals are above tolerance, fall back to full diagonalization for that k-point:

```rust
let result = diagonalize_subspace(h, n_bands, prev_v);
let residuals = eigenvalue_residuals(h, &result.eigenvalues, &result.eigenvectors);
let max_res = residuals.iter().cloned().fold(0.0_f64, f64::max);

if max_res > 1e-6 {
    // Subspace approximation insufficient, do full diag
    diagonalize_lowest(h, n_bands)
} else {
    result
}
```

### Step 4: Integration with SCF loop

Replace line 251 in `src/scf/mod.rs`:

```rust
// Before:
dense::diagonalize_lowest(&h, params.n_bands)

// After:
let prev_v = prev_wavefunctions.as_ref().map(|wfns| &wfns[ik_index]);
dense::diagonalize_subspace(&h, params.n_bands, prev_v)
```

Since the eigensolve is parallelized over k-points (line 242), `prev_wavefunctions` needs to be indexed per k-point. The `par_iter` closure captures the previous vectors by reference.

### Step 5: Subspace alignment (phase consistency)

Between iterations, eigenvectors may acquire arbitrary phase factors or rotate within degenerate subspaces. For subspace diag this is handled automatically (the projection absorbs any rotation). For warm-starting an iterative solver, explicit alignment may be needed:

```rust
/// Align new eigenvectors with previous ones to minimize phase jumps.
/// Computes overlap S = V_old^H * V_new and applies a unitary transform.
fn align_subspace(v_old: &DMatrix<Complex64>, v_new: &mut DMatrix<Complex64>) {
    let overlap = v_old.adjoint() * &*v_new;
    // SVD of overlap: U * Sigma * V^H
    // Optimal alignment: V_new <- V_new * V * U^H
    // For non-degenerate case, just fix phases: multiply each column by conj(sign(S_ii))
    for i in 0..overlap.ncols().min(overlap.nrows()) {
        let s = overlap[(i, i)];
        if s.norm() > 1e-10 {
            let phase = s / s.norm();
            v_new.column_mut(i).iter_mut().for_each(|x| *x *= phase.conj());
        }
    }
}
```

## Performance Impact

For a system with n_pw = 200, n_bands = 8:
- Full diag: O(200^3) = 8M operations per k-point
- Subspace diag: O(200^2 * 8) = 320K operations per k-point
- Speedup: ~25x on the eigensolve step (25-35% of SCF time)
- Net SCF speedup: ~20-30% for small systems

The speedup is even larger for bigger systems where n_pw >> n_bands.

The key caveat is accuracy: subspace diag is only valid when the Hamiltonian changes slowly. During early SCF iterations (large density changes), the subspace may be inaccurate and full diag should be used. The residual check (Step 3) handles this automatically.

## Acceptance Criteria

1. **Subspace diag produces correct eigenvalues:** For a converged SCF (H unchanged), subspace diag with the correct eigenvectors returns the same eigenvalues as full diag (within 1e-12 eV).
2. **SCF convergence is unchanged:** The final converged energy matches full-diag SCF to within 1e-8 eV.
3. **Performance improvement:** Subspace diag is measurably faster than full diag when n_bands << n_pw. Benchmark on Si with n_pw=200, n_bands=8.
4. **Graceful fallback:** When residuals exceed tolerance (early iterations, large H changes), the code falls back to full diag without affecting correctness.
5. **First iteration:** Without previous eigenvectors, full diag is used (no crash, no degradation).
6. **Parallel correctness:** Subspace diag works correctly within rayon parallel k-point loop.
