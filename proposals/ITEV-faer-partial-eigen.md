---
id: ITEV
status: active
priority: high
complexity: medium
risk: medium
depends_on: []
blocks: []
supersedes: [DVSN]
---

# ITEV: Iterative Eigensolver via `faer::matrix_free::eigen::partial_self_adjoint_eigen`

> **Scope.** Replace the full dense Hermitian eigensolve
> (`faer::self_adjoint_eigen`, called per k-point per SCF iteration) with
> faer's built-in Krylov-subspace iterative eigensolver, which computes only
> the lowest `n_bands` eigenpairs. Composable with WFRX warm-start.

## Problem

Profiling (Apple M3 Max, Si SCF, `si_scf_converged.yaml`, ecut=200 Ry, 4×4×4 MP
k-grid → 10 irr. k-points, n_pw=259) pins the eigensolver at roughly
**85-90% of user-code CPU time** in the SCF hot path:

```
=== sample-based profile, Si ecut=200, 5×SCF runs, 4 s sampling ===
Total samples: 31 197
  Idle (cond_wait / yield / spin):  39.7%
  Rayon pool overhead (truncated):  27.6%
  User code (classified):           32.7%
    eigensolver (faer+gemm):        60.7%  of classified user CPU
    V_NL apply (add_to_hamiltonian): 7.1%
    FFT:                             2.1%
    build_H_with_v_eff:              1.2%
    symmetrize_density:              0.1%
    XC grid:                         0.0%
```

The rayon "pool overhead" bucket is almost entirely truncated stacks **inside
faer's internal parallelism** (spindle for apply_block_householder / matmul)
and inside the per-k-point `par_iter` eigensolve block in
`src/scf/mod.rs:299-307`. Re-attributing it by call-site proportion places
eigensolver at roughly 85-90% of SCF user-mode compute.

Microbenchmark numbers (criterion, post-FFTB / post-FMAD, machine locked,
Apple M3 Max):

| n_pw | `faer_eigen` (full) | `vnl_apply` | `vnl_new` (init-only) | `kinetic` |
|------|---------------------|-------------|------------------------|-----------|
| 89   |   1.07 ms           |  0.53 ms    |   7.17 ms              |  2.77 µs  |
| 259  |  56.3  ms           |  4.70 ms    |  22.1  ms              | 22.6  µs  |
| 725  | 835.8  ms           | 28.6  ms    |  47.2  ms              | 166   µs  |

`faer_eigen` scales cubically and dominates absolutely at production `n_pw`.
At n_pw=259, **1 eigensolve ≈ 12× V_NL apply**. At n_pw=725 the ratio is
**29×**. Meanwhile `diagonalize_lowest` discards all but the lowest
`n_bands=8` of those ~725 eigenpairs (wastes ~99%).

(Note: the prior logbook's `n725 → 73.5 ms` entry was noise from a
concurrent process. Re-measured under a held machine lock, the correct
figure is 726-942 ms.)

## Why this is the #1 remaining perf target

FFTB (FFT buffer reuse, −23 to −33 % on `scf_iter_20x`) and FMAD (−3 to
−4 % on `lda_xc_grid`) have both landed. The next bottleneck ranking:

1. **Eigensolve** — 60-90 % of SCF CPU, O(n³), discards 99 % of its output
2. V_NL `add_to_hamiltonian` — 7 % of SCF CPU, O(n² · n_proj²)
3. FFT — 2 % of SCF CPU (FFTB already addressed)
4. Everything else — < 2 % combined

No other single bottleneck even comes close to the eigensolver. Any
serious further SCF speedup has to address this one first.

## References

- Davidson, E.R., J. Comp. Phys. 17, 87 (1975) — canonical Davidson method
- Knyazev, A.V., SIAM J. Sci. Comput. 23, 517 (2001) — LOBPCG
- Lehoucq, R.B. & Sorensen, D.C., SIAM J. Mat. Anal. Appl. 17, 789 (1996) —
  implicitly-restarted Arnoldi (the algorithm faer uses)
- Kresse, G. & Furthmuller, J., Phys. Rev. B 54, 11169 (1996) — VASP / RMM-DIIS
- faer 0.24 source:
  `src/operator/self_adjoint_eigen/mod.rs::partial_self_adjoint_eigen_imp`,
  `src/operator/eigen/mod.rs::partial_self_adjoint_eigen` (public wrapper)

## Proposed approach

**Use the upstream solver. Do NOT hand-roll.**

faer 0.24 (already in our dependency tree at `Cargo.toml`) ships a
production-ready implicitly-restarted Arnoldi partial Hermitian eigensolver:

```rust
faer::matrix_free::eigen::partial_self_adjoint_eigen(
    eigvecs: MatMut<'_, T>,
    eigvals: &mut [T],
    A: &dyn LinOp<T>,
    v0: ColRef<'_, T>,         // starting vector (WFRX warm-start fits here)
    tolerance: T::Real,
    par: Par,
    stack: &mut MemStack,
    params: PartialEigenParams, // { min_dim, max_dim, max_restarts }
) -> PartialEigenInfo;
```

Supporting evidence this is production-quality:
- Matrix-free via `faer::operator::LinOp`: only `H · v` is required, not a
  dense `H`. Aligns with the long-term sparse/matrix-free story (SPRS).
- Used inside faer for its own partial SVD and sparse eigen paths.
- Unit-tested in faer against reference eigendecompositions (see
  `test_arnoldi_real`, `test_arnoldi_cplx`, `test_toeplitz` in that file).
- Interface accepts `v0` directly → native WFRX warm-start (no wrapping
  required).

This eliminates the bulk of the DVSN implementation burden:
- no hand-rolled Davidson / LOBPCG kernel,
- no hand-rolled Gram-Schmidt / restart logic,
- no hand-rolled preconditioner-safety plumbing,
- we only need: (a) a `LinOp` impl for our Hamiltonian, (b) shift-invert
  (or negation) trick to get algebraically-smallest instead of
  largest-magnitude eigenvalues, (c) a tolerance and restart budget picker,
  (d) correctness tests.

### Shift strategy (largest-magnitude → lowest algebraic)

`partial_self_adjoint_eigen` returns eigenpairs sorted by descending
magnitude. For the Kohn-Sham Hamiltonian we want the algebraically-lowest
`n_bands`. Two clean options:

- **Spectrum flip:** define a `LinOp` that returns `−H v`. Run the solver;
  the largest-magnitude eigenvalues of `−H` correspond to the most negative
  eigenvalues of `H` (fine for bound states, which are the ones we want).
  Flip signs on output.
- **Shift-and-flip:** `A' v = (σ I − H) v` for a shift `σ` chosen above the
  band cutoff (say `σ = max(diag(H)) + ε`). Then largest-magnitude
  eigenvalues of `A'` correspond to lowest of `H`, all positive, which
  avoids any symmetry concerns with the sign flip. Eigenvalue of `H` is
  `σ − λ_A'`.

Recommended: shift-and-flip — numerically cleaner and matches standard
practice in ARPACK-style shift-and-invert when inversion is skipped. For
our matrices (`kinetic + V_eff + V_NL`), `max(diag(H)) ≈ k_max² + |V_eff|∞`
is cheap to estimate.

### Integration point

Replace the sole call site in `src/scf/mod.rs`:

```rust
// current:
dense::diagonalize_lowest(&h, ctx.params.n_bands)
```

with a thin wrapper in `src/eigensolver/iterative.rs`:

```rust
pub fn diagonalize_lowest_iterative(
    h: &faer::Mat<Complex64>,
    n_bands: usize,
    v0: Option<&faer::Mat<Complex64>>,
    tol: f64,
) -> Result<EigenResult>;
```

that (a) builds a `LinOp` wrapper implementing `H · v` via faer's mat-vec,
(b) picks `σ`, (c) runs `partial_self_adjoint_eigen`, (d) converts back to
our `EigenResult` type.

Note on parallelism: `Par::Seq` internally; the enclosing k-point
`par_iter` already saturates cores, so we do NOT want `Par::Rayon` inside
the k-point loop (it would oversubscribe). Worth benchmarking both; but
prior: seq is usually right when there's coarser parallelism outside.

### WFRX integration

WFRX (subspace diag with previous eigenvectors, currently status=active,
priority=low) fits naturally: pass the first column of the previous
iteration's eigvec matrix as `v0`. Warm-started Arnoldi typically cuts
iteration count by 2-5× when the Hamiltonian perturbation is small (which
is exactly what happens late in SCF). WFRX then becomes: store
`Vec<faer::Mat<Complex64>>` per k-point between SCF iterations and feed
the appropriate column into `v0`. Combined ITEV+WFRX is the realistic
production configuration.

### Fallback

Keep `diagonalize_lowest` as a selectable backend. YAML setting
`scf.eigensolver: "dense" | "iterative"` (default `dense` for one release,
flip to `iterative` after validation). CFGN already proposes exposing
numerics as settings — this fits its scheme.

### Expected speedup

Three regimes:

- **Small systems (n_pw ≤ 100):** dense may win. ~1 ms full diag is already
  near the floor of Arnoldi's fixed overhead. Expect parity or slight
  regression; this is why we keep the dense backend.
- **Medium (n_pw = 200-500):** expect 3-10× per-k-point. With warm-start
  (WFRX), expect 5-20× late in SCF. Concrete projection for n_pw = 259:
  dense = 56 ms; iterative (cold) ≈ 6-15 ms; iterative (warm) ≈ 3-8 ms.
- **Large (n_pw > 500):** expect 10-50× because O(n² · n_bands · k_iter) vs
  O(n³). At n_pw = 725: dense = 836 ms; iterative (cold) ≈ 30-80 ms;
  iterative (warm) ≈ 10-30 ms.

**Net SCF wall-time impact:** at n_pw = 725 the eigensolver is ~85 % of CPU;
a 10× eigensolver speedup converts a 100-unit SCF into 15 + 10 = 25
units ≈ **4× overall SCF speedup**. At n_pw = 259 the proportion is similar
but the absolute gain is smaller (SCF is already subsecond).

Conservative projected speedups at the n_pw = 259 baseline (ecut = 200,
converged yaml, 9 iter, 10 k-points):
- Current wall: 1.13 s (2.96 s CPU, 2.6× parallel)
- Projected with ITEV cold: 0.45 s (0.40× wall, 2.5× speedup)
- Projected with ITEV+WFRX warm: 0.30 s (0.27× wall, 3.8× speedup)

Order-of-magnitude only — actual numbers come from benchmarking.

## Relationship to existing proposals

- **DVSN** (Davidson/LOBPCG, hand-rolled, large + medium): obsoleted by this
  proposal. Mark as superseded on merge.
- **WFRX** (subspace warm-start, medium + low): still valid; plugs into ITEV
  via the `v0` parameter. Should be re-scoped as "ITEV warm-start" once
  ITEV lands. WFRX's subspace-rotation variant (project H into old
  subspace, dense-solve small matrix) is no longer the recommended path —
  ITEV+WFRX.v0 dominates it.
- **SPRS** (sparse matrices, large + medium): ITEV's `LinOp` wrapper gives us
  the matrix-free API SPRS needs. So ITEV doesn't block SPRS — rather,
  once ITEV lands, SPRS becomes "implement a cheaper `H · v` without
  building the dense `H`," which is the high-value half of SPRS.

So ITEV subsumes part of DVSN and unblocks the valuable half of SPRS while
keeping WFRX relevant.

## Implementation plan

### Phase 1: wrapper + correctness

1. Add `src/eigensolver/iterative.rs` with:
   - `struct HLinOp<'a>(&'a faer::Mat<Complex64>, f64)` holding `H` and shift `σ`.
   - `impl LinOp<c64> for HLinOp` computing `(σ · v − H · v)`.
   - `diagonalize_lowest_iterative(h, n_bands, v0, tol)` calling
     `partial_self_adjoint_eigen`, converting back eigenvalues via
     `λ_H = σ − λ_A'` and re-sorting ascending.
2. Port all `diagonalize_lowest` correctness tests
   (`src/eigensolver/dense.rs::tests`) to run against iterative too.
3. Add an integration test: iterative SCF energy matches dense SCF energy to
   within `1e-8 eV` on Si (ecut = 100, gamma) and Si (ecut = 200, 4×4×4).

### Phase 2: SCF integration with switch

4. Add `ScfParams::eigensolver: EigensolverKind { Dense, Iterative }`
   (default `Dense`) with YAML wiring.
5. In `src/scf/mod.rs`, branch on the kind at the two call sites (non-spin
   + spin).
6. Validate QE-match regression tests pass with `Iterative`.

### Phase 3: bench + tune

7. Extend `benches/scf_benchmarks.rs`:
   - `iterative_cold_n{89,259,725}` — first-call (no warm-start)
   - `iterative_warm_n{89,259,725}` — v0 = first column of prior result
8. Pick defaults for `min_dim`, `max_dim`, `max_restarts`, `tolerance`.
   Starting guesses: `min_dim = max(2 · n_bands, 16)`,
   `max_dim = 4 · n_bands`, `max_restarts = 200`, `tolerance = 1e-8 Ry`.

### Phase 4: WFRX warm-start

9. In `ScfContext`, store `prev_eigvecs: Option<Vec<faer::Mat<Complex64>>>`.
10. On iteration ≥ 2, pass `prev_eigvecs[ik].col(0)` as `v0`.
11. Bench and confirm the expected 2-5× iteration reduction.

### Phase 5: promote + deprecate

12. Flip default to `Iterative`.
13. Archive DVSN; re-scope WFRX to "ITEV warm-start (already done)" and
    archive.

## Acceptance criteria

1. **Correctness:** `diagonalize_lowest_iterative` produces the same lowest
   `n_bands` eigenvalues as `diagonalize_lowest` within `1e-10 Ry` on Si
   test matrices. Eigenvectors are orthonormal to `1e-12` and satisfy
   `‖H v − λ v‖ / ‖v‖ < tolerance`.
2. **SCF convergence unchanged:** Si (ecut = 100 and ecut = 200) converges
   to the same total energy (within `1e-8 eV`) and same iteration count ±1
   as dense.
3. **QE validation:** Tier-1+Tier-2 QE-match tests pass with iterative
   backend.
4. **Benchmarks:** microbenchmark `iterative_cold_n725` ≥ 5× faster than
   `faer_eigen_n725` dense.
5. **Graceful fallback:** on convergence failure within `max_restarts`, fall
   back to dense with a `log::warn!` — do NOT panic.
6. **Coexistence:** `ScfParams::eigensolver = Dense` bit-reproduces current
   behavior. Both backends covered by tests.

## Risks

- **Tolerance selection:** too loose → SCF doesn't converge; too tight →
  Arnoldi never stops. Plan: default `tol = 1e-8 Ry`, benchmark convergence
  iterations; tighten if needed.
- **Degenerate eigenvalues at high-symmetry k-points:** Krylov methods can
  have slow convergence at degeneracies. Test with Si X-point (known
  3-fold degeneracy) explicitly.
- **First-iteration cold-start cost:** without a good `v0`, cold Arnoldi
  may be slower than dense at small n. Benchmark; prefer dense when
  `n_pw ≤ 100` (hybrid strategy).
- **Shift estimate robustness:** if `σ` underestimates `max(eigenvalue)`,
  the spectral flip produces wrong eigenpair ordering. Use
  `σ = max(diag(H)) + ‖H − diag(H)‖_∞` or a safe over-estimate.
