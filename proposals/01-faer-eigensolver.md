# Proposal: Replace nalgebra-lapack with faer

## Motivation

The eigensolver is the single largest bottleneck in the SCF loop, consuming 25-35% of wall time. The current implementation (`src/eigensolver/dense.rs`) calls LAPACK `zheev` through unsafe FFI with manual workspace queries, flat buffer copies, and raw pointers. It is single-threaded and always computes the full eigendecomposition even when only the lowest `n_bands` are needed (`diagonalize_lowest` truncates after the fact).

`faer` is a pure-Rust linear algebra library with performance matching or exceeding vendor LAPACK on modern hardware. It provides:

- Safe, ergonomic Hermitian eigendecomposition (no unsafe FFI)
- Multi-threaded eigensolve by default (using rayon internally)
- Partial eigendecomposition support via Krylov methods (future path)
- No platform-specific feature flags (`lapack-accelerate` goes away)

## Scope of Changes

### Dependencies

Remove:
```toml
nalgebra-lapack = { version = ">=0.27", ... }
lapack = { version = ">=0.20", ... }
```

Add:
```toml
faer = ">=0.21"
```

Keep `nalgebra` for `Vector3<f64>` geometry operations (crystal, k-points, G-vectors). `faer` handles all dense matrix algebra.

### Phase 1: Eigensolver replacement

**File: `src/eigensolver/dense.rs`**

Replace the entire `diagonalize_hermitian` function. Current code (lines 19-89) does:
1. Manual column-major copy into flat `Vec<Complex64>` (lines 32-37)
2. Workspace query via `lapack::zheev` with `lwork = -1` (lines 44-59)
3. Workspace allocation + actual diagonalization (lines 62-80)
4. Reconstruction of `DMatrix` from flat buffer (line 83)

New code with faer:
```rust
use faer::prelude::*;
use faer::complex_native::c64;

pub fn diagonalize_hermitian(h: &faer::Mat<c64>) -> EigenResult {
    let n = h.nrows();
    let decomp = h.selfadjoint_eigendecomposition(faer::Side::Lower);
    let eigenvalues: Vec<f64> = (0..n).map(|i| decomp.s().column_vector()[i]).collect();
    let eigenvectors = decomp.u().to_owned();
    EigenResult { eigenvalues, eigenvectors }
}
```

This eliminates all unsafe code, workspace management, and manual buffer copies.

**File: `src/eigensolver/dense.rs` — `EigenResult` type**

Change the eigenvector type:
```rust
pub struct EigenResult {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: faer::Mat<c64>,  // was: DMatrix<Complex64>
}
```

### Phase 2: Hamiltonian matrix type migration

**File: `src/hamiltonian.rs`**

Change `build_kinetic` and `build_hamiltonian` to return `faer::Mat<c64>` instead of `DMatrix<Complex64>`. The construction logic is identical — just different matrix API:

```rust
// Before (line 13):
let mut h = DMatrix::zeros(n, n);
h[(i, i)] = Complex64::new(ke, 0.0);

// After:
let mut h = faer::Mat::<c64>::zeros(n, n);
h[(i, i)] = c64::new(ke, 0.0);
```

**File: `src/scf/mod.rs` — `build_hamiltonian_with_v_eff` (lines 412-440)**

Same pattern: replace `nalgebra::DMatrix::zeros(n, n)` with `faer::Mat::<c64>::zeros(n, n)`. Indexing syntax is identical.

**File: `src/potential/nonlocal.rs` — `add_to_hamiltonian` (lines 112-202)**

Change the `h` parameter type from `&mut nalgebra::DMatrix<Complex64>` to `&mut faer::Mat<c64>`. The mutation pattern `h[(ig, jg)] += ...` works the same way.

### Phase 3: Density computation

**File: `src/scf/density.rs`**

`compute_density` accesses wavefunctions as `wfn[(ig, ib)]` (line 57). With faer, column access is `wfn.read(ig, ib)` or use column views. Minimal change.

### Phase 4 (future): Iterative eigensolver

With faer's `SparseColMat` and built-in Krylov infrastructure, implement a LOBPCG or Davidson solver that computes only the lowest `n_bands` eigenvalues. This would change the scaling from O(n^3) to O(n^2 * n_bands) for the eigensolve step, which is the single highest-impact optimization possible.

## Migration Strategy

1. Add `faer` alongside existing deps. Implement `diagonalize_hermitian_faer` as a parallel function.
2. Add a test comparing faer vs LAPACK results (eigenvalues within 1e-12, eigenvectors unitary).
3. Run SCF benchmarks comparing wall time.
4. Once validated, swap the Hamiltonian type from `DMatrix<Complex64>` to `faer::Mat<c64>` across the codebase.
5. Remove `nalgebra-lapack` and `lapack` from Cargo.toml.

## Risks

- `faer` uses its own complex type (`c64`) rather than `num_complex::Complex64`. Conversion is needed at the FFT boundary (faer matrices <-> Complex64 vectors for FFT). This is a `transmute`-safe zero-cost conversion since both are `[f64; 2]` layout.
- faer's eigendecomposition ordering may differ from LAPACK (ascending vs arbitrary). Verify in tests.
- nalgebra `Vector3<f64>` remains throughout geometry code — no conflict, just two linear algebra crates for different purposes.

## Expected Impact

- **Performance:** ~2x speedup on the eigensolve step from multi-threading alone. For a 200-basis system, this is ~15-20% total SCF speedup.
- **Code quality:** Eliminates all unsafe LAPACK FFI, workspace management, and buffer copies.
- **Portability:** No more `lapack-accelerate` feature flag. Builds identically on macOS and Linux.
- **Future path:** Enables iterative eigensolvers without external Fortran dependencies.
