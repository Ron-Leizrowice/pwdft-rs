# Proposal 20: Comprehensive Error Handling Cleanup

## Problem

The codebase already uses `thiserror` (v>=2) and has a `PwdftError` enum with 5 variants (`src/error.rs`). However, 25 production-code `unwrap()`/`expect()`/`panic!()` calls remain, concentrated in:

| Location | Count | Risk |
|----------|-------|------|
| GPU module (`gpu/mod.rs`) | 3 | High — in SCF hot path, channel failure = crash |
| Eigensolver (`eigensolver/dense.rs`) | 1 | High — eigendecomposition failure = crash |
| FFT (`fft.rs`) | 4 | Medium — shape mismatch = crash |
| Non-local potential (`potential/local.rs`) | 3 | Medium — PP lookup failure = crash |
| Main entry point (`main.rs`) | 3 | Medium — user-facing panics |
| Symmetry (`symmetry/density.rs`) | 2 | Low — setup code |
| Mixing (`scf/mixing.rs`) | 1 | Low — convergence path |
| Other | 8 | Low — mostly initialization |

The existing `PwdftError` enum lacks variants for GPU, eigensolver, and FFT failures, so these modules resort to panics.

## Implementation

### Step 1: Extend PwdftError with new variants

In `src/error.rs`:

```rust
#[derive(Debug, Error)]
pub enum PwdftError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("SCF did not converge after {iterations} iterations (delta = {delta:.2e})")]
    ConvergenceFailure { iterations: usize, delta: f64 },

    #[error("invalid input: {0}")]
    InvalidInput(String),

    #[error("missing pseudopotential for element {0}")]
    MissingPseudopotential(String),

    #[error("parse error: {0}")]
    Parse(String),

    // --- New variants ---

    #[error("eigendecomposition failed for {size}x{size} matrix: {detail}")]
    Eigensolver { size: usize, detail: String },

    #[error("GPU error: {0}")]
    Gpu(String),

    #[error("FFT error: {0}")]
    Fft(String),

    #[error("numerical error: {0}")]
    Numerical(String),
}
```

### Step 2: Fix GPU module unwraps

In `src/gpu/mod.rs`:

```rust
// Line 227 — buffer pool access:
// Before:
pool.g_squared_buf.as_ref().unwrap()
// After:
pool.g_squared_buf.as_ref()
    .ok_or_else(|| PwdftError::Gpu("g_squared buffer not allocated".into()))?

// Lines 439-442 — staging buffer readback:
// Before:
sender.send(result).unwrap();
receiver.recv().unwrap().unwrap()
// After:
sender.send(result)
    .map_err(|_| PwdftError::Gpu("GPU readback channel closed".into()))?;
receiver.recv()
    .map_err(|_| PwdftError::Gpu("GPU readback channel failed".into()))?
    .map_err(|e| PwdftError::Gpu(format!("GPU buffer mapping failed: {e}")))?
```

This requires the GPU kernel functions (`hartree_potential`, `lda_xc`, `v_eff_assembly`) to return `Result<Vec<...>>` instead of `Vec<...>`. Callers in `run_scf()` add `?`.

### Step 3: Fix eigensolver expect

In `src/eigensolver/dense.rs`:

```rust
// Line 31:
// Before:
.expect("faer eigendecomposition failed")
// After:
.map_err(|_| PwdftError::Eigensolver {
    size: n,
    detail: "faer returned no eigendecomposition".into(),
})?
```

This requires `diagonalize_hermitian` to return `Result<EigenResult>`. Since it's called in the parallel k-point loop inside `run_scf`, the `par_iter().map()` closure must return `Result<EigenResult>` and be collected with error propagation:

```rust
let kpoint_results: Result<Vec<_>> = kpoints
    .par_iter()
    .map(|kp| -> Result<EigenResult> {
        let mut h = build_hamiltonian_with_v_eff(...);
        vnl.add_to_hamiltonian(&mut h, ...);
        dense::diagonalize_lowest(&h, params.n_bands)
    })
    .collect();
let kpoint_results = kpoint_results?;
```

Rayon's `par_iter().collect()` supports `Result<Vec<T>>` natively.

### Step 4: Fix FFT unwraps

In `src/fft.rs`, the `unwrap()` calls on `Array3::from_shape_vec` (lines 52, 68) are guaranteed safe by the preceding length assertion. Convert to `debug_assert` + comment:

```rust
// The assert at line 50 guarantees data.len() == nx*ny*nz,
// so from_shape_vec cannot fail.
let mut a = Array3::from_shape_vec((nx, ny, nz), data.to_vec())
    .expect("BUG: FFT data length mismatch despite assertion");
```

For `as_slice().unwrap()` (lines 59, 75), ndarray guarantees contiguous layout for standard-order Array3. Add a comment:

```rust
// Array3 with default (row-major) layout is always contiguous.
data.copy_from_slice(b.as_slice().expect("BUG: Array3 should be contiguous"));
```

These are internal invariants, not user-facing errors, so `expect` with a BUG prefix is appropriate — it signals a programming error rather than a runtime condition.

### Step 5: Fix main.rs panics

In `src/main.rs`:

```rust
// Line 92 — missing SCF config:
// Before:
let scf_config = config.scf.as_ref()
    .expect("SCF calculation requires [scf] section...");
// After:
let scf_config = config.scf.as_ref()
    .ok_or_else(|| PwdftError::InvalidInput(
        "SCF calculation requires [scf] section with pseudopotential paths".into()
    ))?;

// Line 106 — missing pseudopotential path:
// Before:
.unwrap_or_else(|| panic!("no pseudopotential path for element {sym}"))
// After:
.ok_or_else(|| PwdftError::MissingPseudopotential(sym.clone()))?;
```

### Step 6: Fix potential/local.rs unwraps

In `src/potential/local.rs`, pseudopotential lookup unwraps should return descriptive errors:

```rust
// Before:
let pp = find_pp_for_z(z, pseudopotentials).unwrap();
// After:
let pp = find_pp_for_z(z, pseudopotentials)
    .ok_or_else(|| PwdftError::MissingPseudopotential(
        format!("Z={z} not found in loaded pseudopotentials")
    ))?;
```

## Migration Order

1. Add new `PwdftError` variants (non-breaking)
2. Fix `main.rs` panics (highest user-facing impact)
3. Fix eigensolver `expect` → `Result` (requires signature change, most callers affected)
4. Fix GPU unwraps → `Result` (requires signature change, behind feature flag)
5. Document remaining `expect` calls in FFT as internal invariants
6. Fix `potential/local.rs` unwraps

## Acceptance Criteria

1. **Zero panics on invalid input:** Malformed TOML, unknown elements, missing pseudopotentials all produce `PwdftError` with descriptive messages, not panics.
2. **GPU failures graceful:** If GPU device is lost or buffer allocation fails, the error propagates to `main()` with a clear message, not a channel panic.
3. **Eigensolver failures graceful:** If faer fails (e.g., non-Hermitian matrix due to a bug), the error includes matrix size and propagates cleanly.
4. **No behavior change on valid inputs:** All existing tests pass unchanged.
5. **Remaining `expect` calls documented:** Any `expect` that remains has a `BUG:` prefix explaining why it cannot fail under normal operation.
6. **`cargo clippy` clean:** No new warnings introduced.
