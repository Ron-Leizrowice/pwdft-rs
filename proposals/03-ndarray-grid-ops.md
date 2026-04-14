# Proposal: Adopt ndarray for grid operations

## Motivation

Grid-based data (densities, potentials, FFT buffers) is currently stored as flat `Vec<f64>` or `Vec<Complex64>` with manual 3D indexing throughout the codebase. Examples:

- `src/scf/mod.rs` line 113: `i1 * self.dims[1] * self.dims[2] + i2 * self.dims[2] + i3`
- `src/scf/density.rs` line 62: flat iteration over `psi_g` with separate index tracking
- `src/scf/mixing.rs` lines 64-79: manual dot products via `fn dot(a: &[f64], b: &[f64])`
- `src/potential/xc.rs`: element-wise operations on flat arrays

`ndarray` provides typed multi-dimensional arrays with:
- 3D indexing (`arr[[ix, iy, iz]]`) that catches out-of-bounds at the type level
- BLAS-backed `dot()` for the mixing module
- `Zip` for fused parallel element-wise operations (avoids intermediate allocations)
- `ArrayView` for zero-copy slicing
- Interoperability with FFTW via `ndarray` views (if Proposal 02 is adopted)

## Scope of Changes

### Dependencies

Add:
```toml
ndarray = { version = ">=0.16", features = ["rayon"] }
```

The `rayon` feature enables `par_azip!` for parallel element-wise ops.

### Phase 1: Grid types (incremental, no breaking changes)

**File: `src/scf/mod.rs` — `FftGrid`**

Add helper methods that return `ndarray::Array3` views:

```rust
use ndarray::Array3;

impl FftGrid {
    /// Wrap a flat density vector as a 3D array view.
    fn as_array3<'a>(&self, data: &'a [f64]) -> ndarray::ArrayView3<'a, f64> {
        let [nx, ny, nz] = self.dims;
        ndarray::ArrayView3::from_shape((nx, ny, nz), data).unwrap()
    }

    fn as_array3_mut<'a>(&self, data: &'a mut [f64]) -> ndarray::ArrayViewMut3<'a, f64> {
        let [nx, ny, nz] = self.dims;
        ndarray::ArrayViewMut3::from_shape((nx, ny, nz), data).unwrap()
    }
}
```

This lets existing `Vec<f64>` data be viewed as 3D without allocation, enabling gradual migration.

### Phase 2: Mixing module

**File: `src/scf/mixing.rs`**

Replace the hand-rolled `dot` function (line 102-104) with ndarray's BLAS-backed dot:

```rust
// Before:
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

// After:
use ndarray::ArrayView1;
fn dot(a: &[f64], b: &[f64]) -> f64 {
    ArrayView1::from(a).dot(&ArrayView1::from(b))
}
```

This uses BLAS `ddot` under the hood when `ndarray` is compiled with BLAS support, and is auto-vectorized otherwise.

The history vectors (`history_in`, `history_res`) could also become `Vec<Array1<f64>>` for cleaner residual computation (lines 64-68):

```rust
// Before:
let dr_i: Vec<f64> = self.history_res[i].iter()
    .zip(self.history_res[last].iter())
    .map(|(&a, &b)| a - b).collect();

// After:
let dr_i = &self.history_res[i] - &self.history_res[last];
```

### Phase 3: XC and Hartree grid operations

**File: `src/potential/xc.rs`**

The `lda_xc_grid` function maps over a flat `&[f64]`. With ndarray + rayon:

```rust
use ndarray::parallel::par_azip;

pub fn lda_xc_grid(rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let rho = ArrayView1::from(rho_r);
    let mut exc = Array1::<f64>::zeros(rho.len());
    let mut vxc = Array1::<f64>::zeros(rho.len());

    par_azip!((r in &rho, e in &mut exc, v in &mut vxc) {
        let (ei, vi) = lda_xc(*r);
        *e = ei;
        *v = vi;
    });

    (exc.into_raw_vec(), vxc.into_raw_vec())
}
```

**File: `src/scf/mod.rs` — `assemble_v_eff` (lines 398-409)**

```rust
// Before: three parallel zips producing a new Vec
// After:
fn assemble_v_eff(v_local: &[Complex64], v_h: &[Complex64], v_xc: &[Complex64]) -> Vec<Complex64> {
    let mut result = Array1::from(v_local.to_vec());
    let vh = ArrayView1::from(v_h);
    let vxc = ArrayView1::from(v_xc);
    azip!((r in &mut result, &h in &vh, &x in &vxc) *r += h + x);
    result.into_raw_vec()
}
```

### Phase 4: Density symmetrization

**File: `src/symmetry/density.rs`**

Density symmetrization operates on a flat buffer with 3D index arithmetic. Wrapping as `Array3` makes the rotation/mapping logic clearer and bounds-checked.

## Migration Strategy

ndarray is additive — it wraps existing `Vec` data via views without changing storage. Migration can be done module-by-module:

1. Add the dependency, create array views in `FftGrid`
2. Migrate `mixing.rs` (smallest scope, easy to validate)
3. Migrate XC and Hartree grid ops
4. Migrate density symmetrization
5. Optionally, change storage types from `Vec<f64>` to `Array1<f64>` for the density/potential buffers that persist across SCF iterations

## Risks

- Adds a dependency for convenience rather than capability. All current code works correctly.
- ndarray's `Array` is heap-allocated like `Vec` — no performance difference for element access. The benefit is ergonomic and for BLAS-backed operations.
- If combined with faer (Proposal 01), there are two matrix ecosystems. ndarray is for grid data (1D/3D real arrays); faer is for dense complex matrices (Hamiltonian/eigenvectors). Clean separation.

## Expected Impact

- **Performance:** Modest. BLAS dot products in mixing are faster for large grids. `par_azip!` fuses operations to avoid intermediate allocations.
- **Code quality:** Eliminates manual 3D index arithmetic. Catches dimension mismatches at compile time.
- **Interop:** If FFTW is adopted (Proposal 02), `ndarray` views can wrap FFTW-aligned buffers directly.
