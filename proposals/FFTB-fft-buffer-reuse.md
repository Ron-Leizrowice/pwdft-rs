---
id: FFTB
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# FFTB: FFT Buffer Reuse

## Problem

Every call to `FFT3D::forward()` or `FFT3D::inverse()` allocates two `Array3<Complex64>` temporaries:

```rust
let mut a = Array3::from_shape_vec((nx, ny, nz), data.to_vec()).unwrap();  // clone + reshape
let mut b = Array3::zeros((nx, ny, nz));                                   // zero-fill alloc
```

For a 32^3 grid (32,768 Complex64 = 512 KB per array), each FFT call allocates ~1 MB. FFTs are called many times per SCF iteration:
- `density_r_to_g` (density → G-space)
- `real_to_g_space` (V_xc → G-space)
- `fft.inverse` in `compute_density` (per band, per k-point)
- Kerker preconditioning (2 FFTs per mix)
- Convergence energy computation

For nspin=2 this doubles. A typical Si SCF iteration does ~20 FFT calls, meaning ~20 MB of throwaway allocations per iteration.

## Implementation

### Step 1: Add buffer fields to FFT3D

```rust
pub struct FFT3D {
    dims: [usize; 3],
    fwd_handlers: [FftHandler<f64>; 3],
    inv_handlers: [FftHandler<f64>; 3],
    buf_a: Array3<Complex64>,
    buf_b: Array3<Complex64>,
}
```

### Step 2: Allocate once in `new()`

```rust
pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
    Self {
        dims: [nx, ny, nz],
        fwd_handlers: [ /* unchanged */ ],
        inv_handlers: [ /* unchanged */ ],
        buf_a: Array3::zeros((nx, ny, nz)),
        buf_b: Array3::zeros((nx, ny, nz)),
    }
}
```

### Step 3: Rewrite forward/inverse to reuse buffers

```rust
pub fn forward(&mut self, data: &mut [Complex64]) {
    let [nx, ny, nz] = self.dims;
    assert_eq!(data.len(), nx * ny * nz);

    self.buf_a.as_slice_mut().unwrap().copy_from_slice(data);
    self.buf_b.fill(Complex64::new(0.0, 0.0));

    ndfft(&self.buf_a, &mut self.buf_b, &self.fwd_handlers[0], 0);
    ndfft(&self.buf_b, &mut self.buf_a, &self.fwd_handlers[1], 1);
    ndfft(&self.buf_a, &mut self.buf_b, &self.fwd_handlers[2], 2);

    data.copy_from_slice(self.buf_b.as_slice().unwrap());
}
```

Same pattern for `inverse()`.

### Step 4: Update compute_density

`src/scf/density.rs:58` creates a new `FFT3D` per k-point inside the parallel fold. This is fine with the buffer approach — each clone gets its own buffers. No change needed if `FFT3D` derives `Clone`, which it will since `Array3` is `Clone`.

However, the fold identity now allocates the buffers once per thread-chunk rather than per-FFT-call, which is the desired behavior.

## Verification

```bash
cargo test                         # roundtrip, Parseval tests still pass
cargo test --test qe_validation    # physics unchanged
cargo bench --bench scf_benchmarks # measure allocation reduction
```

## Estimated Effort

Under an hour. The ndrustfft API requires `&Array3` inputs (not slices), so the buffer approach is the natural fit.
