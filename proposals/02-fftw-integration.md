# Proposal: Replace rustfft with FFTW3

## Motivation

The current FFT module (`src/fft.rs`) implements 3D transforms by batching 1D `rustfft` calls across three sequential sweeps (z, y, x). The y and x dimensions require manual gather/scatter with thread-local buffers, and the x-dimension uses unsafe raw pointer arithmetic (lines 121-141) to work around Rust's borrow checker for parallel strided access.

FFTW3 provides native 3D complex-to-complex transforms with:
- Hardware-optimized plans (SIMD: SSE2/AVX2/NEON on Apple Silicon)
- Cache-oblivious algorithms that eliminate the gather/scatter overhead
- Typically 2-5x faster than rustfft for the 32^3 to 128^3 grids used here
- Well-validated numerics (the reference FFT implementation in computational science)

The density computation (`src/scf/density.rs`) performs `n_bands * n_kpts` inverse FFTs per SCF iteration (line 59), making FFT throughput a direct multiplier on 15-20% of total wall time.

## Scope of Changes

### Dependencies

Remove:
```toml
rustfft = ">=6.2"
```

Add:
```toml
fftw = ">=0.8"
```

System requirement: `brew install fftw` (macOS) or `apt install libfftw3-dev` (Linux). FFTW3 is a near-universal dependency in scientific computing.

### Option A: Direct replacement (recommended)

**File: `src/fft.rs` — full rewrite**

Replace the entire `FFT3D` struct. Current implementation is ~143 lines of manual batching. New implementation:

```rust
use fftw::plan::*;
use fftw::types::*;
use fftw::array::AlignedVec;

pub struct FFT3D {
    dims: [usize; 3],
    plan_fwd: C2CPlan64,
    plan_inv: C2CPlan64,
}

impl FFT3D {
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        // FFTW uses row-major (C) order: dims are [nx, ny, nz]
        let plan_fwd = C2CPlan64::aligned(
            &[nx, ny, nz],
            Sign::Forward,
            Flag::MEASURE,
        ).unwrap();
        let plan_inv = C2CPlan64::aligned(
            &[nx, ny, nz],
            Sign::Backward,
            Flag::MEASURE,
        ).unwrap();
        Self { dims: [nx, ny, nz], plan_fwd, plan_inv }
    }

    pub fn forward(&self, data: &mut [c64]) {
        // FFTW operates in-place on aligned memory
        self.plan_fwd.c2c(data, data).unwrap();
    }

    pub fn inverse(&self, data: &mut [c64]) {
        self.plan_inv.c2c(data, data).unwrap();
    }
}
```

This eliminates:
- The three-pass dimension sweep (lines 68-143)
- All gather/scatter buffers (lines 96-110, 126-141)
- The unsafe raw pointer arithmetic for x-dimension parallelism (lines 121-141)
- Per-worker scratch allocations (lines 78, 96, 126)

**FFTW plan flags:**
- `MEASURE`: spends a few seconds at startup profiling different algorithms for the given grid size. Amortized over many SCF iterations.
- `PATIENT`: more thorough search, ~10s startup. Worth it for production runs.
- `ESTIMATE`: instant startup, ~30% slower execution. Good for tests.

### Option B: Feature-gated (if pure-Rust fallback is desired)

Keep `rustfft` as the default, add `fftw` behind a feature flag:

```toml
[features]
fftw = ["dep:fftw"]
```

Use a trait to abstract the FFT interface. This adds complexity but preserves zero-dependency builds.

### Callers to update

The `FFT3D` API (`forward`, `inverse`, `inverse_normalized`, `dims`, `total_size`) stays identical. Callers don't change. The only type-level difference is if FFTW uses its own aligned buffer type — in that case, add conversion at the boundary.

**Files that call FFT:**
- `src/scf/mod.rs` — `density_r_to_g` (line 443), `real_to_g_space` (line 386)
- `src/scf/density.rs` — `compute_density` (line 59)
- Both construct `FFT3D::new()` — no change needed

### FFTW thread safety

FFTW plan creation is not thread-safe. The `FFT3D::new` call in `compute_density` (line 48, inside rayon fold) creates one planner per worker thread. With FFTW, either:
1. Create the plan once outside the parallel region and share via `Arc` (preferred)
2. Use `fftw_make_planner_thread_safe()` (available in FFTW 3.3.6+)

Option 1 is cleaner: pass `&FFT3D` into the closure and use FFTW's thread-safe execution (only planning is non-thread-safe; execution is safe for distinct buffers).

## Testing Strategy

1. Port all existing FFT tests (`test_fft_roundtrip`, `test_fft_parseval`, grid size tests)
2. Add a cross-validation test: generate random data, FFT with both rustfft and FFTW, compare results to machine epsilon
3. Run full SCF and compare converged energies to ensure numerical equivalence

## Risks

- System dependency: requires FFTW3 installed. Manageable via `brew`/`apt` but is a departure from pure-Rust.
- FFTW planning can be slow with `MEASURE`/`PATIENT` flags. Use `ESTIMATE` in tests, `MEASURE` in production.
- Memory alignment: FFTW prefers 16-byte aligned buffers. The `fftw` crate handles this with `AlignedVec`, but conversion from `Vec<Complex64>` may require a copy.

## Expected Impact

- **Performance:** 2-5x faster FFTs. For a 64^3 grid with 8 bands and 4 k-points, this saves ~32 FFTs per SCF iteration, each ~2-4x faster.
- **Code quality:** Eliminates ~80 lines of manual batching, gather/scatter, and unsafe pointer math.
- **Correctness:** FFTW is the gold standard for FFT correctness in scientific computing.
