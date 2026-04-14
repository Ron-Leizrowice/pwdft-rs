# Proposal: Replace nalgebra-lapack with faer

**Status: COMPLETED**

## Result

Replaced LAPACK zheev (unsafe FFI via nalgebra-lapack + Apple Accelerate) with faer (pure Rust, safe, multi-threaded). Full migration: `DMatrix<Complex64>` → `faer::Mat<Complex64>` across 5 source files and 4 test files. Removed `nalgebra-lapack` and `lapack` dependencies entirely. nalgebra retained for `Vector3<f64>` geometry only.

## Benchmark (M2 Mac, criterion, release)

| Matrix size | LAPACK/Accelerate | faer | Speedup |
|-------------|-------------------|------|---------|
| n=89 (ecut=100) | 949 µs | 906 µs | 1.05× |
| n=259 (ecut=200) | 12.7 ms | 6.9 ms | 1.84× |
| n=725 (ecut=400) | 206 ms | 74 ms | 2.8× |
| n=1363 (ecut=600) | **CRASHES** | works | — |

## Key findings

- faer's advantage grows with matrix size due to better multi-threading
- Apple Accelerate crashes at n=1363 due to a libc++ TMO bug in release mode
- `faer::c64` is literally `num_complex::Complex64` (same type alias) — zero conversion needed for individual values
- The `profile.dev.package` optimization for faer/gemm/pulp is essential — 230× debug-mode slowdown without it
- Eigenvalues match LAPACK to 1e-10 on realistic Si Hamiltonians; eigenvectors verified unitary

## Commits

- `1e72c50` Add faer eigensolver with benchmark: 1.8-2.8× faster than Accelerate
- `e68a9c1` Replace LAPACK eigensolver with faer, remove Accelerate dependency
