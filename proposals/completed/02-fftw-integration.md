# Proposal: Replace rustfft with FFTW3

**Status: COMPLETED (different outcome than proposed)**

## Original proposal

Replace the hand-rolled rustfft 3D FFT with FFTW3 for 2-5× speedup.

## What actually happened

Three iterations:

1. **FFTW with copy**: Added FFTW behind `--features fftw`. Fast at small grids (11× at 20³) but slower at 32³+ due to aligned buffer copy overhead.

2. **FFTW zero-copy + threading**: Eliminated copies via in-place `fftw_sys` FFI. Added adaptive threading (1 thread below 32³, all cores above). FFTW won at every size: 3-19×.

3. **Safety audit → ndrustfft**: Audit found 12 unsafe sites across both backends — strict provenance violations, planning buffer use-after-free risk, raw pointer arithmetic. Replaced the entire hand-rolled rustfft backend with **ndrustfft** (safe ndarray-based 3D FFT wrapper), achieving zero unsafe in the default backend. Then removed FFTW entirely to go pure Rust.

## Final result

Pure Rust, zero unsafe, zero system dependencies. Uses ndrustfft (ndarray + rustfft).

## Benchmark (M2 Mac, criterion, release)

| Grid | Old rustfft (unsafe) | ndrustfft (safe) | Notes |
|------|---------------------|-----------------|-------|
| 16³ | 174 µs | 27 µs | 6× faster |
| 20³ | 236 µs | 108 µs | 2× faster |
| 24³ | 231 µs | 94 µs | 2.5× faster |
| 32³ | 278 µs | 329 µs | ~same |
| 48³ | 457 µs | 1600 µs | slower (no rayon) |

ndrustfft is faster at typical DFT grid sizes (16³-24³) and comparable at 32³. The 48³ regression is single-threaded vs the old rayon-parallel approach, but irrelevant since the eigensolve dominates wall time.

## Key findings

- FFTW is genuinely faster (3-19×) but adds system dependencies and unavoidable FFI unsafe
- ndrustfft provides safe 3D FFT with good performance at typical DFT sizes
- The strict provenance violation in the old rustfft x-dimension transform (`ptr as usize → usize as ptr`) was a real correctness risk
- The FFTW planning buffer bug (Vec deallocated while FFTW may reference it) was a latent memory safety issue

## Commits

- `e8718b1` Add FFTW3 backend behind feature flag, keep rustfft as default
- `f61dbe8` FFTW: zero-copy in-place transforms via fftw_sys
- `2babbea` FFTW: add adaptive threading, wins at all grid sizes
- `13a170b` Replace unsafe rustfft 3D with safe ndrustfft backend
- `4ba7314` Remove FFTW backend, go pure Rust with ndrustfft
