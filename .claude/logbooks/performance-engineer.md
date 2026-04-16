# Performance Engineer Logbook

Entries: date, measurements (actual numbers), bottleneck findings, proposals assessed. Always include hardware context.

## 2026-04-16 — Orientation (no benchmarks run yet)

**Proposal audit (code inspection only, not profiled):**
- FFTB: well-motivated, FFT allocates fresh Array3 every call. Impact unknown without profiling.
- XCPR: XC parallelization straightforward but may not help at small grid sizes (4096 points). Spin-channel join is a clean win.
- WFRX: theoretically sound but 25x claim is optimistic. Real gain likely 5-15% net.
- DVSN: strategically correct long-term target, not highest near-term impact.

**Gap found:** No proposal for allocator selection (jemalloc/mimalloc). macOS libmalloc underperforms under multithreaded pressure. Trivial integration via feature flag. Worth proposing.

**TODO next session:** Run `cargo bench`, profile SCF with Instruments.app, establish actual bottleneck distribution before proposing anything.

## 2026-04-16 — Baseline Profiling (machine-lock held, no contention)

**Hardware:** Apple M2, macOS Darwin 25.3.0. All runs with machine lock acquired.

### End-to-end SCF — Si FCC 2 atoms (pre-compiled binary)

| Config | ecut (Ry) | k-grid (irr.) | Iters | Wall | User | Peak RSS | Peak Footprint |
|--------|-----------|---------------|-------|------|------|----------|----------------|
| si_scf.yaml | 100 | 2×2×2 (3) | 9 | 0.11s | 0.19s | 12 MB | 10 MB |
| si_scf_converged.yaml | 200 | 4×4×4 (10) | 9 | 0.34s | 2.09s | 75 MB | 73 MB |

Memory scales ~6x (3→10 k-points, 3x basis). Dominated by per-k-point Hamiltonian matrices (n_pw² complex).
Rayon parallelism effective: 2.09s user / 0.34s wall ≈ 6x speedup across 10 k-points.

### Criterion Microbenchmarks

**Eigensolver (faer Hermitian, SCF bottleneck, O(n³)):**

| n_pw | Time |
|------|------|
| 89 | ~1.2–1.9 ms (high variance at 20 samples) |
| 259 | 7.36 ms |
| 725 | 73.5 ms |

Scaling: 259→725 is 10x for 2.8x n — consistent with O(n³).

**Hamiltonian construction:**

| n_pw | Kinetic | V_NL build | V_NL apply |
|------|---------|------------|------------|
| 89 | 2.1 µs | 5.27 ms | 388 µs |
| 259 | 16.9 µs | 15.1 ms | 3.22 ms |
| 725 | 109 µs | 42.3 ms | 25.9 ms |

**FFT (3D complex, ndrustfft):**

| Grid | Forward | Roundtrip |
|------|---------|-----------|
| 16³ | 28.7 µs | 57.7 µs |
| 20³ | 112 µs | 225 µs |
| 24³ | 94 µs | 190 µs |
| 32³ | 349 µs | 704 µs |
| 48³ | 1.65 ms | 3.32 ms |

Note: 24³ faster than 20³ — likely favorable FFT factorization (24 = 2³×3).

### Bottleneck Ranking (at n_pw=725, production size)

1. **Eigensolver: 73.5 ms** — dominant cost, called once per k-point per SCF step
2. **V_NL build: 42.3 ms** — #2, once per k-point per SCF step
3. **V_NL apply: 25.9 ms** — #3, adds KB projectors to Hamiltonian
4. **FFT: ~3.3 ms** — cheap, <5% of one eigensolve even at 48³
5. **Kinetic/basis: negligible** — µs-scale

Per SCF iteration at converged settings (10 irr. k-points): ~730 ms eigensolver + 423 ms V_NL build + 259 ms V_NL apply ≈ 1.4s compute, 9 iterations ≈ 12.7s serial, wall 0.34s with rayon.

### Implications for Proposals

- **FFTB** (FFT buffer reuse): confirmed low-impact. FFT is <5% of cost. Worth doing for cleanliness but not a performance win.
- **XCPR** (XC parallelization): grid ops are cheap; marginal gain expected.
- **WFRX** (wavefunction reuse): could skip V_NL rebuild — would save 42 ms/k-point/iter if applicable.
- **Eigensolver optimization** is the highest-impact target (not yet proposed). Iterative solvers (Davidson/LOBPCG) with warm-start from previous SCF step would slash the dominant cost.
- **Allocator proposal** still relevant — 75 MB peak with rayon threads means malloc pressure.

**TODO next session:** Profile with Instruments.app for per-function breakdown within eigensolver and V_NL build. Consider proposing iterative eigensolver.
