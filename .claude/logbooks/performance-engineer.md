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

## 2026-04-17 — XCPR Step 1 implemented (PR #19)

**Scope:** Parallelized `lda_xc_grid` + `lda_xc_spin_grid` with rayon `par_iter`, gated by empirical size threshold `XC_PARALLEL_THRESHOLD = 16 384`. Step 2 (spin-channel diag in `src/scf/mod.rs`) deferred — SPXC branch is active there.

**Calibration (Apple M2, machine lock held):**
- Rayon fork/join/unzip overhead on this machine: ~70 µs per region.
- Sequential `lda_xc_grid`: ~11 ns/point. Sequential `lda_xc_spin_grid`: ~26 ns/point.
- n=4096 naive parallel: 113 µs (vs 43 µs serial) — 161% regression. Confirmed threshold must be above 4 k.
- n=16 384 parallel: 160 µs unpolarized, 197 µs spin — net positive, modest.
- n=32 768: 2.09x unpolarized, 3.45x spin.
- n=262 144 (64³): 6.86x unpolarized, 10.07x spin (3.0 ms → 442 µs; 7.8 ms → 777 µs).

**Observations worth flagging:**
- Revised prior claim: XCPR is not "marginal" — for 48³ and 64³ grids, spin XC cost drops by ~10x. For 32³ (typical production) it's still 2-3x. Small/gate grids are guarded by the threshold.
- Grid sizes 24³ (13 824 pts) fall *just below* the threshold. If we find real SCF calls spending time at 24³ XC, consider lowering threshold to 8 192 — but `lda_xc_grid` would likely regress.
- Competing `spin_polarization` integration test runs ~3-5 min at 1400% CPU — a background agent ran it during my first benchmark pass and corrupted the timings. Always verify load avg < 5 before trusting criterion numbers. Re-ran bench clean after it exited.

**Tangential ideas:**
- The n=4096 pure-sequential cost is 43 µs; rayon costs ~70 µs. A thread-pool-reuse scheme (pin workers, skip fork/join) could lower the breakeven to ~1 000 pts. Not worth it for XC alone — but same pattern applies to Hartree, V_eff assembly, density reconstruction — might be a cross-cutting proposal ("hot-loop rayon pool").
- XC still dominated by `cbrt` and `ln`. CBRT proposal (replace cbrt with Newton-refined f64::from_bits hack) could stack with this.

**Next session TODO:**
- Step 2 after SPXC merges.
- Profile full SCF end-to-end to confirm XCPR net impact on wall time (XC is maybe 5-10% of SCF, so expected overall SCF speedup at 32³: ~3-5%, at 64³: ~8-15%).

## 2026-04-17 — FMAD applied (PR pending)

**Scope:** 24 `mul_add` substitutions in hot-path kernels. Test code + out-of-scope files (mod.rs, fft.rs, symmetry/kpoints.rs) excluded per task brief.

- `src/potential/xc.rs`: 13 sites (PZ correlation, spin interpolation, Slater exchange)
- `src/numerics.rs`: 4 sites (Simpson boundary terms)
- `src/scf/smearing.rs`: 3 sites (MP, cold, entropy)
- `src/potential/nonlocal.rs`: 2 sites (spherical Bessel + Legendre recurrences)
- `src/scf/mixing.rs`: 2 sites (Broyden linear + corrected mix)

**Gotcha:** `ec_ha = (c * rs).mul_add(ln_rs, a.mul_add(ln_rs, b))` forced type annotation on `a, b, c, d: f64` (untyped literals no longer deducible through mul_add).

**Measurements (Apple M2, machine lock held):**
- `lda_xc_grid_n256`: 2.89 → 2.80 µs (-2.9%)
- `lda_xc_grid_n512`: 5.74 → 5.50 µs (-4.2%)
- `lda_xc_grid_n4096`: 45.97 → 44.16 µs (-3.9%)
- `lda_xc_grid_n16384`: 160 → 154 µs (-3.9%)
- `lda_xc_spin_grid_*`: within noise (~±2%). Expected — spin kernel dominated by `cbrt`, not additive ops.
- Parallel sizes (n≥32k) noisy due to rayon fork/join variance; another agent's debug tests occasionally pollute even under lock (compilation happens outside lock).

**Takeaways:**
- FMA win is measurable but modest (~3-4%) on unpolarized XC, consistent with the fact that cbrt/ln dominate. CBRT proposal would stack multiplicatively.
- All 81 suboptimal_flops warnings reduced to 54, with all remaining in test code, mod.rs orchestration, or symmetry/pseudopot (out of scope).
- No tests broke; relative_eq tolerances absorb the ULP-level bit diff from FMA.
- For future bench runs: watch for concurrent debug test builds at peak (spin_polarization at 1400% CPU masked by lock scheduler). Sequential benchmarks (n ≤ 16 384) are robust; parallel benchmarks (n ≥ 32 768) require quiet machine for reliable criterion stats.

## 2026-04-17 — Post-FFTB/FMAD profiling pass; new ITEV proposal

Task: identify next perf target after FFTB/FMAD landed. Lock held 702 s, released before writing.

**Baseline re-measured (criterion, machine-locked, Apple M2):**

| n_pw | `faer_eigen` | `vnl_apply` | `vnl_new` | `kinetic` |
|------|--------------|-------------|-----------|-----------|
| 89   |   1.07 ms    |  0.53 ms    |   7.17 ms |  2.77 µs  |
| 259  |  56.3  ms    |  4.70 ms    |  22.1  ms | 22.6  µs  |
| 725  | 835.8  ms    | 28.6  ms    |  47.2  ms | 166   µs  |

Prior logbook's `n725 = 73.5 ms` was noise. O(n³) scaling confirmed.

**End-to-end SCF profile** (`sample` 4 s, 5×Si ecut=200 runs, n_pw=259, 31 197 samples):
- Idle: 39.7 %
- Rayon pool overhead (truncated stacks): 27.6 %
- Classified user code: 32.7 %
  - Eigensolver (faer+gemm): **60.7 % of user CPU**
  - V_NL apply: 7.1 %
  - FFT: 2.1 %; build_H: 1.2 %; symm: 0.1 %; XC: 0.0 %
- Re-attributing rayon pool overhead by call-site → **eigensolver ≈ 85-90 % of SCF user CPU**.

**Decisive finding:** `faer 0.24` ships
`faer::matrix_free::eigen::partial_self_adjoint_eigen` — an upstream
implicitly-restarted Arnoldi partial Hermitian eigensolver with `LinOp`
and `v0` warm-start. This obsoletes DVSN's hand-rolled Davidson plan and
drops the complexity from `large` to `medium`.

**Proposal opened:** `ITEV-faer-partial-eigen.md` (new). Supersedes
DVSN; WFRX remains as warm-start for ITEV. Projected 2.5-4× SCF speedup
at production sizes (n_pw ≥ 259), 10-50× on the eigensolve step at
n_pw = 725. Shift-and-flip (`A' = σI − H`) to map largest-magnitude →
algebraically-lowest. Keep dense backend as selectable fallback.

**Tangential ideas:**
- V_NL apply @ n_pw=725 is 29 ms = #2 bottleneck after eigensolve. Inner
  double-loop `for ig × for jg × for i_proj × for j_proj` is O(n_pw² · n_proj²).
  Worth a proposal *after* ITEV lands: projector-sum lifted to O(n_pw · n_proj²)
  blocked-matmul (V_NL = Σ_α β_α D_α β_α^H can be GEMM'd).
- `vnl_new` @ n_pw=725 is 47 ms but only runs ONCE per k-point per SCF
  (cached in ScfContext) — not a hot target.
- Rayon pool overhead 27.6 % of samples is high; likely from faer's inner
  spindle + outer k-point par_iter double-nesting. ITEV should drop this
  because each k-point call becomes much shorter → less in-faer parallelism.

**Next session TODO:**
- Wait for EM decision on ITEV priority/approval.
- If ITEV approved: Phase 1 wrapper + correctness tests.
- If blocked: consider V_NL blocked-matmul as a parallel-track proposal.
