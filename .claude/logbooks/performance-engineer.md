# Performance Engineer Logbook

Entries: date, measurements (actual numbers), bottleneck findings, proposals assessed. Always include hardware context.

## 2026-04-16 — Orientation

Proposal audit (inspection only): FFTB well-motivated (fresh `Array3` per call); XCPR straightforward, may not help small grids (4096 pts); WFRX sound but 25× claim optimistic (likely 5–15% net); DVSN long-term. Gap: no allocator proposal (jemalloc/mimalloc) — macOS libmalloc underperforms under multithreaded pressure.

## 2026-04-16 — Baseline profiling (Apple M2, lock held)

End-to-end SCF on Si FCC 2 atoms: `si_scf.yaml` (ecut=100 Ry, 2×2×2) 0.11 s wall, 12 MB peak; `si_scf_converged.yaml` (ecut=200 Ry, 4×4×4) 0.34 s wall, 75 MB peak. Rayon ≈6× across 10 k-points. Memory dominated by per-k Hamiltonian (n_pw² complex).

**Bottleneck ranking at n_pw=725:** eigensolver 73.5 ms (dominant); V_NL build 42.3 ms; V_NL apply 25.9 ms; FFT ~3 ms; kinetic/basis µs. Per-iter at converged settings: ~1.4 s compute, ~12.7 s serial, 0.34 s wall with rayon. Note: n_pw=725 eigensolver later re-measured at 836 ms — the 73.5 ms figure was criterion noise (see 2026-04-17 ITEV entry).

**Implications:** FFTB is <5% of cost (do for cleanliness); XCPR marginal; WFRX could skip V_NL rebuild (~42 ms/k/iter); eigensolver optimisation is the highest-impact unproposed target; allocator still relevant at 75 MB peak.

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

## 2026-04-18 — ITEV Phase 1+2 landed (PR #45)

**Scope:** Wrapper (`src/eigensolver/iterative.rs`) + SCF dispatch + YAML switch + 8 correctness tests. Default still Dense. Iterative is opt-in.

**Correctness proven (unit-level):**
- 8 tests pass, incl. real Si diamond Γ-point H (n_pw=89, triple-degenerate valence) vs dense to < 1e-6 eV.
- Degeneracy collapse fixed via `n_request = max(n_bands+4, 1.5*n_bands)`, keep lowest n_bands after over-request.
- Shift strategy: tight per-row Gershgorin bound (`max_i (diag_i + Σ_{j≠i} |H_{ij}|) + ε`). Initial loose shift caused Lanczos reorthogonalization hangs; tightening resolved them for single calls.

**BLOCKER: upstream faer 0.24 bug.** `operator/self_adjoint_eigen/iterate_lanczos` (lines 42-59) inner Gram-Schmidt loop spins indefinitely when a Krylov vector becomes numerically null. Reproduces non-deterministically on multi-iteration Si SCF at n_pw=89 — single calls are fine, SCF hangs after ~3-10 iters. Confirmed via `sample`: stuck in `norm_l2_simd_pairwise_rows`. No user-code workaround; we need upstream fix or a replacement (SPRS matrix-free).

**What this means for ITEV:**
- Bench `iterative_cold_n{N}` would hang; disabled.
- End-to-end SCF integration test is `#[ignore]`d.
- Cannot flip default to Iterative until upstream fixed.
- Correctness scaffolding is in place — the shift/over-request/fallback logic is validated and ready.

**Measurements (sample-based, n=89 Si H, single call):**
- Iterative: 2-5 ms per diagonalize (consistent).
- Dense: similar at this size (n < 100 dense Hessenberg dominates).
- True speedup only visible at n_pw ≥ 259 — where the SCF hang also manifests, so benchmarking blocked.

**Tangential ideas:**
- File issue upstream against faer re: `iterate_lanczos` reorthogonalization fragility; link to Horst & Meurant (2000) on selective reorthogonalization.
- Alternative: implement our own matrix-free partial eigensolver (SPRS could revive DVSN's Davidson-style plan on top of our own LinOp). Bigger scope.
- Alternative: investigate LOBPCG via faer or hand-rolled — avoids Lanczos entirely.
- V_NL blocked-matmul remains a parallel-track target (29 ms at n_pw=725, #2 bottleneck after eigensolve).

**Files touched:** src/eigensolver/{iterative.rs (new), mod.rs}, src/scf/mod.rs, src/settings.rs, benches/scf_benchmarks.rs, tests/itev_iterative_eigensolver.rs (new), Cargo.toml (+ dyn-stack).

## 2026-04-18 — VNLM landed (PR #49)

**Scope:** `src/potential/nonlocal.rs` only. Replaced nested (ig,jg,i,j)
scalar loop with `H += B · D · B^H` single-GEMM. Expanded KB projector
matrix B has one channel per (atom, projector, m); D is block-diagonal;
addition theorem `Σ_m Y_lm Y*_lm = (2l+1)/(4π) P_l(cos θ)` makes it exact.

**Measured (Apple M2, lock held):**

| n_pw | vnl_apply | vnl_apply | speedup | vnl_new | vnl_new |
|------|-----------|-----------|---------|---------|---------|
|      | before    | after     |         | before  | after   |
| 89   | 398 µs    | 80 µs     | 5.0×    | 5.41 ms | 7.08 ms |
| 259  | 3.18 ms   | 776 µs    | 4.1×    | 15.5 ms | 21.0 ms |
| 725  | 27.4 ms   | 5.24 ms   | 5.2×    | 43.2 ms | 78.5 ms |

Profile matched logbook claim exactly (27 ms ≈ 29 ms). Ratio logbook
quoted was slightly higher at measurement time; current n_pw=725 is 27 ms.

**Surprise:** `vnl_new` regressed ~1.5× because we now build the full B
and D·B^H at construction. One-time cost per k-point (cached in
`ScfContext.vnl_cache`). Break-even at ~2 SCF iters; 15-iter SCF net
2.9× faster V_NL per k-point. Worth it.

**Correctness:** 265 tests pass both feature flags. Analytic G=0 check
(tol 1e-8), hermiticity (1e-10), shell degeneracy (1e-8), Python
reference cross-check all pass. Added `test_ylm_addition_theorem` pinning
the algebraic identity to 1e-12.

**Implementation note:** real Y_lm table built per k-point via QE's
ylmr2 recurrence. Only lmax = max-over-PPs l. Memory ~n_pw × (lmax+1)²
f64 + 2 × n_pw × n_channels complex — ~200 KB at n_pw=725 — negligible.

**Tangential ideas:**
- vnl_new regression is from B-fill + D·B^H construction loops. Could
  parallelize with rayon over atoms or use a blocked matmul for D·B^H.
  Not worth it — we're below 80 ms at n_pw=725 and it's already out of
  the SCF hot path.
- Next perf target after eigensolver: profile whole-SCF again post-ITEV
  (when iterative fixes land) to see if V_NL moved off the podium.
- ITEV upstream faer issue still blocking iterative eigen — independent.

## 2026-04-18 — Post-MXBA headline benchmark pass

Bench-only session, read-only source. Lock held 546 s. Full report at
`proposals/completed/PERF-2026-04-18-benchmark-pass.md`.

**Headline:** Si n_pw=725 hot-path (eig + vnl_apply) **1.47x** per iter
(110.9 ms -> 75.3 ms). V_NL-only, 15-iter amortized: **4.32x**
(430.8 ms -> 99.7 ms per k-point). Entire measurable speedup is VNLM.

**Measurements (Apple M2, lock held, criterion 100-sample):**

| Bench              | 2026-04-16 | 2026-04-18 | Delta |
|--------------------|------------|------------|-------|
| vnl_apply_n89      |  391 µs    |  67.7 µs   | 5.8x  |
| vnl_apply_n259     | 3.22 ms    |  485 µs    | 6.6x  |
| vnl_apply_n725     | 25.9 ms    | 3.75 ms    | 6.9x  |
| vnl_new_n725       | 42.3 ms    | 43.4 ms    | +2.5% |
| faer_eigen_n89     | 1.52 ms    |  931 µs    | -39%  |
| faer_eigen_n259    | 7.36 ms    | 10.04 ms   | +36%* |
| faer_eigen_n725    | 85.0 ms    | 71.6 ms    | -16%* |
| lda_xc_grid_n4096  |  45.97 µs  | 42.7 µs    | -7%   |
| scf_iter_20x_32^3  |  n/a       | 13.85 ms   | new   |

*= unexplained; see anomalies.

**Reproducibility notes:**
- VNLM reproduced above PR #49's claim (6.9x actual vs 5.2x at n=725).
- VNLM's reported 1.5x vnl_new regression at n=725 (43.2->78.5 ms)
  **did not reproduce**; today reads 43.4 ms, essentially pre-VNLM.
  Either MODR-B hid the cost or original regression was variance.
- FMAD (3-4% on lda_xc_grid_*) survived MODR refactors.
- XCPR parallel regime (n>=16384) still active.
- FFTB's benefit invisible on standalone forward_* (within 1-2% of
  2026-04-16); use `scf_iter_20x_*` for FFTB impact — no pre-FFTB
  reference bench exists in the repo.

**Current bottleneck ranking at n_pw=725 (post-VNLM):**
1. Eigensolve (dense faer) — 71.6 ms/iter, ~95% of per-iter cost.
   Only ITEV moves this, blocked on faer 0.24 Lanczos upstream bug.
2. V_NL new — 43.4 ms one-time per k-point, amortizable.
3. V_NL apply — 3.75 ms/iter, no longer a target.
4. FFT — ~14 ms for 20 calls at 32^3.
5. XC grid — sub-ms at all realistic sizes.

**Anomalies flagged for FLUP:**
- **ANOM-1 / faer_eigen_n259 +36%:** wide CI (+/-15%); re-bench with
  --measurement-time 15 to confirm. If real, bisect today's 13 PRs.
- **ANOM-2 / faer_eigen_n725 -16%:** unclaimed free win; verify.
- **ANOM-3 / vnl_new_n725 no regression:** VNLM's break-even caveat
  may be stale. Proposed FLUP ID **VNLT**.

**Tangential ideas:**
- Per-iter ratio eig:vnl_apply jumped ~3:1 -> ~19:1. Anything that
  moves eigensolver off top (ITEV, WFRX subspace reuse) is now
  overwhelmingly the highest-leverage perf work.
- n_pw=259 faer_eigen wide CI suggests rayon global pool + inner faer
  parallel contention; worth a targeted profile with ANOM-1.

**Next session TODO:**
- Revisit when ITEV unblocks (faer upstream fix).
- If EM wants ANOM-1 investigated, re-bench longer.
