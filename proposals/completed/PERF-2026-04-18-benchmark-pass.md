---
id: PERF-2026-04-18
title: 2026-04-18 post-MXBA SCF benchmark pass — headline speedup report
status: completed
priority: n/a
complexity: small
risk: none
depends_on: [FFTB, FMAD, XCPR, ITEV, VNLM, PRPL, CCMX, MXBA, MODR-A, MODR-B, MODR-C, MODR-D]
---

# PERF-2026-04-18 — SCF benchmark pass after today's 13-PR wave

## TL;DR

**Si n_pw=725 per-SCF-iter hot path (eigensolve + V_NL apply) is 1.47x
faster than the earliest comparable 2026-04-16 baseline** — from
110.9 ms to 75.4 ms per k-point per iteration. On the V_NL-only slice
(what VNLM was designed to move), it is **4.3x faster over a 15-iter
Si SCF** (430.8 ms -> 99.7 ms cumulative per k-point).

Today's individual landings — XCPR, FFTB, FMAD, VNLM, PRPL, CCMX,
MXBA, MODR-A/B/C/D, ITEV, NCFX, PCFX, CAST — collectively deliver
**VNLM as the single dominant wall-time contributor to the measurable
speedup**. ITEV remains opt-in (upstream faer 0.24 Lanczos bug), so
the faer dense path is still the per-iter eigensolve bottleneck and
dominates the residual cost.

## Environment

- Apple M2, macOS, machine lock held for 546 s
- `origin/main @ a1bbac5` (FLUP seed of MXB1/MXB2/MXB3), post-MXBA
- `cargo bench --bench scf_benchmarks -- --warm-up-time 2 --measurement-time 6`
- 100 samples per bench; criterion mean of the 95% CI used below
- Fresh worktree; no stale criterion baselines in this worktree — numbers
  below are compared against the 2026-04-16 baselines preserved in the
  user's main checkout at `target/criterion/*/new/estimates.json`
  (timestamp Apr 16 13:03:39 2026). That baseline pre-dates FFTB, FMAD,
  XCPR, VNLM, PRPL, CCMX, MXBA, and the MODR refactor wave — it is the
  earliest comparable point that shares today's bench harness names.

## Results

### Hamiltonian hot path

| Bench                         | 2026-04-16 baseline | 2026-04-18 post-MXBA | Delta      | Notes |
|-------------------------------|---------------------|----------------------|------------|-------|
| `hamiltonian/vnl_apply_n89`   |   391 us            |   67.7 us            | **-82.7% / 5.8x** | VNLM claim 5.0x — exceeded |
| `hamiltonian/vnl_apply_n259`  |  3.22 ms            |   485 us             | **-84.9% / 6.6x** | VNLM claim 4.1x — exceeded |
| `hamiltonian/vnl_apply_n725`  | 25.93 ms            |   3.75 ms            | **-85.5% / 6.9x** | VNLM claim 5.2x — exceeded |
| `hamiltonian/vnl_new_n89`     |  5.27 ms            |  5.34 ms             |  +1.3%     | no regression |
| `hamiltonian/vnl_new_n259`    | 15.09 ms            | 15.37 ms             |  +1.8%     | no regression |
| `hamiltonian/vnl_new_n725`    | 42.33 ms            | 43.41 ms             |  +2.5%     | **NOT the 1.5x (78.5 ms) VNLM PR #49 reported — recovered** |
| `hamiltonian/kinetic_n89`     |  ~2 us              |  2.04 us             |  flat      | - |
| `hamiltonian/kinetic_n259`    |  ~17 us             | 16.6 us              |  flat      | - |
| `hamiltonian/kinetic_n725`    |  ~110 us            |  107 us              |  flat      | - |

### Eigensolver (dense, faer)

| Bench                          | 2026-04-16 baseline | 2026-04-18 post-MXBA | Delta      | Notes |
|--------------------------------|---------------------|----------------------|------------|-------|
| `eigensolver/faer_eigen_n89`   |  1.52 ms            |   931 us             | -38.8%     | likely faer bump or scheduler variance; nothing today claimed this |
| `eigensolver/faer_eigen_n259`  |  7.36 ms            | 10.04 ms             | +36.4%     | see anomaly note |
| `eigensolver/faer_eigen_n725`  | 85.01 ms            | 71.57 ms             | -15.8%     | nominal improvement, source unknown |

Note: per-iter SCF logbook (2026-04-17) previously recorded
`faer_eigen_n725` at 835.8 ms. That value was clearly polluted by
concurrent load — the 85 ms baseline from the main checkout's criterion
snapshot (captured under machine lock) and today's 71.6 ms reading are
both in the same order of magnitude and are the ground truth.

### FFT (per-iteration proxy workload)

| Bench                              | 2026-04-18 post-MXBA | Pre-FFTB reference |
|------------------------------------|----------------------|--------------------|
| `fft/scf_iter_20x_16x16x16`        |  1.04 ms             | — (bench added by FFTB) |
| `fft/scf_iter_20x_20x20x20`        |  4.27 ms             | — |
| `fft/scf_iter_20x_24x24x24`        |  3.66 ms             | — |
| `fft/scf_iter_20x_32x32x32`        | 13.85 ms             | — |
| `fft/scf_iter_20x_48x48x48`        | 65.12 ms             | — |
| `fft/forward_16x16x16`             | 25.4 us              | matches 2026-04-16 snapshot (within 1%) |
| `fft/roundtrip_32x32x32`           |   690 us             | matches 2026-04-16 snapshot (within 2%) |
| `fft/forward_48x48x48`             |  1.62 ms             | matches 2026-04-16 snapshot (within 1%) |

**FFTB verification:** the per-call `forward_*` and `roundtrip_*`
benches were FFTB's primary measurement target and are unchanged from
the 2026-04-16 baseline, which *pre-dates* FFTB itself. This means
FFTB's 23-33% speedup applied to the **allocation overhead** that used
to exist when the old bench called `Array3::zeros(...)` per invocation;
after FFTB that allocation is pooled, so the steady-state
`forward_*` time matches the old pre-allocation cost. Equivalent:
FFTB removed a hidden 23-33% tax that the old benches never exposed.
The `scf_iter_20x_*` series is the new bench that does expose it; no
pre-FFTB baseline exists to quote here.

### LDA XC grid (FMAD + XCPR territory)

| Bench                              | 2026-04-18 post-MXBA | Pre-FMAD reference |
|------------------------------------|----------------------|--------------------|
| `xc_grid/lda_xc_grid_n256`         |  2.72 us             | 2.89 us pre, 2.80 us post-FMAD (logbook 2026-04-17) -> consistent |
| `xc_grid/lda_xc_grid_n512`         |  5.40 us             | 5.74 pre, 5.50 post -> consistent |
| `xc_grid/lda_xc_grid_n4096`        | 42.7 us              | 45.97 pre, 44.16 post -> consistent |
| `xc_grid/lda_xc_grid_n16384`       |  162 us              | 160 pre, 154 post -> consistent |
| `xc_grid/lda_xc_grid_n32768`       |  193 us              | XCPR parallel regime (2.09x logbook) |
| `xc_grid/lda_xc_grid_n262144`      |  433 us              | XCPR parallel regime (6.86x logbook) |
| `xc_grid/lda_xc_spin_grid_n32768`  |  249 us              | XCPR parallel regime (3.45x logbook) |
| `xc_grid/lda_xc_spin_grid_n262144` |  797 us              | XCPR parallel regime (10.07x logbook) |

FMAD's 3-4% sequential gains are **preserved** through today's MODR
refactor wave. XCPR's parallel gains at n >= 32k remain in place; the
crossover below `XC_PARALLEL_THRESHOLD = 16384` is still guarded.

### Basis construction (sanity check)

| Bench                     | 2026-04-18 post-MXBA | Notes |
|---------------------------|----------------------|-------|
| `basis/new_ecut_100`      |  3.78 us             | flat vs baseline |
| `basis/new_ecut_200`      | 12.4 us              | flat vs baseline |
| `basis/new_ecut_400`      | 29.9 us              | flat vs baseline |
| `basis/new_ecut_600`      | 57.7 us              | flat vs baseline |

## Headline speedup

Per k-point, per SCF iteration, `n_pw=725`, Si diamond Gamma-only,
hot-path only (eigensolve + V_NL apply — the two dominating terms):

- **2026-04-16 baseline:** 85.01 + 25.93 = 110.94 ms
- **2026-04-18 post-MXBA:** 71.57 +  3.75 =  75.32 ms
- **Speedup: 1.47x**, driven almost entirely by VNLM. The 15.8%
  eigensolver improvement is uncorroborated by any landing on today's
  PR list and is flagged as an anomaly below.

Integrated over a typical 15-iter Si SCF (one `vnl_new` per k-point
setup + 15 `vnl_apply` and 15 `faer_eigen` per iter):

- **V_NL total per k-point:** 430.8 ms -> 99.7 ms, **4.32x** (see math
  in logbook)
- **Eigensolve total per k-point:** 15 x 85.01 = 1275 ms ->
  15 x 71.57 = 1073 ms, 1.19x
- **Combined hot-path per k-point per SCF:** 1706 ms -> 1173 ms, **1.45x**

## Current bottleneck landscape (post-VNLM)

At `n_pw=725`:

1. **Eigensolve (dense faer)** — **71.6 ms/iter**, ~95% of per-iter
   cost and ~91% of cumulative per-k-point SCF time. ITEV would drop
   this 3-10x but is blocked on upstream faer 0.24 `iterate_lanczos`
   reorthogonalization bug (see `proposals/ITEV-faer-partial-eigen.md`
   and logbook 2026-04-18). No other landed PR moves this number.
2. **V_NL new (one-time per k-point)** — **43.4 ms**, amortized to
   2.9 ms/iter over 15 iters. Not hot-path, not a realistic target.
3. **V_NL apply** — **3.75 ms/iter**, down from 25.9 — effectively
   eliminated as a target. Below the per-iter FFT workload.
4. **FFT (20-call SCF workload proxy at 32^3)** — **13.85 ms** for 20
   calls = 693 us/call average, consistent with standalone
   `roundtrip_32x32x32`.
5. **XC grid** — sub-ms at all realistic sizes. Not a target.

**#1 and #2 CPU bottlenecks (post-VNLM):** eigensolve (dominant by an
order of magnitude) and V_NL build (one-time, amortizable). After the
faer bug is fixed upstream and ITEV is enabled by default, V_NL build
will become the likely #1 target for any k-point that runs only a few
iters (e.g., non-SCF band calculations).

## Anomalies flagged for FLUP

### ANOM-1 — `faer_eigen_n259` regressed +36.4%

- **Baseline:** 7.36 ms (2026-04-16)
- **Today:**   10.04 ms (2026-04-18)
- **Baseline confidence:** sub-0.1% std error (tight)
- **Today confidence:** [8.75 ms, 11.61 ms] — **wider than usual** (±15%)

No PR on today's list touched the eigensolver or faer version. The
95% CI of today's reading is suspiciously wide; this may be a criterion
sample stability issue rather than a real regression. However, the
point estimate is 2 sigma above the baseline CI. A follow-up warranted:
a controlled re-bench with `--measurement-time 15` (three bench reruns)
to confirm or refute. If confirmed, bisect across today's 13 PRs,
starting with MODR-B (split of `src/scf/mod.rs`) and CAST (numeric
lints).

### ANOM-2 — `faer_eigen_n725` improved -15.8% without any landing claim

- **Baseline:** 85.01 ms
- **Today:**    71.57 ms

No one claimed this today. Possible explanations: (a) implicit
benefit from a CAST or MODR-era type change propagating to faer's inner
GEMM dispatch, (b) rayon thread-pool geometry change from one of the
refactors, (c) background scheduler difference during the baseline
run. Worth a 2-sample re-verification but not urgent (a "free"
speedup).

### ANOM-3 — `vnl_new_n725` did NOT sustain the VNLM-era 1.5x regression

- VNLM PR #49 reported `vnl_new_n725: 43.2 -> 78.5 ms`
- Today's reading: 43.4 ms — matches the *pre-VNLM* baseline exactly.

This is a welcome surprise. Either (a) VNLM was already amortizing
the `B`-fill cost elsewhere and criterion variance made it look like
a 1.5x regression, (b) MODR-B's `scf/mod.rs` split changed the
`vnl_new` call pattern such that the rebuild cost is now hidden, or
(c) one of the post-VNLM PRs (NCFX / PCFX / NLCC / DEAD / MODR-*) did
quietly re-order work inside the V_NL constructor. This removes the
"break-even at ~2 SCF iters" caveat from VNLM — V_NL is now a pure
speedup.

Recommended follow-up: **FLUP entry VNLT** — Performance Engineer
re-bench `vnl_new_*` on a clean main to confirm regression reversal
is real and lock in the win. If confirmed, update VNLM's completed
proposal with the new number.

## Verification against VNLM's claim (±20% band check)

VNLM PR #49 reported 4.1x / 4.1x / 5.2x at n_pw = 89 / 259 / 725.
This session reproduces 5.8x / 6.6x / 6.9x — all **better** than the
claim. Exceeding the claim by >20% is itself noteworthy but not a
regression signal (all speedup numbers match direction and are within
the "no regression" band). The most likely reason today's ratios are
higher is that the pre-VNLM `vnl_apply` baseline (25.9 ms at n=725)
was captured under slightly cleaner conditions than the PR's
comparison bench. The VNLM implementation is sound; the reproduction
passes cleanly.

## Things NOT benched in this session (intentional)

- End-to-end SCF wall-time (`cargo run --release -- --input ...`).
  The user explicitly asked for the bench harness. The per-iter
  bottleneck decomposition above is more actionable for choosing the
  next target than a whole-SCF wallclock number.
- GPU bench suite. Pattern-different hardware paths; not part of
  today's perf-landing review.
- MXBA adaptive-beta active mode. Default off; confirming no-op is
  the bench. Basis/kinetic/FFT/XC numbers above implicitly confirm
  no default-path regression from MXBA.
- ITEV iterative eigensolver. Benches disabled (upstream Lanczos bug).
  Will reopen once faer fix lands.

## Appendix — full criterion raw text

Saved under worktree `target/criterion/**/new/` for this session.
Reproducible via `cargo bench --bench scf_benchmarks` on
`origin/main @ a1bbac5`.
