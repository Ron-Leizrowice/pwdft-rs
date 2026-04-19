# Performance Engineer Logbook

Entries: date, measurements (actual numbers), bottleneck findings, proposals assessed. Always include hardware context.

## 2026-04-19 — PROF landed: samply is the canonical profiler

`samply` is now named as the single canonical profiler in CLAUDE.md § Observability, profiling, benchmarking. Install via `cargo install samply`; always hold the machine lock while recording. `cargo flamegraph`, `tracing-flame`, and hand-rolled `Instant::now()` timers are off-menu for new work. Instruments.app stays as a fallback only for Metal GPU timeline questions. See CLAUDE.md § Profiling recipe (samply) for the invocation; § When to reconsider `tracing` captures the triggers that would reopen the decision.

## 2026-04-19 — GOPT PR-B landed (BufferPool extended to XC + V_eff); PR #106

F-4 + F-11 (§4 PR-B). Pool now covers all 3 kernels: 5 complex + 3 real-scalar + 2 staging + 2 uniform + 3 cached bind groups. Steady-state SCF iter now issues zero `create_buffer` and zero `create_bind_group` — only `queue.write_buffer` uploads. Fresh-alloc fallback preserved.

Apple M3 Max, lock held, `cargo bench --features gpu --bench gpu_benchmarks -- --quick`:

| Kernel  | 32³ fresh→pool | 64³ fresh→pool | 128³ fresh→pool |
|---------|:--------------:|:--------------:|:---------------:|
| hartree | 1.44→1.37 ms   | 2.22→1.99 ms   | 10.05→6.54 ms   |
| v_eff   | 1.55→1.50 ms   | 3.13→2.96 ms   | 17.87→14.35 ms  |
| lda_xc  | 1.38→1.35 ms   | 1.63→1.52 ms   | 5.24→3.18 ms    |

64³ (production scale): **5–11%** across all 3 kernels — matches proposal §2 5–10% estimate. 128³: 24–65% (per-call `create_buffer` cost grows with grid). Alloc count before→after: **11 + 3 → 0 + 0** per SCF iter.

**Correctness:** 7/7 `tests/gpu_consistency.rs` pass at existing tolerances. New `test_gpu_buffer_pool_xc_and_v_eff_match_fresh` asserts bit-identical pool↔fresh for XC + V_eff.

**Note:** GPU kernel per-call time still dominated by wgpu submit + poll overhead (~1.3 ms floor). Pool removes the allocation tax but GPU_MIN_GRID (F-3, PR-A) still has to flip XC back to CPU at n < 64³ — CPU is 8× faster there (176 µs XC CPU vs 1.35 ms GPU at 32³). Makes F-1 chain fusion (PR-C) the next unlock since it amortizes one submit across all 3 kernels.

**Next:** PR-A (F-3 + F-10) still open, then PR-C (F-1 chain fusion) on top of PR-B. PR-D (F-6 vec2 + F-8 cbrt NR) shader-only, independent.

**Gotcha:** machine-lock had races today — `agent-a95465b5` force-acquired mid-test of agent-af590306 because my test PID looked "stale" briefly. Tests completed anyway (hook only blocks new commands). MLFX claims PID-liveness check fixed this — it didn't fully. Worth a second look.

## 2026-04-19 — ALOC F-5 landed (Hamiltonian Mat cache)

Per-k `faer::Mat<Complex64>` scratch now lives in `ScfContext::h_scratch` (length `nspin * n_k`), fully overwritten by new `fill_hamiltonian_with_v_eff`. Zero per-iter `Mat::zeros(n_pw, n_pw)`.

Apple M3 Max, per-k assembly bench (inline old vs new, VNL included, lock held):

| n_pw | alloc_and_fill (old) | fill_into_cached (new) | Δ/call         |
|------|---------------------:|-----------------------:|:---------------|
|  89  |    88.4 µs           |   103 µs (noisy)       | within noise   |
| 259  |   723.6 µs           |   721.1 µs             | −0.3%          |
| 725  |  5.585 ms            |  5.494 ms              | **−91 µs (−1.6%)** |

Saving per iter scales `n_k · 91 µs` at n_pw=725 → ~0.9 ms/iter @ n_k=10. **Allocator traffic saved: ~84 MB × n_iter transient → zero** (the real structural win; 1.68 GB → 0 on a 20-iter Si 4×4×4 ecut=400 run).

**End-to-end SCF wall-times, post-ALOC-F5:** Si Γ ecut=100 523 ms, Γ ecut=200 1.18 s, 2×2×2 ecut=200 1.14 s, 4×4×4 ecut=200 1.58 s.

**Gotcha:** `VNL::add_to_hamiltonian` uses `matmul(Accum::Add, ...)` — accumulates into H. The assembly contract is therefore "fill fully first, then accumulate." Zero-fill skippable only because we rewrote kinetic+V_eff loop to write each `(i,j)` with `=` (diagonal) or `=` (off-diag) rather than `+=`.

**Bit-identity verified:** `aloc_f5_si_scf_is_deterministic` asserts `to_bits()` equality across two back-to-back SCFs on identical inputs. Eigenvalues match exactly at every (k, band). Si SCF pin (50 meV tol): E_total = −229.0566 eV.

**Flagged for Core Engineer:** `tests/vgc5_per_component_si.rs::test_madoc_band_sum_identity_si` fails under `--features gpu` on origin/main (216e050 and current). Residual 5.58e-6 eV vs tol 1e-7 eV — CPU version passes. Pre-existing MADOC vs GPU f32 precision issue, not an ALOC-F5 regression.

**Next-highest perf wins (from ALOC §2.5):** F-7 (psi_g in band loop, 200–1000 µs/iter), F-12 (FFT3D twiddle rebuild in fold init, 100–500 µs/iter), F-2 (XC grid in-place, 50–200 µs/iter).

## 2026-04-19 — WFRX Phase 1 landed (PR #99)

Subspace warm-start on dense eigensolver. **Default OFF** per spec — small-n regression rules out default-on. Numerical equivalence: Si |ΔE| = 1.39e-10 eV (proposal gate = 1e-8; 2 orders tighter).

Apple M3 Max, isolated-kernel bench (lock held):

| n_pw | full_dense | subspace_warm | warm / full |
|------|-----------:|--------------:|------------:|
|   89 |   0.99 ms  |     1.44 ms   |    1.45×    |
|  259 |   7.65 ms  |     7.92 ms   |    1.03×    |
|  725 |  79.6  ms  |    74.2  ms   |    0.93× (7% win) |

**Key takeaway:** faer 0.24's full Hermitian `self_adjoint_eigen` is already highly optimized. The Rayleigh-Ritz project + residual + rotate overhead (~5-10 ms at n=725) dominates the savings for all but production-scale runs. **Prior 25× claim is obsolete** for this faer configuration — the real number is ≤10% at n≥725, and negative at n<200.

**Gotchas captured:**
- Convergence test must use tight `conv_threshold` (1e-8) + good mixer (Broyden β=0.7) to hit 1e-8 eV agreement. At conv_threshold=1e-6 the SCF-level noise floor is ~1e-6 eV, which swallows the WFRX signal. Plain mixing at tight threshold won't converge in reasonable iterations.
- Machine-lock bench contamination persists — 3 other agents ran `spin_polarization` / `vgc5` at 500-1500% CPU while I held the lock. First bench pass was unusable (n=89 full_dense at 40 ms instead of 1 ms). Had to wait ~15 min for load avg to drop from 49 → 29 before numbers stabilized.

**Next-step lever for eigensolve reduction:**
- Per ALOC F-5 projection, the eigensolver is still ~95% of iter cost at n=725. Only **ITEV** (blocked on faer 0.24 Lanczos upstream bug) moves the needle meaningfully here. WFRX Phase 2 (ITEV warm-start) becomes free once ITEV lands — the `prev_wavefunctions` plumbing is already in place.

## 2026-04-19 — ALOC: per-iteration allocation audit (PR #93)

Static read-only audit of the SCF iteration body. 17 findings (F-1 to F-17) across `scf/driver.rs`, `scf/driver_spin.rs`, `scf/energy.rs`, `scf/density.rs`, `scf/potentials.rs`, `scf/mixing/*`, `symmetry/density/g_space.rs`, `potential/xc.rs`. Top wins (ranked µs/iter at Si 4×4×4, n_pw=725, 32³):

| # | Site | Est. µs/iter | Fix |
|---|---|---|---|
| F-5  | `scf/potentials.rs:147` — `faer::Mat::zeros(n_pw, n_pw)` per k | **200 – 2 000** | cache Mat per-k in ScfContext; `h.fill(0)` + kinetic-fill |
| F-7  | `scf/density.rs:65` — `psi_g` per band per k | **200 – 1 000** | hoist into rayon `fold` init closure |
| F-12 | `scf/density.rs:57` — `FFT3D::new` per worker per iter | **100 – 500**   | preallocate FFT3D pool sized to `rayon::current_num_threads()` |
| F-2  | `potential/xc.rs` + `scf/energy.rs` (XC + vxc_g) | **50 – 200** | `lda_xc_grid_into` + workspace |
| F-11 | mixer — 5–6 full-grid allocs/iter | **30 – 100** | `MixerWorkspace` + slice-based `precondition_residual` |

**Aggregate projected savings at production scale (non-spin, Si 4×4×4, n_pw=725, n_grid=32³):**
~640 – 4 060 µs/iter. Spin ≈ 2× on grid findings + F-13 (6× `Vec<f64>` in CCMX basis change). Context: current eigensolve is ~72 ms/iter at n_pw=725 → savings are **1.5 – 5% of iter wall time now, proportionally larger post-WFRX**.

**Gotcha for implementation:** F-5's Mat cache saves the allocator call, but still pays 8.4 MB of zero-writes per k per iter. Verify with a microbench that `fill(Complex64::zero()) + kinetic_fill` really beats `Mat::zeros() + kinetic_fill` — at 725² it might be a wash because the kinetic-fill already touches every slot. If the zero pass dominates, skip it and overwrite every entry unconditionally (kinetic already fills diagonal; V_eff off-diagonal loop fills the rest; VNLM GEMM is an additive pass).

**Out-of-scope flags for next session:**
- `NonlocalPotential::new` allocates heavily (audit shows `form_factor_by_atom` + `proj_l_by_atom` + 3 more per-atom Vecs are each cloned from the `type_ff` cache — O(n_atoms × n_proj × n_pw) total). All one-time at SCF entry → not ALOC's concern, but a separate optimization if Fe supercells come up.
- `compute_density`'s final re-normalization loop (`density.rs:88-97`) mutates in place; no allocations. Fine.
- `density_r_to_g` takes `&mut [Complex64]` and is already zero-alloc internally; caller-side staging was the allocator.

**Fix sequence** (§3 of proposal): Step 1 = skeleton `ScfWorkspace` struct (no savings, CI bit-identical); Step 2 = F-5 + F-7 (headline PR); Step 3 = F-1/2/4/8 together; Step 4 = F-10; Step 5 = F-11; Step 6 = F-12; Steps 7-8 = spin + cleanup.

## 2026-04-19 — VNLT: VNLM vnl_new regression was bench noise

Three clean runs on current main (Apple M3 Max, lock held ~134 s): `hamiltonian/vnl_new_n725` = 44.37 / 44.36 / 44.23 ms (CI < 0.5%). PR #49's 78.5 ms was a single-run criterion outlier. Only two commits touched `src/potential/nonlocal.rs` since VNLM — CAST added `#[allow]` + asserts (zero runtime cost), RDOC a docstring edit; neither can explain a ~35 ms swing.

VNLM now credited as 4.5× per-k-point over 15 iters with zero break-even. Struck VNLB from FLUP (no regression to recover). EIGV/EIGW anomalies still open — same re-bench protocol would clear them.

**Warning for next session:** criterion single-run outliers at n_pw=725 are ~±40% for vnl_new; always do 3 clean runs before believing a regression.

## 2026-04-18 — SYMP landed (PR #64)

**Gotcha preserved:** naive `par_iter_mut` on `symmetrize_density_g`'s per-G loop regressed 72³·8 by 1.5× (26 → 39 ms) — per-item scheduling overhead eats wins when inner work is a few sin_cos + complex MAC. Fix: `par_chunks_mut(ny·nz)` over xy-slabs. Coarsens scheduling; output-write locality better.

Best speedups (Apple M3 Max, lock held): 36³·48 ops **5.24×**, 72³·48 ops **6.00×**. All PCFX tests bit-identical (per-slot reduction stays serial).

**Tangential:** xy-slab chunking pattern may fit `scf/density.rs` accumulation + V_eff assembly. Audit anywhere `par_iter_mut` was tried naively.

## 2026-04-18 — Post-MXBA headline benchmark pass

| Bench              | 2026-04-16 | 2026-04-18 | Delta |
|--------------------|------------|------------|-------|
| vnl_apply_n725     | 25.9 ms    | 3.75 ms    | 6.9×  |
| vnl_new_n725       | 42.3 ms    | 43.4 ms    | +2.5% |
| faer_eigen_n725    | 85.0 ms    | 71.6 ms    | -16%* |
| faer_eigen_n259    | 7.36 ms    | 10.04 ms   | +36%* |

VNLM reproduced above PR #49's claim (6.9× actual vs 5.2× reported). *=unexplained anomalies, logged to FLUP as ANOM-1/2.

**Current bottleneck ranking at n_pw=725 (post-VNLM):**
1. Eigensolve (dense faer) — 71.6 ms/iter, ~95% of per-iter cost. Only ITEV moves this; blocked on faer 0.24 Lanczos upstream bug.
2. V_NL new — 43.4 ms one-time per k-point, amortizable.
3. FFT — ~14 ms for 20 calls at 32³.

Per-iter ratio eig:vnl_apply jumped ~3:1 → ~19:1. Anything that moves eigensolver off top (ITEV, WFRX) is now overwhelmingly highest-leverage.

## 2026-04-18 — ITEV Phase 1+2 landed (PR #45)

Wrapper + SCF dispatch + YAML switch + 8 correctness tests. Default still Dense. Iterative opt-in.

**BLOCKER: upstream faer 0.24 bug.** `operator/self_adjoint_eigen/iterate_lanczos` (lines 42-59) inner Gram-Schmidt loop spins indefinitely when a Krylov vector becomes numerically null. Reproduces non-deterministically on multi-iteration Si SCF at n_pw=89 — single calls are fine, SCF hangs after ~3-10 iters (stuck in `norm_l2_simd_pairwise_rows`). Cannot flip default to Iterative until upstream fixed.

**Gotchas captured for next session:**
- Degeneracy collapse needs `n_request = max(n_bands+4, 1.5*n_bands)`.
- Shift strategy: tight per-row Gershgorin bound. Loose shifts cause Lanczos reorth hangs even on single calls.

**Alternative paths:** hand-rolled LOBPCG or SPRS-based Davidson avoid Lanczos entirely — revisit if upstream fix doesn't land.

## 2026-04-17 — Post-FFTB/FMAD profiling pass; ITEV proposal opened

Baseline re-measured (Apple M3 Max, criterion, lock held):

| n_pw | `faer_eigen` | `vnl_apply` | `vnl_new` |
|------|--------------|-------------|-----------|
| 89   | 1.07 ms      | 0.53 ms     | 7.17 ms   |
| 259  | 56.3 ms      | 4.70 ms     | 22.1 ms   |
| 725  | 835.8 ms     | 28.6 ms     | 47.2 ms   |

**Prior logbook's n725 = 73.5 ms was criterion noise.** O(n³) scaling confirmed — the correct number is ~836 ms.

End-to-end SCF sample profile: eigensolver 60.7% of classified user CPU; after re-attributing rayon pool overhead → **~85-90% of SCF user CPU**. ITEV proposal obsoletes DVSN (faer 0.24 ships `partial_self_adjoint_eigen`).

**Tangential:** Rayon pool overhead 27.6% of samples is high — likely faer's inner spindle + outer k-point par_iter double-nesting. ITEV should drop this because each k-point call becomes much shorter.

## 2026-04-17 — XCPR Step 1 (PR #19) + FMAD

**XCPR threshold calibration:** rayon fork/join overhead ~70 µs/region on Apple M3 Max; sequential `lda_xc_grid` ~11 ns/point, spin variant ~26 ns/point. n=4096 naive parallel is 161% regression. `XC_PARALLEL_THRESHOLD = 16384`. At 64³, spin XC drops 10× — prior "marginal" label was wrong.

**FMAD:** 24 `mul_add` substitutions; measured -3 to -4% on `lda_xc_grid_*`. Spin kernel within noise (dominated by cbrt, not additive).

**Warning:** during first XCPR bench pass, background `spin_polarization` test at 1400% CPU corrupted criterion timings. Always verify load avg < 5 before trusting criterion numbers — compilation happens outside the machine lock.

**Tangential:** rayon pool reuse (pin workers, skip fork/join) could lower breakeven to ~1000 pts — cross-cutting cross-proposal potential.

## 2026-04-16 — Orientation + baseline

Apple M3 Max, end-to-end SCF. `si_scf.yaml` (ecut=100, 2×2×2): 0.11 s wall, 12 MB peak. `si_scf_converged.yaml` (ecut=200, 4×4×4): 0.34 s wall, 75 MB peak. Rayon ~6× across 10 k-points.

FFTB well-motivated (fresh Array3 per call); XCPR marginal on small grids (later revised — see 2026-04-17 entry); WFRX 25× claim optimistic (likely 5-15% net). Gap: no allocator proposal (jemalloc/mimalloc) — macOS libmalloc underperforms under multithreaded pressure.
