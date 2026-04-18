# Performance Engineer Logbook

Entries: date, measurements (actual numbers), bottleneck findings, proposals assessed. Always include hardware context.

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

Three clean runs on current main (Apple M2, lock held ~134 s): `hamiltonian/vnl_new_n725` = 44.37 / 44.36 / 44.23 ms (CI < 0.5%). PR #49's 78.5 ms was a single-run criterion outlier. Only two commits touched `src/potential/nonlocal.rs` since VNLM — CAST added `#[allow]` + asserts (zero runtime cost), RDOC a docstring edit; neither can explain a ~35 ms swing.

VNLM now credited as 4.5× per-k-point over 15 iters with zero break-even. Struck VNLB from FLUP (no regression to recover). EIGV/EIGW anomalies still open — same re-bench protocol would clear them.

**Warning for next session:** criterion single-run outliers at n_pw=725 are ~±40% for vnl_new; always do 3 clean runs before believing a regression.

## 2026-04-18 — SYMP landed (PR #64)

**Gotcha preserved:** naive `par_iter_mut` on `symmetrize_density_g`'s per-G loop regressed 72³·8 by 1.5× (26 → 39 ms) — per-item scheduling overhead eats wins when inner work is a few sin_cos + complex MAC. Fix: `par_chunks_mut(ny·nz)` over xy-slabs. Coarsens scheduling; output-write locality better.

Best speedups (Apple M2, lock held): 36³·48 ops **5.24×**, 72³·48 ops **6.00×**. All PCFX tests bit-identical (per-slot reduction stays serial).

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

Baseline re-measured (Apple M2, criterion, lock held):

| n_pw | `faer_eigen` | `vnl_apply` | `vnl_new` |
|------|--------------|-------------|-----------|
| 89   | 1.07 ms      | 0.53 ms     | 7.17 ms   |
| 259  | 56.3 ms      | 4.70 ms     | 22.1 ms   |
| 725  | 835.8 ms     | 28.6 ms     | 47.2 ms   |

**Prior logbook's n725 = 73.5 ms was criterion noise.** O(n³) scaling confirmed — the correct number is ~836 ms.

End-to-end SCF sample profile: eigensolver 60.7% of classified user CPU; after re-attributing rayon pool overhead → **~85-90% of SCF user CPU**. ITEV proposal obsoletes DVSN (faer 0.24 ships `partial_self_adjoint_eigen`).

**Tangential:** Rayon pool overhead 27.6% of samples is high — likely faer's inner spindle + outer k-point par_iter double-nesting. ITEV should drop this because each k-point call becomes much shorter.

## 2026-04-17 — XCPR Step 1 (PR #19) + FMAD

**XCPR threshold calibration:** rayon fork/join overhead ~70 µs/region on Apple M2; sequential `lda_xc_grid` ~11 ns/point, spin variant ~26 ns/point. n=4096 naive parallel is 161% regression. `XC_PARALLEL_THRESHOLD = 16384`. At 64³, spin XC drops 10× — prior "marginal" label was wrong.

**FMAD:** 24 `mul_add` substitutions; measured -3 to -4% on `lda_xc_grid_*`. Spin kernel within noise (dominated by cbrt, not additive).

**Warning:** during first XCPR bench pass, background `spin_polarization` test at 1400% CPU corrupted criterion timings. Always verify load avg < 5 before trusting criterion numbers — compilation happens outside the machine lock.

**Tangential:** rayon pool reuse (pin workers, skip fork/join) could lower breakeven to ~1000 pts — cross-cutting cross-proposal potential.

## 2026-04-16 — Orientation + baseline

Apple M2, end-to-end SCF. `si_scf.yaml` (ecut=100, 2×2×2): 0.11 s wall, 12 MB peak. `si_scf_converged.yaml` (ecut=200, 4×4×4): 0.34 s wall, 75 MB peak. Rayon ~6× across 10 k-points.

FFTB well-motivated (fresh Array3 per call); XCPR marginal on small grids (later revised — see 2026-04-17 entry); WFRX 25× claim optimistic (likely 5-15% net). Gap: no allocator proposal (jemalloc/mimalloc) — macOS libmalloc underperforms under multithreaded pressure.
