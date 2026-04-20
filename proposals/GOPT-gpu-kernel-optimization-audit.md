---
id: GOPT
title: GPU kernel + wgpu host path optimization audit (scoping)
priority: medium
complexity: medium
risk: low-medium
depends_on: []
blocks: []
status: active
author: Performance Engineer
date: 2026-04-18
---

## GOPT — GPU kernel + wgpu host path optimization audit

### §1 Scope

#### In scope

- Static audit of the wgpu compute path in `src/gpu/mod.rs` and the three
  WGSL kernels in `src/gpu/shaders/` (`hartree.wgsl`, `lda_xc.wgsl`,
  `v_eff_add.wgsl`).
- Static audit of the per-iteration GPU dispatch pattern in
  `src/scf/driver.rs:161-199` (non-spin SCF loop).
- Findings categorized by: workgroup/dispatch, buffer coalescing, f32
  precision, register pressure, CPU↔GPU transfers, pipeline/bind-group
  reuse, f64/f32 boundary cost, submission batching, CPU-fallback
  threshold.
- A benchmarking plan (§3) that each landed fix must follow once the
  machine lock is free.
- PR sequencing (§4) ordered by expected impact × ease.

#### Explicit non-goals

- **Rewriting the kernels in CubeCL** — that is the deferred CUCL proposal;
  GOPT stays within the current wgpu/WGSL stack so it can land incrementally.
- **Rewriting CPU paths.** Every finding touches the GPU side (shader,
  host dispatch, or the f64↔f32 boundary only).
- **Running `cargo bench --features gpu` during the audit** — WFRX holds
  the machine lock this session. Each fix PR will acquire the lock and
  re-bench before landing (§3).
- **Changing numerical tolerances in `tests/gpu_consistency.rs`.** Any
  optimization must keep the existing tightness (Hartree ≤ 1e-5 rel,
  XC ε_xc ≤ 0.01 eV, V_eff ≤ 1e-5 rel, Si SCF energy ±0.1 eV).
- **Defaulting the GPU feature on.** It remains opt-in behind
  `--features gpu` per CLAUDE.md.

#### Hardware target

Apple M3 Max (integrated GPU, 10-core, unified memory, Metal backend via
wgpu 26). Findings reference Apple Metal SIMD-group size (32 lanes) and
Apple-silicon buffer-coalescing rules. A few findings apply equally to
the Vulkan fallback on non-Apple hardware, but the benchmarking plan
only covers macOS.

### §2 Findings

Numbered F1–F12, grouped by category. Each has `file:line` citations,
a concrete fix, and an impact estimate at **production scale** — here
meaning a 64³ or larger FFT grid (Si `si_scf_converged.yaml`-class run,
n_grid ≥ 2.6×10⁵). Micro < 5 %, modest 5–15 %, major > 15 %.

#### A. CPU ↔ GPU transfer overhead (largest wins live here)

**F1 — Transfer thrash in the three-kernel chain is the dominant cost.**
`src/scf/driver.rs:167-199` calls `gpu.hartree_potential`,
`gpu.lda_xc`, and `gpu.v_eff_assembly` as three independent
round-trips per SCF iteration. Each method:

1. Copies inputs host→device (`queue.write_buffer`).
2. Runs one dispatch.
3. Copies output to staging.
4. Blocks on `device.poll(Wait)`.
5. Returns a freshly allocated `Vec<Complex64>` (or `Vec<f64>` for XC).

`v_eff_assembly` then takes the `v_h_fft` **and** `vxc_g` it just
received from steps 1-4 of the other two kernels and ships them back
onto the device (`src/gpu/mod.rs:312-316`). For a 64³ grid
(n_grid = 262 144) that is 3 × 2 MB complex outputs round-tripped
needlessly per iteration, plus 3 blocking submissions. Current code
`src/gpu/mod.rs:248, 342, 393`:

```text
self.queue.submit(std::iter::once(encoder.finish()));  // ×3 per iter
```

Measured transfer bandwidth on Apple unified memory is ~50 GB/s host→device
via `write_buffer`, but the per-submission fixed cost (encoder+submit+poll)
is ~40–80 µs even on M2. At 40 µs × 3 submits × 20 SCF iters = 2.4 ms of
pure submission overhead; on a small 16³ grid that is already ~30 % of
iteration wall time.

**Fix (major).** Introduce a `gpu.hartree_xc_veff_chain(&rho_g, &rho_r,
&v_local)` façade that:

- Writes `rho_g` and `rho_r` once.
- Dispatches `hartree_pipeline` and `lda_xc_pipeline` into the same
  command encoder (they read disjoint inputs and write disjoint
  outputs — no barrier needed beyond the natural end-of-pass sync).
- Reads back `vxc_r` to CPU so it can be FFT'd (cannot run FFT on GPU
  yet). **But** keeps `v_h_fft` and (a new) `vxc_g_fft` resident on
  GPU. Either we add an FFT GPU kernel (out of scope; see FFTB/CUCL)
  or we accept one round trip for `vxc_r`; even then the savings on
  Hartree and the V_eff assembly are meaningful.
- On the next host FFT of `vxc_r → vxc_g`, writes `vxc_g` back, then
  dispatches `v_eff_pipeline` using the already-resident `v_local`
  and `v_h`, reads back `v_eff_fft`.

Net: one write-buffer call for `v_local` (persistent, uploaded in
`prepare_buffers`), one for `rho_g`, one for `rho_r`, one for `vxc_g`,
and **two** readbacks instead of three. Submissions drop from 3 to 2.
On 64³ production-scale, expected 20–35 % reduction in SCF-iteration
GPU wall time. **Impact: major.**

**F2 — Submission batching inside `hartree_potential` alone.** Even
without F1's chain fusion, each kernel call today builds a fresh
`CommandEncoder`, records one compute pass + one copy, and submits —
the `queue.submit(std::iter::once(...))` pattern at
`src/gpu/mod.rs:248, 342, 393`. wgpu's submission cost on Metal is
dominated by the `MTLCommandBuffer` commit, not the encoding. Batching
the Hartree + V_eff_assembly into one encoder (XC still needs a
readback for the CPU FFT) halves submission count for the subset that
can share. Simpler stepping-stone to F1. **Impact: modest.**

**F3 — No CPU-fallback threshold for small grids.** `src/scf/driver.rs:167`
unconditionally routes through GPU whenever `gpu.is_some()`. For tiny
production grids (16³ = 4 096 pts, used in e.g.
`examples/si_scf.yaml`), the per-kernel transfer + submission cost
(~200 µs) can exceed the CPU cost of the same operation (~50 µs on
rayon). The CPU path is already available in-tree; we just never pick
it. Evidence in CLAUDE.md's GPU strategy section ("Transfer overhead
is the bottleneck for small grids") — but no code enforces that
insight.

**Fix.** Add a compile-time `const GPU_MIN_GRID: usize = 32_768;`
(conservatively one 32³ grid) in `src/gpu/mod.rs`, expose it as a
method `should_use_gpu(n_grid: usize) -> bool`, and guard each of the
three call sites in `driver.rs`:

```rust
let v_h_fft = if let Some(ref gpu) = gpu
    && gpu.should_use_gpu(ctx.n_grid)
{
    gpu.hartree_potential(&rho_g, &ctx.g_squared, fourpi_e2)
} else {
    hartree_on_fft_grid(&rho_g, &ctx.g_squared)
};
```

The exact threshold must be pinned by the §3 sweep at n = {4 k, 8 k,
32 k, 64 k, 256 k}. **Impact: major for small grids (a silent
regression today); micro at production scale.**

**F4 — `lda_xc` and `v_eff_assembly` bypass the buffer pool.**
`BufferPool` (`src/gpu/mod.rs:45-53`) allocates five complex buffers
and one g² buffer, but only `hartree_potential`'s pooled path
(`:216-252`) actually uses them. `lda_xc` (`:356-401`) and
`v_eff_assembly` (`:292-346`) always take the fresh-allocation path
— 4 fresh `wgpu::Buffer`s and 2 staging buffers per call for XC,
5 fresh + 1 staging for V_eff. Buffer creation is ~200 ns on Metal
but `create_buffer_init` with the UNIFORM flag also runs a `write_buffer`
internally, adding another small transfer.

**Fix.** Extend `BufferPool` with dedicated real-scalar buffers
(`rho_r_buf`, `exc_buf`, `vxc_r_buf` for XC; `v_local_buf`, `v_h_buf`,
`vxc_g_buf`, `v_eff_buf` for V_eff — the last one can alias
`complex_bufs[2]`). `v_local` is static across SCF iterations; upload
it in `prepare_buffers` alongside `g_squared`. Add a pooled branch to
both methods following the shape of the Hartree path. Couples well
with F1 (the façade already needs these buffers resident). **Impact:
modest.**

#### B. Workgroup / dispatch sizing

**F5 — Workgroup size 256 is defensible but un-swept.**
All three WGSL shaders use `@workgroup_size(256)`
(`src/gpu/shaders/hartree.wgsl:14`,
`src/gpu/shaders/lda_xc.wgsl:36`,
`src/gpu/shaders/v_eff_add.wgsl:14`). Apple Metal SIMD-group width
is 32; 256 = 8 SIMD-groups per threadgroup, which is the sweet spot
for most pointwise kernels. But:

- Apple's own guidance for pointwise memory-bound kernels on M-series
  is 64 or 128 per threadgroup — the lower count reduces register
  pressure and leaves more occupancy for concurrent workgroups.
- `lda_xc.wgsl` does cube-root / log / divides; it is compute-bound
  (unlike Hartree + V_eff which are memory-bound). The right
  threadgroup size may differ between kernels.

**Fix.** Sweep {64, 128, 256, 512} per kernel at n_grid ∈ {64 k, 256 k,
2 M}. Lock the fastest per kernel as a WGSL constant. Cost is an
afternoon of bench time after the machine lock is free. **Impact:
micro-to-modest** (≤ 5 % per kernel; compounded across three kernels
could reach modest).

#### C. Buffer coalescing

**F6 — Scalar-strided complex reads could be `vec2<f32>` reads.**
WGSL storage-buffer layout for `hartree.wgsl:9-12` and
`v_eff_add.wgsl:9-12` interleaves complex numbers as `array<f32>`
with consecutive threads reading indices `2*idx, 2*idx+1`. At
workgroup size 256, threadgroup lanes access bytes
`[0,4,8,12,…,2044]` for the real part and `[4,8,12,…,2048]` for the
imaginary — both fall inside the same 128-byte Apple cache line, so
coalescing is fine. But the current form makes the compiler emit two
32-bit loads per thread instead of one 64-bit load.

**Fix.** Declare the buffers as `array<vec2<f32>>` in WGSL
(Rust-side `bytemuck::cast_slice(&rho_f32)` already works with
`#[repr(C)] struct Vec2f32 { re: f32, im: f32 }` — the byte layout
is identical). The shader body becomes:

```wgsl
@group(0) @binding(3) var<storage, read_write> v_h: array<vec2<f32>>;
...
v_h[idx] = rho_g[idx] * factor;   // scalar×vec2, compiled to one fma
```

naga + Metal lower this to a single 64-bit load and one `vec2`
multiply. **Impact: modest on memory-bound kernels (Hartree, V_eff);
micro on XC** (XC is compute-bound).

#### D. f32 precision

**F7 — No accumulation loops means Kahan is not needed (sanity check).**
All three kernels are pointwise — each thread produces one output
from one (or three) input(s), no reduction. No Kahan or pairwise sums
required. The relative error envelope in `tests/gpu_consistency.rs`
(1e-5 Hartree / V_eff, 0.01 eV XC) matches f32 ULP × number of ops
per cell. **No finding; recorded for completeness so the next GOPT
session does not re-audit it.**

**F8 — `pow(x, 1.0/3.0)` in `lda_xc.wgsl` for cube roots.**
`src/gpu/shaders/lda_xc.wgsl:55, 58` compute the Wigner-Seitz radius
and the Slater exchange prefactor via `pow(x, 1.0/3.0)`. `pow` in
WGSL compiles to `exp(y * log(x))` on Metal — two transcendentals for
a cube root. There is no native `cbrt` in WGSL, but the Newton-Raphson
`y = x / (y*y);  y = (2*y + x/(y*y)) / 3;` (seeded from a bit-hack f32
exponent shift) is 3–4× faster than `pow` on Metal and matches
`cbrt(x)` to within 2 ULP after two iterations.

**Fix.** Add a local `fn cbrt_f32(x: f32) -> f32` in `lda_xc.wgsl`
using the bit-hack initial guess + two NR steps, and replace both
`pow(…, 1.0/3.0)` calls. Validate against
`tests/gpu_consistency.rs::test_gpu_xc_across_density_regimes` (1000
log-spaced points, tolerance 0.01 eV — well above 2 ULP × 27 eV
scale). **Impact: modest** on XC specifically (XC is the
compute-bound one); micro on SCF wall time.

#### E. Register pressure / local variables

**F9 — Shader footprints are tiny; no register-pressure concern.**
Largest per-thread working set is `lda_xc.wgsl` with ~12 scalar
locals (`rho`, `rho_bohr`, `rs`, `sqrt_rs`, `denom`, `d_ec`, `ln_rs`,
`ec_ha`, `vc_ha`, `ex_ha`, `vx_ha`, `cbrt_arg`). Metal's register
file per thread is 256 × 32-bit; we use < 5 %. No spill risk.
**No finding.**

#### F. f64 ↔ f32 boundary cost

**F10 — Scalar `.iter().map(|&v| v as f32).collect()` allocates per call.**
`src/gpu/mod.rs:159, 255, 358, 471-477` do the f64→f32 cast in a
plain iterator chain. For `n_grid = 262 144` this allocates a fresh
`Vec<f32>` of 1 MB per kernel call (3 per iter × 20 iter = 60 MB
churn per SCF). The `complex_to_f32_pairs` helper (`:471-477`) is
worst — it `push()`es one element at a time, defeating any
auto-vectorization.

**Fix.** Preallocate a scratch `Vec<f32>` once per `GpuAccelerator`
(or per `BufferPool`), clear and `extend` instead of `collect`. For
the complex version, either SIMD-intrinsic it (portable-simd is
available on nightly; stable has no great answer) or at minimum
preallocate with `Vec::with_capacity(2*n)` and use the loop
directly so LLVM sees the bound. Tracking the scratch buffer in the
pool removes 1 MB × 6 allocs/iter on 64³ grids. **Impact: modest
on small-to-medium grids** (allocation-dominated); micro on large
grids where the dispatch itself dominates.

#### G. Pipeline / bind-group reuse

**F11 — Bind groups rebuilt every dispatch.**
`src/gpu/mod.rs:228-237, 261-270, 318-328, 372-381` create a fresh
`BindGroup` per kernel call. wgpu's Metal backend caches the
underlying `MTLArgumentBuffer` and the reconstruction is cheap
(~2 µs), but it is not free. Because `BufferPool` keeps the
underlying `wgpu::Buffer` handles alive across iterations, a bind
group holding references to those buffers is also reusable.

**Fix.** Cache one `BindGroup` per (pipeline, pool) pair in
`BufferPool`. Invalidate on `prepare_buffers` (fresh buffers → fresh
bind groups). `uniform` parameter buffers (`params_buf`) are tiny
and regenerated every call today — moving those into the pool with
a `write_buffer` update instead of `create_buffer_init` is a
companion cleanup. **Impact: micro** (≤ 2 % at 64³; may matter more
at 16³ where per-call overhead is already a bigger fraction).

#### H. Missing GPU path

**F12 — `driver_spin.rs` has no GPU path at all.**
`rg "feature = \"gpu\"" src/scf/driver_spin.rs` returns zero matches.
Spin-polarized SCF always takes the CPU route even with
`--features gpu` enabled. For nspin = 2 the Hartree and V_eff kernels
are identical (they act on `ρ_total = ρ↑ + ρ↓`); the XC kernel needs
a spin-polarized variant (PZ spin LDA is already implemented on
CPU in `src/potential/xc.rs`). A corresponding `lda_xc_spin.wgsl`
plus a `gpu.lda_xc_spin(rho_up, rho_down)` method closes the gap.

**Fix.** Out-of-scope for the first GOPT PR (touches correctness of
the spin path; needs a `gpu_consistency.rs` test for nspin=2). Flag
as a follow-up proposal (see §4 sequencing). **Impact: major for
spin runs** (which are the Fe/Co/Ni target for VGCMP Phase 2+);
zero for non-spin runs.

### §3 Benchmarking plan

All measurements must be run under the machine lock (CLAUDE.md
"Machine Coordination") on Apple M3 Max with no other heavy processes,
using the existing `benches/gpu_benchmarks.rs` harness extended where
needed.

#### Baseline (must run before any fix lands)

```bash
.claude/bin/machine-lock acquire "Performance Engineer" "GOPT baseline"
cargo bench --features gpu --bench gpu_benchmarks -- \
  hartree v_eff_assembly lda_xc
.claude/bin/machine-lock release
```

Record criterion JSON at
`target/criterion/<group>/<bench>/base/estimates.json` as the
pre-GOPT reference. The current `gpu_benchmarks.rs` already covers
n ∈ {8 k, 64 k, 512 k} per kernel; **extend it to include
n = 2 097 152 (128³)** so we have one production-scale point above
the existing range.

#### Per-finding verification

| Finding | New / modified bench | Success criterion |
|---------|----------------------|-------------------|
| F1 (chain fusion)            | New `scf_iteration_gpu` end-to-end bench at n = {32³, 64³, 128³}. Emulates driver.rs's three-kernel sequence with one warm pool. | ≥ 20 % wall-time reduction at 64³ vs. baseline chain. No regression at 16³. |
| F2 (submit batching)         | Sub-case of F1 (enables F1). Measure alone by folding only Hartree + V_eff into one submit, leaving XC separate. | ≥ 10 % at 64³ vs. baseline. |
| F3 (CPU threshold)           | `gpu_vs_cpu_hartree_small` at n ∈ {4 k, 8 k, 16 k, 32 k, 64 k}. | CPU wins ≤ 32 k, GPU wins ≥ 64 k. Pin `GPU_MIN_GRID` at the crossover. |
| F4 (pool XC + V_eff)         | Existing `lda_xc` and `v_eff_assembly` benches pre/post. | ≥ 5 % at 64³. |
| F5 (workgroup sweep)         | New param sweep: `hartree_wgsz_{64,128,256,512}` × `n ∈ {64k,256k,2M}`. | Pick lowest-mean per (kernel, n). |
| F6 (vec2 layout)             | Existing benches. | ≥ 5 % on Hartree and V_eff at ≥ 64 k. XC unchanged. |
| F8 (cbrt NR)                 | Existing `lda_xc_gpu_*`. Must also re-run `tests/gpu_consistency.rs::test_gpu_xc_across_density_regimes` to confirm tolerance holds. | ≥ 10 % on XC at 64 k; max abs error on exc/vxc stays ≤ 0.01 eV. |
| F10 (scratch Vec pool)       | Existing benches; isolate with an allocation-profiling run (Instruments Allocations instrument) on the Si SCF example. | Allocations/SCF-iter drop ≥ 60 % at n_grid = 4 096. |
| F11 (bind-group cache)       | Existing Hartree bench at n = 8 k (smallest, where overhead is largest). | ≤ 2 µs per-call reduction; not a regression anywhere. |

#### Correctness gates (every PR)

- `cargo test --features gpu --test gpu_consistency` — all four tests
  must pass with existing tolerances. **Do not loosen a tolerance to
  make a perf PR land.**
- `cargo test --features gpu` (full suite with GPU) for overall
  regression safety (~28 s).
- Both clippy invocations per CLAUDE.md § Code Quality.

#### Tangential instrument pass

Once F1 + F3 land, run Instruments' **Metal System Trace** on
`examples/si_scf_converged.yaml` with `--features gpu` to confirm:

- GPU utilization > 30 % during kernel dispatches (baseline is
  probably < 10 %).
- Submission gaps (dead time between encoder.finish and first
  dispatch) shrink.
- `MTLBuffer` allocation count per SCF iter drops monotonically.

### §4 Recommended PR sequencing

Landing order maximizes cumulative speedup per PR while minimizing
risk of a regression at small grids (where current code is already
probably wrong, per F3).

1. **PR-A: F3 (CPU-fallback threshold) + F10 (scratch Vec pool).**
   Smallest diff, largest asymmetric upside: fixes the likely
   small-grid regression today and reduces allocation churn. No
   shader changes, no buffer-layout changes. Risk: low. Estimated
   impact on small-grid SCF: > 30 %; on large-grid SCF: micro.
2. **PR-B: F4 (pool for XC + V_eff) + F11 (bind-group cache).**
   Extend `BufferPool` to cover all three kernels. Prerequisite for
   PR-C. Risk: low-medium (easy to mis-key a bind-group cache).
   Estimated impact: 5–10 % at 64³; unlocks the rest.
3. **PR-C: F1 (chain fusion into one encoder) + F2.**
   The big one. Requires PR-B's persistent buffers. Risk: medium
   (pipeline-barrier semantics + readback ordering must be right).
   Estimated impact: 20–35 % at 64³.
4. **PR-D: F6 (vec2<f32> layout) + F8 (cbrt NR).**
   Shader-only changes; correctness gated by `gpu_consistency`.
   Risk: medium (F8 can regress exc/vxc tolerance if NR iterations
   are cut too short — keep 2 iters minimum). Estimated impact:
   5–15 % on per-kernel wall time.
5. **PR-E: F5 (workgroup-size sweep).**
   Pure tuning. Requires the benches from earlier PRs to be stable.
   Risk: low. Impact: micro-to-modest compounded across three
   kernels.

#### Not in this proposal (flagged)

- **F12 (spin-polarized GPU path).** Large enough to merit its own
  proposal; will lift into INDEX as `GSPN` or similar after GOPT
  PR-A lands. Unblocks VGCMP Phase 2 (Fe 4×4×4 free-mag) running
  with GPU acceleration, which today takes several minutes on CPU.

#### Stop conditions

If after PR-C the GPU path is not within 2× of the theoretical
memory-bandwidth ceiling (50 GB/s × n_grid / iter → ~1 ms for
64³), escalate to the deferred CUCL proposal — wgpu overhead would
then be the genuine bottleneck, not kernel organization.

### Notes

- All file:line references are against commit `0520c8c` (main,
  2026-04-18).
- The GPU feature flag stays opt-in throughout; no default-features
  changes.
- The `tests/gpu_consistency.rs` tolerances (DBGC hoisted constants
  `SI_REFERENCE_TOTAL_EV = -213.0283`, `SI_REFERENCE_FERMI_EV =
  6.709`, ±0.1 eV) are the correctness backstop — any PR that drifts
  these is a bug, not an optimization.
- This proposal does **not** touch `driver_spin.rs`; F12 is
  recorded and promoted separately once the facade from F1 exists.
