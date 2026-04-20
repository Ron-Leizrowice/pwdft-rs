---
id: PMTL
status: active
priority: medium
complexity: medium
risk: medium
depends_on: []
blocks: [CUCL]
---

# PMTL: Split GPU backend into its own workspace package `pwdft-metal`

## Problem

GPU acceleration currently lives **inside** `pwdft-core` behind a `gpu` feature flag. The split gives us the worst of both worlds: GPU code lives in the hot Rust crate (slowing every compile when enabled), while the split across `#[cfg(feature = "gpu")]` and `#[cfg(not(feature = "gpu"))]` branches makes the SCF driver hard to read.

**Concrete evidence.**

| Metric | Count | Source |
|---|---|---|
| GPU Rust LOC | 980 | `pwdft/pwdft-core/src/gpu/mod.rs` |
| WGSL shader LOC | 140 | `pwdft/pwdft-core/src/gpu/shaders/{hartree,lda_xc,v_eff_add}.wgsl` |
| `#[cfg(feature = "gpu")]` branches in the SCF driver | 11 | `pwdft/pwdft-core/src/scf/driver.rs:84..450` |
| GPU-specific deps (`wgpu`, `pollster`, `bytemuck`) | 3 | `pwdft/pwdft-core/Cargo.toml` — all `optional = true` |
| Per-feature clippy invocations required by CI | 2 | `.github/workflows/ci.yml` — clippy default + clippy `--features gpu` |
| GPU integration tests | 1 | `tests/gpu_consistency.rs` |
| GPU benches | 1 | `benches/gpu_benchmarks.rs` |

**Why the status quo hurts.**

1. **Build-time tax.** Every CPU-only build still compiles the `wgpu`/`pollster`/`bytemuck` dependency graph when `cargo build --features gpu` is invoked — and CI always runs both variants (`clippy` default + `clippy --features gpu`). wgpu pulls in ~40 transitive crates. Splitting lets CPU consumers skip that cost entirely; GPU consumers still pay it but in a crate that's explicitly about GPU.
2. **Driver-file branching.** `scf/driver.rs` has 11 `#[cfg(feature = "gpu")]` blocks, most paired with a `#[cfg(not(feature = "gpu"))]` fallback. Each is a micro-fork of the SCF iteration. Readers have to mentally apply the flag at every step.
3. **Orthogonality signal.** `pwdft-core` is the physics engine. A GPU backend is an optimization strategy — conceptually a *plugin*. Burying it inside the core crate signals the opposite, and RWHK C1 already flagged one GPU-related correctness gap (GPU-vs-GPU comparison) that would have been harder to sneak in with a dependency-inverted boundary.
4. **Blocks CubeCL.** The deferred CUCL proposal (CubeCL GPU kernels) requires a clean backend boundary to allow switching wgpu ↔ CubeCL ↔ future accelerators without pinning every call-site. Today every potential second backend would have to rewrite the same cfg-dance.
5. **Name mismatch.** The crate feature is `"gpu"` but the current code hardcodes the wgpu/Metal pipeline. A second backend would need either a `"gpu-cubecl"` feature in the same crate (worse fragmentation) or — cleanly — its own crate. The name `pwdft-metal` declares the backend explicitly: Apple-silicon Metal via wgpu's Metal dispatch.

## Research

### Proposed workspace layout

```text
Cargo.toml                               # [workspace] members += "pwdft/pwdft-metal"
pwdft/
├── pwdft-core/                          # physics engine; NO gpu feature, NO wgpu dep
├── pwdft-metal/                         # new — GPU backend plugin
│   ├── Cargo.toml
│   ├── src/
│   │   ├── lib.rs                       # `GpuAccelerator`, `BufferPool`, public API
│   │   ├── hartree.rs
│   │   ├── xc.rs
│   │   ├── v_eff.rs
│   │   └── shaders/
│   │       ├── hartree.wgsl
│   │       ├── lda_xc.wgsl
│   │       └── v_eff_add.wgsl
│   ├── tests/                           # was tests/gpu_consistency.rs
│   └── benches/                         # was benches/gpu_benchmarks.rs
├── pwdft-validation/
└── faer/
```

### Dependency inversion contract

`pwdft-core` defines a small **trait** covering the grid operations that GPU offloads (Hartree, LDA XC, V_eff assembly):

```rust
// pwdft-core/src/gpu_backend.rs  — new, no GPU code
pub trait GridAccelerator: Send + Sync {
    fn hartree(&self, rho_g: &[Complex64], g_shell: &[f64], omega: f64) -> Vec<Complex64>;
    fn lda_xc(&self, rho: &[f64]) -> (Vec<f64>, f64);
    fn v_eff(&self, v_local: &[f64], v_h: &[f64], v_xc: &[f64]) -> Vec<f64>;
    fn buffer_pool_stats(&self) -> BufferPoolStats;   // for observability
}
```

`pwdft-metal` implements the trait (`impl GridAccelerator for MetalAccelerator`). SCF driver takes `Option<&dyn GridAccelerator>` (or `Option<Arc<dyn GridAccelerator>>`); on `None` it runs the existing CPU path. No `cfg` branches.

The binary crate (`pwdft`) depends on both `pwdft-core` and (optionally) `pwdft-metal` and wires them at startup:

```rust
// pwdft/src/main.rs
#[cfg(feature = "metal")]
let accel = Some(Arc::new(pwdft_metal::MetalAccelerator::new()?) as Arc<dyn GridAccelerator>);
#[cfg(not(feature = "metal"))]
let accel: Option<Arc<dyn GridAccelerator>> = None;

scf::run_scf(&params, accel.as_deref())?;
```

The feature flag **moves** from `pwdft-core` to the binary target. Nothing in the physics crate cares about backends.

### Naming: why `pwdft-metal` and not `pwdft-gpu`

| Option | Pros | Cons |
|---|---|---|
| `pwdft-gpu` | generic; survives future backends | lies — current impl is wgpu/Metal only; invites `pwdft-gpu/src/cubecl/` growth |
| `pwdft-metal` | honest about what it is today | doesn't cover a CUDA/ROCm/CubeCL port |
| `pwdft-wgpu` | matches the dep | wgpu is the dispatcher, not the target |
| `pwdft-backend-metal` + `pwdft-backend-cubecl` | scales | bike-shedding; premature |

Pick **`pwdft-metal`**. When CUCL lands, add a sibling crate `pwdft-cubecl`; both implement the same `GridAccelerator` trait and the binary picks one at build time. No shared-crate feature fragmentation.

### What stays in `pwdft-core`

- The `GridAccelerator` trait definition.
- The CPU fallbacks in `scf/potentials.rs` / `potential/xc.rs` — unchanged.
- All integration tests that don't touch GPU.
- Nothing else.

### What moves to `pwdft-metal`

- `src/gpu/mod.rs` → `pwdft-metal/src/lib.rs` (split by responsibility: `hartree.rs`, `xc.rs`, `v_eff.rs`, `buffer_pool.rs`).
- `src/gpu/shaders/*.wgsl` → `pwdft-metal/src/shaders/`.
- `tests/gpu_consistency.rs` → `pwdft-metal/tests/`.
- `benches/gpu_benchmarks.rs` → `pwdft-metal/benches/`.
- The three optional deps (`wgpu`, `pollster`, `bytemuck`) become mandatory deps of `pwdft-metal` (no more `optional = true`).

## Implementation

### Phase A — Define the trait in `pwdft-core`

1. Create `pwdft-core/src/gpu_backend.rs` with the `GridAccelerator` trait + `BufferPoolStats` struct. No behaviour, just the contract.
2. Re-export from `lib.rs` under a `gpu_backend` module.
3. Thread `Option<&dyn GridAccelerator>` through `scf::run_scf` → `scf::driver::run_scf_unpolarized` → the inner iteration helpers. Delete the 11 `#[cfg(feature = "gpu")]` branches — replace with `match accel { Some(a) => a.hartree(...), None => cpu_hartree(...) }`.
4. Unit test: a no-op `GridAccelerator` that forwards to the CPU path; assert byte-identical output vs the no-accelerator branch.

### Phase B — Create `pwdft-metal`

1. `cargo new --lib pwdft/pwdft-metal`.
2. Add to workspace `members` in root `Cargo.toml`.
3. Move `src/gpu/**` verbatim into `pwdft-metal/src/` (keep structure, just relocate). Split `mod.rs` (980 LOC) into per-operation files.
4. Move shaders to `pwdft-metal/src/shaders/` and update `include_str!` paths.
5. Add `pwdft-metal` dep on `pwdft-core` for the trait.
6. Implement `GridAccelerator for MetalAccelerator` — a ~30-line `impl` that delegates to the existing methods (which already do the work; just type-thread them through the trait).
7. Move `tests/gpu_consistency.rs` → `pwdft-metal/tests/`.
8. Move `benches/gpu_benchmarks.rs` → `pwdft-metal/benches/`.

### Phase C — Rewire the binary

1. `pwdft-core/Cargo.toml`: delete the `[features] gpu = [...]` block, delete the three optional deps, delete `src/gpu/**`.
2. `pwdft/src/main.rs` (binary under `pwdft-core`): add a `metal` feature that depends on `pwdft-metal`.
3. `Cargo.toml` binary: `[target.'cfg(target_os = "macos")'.dependencies] pwdft-metal = { path = "../pwdft-metal", optional = true }`.
4. At startup, construct the accelerator if the feature is on; pass `None` otherwise.
5. Update `CLAUDE.md § Build & Run`: `cargo build --features gpu` → `cargo build --features metal`.

### Phase D — CI + tooling

1. `.github/workflows/ci.yml`: replace the two clippy invocations with a single `cargo clippy --workspace --all-targets`, then add a second step `cargo clippy -p pwdft-metal --all-targets` (mirrors the pattern we already use for the Python validation job).
2. Drop the `--features gpu` step — it's now just another workspace member, linted by the workspace clippy.
3. Test step: `cargo test --workspace` (Tier 1 default). GPU tests inside `pwdft-metal` are gated by a runtime "does an adapter exist?" check that skips in CI (as today).

### Phase E — Cleanup sweep

1. `CLAUDE.md § Architecture` GPU section: update paths, rename `gpu` feature → `metal`.
2. `.claude/agents/performance-engineer.md` and other agent defs: `gpu feature` → `metal feature` + path `pwdft-metal/src/...`.
3. Grep for `src/gpu/` and `feature = "gpu"` in proposals/ — leave completed proposals as historical; update active ones (CUCL, GOPT, CFGN) to reference the new layout.

### Ordering

| Phase | Depends on | Unlocks |
|---|---|---|
| A | — | B |
| B | A | C |
| C | B | working `cargo build --features metal` |
| D | C | CI parity |
| E | C | docs coherence |

Phases A+B+C land together (one atomic PR — the trait move is not useful without the crate split). D+E are separate cleanup PRs.

## Verification

- `cargo build -p pwdft-core` succeeds with *no* `wgpu`/`pollster`/`bytemuck` in the compiled artifact (verify via `cargo tree -p pwdft-core -e normal`).
- `cargo build --features metal` (on macOS) produces a binary with the wgpu backend linked.
- `cargo test --workspace` — all existing tests pass, `gpu_consistency` runs under `pwdft-metal`.
- Cold-build time for `cargo build -p pwdft-core` drops (rough estimate: −30 to −60 s on a fresh target dir, from removing ~40 wgpu-transitive crates).
- `rg '#\[cfg\(feature\s*=\s*"gpu"\)\]' pwdft/pwdft-core/src/` returns zero hits.
- SCF driver reads linearly top-to-bottom — no more paired cfg forks.
- The trait indirection cost is measured: on the existing `scf_benchmarks`, CPU-only path should be within noise (±1%) of pre-refactor baseline — dyn-dispatch on `Option<&dyn GridAccelerator>` evaluates one branch per iteration, not per grid-point.

**Risk register.**

- **Trait-object overhead.** `Option<&dyn GridAccelerator>` introduces one indirect call per grid-level op (3× per SCF iteration — Hartree, XC, V_eff). Each op is a single bulk call, not a per-point invocation, so the overhead is amortized over a full grid. Bench before merging to confirm <1% regression.
- **Metal-only boundaries.** wgpu on macOS targets Metal natively; on Linux it maps to Vulkan. `pwdft-metal` should compile and run on both — the name is about Apple-silicon priority, not an OS gate. Keep the crate portable.
- **Future CubeCL crate.** When CUCL un-defers, `pwdft-cubecl` slots in as a parallel sibling with no changes to `pwdft-core`. Having to pick one at build time is acceptable; both-at-once is a later story.
