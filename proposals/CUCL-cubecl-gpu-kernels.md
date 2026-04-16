---
id: CUCL
status: deferred
priority: low
complexity: large
risk: high
depends_on: []
blocks: []
---

# CUCL: CubeCL GPU Kernels

> **Note:** Line numbers reference the pre-ScfContext codebase. Verify locations before implementing.

## Motivation

The current GPU module (`src/gpu/mod.rs`) uses wgpu with hand-written WGSL shaders (`src/gpu/shaders/`). This works for the three current kernels (Hartree, LDA XC, V_eff assembly), but has limitations:

1. **WGSL is verbose and untyped:** No Rust type checking, no IDE support, string-embedded via `include_str!`. Bugs in shader code are only caught at GPU runtime.

2. **f32 only:** WGSL has limited f64 support. The current approach converts f64->f32 at the CPU-GPU boundary, which is fine for grid operations but limits what can move to GPU.

3. **Expanding GPU scope:** The highest-value GPU targets — non-local potential application (O(n_pw^2) pair iteration) and batched FFTs — require more complex kernels than the current element-wise operations. Writing these in WGSL is error-prone.

**CubeCL** is a Rust GPU compute framework that compiles Rust code to WGSL, Metal Shading Language, CUDA PTX, or ROCm. It would let you write GPU kernels in Rust with full type checking, then compile to the appropriate backend.

## Dependencies

```toml
cubecl = { version = ">=0.6", optional = true, features = ["wgpu"] }
```

Replace or supplement:
```toml
# Current GPU deps become optional/removable:
wgpu = { version = ">=24.0", optional = true }
pollster = { version = ">=0.4", optional = true }
bytemuck = { version = ">=1.14", features = ["derive"], optional = true }
```

CubeCL uses wgpu internally but provides a higher-level API.

## Scope of Changes

### Phase 1: Port existing kernels

**Replace `src/gpu/shaders/hartree.wgsl` with Rust:**

Current WGSL (Hartree potential):
```wgsl
@compute @workgroup_size(256)
fn hartree(@builtin(global_invocation_id) id: vec3<u32>) {
    let idx = id.x;
    if (idx >= arrayLength(&rho_g)) { return; }
    let g2 = g_squared[idx];
    if (g2 < 1e-20) { v_h[idx] = vec2<f32>(0.0, 0.0); return; }
    let rho = rho_g[idx];
    let scale = fourpi_e2 / g2;
    v_h[idx] = vec2<f32>(rho.x * scale, rho.y * scale);
}
```

CubeCL equivalent in Rust:
```rust
#[cube(launch)]
fn hartree_kernel(
    rho_g: &Tensor<f32>,
    g_squared: &Tensor<f32>,
    v_h: &mut Tensor<f32>,
    fourpi_e2: f32,
) {
    let idx = ABSOLUTE_POS;
    let g2 = g_squared[idx];
    if g2 < 1e-20 {
        v_h[idx * 2] = 0.0;
        v_h[idx * 2 + 1] = 0.0;
        return;
    }
    let scale = fourpi_e2 / g2;
    v_h[idx * 2] = rho_g[idx * 2] * scale;
    v_h[idx * 2 + 1] = rho_g[idx * 2 + 1] * scale;
}
```

Benefits: Rust type checking, no string embedding, IDE support, unit-testable on CPU.

### Phase 2: GPU-accelerated non-local potential

The non-local potential (`src/potential/nonlocal.rs`, lines 154-200) has an O(n_pw^2) pair loop that is the secondary computational bottleneck. This is a natural GPU target:

```rust
#[cube(launch)]
fn nonlocal_kernel(
    form_factors: &Tensor<f32>,     // [n_proj, n_pw]
    dij: &Tensor<f32>,               // [n_proj, n_proj]
    q_vecs: &Tensor<f32>,            // [n_pw, 3]
    q_norms: &Tensor<f32>,           // [n_pw]
    atom_phases: &Tensor<f32>,        // [n_atoms, n_pw, 2] (complex)
    h_out: &mut Tensor<f32>,          // [n_pw, n_pw, 2] (complex)
    // ...
) {
    // Each workgroup handles one (ig, jg) pair
    let ig = ABSOLUTE_POS_X;
    let jg = ABSOLUTE_POS_Y;
    // ... compute projector sum, structure factor, accumulate into h_out
}
```

This moves the O(n_pw^2) computation to GPU where it can exploit massive parallelism. For n_pw=200, that's 40,000 independent pair computations — ideal GPU workload.

### Phase 3: GPU batched FFT (exploratory)

CubeCL does not provide FFT primitives, but it could be used to write a GPU-side FFT for the density computation (n_bands * n_kpts inverse FFTs per SCF iteration). This is more complex and may be better served by a dedicated GPU FFT library.

### Restructured GPU module

```
src/gpu/
├── mod.rs              # GpuAccelerator: device init, buffer management
├── kernels/
│   ├── hartree.rs      # Hartree kernel (CubeCL)
│   ├── lda_xc.rs       # LDA XC kernel (CubeCL)
│   ├── v_eff.rs        # V_eff assembly kernel (CubeCL)
│   └── nonlocal.rs     # Non-local potential kernel (CubeCL, Phase 2)
└── shaders/            # (removed, replaced by CubeCL kernels)
```

## Trade-offs vs Current Approach

| Aspect | wgpu + WGSL (current) | CubeCL |
|--------|----------------------|--------|
| Kernel language | WGSL strings | Rust |
| Type safety | Runtime errors | Compile-time |
| Debugging | printf in shader | Rust tests on CPU |
| Backend support | Metal, Vulkan | Metal, Vulkan, CUDA, ROCm |
| Maturity | Stable, well-tested | Newer, evolving API |
| Learning curve | WGSL syntax | CubeCL macros |
| Complex kernels | Painful | Natural |
| f64 support | Limited | Backend-dependent |

## Risks

- **Maturity:** CubeCL is younger than wgpu. API may change between versions.
- **Compile times:** Procedural macros (`#[cube]`) add to compilation time.
- **Performance:** CubeCL-generated WGSL may not be as optimized as hand-written WGSL for simple kernels. Benchmark before committing.
- **Complexity:** For the three existing simple kernels, CubeCL adds framework overhead without much benefit. The value is in enabling Phase 2 (non-local on GPU) and Phase 3 (GPU FFT), which are genuinely complex.

## Recommendation

**Wait on this until the non-local potential GPU kernel is needed.** The current three WGSL shaders are simple, correct, and fast. CubeCL's value proposition is for the more complex kernels (non-local, FFT) that aren't implemented yet. When those become a priority, CubeCL avoids writing 200+ lines of WGSL for the non-local pair loop.

If CUDA support becomes important (for non-Apple hardware), CubeCL's multi-backend compilation is a stronger argument for adoption.

## Expected Impact

- **Developer experience:** Write and test GPU kernels in Rust instead of WGSL.
- **Performance (Phase 2):** Moving non-local potential to GPU could yield 5-10x speedup on that component (15-25% of SCF time), for ~10-20% total SCF improvement.
- **Portability:** Single codebase compiles to Metal, Vulkan, and CUDA.
- **Priority:** Medium-low. Evaluate when non-local GPU acceleration is on the roadmap.
