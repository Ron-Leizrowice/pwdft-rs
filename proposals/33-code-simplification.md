# Proposal 33: Code Simplification and Deduplication

## Problem

A full-codebase review across reuse, quality, and efficiency identified 25+ concrete issues. The most impactful fall into four categories:

1. **Duplicated constants and functions** — the same physics (E2 Coulomb, PZ correlation, j0 sinc, structure factors) is implemented multiple times, creating drift risk
2. **Hot-path inefficiencies** — `powf(1.0/3.0)` instead of `cbrt()`, per-iteration allocations for immutable data, sequential XC evaluation
3. **Copy-pasted SCF logic** — `run_scf` and `run_scf_spin` share ~70% of their code, with the spin path missing optimizations the unpolarized path has
4. **Bugs** — convergence failure always reports `delta: 0.0`, hiding the actual convergence state

## Findings

### Tier 1 — High impact, low risk

#### 1.1 E2 Coulomb constant defined in 4 places

| Location | Name | Scope |
|----------|------|-------|
| `src/consts.rs:13` | `E2_COULOMB` | `pub` (canonical) |
| `src/potential/hartree.rs:11` | `E2` | `pub` |
| `src/pseudopotential/mod.rs:101` | `E2` | `const` (function-local) |
| `src/pseudopotential/upf.rs:236` (test) | `e2` | `let` binding |

All four have the same value `14.399_645_351_950_548`. The `hartree::E2` is imported by `ewald.rs`, `scf/energy.rs`, `scf/mod.rs`, `gpu/mod.rs`, `tests/gpu_consistency.rs`, and `benches/gpu_benchmarks.rs`.

**Fix:** Delete `hartree::E2` and `pseudopotential::E2`. Use `crate::consts::E2_COULOMB` everywhere. Update imports in ewald.rs, scf/energy.rs, scf/mod.rs, gpu/mod.rs, and test/bench files.

#### 1.2 `powf(1.0/3.0)` → `cbrt()` in hot XC loops

`f64::cbrt()` is a dedicated hardware instruction, 3-5x faster than `powf()` (which routes through `exp(ln(x)/3)`). Found at:

| File | Line(s) | Expression |
|------|---------|------------|
| `src/potential/xc.rs` | 90 | `(3.0 * rho_bohr / PI).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 111 | `(3.0 / (4.0 * PI * rho_bohr)).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 237 | `(6.0 * rho_up_bohr / PI).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 242 | `(6.0 * rho_down_bohr / PI).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 274 | `(3.0 / (4.0 * PI * rho_bohr)).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 283 | `2.0_f64.powf(4.0 / 3.0)` |
| `src/potential/xc.rs` | 286-287 | `(1.0 + zeta).max(0.0).powf(4.0 / 3.0)` |
| `src/potential/xc.rs` | 292-293 | `(1.0 + zeta).max(0.0).powf(1.0 / 3.0)` |
| `src/ewald.rs` | 52 | `(n_atoms as f64 * PI / omega).powf(1.0 / 3.0)` |
| `src/scf/mixing.rs` | ~195 | Kerker G² threshold uses powf |

**Fix:** Replace `x.powf(1.0 / 3.0)` → `x.cbrt()`. Replace `x.powf(4.0 / 3.0)` → `x.cbrt() * x`. The constant `2.0_f64.powf(4.0 / 3.0)` can become `2.0_f64.cbrt() * 2.0`.

#### 1.3 Precompute `rho_core_half` in spin SCF

`src/scf/mod.rs:416-417`:
```rust
let rho_up_xc = add_core_density(&rho_up_r, &ctx.rho_core_r.iter().map(|&c| c / 2.0).collect::<Vec<_>>());
let rho_down_xc = add_core_density(&rho_down_r, &ctx.rho_core_r.iter().map(|&c| c / 2.0).collect::<Vec<_>>());
```

Every SCF iteration allocates two `Vec<f64>` of size `n_grid` to hold `rho_core / 2.0`. This data is constant.

**Fix:** Compute once before the loop:
```rust
let rho_core_half: Vec<f64> = ctx.rho_core_r.iter().map(|&c| c / 2.0).collect();
```
Then use `&rho_core_half` inside the loop.

#### 1.4 Use `assemble_v_eff` in spin-polarized path

`src/scf/mod.rs:424-427` assembles V_eff inline with sequential iteration:
```rust
let v_eff_up: Vec<Complex64> = ctx.v_local_fft.iter().zip(v_h_fft.iter()).zip(vxc_up_g.iter())
    .map(|((&vl, &vh), &vxc)| vl + vh + vxc).collect();
```

The function `assemble_v_eff` in `scf/energy.rs:125-137` does exactly this but with `par_iter()`.

**Fix:** Replace both inline blocks with:
```rust
let v_eff_up = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_up_g);
let v_eff_down = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_down_g);
```

#### 1.5 Replace inline FFT normalization with `density_r_to_g`

`src/scf/mod.rs:265-273` manually reimplements `density_r_to_g`:
```rust
let mut rho_g_new = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
for (i, &r) in rho_r_new.iter().enumerate() {
    rho_g_new[i] = Complex64::new(r, 0.0);
}
ctx.grid.fft.forward(&mut rho_g_new);
let fft_norm = 1.0 / ctx.n_grid as f64;
for v in &mut rho_g_new {
    *v *= fft_norm;
}
```

This is identical to `density_r_to_g(&mut ctx.grid.fft, &rho_r_new, &mut rho_g_new)`, which is already used elsewhere in the same function (lines 160, 337).

**Fix:** Replace with `density_r_to_g(&mut ctx.grid.fft, &rho_r_new, &mut rho_g_new);`

#### 1.6 Implement `real_to_g_space` via `density_r_to_g`

`src/scf/energy.rs` has two near-identical functions:
- `density_r_to_g` (lines 101-110): writes into pre-allocated buffer
- `real_to_g_space` (lines 113-122): allocates and returns a new `Vec`

**Fix:**
```rust
pub(crate) fn real_to_g_space(data_r: &[f64], fft: &mut FFT3D) -> Vec<Complex64> {
    let mut data_g = vec![Complex64::new(0.0, 0.0); data_r.len()];
    density_r_to_g(fft, data_r, &mut data_g);
    data_g
}
```

### Tier 2 — Bug fixes

#### 2.1 Convergence failure reports wrong delta

`src/scf/mod.rs:340-343` (unpolarized):
```rust
Err(PwdftError::ConvergenceFailure {
    iterations: ctx.params.max_iter,
    delta: density_diff(&rho_r, &rho_r, ctx.omega, ctx.n_grid), // always 0!
})
```

`src/scf/mod.rs:597-600` (spin):
```rust
Err(PwdftError::ConvergenceFailure {
    iterations: ctx.params.max_iter,
    delta: 0.0, // hardcoded!
})
```

**Fix:** Track the last computed `delta` in a variable before the loop, update it each iteration, and use it in the error:
```rust
let mut last_delta = f64::INFINITY;
for iter in 0..ctx.params.max_iter {
    // ...
    last_delta = delta;
    // ...
}
Err(PwdftError::ConvergenceFailure { iterations: ctx.params.max_iter, delta: last_delta })
```

### Tier 3 — Performance

#### 3.1 FFT allocates 3 temporary arrays per call

`src/fft.rs:53-65` (forward) and `69-81` (inverse) both do:
```rust
let mut a = Array3::from_shape_vec((nx, ny, nz), data.to_vec()).unwrap();
let mut b = Array3::zeros((nx, ny, nz));
```

For a 32^3 grid, each FFT call allocates ~1 MB. FFTs are called many times per SCF iteration (density, vxc, Kerker, convergence energy).

**Fix:** Store `a` and `b` as fields of `FFT3D`. On each call, copy data into `a` with `a.as_slice_mut().unwrap().copy_from_slice(data)` and zero `b` with `b.fill(Complex64::new(0.0, 0.0))`.

```rust
pub struct FFT3D {
    dims: [usize; 3],
    fwd_handlers: [FftHandler<f64>; 3],
    inv_handlers: [FftHandler<f64>; 3],
    buf_a: Array3<Complex64>,
    buf_b: Array3<Complex64>,
}
```

#### 3.2 XC grid evaluation is sequential

`src/potential/xc.rs:48-59` (`lda_xc_grid`) and `194-211` (`lda_xc_spin_grid`) process each grid point sequentially. Each point involves `cbrt()` and `ln()` calls — embarrassingly parallel.

**Fix:** Use rayon `par_iter` for both functions. For the unpolarized case:
```rust
pub fn lda_xc_grid(rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    use rayon::prelude::*;
    let results: Vec<_> = rho_r.par_iter().map(|&rho| lda_xc(rho)).collect();
    let exc = results.iter().map(|xc| xc.exc).collect();
    let vxc = results.iter().map(|xc| xc.vxc).collect();
    (exc, vxc)
}
```

#### 3.3 Spin-up and spin-down diagonalization run sequentially

`src/scf/mod.rs:430-440`: both spin channels use `par_iter` internally but the two blocks run sequentially. They share no mutable state.

**Fix:** Use `rayon::join` to run both spin channels concurrently:
```rust
let (kpoint_results_up, kpoint_results_down) = rayon::join(
    || ctx.kpoints.par_iter().enumerate().map(|(ik, kp)| { /* up */ }).collect::<Vec<_>>(),
    || ctx.kpoints.par_iter().enumerate().map(|(ik, kp)| { /* down */ }).collect::<Vec<_>>(),
);
```

### Tier 4 — Code quality cleanup

#### 4.1 `RHO_FLOOR` constant defined but unused

`consts.rs:19` defines `RHO_FLOOR = 1e-20` but xc.rs uses hardcoded `1e-30` in 4 places (lines 27, 177, 226, 269).

**Fix:** Update `RHO_FLOOR` to `1e-30` and use it in all 4 locations. Or keep `1e-20` and update the xc code — but the values should be consistent.

#### 4.2 `consts::PI` re-exports `std::f64::consts::PI`

Only used in `crystal.rs`. Every other file imports from `std::f64::consts::PI` directly.

**Fix:** Remove `consts::PI`, use `std::f64::consts::PI` in crystal.rs.

#### 4.3 GPU `read_staging_buffer` double-copies data

`src/gpu/mod.rs:445-449`:
```rust
let result: Vec<f32> = bytemuck::cast_slice(&data).to_vec();
drop(data);
buffer.unmap();
result[..n_floats].to_vec() // second copy!
```

**Fix:**
```rust
let result = bytemuck::cast_slice(&data)[..n_floats].to_vec();
drop(data);
buffer.unmap();
result
```

#### 4.4 Ewald structure factor uses manual cos/sin

`src/ewald.rs:72-80` accumulates `s_re` and `s_im` manually instead of using `Complex64::cis()`.

**Fix:**
```rust
let s: Complex64 = positions.iter().enumerate()
    .map(|(i, pos)| charges[i] * Complex64::cis(g.dot(pos)))
    .sum();
let s_sq = s.norm_sqr();
```

#### 4.5 `PseudopotentialData::n_projectors` and `has_nlcc` are redundant

- `n_projectors` always equals `beta_projectors.len()`
- `has_nlcc` always equals `!core_charge.is_empty()`

**Fix:** Remove the fields. Add methods:
```rust
pub fn n_projectors(&self) -> usize { self.beta_projectors.len() }
pub fn has_nlcc(&self) -> bool { !self.core_charge.is_empty() }
```
Update all usages (including UPF parser construction and test mock objects).

#### 4.6 `_z_val` unused parameter

`src/scf/initial_density.rs:149`: `_z_val: f64` is passed but never used.

**Fix:** Remove the parameter and update the call site.

### Tier 5 — Deferred (larger refactors)

These are real issues but require more design work:

- **Unify `run_scf` and `run_scf_spin`**: ~70% shared code. Would need a spin-generic SCF body.
- **Deduplicate PZ correlation**: `perdew_zunger_correlation` duplicates `pz_correlation_rs(rs, false)`.
- **Test fixture consolidation**: `si_crystal()` defined ~15 times across test modules.
- **Separate `hartree.rs` legacy functions**: `hartree_potential` and `hartree_energy` using G-vectors are superseded by `g_squared`-based versions in `scf/energy.rs`.
- **KB non-local separable form**: O(n_pw^2 * n_proj^2) could become O(n_pw * n_proj + n_proj^2).

## Implementation

Steps ordered by dependency and risk:

1. **Constants** (1.1, 4.1, 4.2): Consolidate E2, use RHO_FLOOR, remove PI re-export
2. **Easy dedup** (1.4, 1.5, 1.6, 4.3, 4.4): Use existing functions, fix double-copy, use cis()
3. **Performance** (1.2, 1.3): cbrt(), precompute rho_core_half
4. **Bug fix** (2.1): Track last delta for convergence failure
5. **FFT buffers** (3.1): Requires changing FFT3D struct — most invasive
6. **Parallelization** (3.2, 3.3): rayon in XC grid, join for spin channels
7. **Struct cleanup** (4.5, 4.6): Remove redundant fields — touches many files

## Verification

```bash
cargo clippy -q --all-targets          # no new warnings
cargo test                             # all tests pass
cargo test --test qe_validation        # physics unchanged
cargo bench --bench scf_benchmarks     # performance comparison (before/after)
```

## Estimated Effort

Tiers 1-4: one focused session. Tier 5 (larger refactors): separate proposals when ready.
