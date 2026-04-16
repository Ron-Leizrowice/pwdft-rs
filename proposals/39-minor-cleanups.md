# Proposal 39: Minor Code Quality Cleanups

## Problem

Six small, independent code quality issues found during codebase review. Each is low-risk and self-contained.

## Findings and Fixes

### 1. GPU `read_staging_buffer` double-copies data

`src/gpu/mod.rs:445-449`:
```rust
let result: Vec<f32> = bytemuck::cast_slice(&data).to_vec();  // copy 1
drop(data);
buffer.unmap();
result[..n_floats].to_vec()  // copy 2
```

**Fix:**
```rust
let result = bytemuck::cast_slice(&data)[..n_floats].to_vec();  // single copy
drop(data);
buffer.unmap();
result
```

### 2. Ewald structure factor uses manual cos/sin

`src/ewald.rs:72-80` manually accumulates `s_re` and `s_im`:
```rust
let mut s_re = 0.0;
let mut s_im = 0.0;
for (i, pos) in positions.iter().enumerate() {
    let phase = g.dot(pos);
    s_re += charges[i] * phase.cos();
    s_im += charges[i] * phase.sin();
}
let s_sq = s_re * s_re + s_im * s_im;
```

**Fix:** Use `Complex64::cis()` (already used elsewhere in the codebase):
```rust
let s: Complex64 = positions.iter().enumerate()
    .map(|(i, pos)| charges[i] * Complex64::cis(g.dot(pos)))
    .sum();
let s_sq = s.norm_sqr();
```

This requires adding `use num_complex::Complex64;` to the imports in `ewald.rs`.

### 3. `PseudopotentialData::n_projectors` is redundant

`src/pseudopotential/mod.rs:34`: `pub n_projectors: usize` always equals `beta_projectors.len()`. Used in 9 places across src/ (mostly test assertions like `assert_eq!(pp.beta_projectors.len(), pp.n_projectors)`).

**Fix:** Remove the field. Add a method:
```rust
pub fn n_projectors(&self) -> usize { self.beta_projectors.len() }
```

Update usages:
- `src/pseudopotential/upf.rs:122`: remove from struct construction
- `src/potential/nonlocal.rs:93`: `pp.n_projectors` → `pp.n_projectors()` (or `pp.beta_projectors.len()`)
- `src/ewald.rs:194,208`: remove from mock PP construction
- Tests: remove redundant `assert_eq!(pp.beta_projectors.len(), pp.n_projectors)` assertions

### 4. `PseudopotentialData::has_nlcc` is redundant

`src/pseudopotential/mod.rs:42`: `pub has_nlcc: bool` always equals `!core_charge.is_empty()`. Used in 4 places in src/.

**Fix:** Remove the field. Add a method:
```rust
pub fn has_nlcc(&self) -> bool { !self.core_charge.is_empty() }
```

Update:
- `src/pseudopotential/upf.rs:125`: remove from struct construction
- `src/scf/potentials.rs:62,73`: `pp.has_nlcc` → `pp.has_nlcc()` (unchanged call syntax)
- `src/ewald.rs:197,211`: remove from mock PP construction

### 5. `_z_val` unused parameter

`src/scf/initial_density.rs:148`: function `add_atomic_density_from_pp` takes `_z_val: f64` but never uses it (underscore prefix suppresses warning).

**Fix:** Remove the parameter. Update the call site (should be one location in the same file).

### 6. `apply_rotation` duplicates `SymmOp::apply`

`src/symmetry/detect.rs:190-196` defines `fn apply_rotation(r: &[[i32; 3]; 3], f: &[f64; 3]) -> [f64; 3]` which is identical to `SymmOp::apply` in `operations.rs:121-128`.

**Fix:** Remove `apply_rotation` from `detect.rs`. Replace calls with:
```rust
SymmOp { rotation: *rotation }.apply(&pos)
```
Or extract the shared logic into a free function in `operations.rs` used by both.

## Verification

```bash
cargo clippy -q --all-targets
cargo test
cargo test --features gpu  # for GPU change
```

## Estimated Effort

Under an hour total. Each fix is independent — can be applied in any order.
