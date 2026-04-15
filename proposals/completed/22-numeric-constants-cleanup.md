# Proposal 22: Numeric Constants and Code Hygiene Cleanup

## Problem

Several recurring code patterns reduce maintainability and risk silent inconsistencies:

### 1. Magic number `1e-20` for G² threshold

Used in 5+ locations with the same value but no shared constant:

- `src/scf/mod.rs:441` — `hartree_on_fft_grid()`
- `src/scf/mod.rs:568` — `compute_total_energy_from_components()`
- `src/potential/hartree.rs:33,59` — `hartree_potential()` and energy
- `src/gpu/shaders/hartree.wgsl:25` — GPU Hartree kernel
- `src/gpu/shaders/lda_xc.wgsl:37` — GPU XC kernel

If one location is updated, the others silently diverge. Additionally, `1e-20` is extremely small. Floating-point computation of `|G|²` near zero can produce values like `1e-28` (exact zero at Γ) or `1e-16` (rounding noise). A threshold of `1e-12` is physically and numerically more appropriate — it corresponds to `|G| ≈ 1e-6 Å⁻¹`, far below any meaningful reciprocal lattice vector.

### 2. `HA_TO_EV` and `BOHR3` defined locally

`src/potential/xc.rs` defines `const HA_TO_EV: f64 = 27.211386245988` twice (lines 83, 106) and computes `let bohr3 = 0.529177210903_f64.powi(3)` twice (lines 86, 107). These should live in `src/consts.rs`.

### 3. `Complex64::new(phase.cos(), phase.sin())` instead of `cis()`

Used in `nonlocal.rs:148`, `initial_density.rs:176`, `ewald.rs:63-64`, `local.rs`. The `Complex64::cis(phase)` function does the same thing with clearer intent and potentially a single `sincos` call.

### 4. Sequential `lda_xc_grid()`

`src/potential/xc.rs:48-58` evaluates XC point-by-point in a serial loop. For grids of 100³+ points, this is an easy parallelization win.

### 5. XC density floor too small

`src/potential/xc.rs:27` uses `rho < 1e-30` as the floor. In f64 arithmetic, densities in vacuum regions are typically ~1e-8 to 1e-15 e/ų. A floor of `1e-30` means the XC functional is evaluated on noise. The GPU shader uses `1e-20` (line 37). These should be consistent and physically motivated.

## Implementation

### Step 1: Add constants to `consts.rs`

```rust
// src/consts.rs

/// Hartree to electronvolt conversion.
pub const HA_TO_EV: f64 = 27.211386245988;

/// Rydberg to electronvolt conversion.
pub const RY_TO_EV: f64 = HA_TO_EV / 2.0;

/// Bohr radius in Ångströms.
pub const BOHR_TO_ANG: f64 = 0.529177210903;

/// Bohr³ in ų (volume conversion factor).
pub const BOHR3_TO_ANG3: f64 = BOHR_TO_ANG * BOHR_TO_ANG * BOHR_TO_ANG;

/// Threshold for treating |G|² as zero (skip G=0 in Coulomb sums).
/// Corresponds to |G| ≈ 1e-6 Å⁻¹, far below any physical reciprocal vector.
pub const G2_ZERO_THRESHOLD: f64 = 1e-12;

/// Minimum electron density for XC evaluation (e/ų).
/// Below this, XC energy and potential are set to zero.
/// Physical vacuum densities are ~1e-8 e/ų; this floor avoids evaluating
/// XC on numerical noise.
pub const RHO_FLOOR: f64 = 1e-20;
```

### Step 2: Replace magic numbers

In all Rust files using `1e-20` for G² checks:

```rust
use crate::consts::G2_ZERO_THRESHOLD;

// Before:
if g2 > 1e-20 { rho.norm_sqr() * fourpi_e2 / g2 } else { 0.0 }

// After:
if g2 > G2_ZERO_THRESHOLD { rho.norm_sqr() * fourpi_e2 / g2 } else { 0.0 }
```

In `potential/xc.rs`:

```rust
use crate::consts::{HA_TO_EV, BOHR3_TO_ANG3, RHO_FLOOR};

fn slater_exchange(rho: f64) -> (f64, f64) {
    let rho_bohr = rho * BOHR3_TO_ANG3;
    // ...
    (ex_ha * HA_TO_EV, vx_ha * HA_TO_EV)
}
```

For the WGSL shaders, the constants are embedded in shader source. Add a comment referencing `consts.rs`:

```wgsl
// Must match G2_ZERO_THRESHOLD in src/consts.rs
if (g2 < 1e-12) { ... }
```

Or, pass the threshold as a uniform parameter to the shader from the CPU side.

### Step 3: Use `Complex64::cis()`

In all locations computing `exp(iφ)`:

```rust
// Before:
Complex64::new(phase.cos(), phase.sin())

// After:
Complex64::cis(phase)
```

Affected files: `nonlocal.rs`, `initial_density.rs`, `ewald.rs`, `local.rs`.

### Step 4: Parallelize `lda_xc_grid`

```rust
use rayon::prelude::*;

pub fn lda_xc_grid(rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let results: Vec<XcPoint> = rho_r.par_iter()
        .map(|&rho| lda_xc(rho))
        .collect();

    let exc = results.iter().map(|xc| xc.exc).collect();
    let vxc = results.iter().map(|xc| xc.vxc).collect();
    (exc, vxc)
}
```

Or, to avoid the intermediate `Vec<XcPoint>`:

```rust
pub fn lda_xc_grid(rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut exc = vec![0.0; rho_r.len()];
    let mut vxc = vec![0.0; rho_r.len()];

    rho_r.par_iter().zip(exc.par_iter_mut()).zip(vxc.par_iter_mut())
        .for_each(|((&rho, exc_i), vxc_i)| {
            let xc = lda_xc(rho);
            *exc_i = xc.exc;
            *vxc_i = xc.vxc;
        });

    (exc, vxc)
}
```

### Step 5: Remove unused `_z_val` parameter

In `src/scf/initial_density.rs:148`:

```rust
// Before:
fn add_atomic_density_from_pp(
    pp: &PseudopotentialData,
    tau: &nalgebra::Vector3<f64>,
    _z_val: f64,
    ...
)

// After:
fn add_atomic_density_from_pp(
    pp: &PseudopotentialData,
    tau: &nalgebra::Vector3<f64>,
    ...
)
```

Update callers accordingly.

## Acceptance Criteria

1. **No magic numbers:** `grep -rn '1e-20' src/` returns only the WGSL shaders (which have a comment referencing `consts.rs`) and `consts.rs` itself.
2. **Single source of truth:** `HA_TO_EV`, `BOHR3_TO_ANG3`, `G2_ZERO_THRESHOLD`, and `RHO_FLOOR` each defined once in `consts.rs`.
3. **`Complex64::cis`:** No remaining `Complex64::new(_.cos(), _.sin())` patterns.
4. **XC parallelized:** `lda_xc_grid` uses rayon. Benchmark shows improvement for grids > 50³.
5. **All tests pass:** No numerical regressions from threshold changes.
6. **GPU/CPU consistency:** WGSL shader thresholds match `consts.rs` values (documented in comments or passed as uniforms).
