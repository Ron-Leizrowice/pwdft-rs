# Proposal 23: XC GPU Shader Precision Audit and CPU/GPU Consistency

**Status:** Partially done. Constants unified (`G2_ZERO_THRESHOLD`, `RHO_FLOOR` in `consts.rs`), density floor consistent. CPU/GPU consistency test and error budget documentation NOT yet added.

## Problem

The LDA XC GPU shader (`src/gpu/shaders/lda_xc.wgsl`) runs in f32 while the CPU path (`src/potential/xc.rs`) runs in f64. Beyond the expected precision difference, there are specific issues:

### 1. Slater exchange formula

CPU (`xc.rs:90-91`):

```rust
let cbrt = (3.0 * rho_bohr / PI).powf(1.0 / 3.0);
let ex_ha = -0.75 * cbrt;
```

GPU (`lda_xc.wgsl:50-51`):

```wgsl
let cbrt_arg = pow(3.0 * rho_bohr / PI, 1.0 / 3.0);
let ex_ha = -0.75 * cbrt_arg;
```

Both use `ε_x = -0.75 × (3ρ/π)^{1/3}` which is the Slater exchange formula. This is correct — the full expression is `ε_x = -(3/4)(3/π)^{1/3} ρ^{1/3}`, where `-(3/4)(3/π)^{1/3} = -0.75 × (3/π)^{1/3}`. However, some references write it as `ε_x = -C_x × ρ^{1/3}` with `C_x = (3/4)(3/π)^{1/3} ≈ 0.7386`. The code factoring is equivalent but could benefit from a comment citing the formula source.

### 2. BOHR3 constant precision

GPU (`lda_xc.wgsl:16`):

```wgsl
const BOHR3: f32 = 0.14818471;
```

Exact value: `0.529177210903³ = 0.148184706...`. The shader value `0.14818471` differs by `4e-8` (relative), which is within f32 precision (~7 digits), but could be more precise:

```wgsl
const BOHR3: f32 = 0.1481847;  // 0.529177^3, 7 significant digits
```

### 3. Density floor inconsistency

CPU (`xc.rs:27`): `rho < 1e-30`
GPU (`lda_xc.wgsl:37`): `rho < 1e-20`

The f32 minimum normal is ~1.2e-38, so `1e-30` is representable in f32 but `1e-30` is unnecessarily small for both. The GPU uses `1e-20`, which is more reasonable. These should match (see Proposal 22).

### 4. PZ correlation parameters

Both CPU and GPU use identical Perdew-Zunger parameters. Spot-checked against the original paper (PRB 23, 5048, 1981) and NIST:

| Parameter | Code | Reference | Match? |
|-----------|------|-----------|--------|
| γ (rs≥1) | -0.1423 | -0.1423 | ✓ |
| β₁ (rs≥1) | 1.0529 | 1.0529 | ✓ |
| β₂ (rs≥1) | 0.3334 | 0.3334 | ✓ |
| A (rs<1) | 0.0311 | 0.0311 | ✓ |
| B (rs<1) | -0.048 | -0.048 | ✓ |
| C (rs<1) | 0.0020 | 0.002 | ✓ |
| D (rs<1) | -0.0116 | -0.0116 | ✓ |

Parameters are correct.

### 5. HA_TO_EV precision

GPU (`lda_xc.wgsl:15`): `const HA_TO_EV: f32 = 27.211386;`
CPU (`xc.rs:83`): `const HA_TO_EV: f64 = 27.211386245988;`

The GPU value truncates to 8 significant digits, which is at the limit of f32 (~7 digits). The value is correct within f32 precision.

## Implementation

### Step 1: Add formula citation comments

In both `xc.rs` and `lda_xc.wgsl`, add clear reference:

```rust
// Slater exchange: ε_x = -(3/4)(3ρ/π)^{1/3} in Hartree
// Reference: Slater, Phys. Rev. 81, 385 (1951)
// Equivalent to ε_x = -C_x ρ^{1/3} with C_x = (3/4)(3/π)^{1/3} ≈ 0.7386
```

```rust
// Perdew-Zunger parametrization of Ceperley-Alder correlation
// Reference: Perdew & Zunger, Phys. Rev. B 23, 5048 (1981), Table I
```

### Step 2: Unify density floor

Use `RHO_FLOOR` from Proposal 22 in both CPU and GPU paths. For the shader, either pass as a uniform or add a comment:

```wgsl
// Must match RHO_FLOOR in src/consts.rs
if (rho < 1e-20) { ... }
```

### Step 3: Add CPU/GPU consistency test

The existing test suite should verify CPU and GPU XC agree within f32 tolerance:

```rust
#[test]
#[cfg(feature = "gpu")]
fn test_gpu_cpu_xc_consistency() {
    let gpu = GpuAccelerator::try_new().unwrap();
    // Test across density range: vacuum to core
    let rho_r: Vec<f64> = (0..1000)
        .map(|i| 10.0_f64.powf(-8.0 + 8.0 * i as f64 / 999.0)) // 1e-8 to 1.0
        .collect();

    let (exc_cpu, vxc_cpu) = xc::lda_xc_grid(&rho_r);
    let (exc_gpu, vxc_gpu) = gpu.lda_xc(&rho_r);

    for i in 0..rho_r.len() {
        let rel_exc = (exc_cpu[i] - exc_gpu[i]).abs() / exc_cpu[i].abs().max(1e-10);
        let rel_vxc = (vxc_cpu[i] - vxc_gpu[i]).abs() / vxc_cpu[i].abs().max(1e-10);
        assert!(rel_exc < 1e-4, "exc mismatch at rho={:.2e}: cpu={:.6e} gpu={:.6e}",
            rho_r[i], exc_cpu[i], exc_gpu[i]);
        assert!(rel_vxc < 1e-4, "vxc mismatch at rho={:.2e}: cpu={:.6e} gpu={:.6e}",
            rho_r[i], vxc_cpu[i], vxc_gpu[i]);
    }
}
```

### Step 4: Document f32/f64 error budget

Add a comment in `gpu/mod.rs` documenting the expected error:

```rust
/// GPU kernels run in f32 for throughput. Expected relative errors vs f64 CPU:
/// - Hartree potential: ~1e-7 (linear operations, dominated by f32 precision)
/// - LDA XC: ~1e-5 (cube root and log introduce amplified rounding)
/// - V_eff assembly: ~1e-7 (linear addition)
///
/// For SCF convergence to 1e-6 eV, f32 error in XC (~1e-5 relative) can cause
/// the GPU path to converge to a slightly different energy (~0.01 meV difference).
/// This is acceptable for production use.
```

## Acceptance Criteria

1. **Formula citations:** Both CPU and GPU XC code cite Slater (1951) and Perdew-Zunger (1981).
2. **Density floor consistent:** CPU and GPU use the same threshold value.
3. **CPU/GPU test passes:** XC values agree within 1e-4 relative error across 5 orders of magnitude in density.
4. **Error budget documented:** Expected f32 vs f64 differences documented in `gpu/mod.rs`.
5. **No numerical regression:** Existing GPU tests pass unchanged.
