---
id: DDUP
status: active
priority: high
complexity: small
risk: low
depends_on: []
blocks: [CFGN]
---

# DDUP: SCF Code Deduplication

## Problem

The spin-polarized SCF path (`run_scf_spin`) contains inline reimplementations of helper functions that already exist in `scf/energy.rs`. It also has per-iteration allocations for constant data. Five specific issues:

1. **Inline V_eff assembly** (lines 424-427): sequential `.iter().zip().map()` instead of `assemble_v_eff()` which uses `par_iter()`
2. **Inline FFT normalization** (lines 265-273 in `run_scf`): 8 lines of manual real-to-G-space conversion that duplicates `density_r_to_g()`
3. **`real_to_g_space` duplicates `density_r_to_g`**: two near-identical functions in `energy.rs`
4. **Per-iteration `rho_core/2` allocation** (lines 416-417): allocates 2 × n_grid `Vec<f64>` every iteration for constant data
5. **Repeated occupation computation** (lines 466-470, 482-484): identical pattern 4 times

## Implementation

### Step 1: Replace inline V_eff assembly with `assemble_v_eff`

In `src/scf/mod.rs`, replace lines 424-427:

```rust
// Before:
let v_eff_up: Vec<Complex64> = ctx.v_local_fft.iter().zip(v_h_fft.iter()).zip(vxc_up_g.iter())
    .map(|((&vl, &vh), &vxc)| vl + vh + vxc).collect();
let v_eff_down: Vec<Complex64> = ctx.v_local_fft.iter().zip(v_h_fft.iter()).zip(vxc_down_g.iter())
    .map(|((&vl, &vh), &vxc)| vl + vh + vxc).collect();

// After:
let v_eff_up = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_up_g);
let v_eff_down = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_down_g);
```

This also gains `par_iter()` parallelism.

### Step 2: Replace inline FFT normalization with `density_r_to_g`

In `src/scf/mod.rs`, replace lines 265-273:

```rust
// Before (8 lines):
let mut rho_g_new = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
for (i, &r) in rho_r_new.iter().enumerate() {
    rho_g_new[i] = Complex64::new(r, 0.0);
}
ctx.grid.fft.forward(&mut rho_g_new);
let fft_norm = 1.0 / ctx.n_grid as f64;
for v in &mut rho_g_new { *v *= fft_norm; }

// After (2 lines):
let mut rho_g_new = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
density_r_to_g(&mut ctx.grid.fft, &rho_r_new, &mut rho_g_new);
```

### Step 3: Implement `real_to_g_space` via `density_r_to_g`

In `src/scf/energy.rs`, replace `real_to_g_space` body (lines 113-122):

```rust
pub(crate) fn real_to_g_space(data_r: &[f64], fft: &mut FFT3D) -> Vec<Complex64> {
    let mut data_g = vec![Complex64::new(0.0, 0.0); data_r.len()];
    density_r_to_g(fft, data_r, &mut data_g);
    data_g
}
```

### Step 4: Precompute `rho_core_half`

In `run_scf_spin`, before the SCF loop (after line 402):

```rust
let rho_core_half: Vec<f64> = ctx.rho_core_r.iter().map(|&c| c / 2.0).collect();
```

Then replace lines 416-417:
```rust
let rho_up_xc = add_core_density(&rho_up_r, &rho_core_half);
let rho_down_xc = add_core_density(&rho_down_r, &rho_core_half);
```

### Step 5: Extract occupation computation helper

Add to `scf/mod.rs` or `smearing.rs`:

```rust
fn compute_occupations(
    eigenvalues: &[Vec<f64>],
    scheme: smearing::SmearingScheme,
    fermi_energy: f64,
    sigma: f64,
    spin_factor: f64,
) -> Vec<Vec<f64>> {
    eigenvalues.iter().map(|evs| {
        evs.iter().map(|&e| smearing::occupation(scheme, e, fermi_energy, sigma, spin_factor)).collect()
    }).collect()
}
```

Replace the 4 identical occupation blocks (lines 238-245, 466-468, 469-471, 482-484, 485-487) with calls to this function.

## Verification

```bash
cargo test
cargo test --test qe_validation    # physics unchanged
cargo test --test spin_polarization # spin behavior unchanged
```

## Estimated Effort

Under an hour. All replacements use existing functions.
