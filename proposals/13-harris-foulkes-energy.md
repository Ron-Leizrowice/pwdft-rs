# Proposal 13: Harris-Foulkes Energy

## Problem

The standard Kohn-Sham total energy uses the output density for double-counting corrections. Before self-consistency is reached, this introduces first-order errors in the density into the energy estimate, causing the reported energy to oscillate during early SCF iterations.

The Harris-Foulkes (HF) energy uses the **input** density for all double-counting corrections but the **output** eigenvalues from the diagonalization:

```
E^HF = sum_i f_i epsilon_i(V[rho_in]) + E_H[rho_in] + E_xc[rho_in] - int V_eff[rho_in] rho_in dr + E_Ewald
```

Equivalently, using the same double-counting correction structure:

```
E^HF = E_band[rho_out] - E_H[rho_in] + E_xc[rho_in] - int V_xc[rho_in] rho_in dr + E_Ewald
```

This is **stationary** at self-consistency: first-order density errors cancel, so E^HF converges to E^KS from above with quadratic error. It provides a better energy estimate during early SCF iterations and serves as a convergence quality indicator:

```
|E^HF - E^KS| -> 0   as   rho_in -> rho_out
```

VASP reports both energies at every iteration. The difference gives a direct, interpretable measure of how far from self-consistency the calculation is.

## References

- Harris, J., Phys. Rev. B 31, 1770 (1985)
- Foulkes, W.M.C. & Haydock, R., Phys. Rev. B 39, 12520 (1989)
- VASP wiki: [Harris-Foulkes functional](https://www.vasp.at/wiki/index.php/Harris-Foulkes_functional)

## Implementation

### Step 1: Compute E^HF alongside E^KS

All ingredients for E^HF are already computed during each SCF iteration — they just aren't assembled into an energy. The key insight is that E^HF uses `rho_in` (the input density) for Hartree and XC, while E^KS uses `rho_out` (the density from new wavefunctions).

In `src/scf/mod.rs`, inside the SCF loop, after the eigensolve and before mixing:

```rust
// E^KS: uses rho_out for corrections (computed after density step)
// E^HF: uses rho_in for corrections (already computed in steps 1-2)

// E_H[rho_in] is already computed as part of v_h_fft
let e_hartree_in: f64 = rho_g.iter().zip(g_squared.iter())
    .map(|(rho, &g2)| {
        if g2 > 1e-20 { rho.norm_sqr() * fourpi_e2 / g2 } else { 0.0 }
    })
    .sum::<f64>() * 0.5 * omega;

// E_xc[rho_in] and int V_xc rho_in dr — from the XC evaluation in step 2
// _exc_r and vxc_r are computed from rho_r (the INPUT density)
let e_xc_in = xc::lda_xc_energy(&rho_r, &_exc_r, omega);
let e_vxc_in: f64 = rho_r.iter().zip(vxc_r.iter())
    .map(|(&rho, &vxc)| rho * vxc * dvol).sum();

// E_band comes from the OUTPUT eigenvalues
let e_band: f64 = eigenvalues_all.iter().zip(occupations.iter()).zip(kpoints.iter())
    .map(|((evs, occs), kp)| {
        evs.iter().zip(occs.iter()).map(|(&e, &f)| f * kp.weight * e).sum::<f64>()
    }).sum();

let e_harris = e_band - e_hartree_in + e_xc_in - e_vxc_in + e_ewald;
```

Note: `rho_r`, `_exc_r`, `vxc_r` in the current SCF loop (lines 220-227) are derived from the **input** density, not the output. So E^HF is essentially free — no extra computation needed.

### Step 2: Log both energies

```rust
info!(
    "SCF iter {}: E^KS={:.6}  E^HF={:.6}  |HF-KS|={:.2e}  Δρ={:.2e}",
    iter + 1, e_ks, e_harris, (e_harris - e_ks).abs(), delta
);
```

### Step 3: Use as convergence indicator

The Harris-Foulkes energy difference can serve as an auxiliary convergence check or diagnostic:

```rust
let hf_diff = (e_harris - e_ks).abs();
if hf_diff > 0.01 && delta < params.conv_threshold {
    log::warn!(
        "Density converged but Harris-Foulkes difference is {:.2e} eV — \
         energy may not be reliable. Consider tightening conv_threshold.",
        hf_diff
    );
}
```

### Step 4: Include in ScfResult

```rust
pub struct ScfResult {
    pub total_energy: f64,
    pub harris_foulkes_energy: f64,  // NEW
    // ... existing fields ...
}
```

## Cost

Zero extra computation. E^HF reuses:
- `e_band` from the eigenvalue sum (already computed for E^KS or will be with Proposal 10)
- `e_hartree_in` from the Hartree potential step (one extra reduction over `rho_g`)
- `e_xc_in`, `e_vxc_in` from the XC step (already evaluated on `rho_r`)
- `e_ewald` cached once per SCF

The only new work is one dot product over `rho_g` for the Hartree energy of the input density, which is O(n_grid).

## Acceptance Criteria

1. **E^HF equals E^KS at convergence:** After SCF converges, `|E^HF - E^KS| < 1e-8 eV`.
2. **E^HF is more stable than E^KS during early iterations:** Plot both vs iteration number; E^HF should vary less than E^KS in the first 5-10 iterations.
3. **E^HF approaches from above:** For a well-behaved system (Si), E^HF >= E^KS at each iteration (Harris-Foulkes is an upper bound when the starting density is reasonable).
4. **`|E^HF - E^KS|` decreases monotonically** in a well-converging calculation.
5. **Both energies logged** at every iteration for monitoring.
6. **No performance regression** — verified by benchmarks showing < 1% time increase.
