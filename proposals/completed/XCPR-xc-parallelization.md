---
id: XCPR
status: completed
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# XCPR: XC and Spin Diagonalization Parallelization

## Problem

Two embarrassingly parallel operations run sequentially:

### 1. XC grid evaluation

`lda_xc_grid` (`src/potential/xc.rs:48-59`) and `lda_xc_spin_grid` (lines 194-211) process each grid point in a sequential `for` loop. Each point requires `cbrt()`, `ln()`, and division — significant per-element work that benefits from parallelism.

### 2. Spin-channel diagonalization

In `run_scf_spin` (`src/scf/mod.rs:430-440`), spin-up k-points are diagonalized in parallel, then spin-down runs in parallel — but the two blocks run sequentially. The eigensolve is the most expensive step per iteration and the two channels share no mutable state.

## Implementation

### Step 1: Parallelize `lda_xc_grid`

```rust
pub fn lda_xc_grid(rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    use rayon::prelude::*;
    let (exc, vxc): (Vec<f64>, Vec<f64>) = rho_r
        .par_iter()
        .map(|&rho| {
            let xc = lda_xc(rho);
            (xc.exc, xc.vxc)
        })
        .unzip();
    (exc, vxc)
}
```

### Step 2: Parallelize `lda_xc_spin_grid`

```rust
pub fn lda_xc_spin_grid(
    rho_up_r: &[f64],
    rho_down_r: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    use rayon::prelude::*;
    let results: Vec<_> = rho_up_r
        .par_iter()
        .zip(rho_down_r.par_iter())
        .map(|(&ru, &rd)| {
            let xc = lda_xc_spin(ru, rd);
            (xc.exc, xc.vxc_up, xc.vxc_down)
        })
        .collect();
    let mut exc = Vec::with_capacity(results.len());
    let mut vxc_up = Vec::with_capacity(results.len());
    let mut vxc_down = Vec::with_capacity(results.len());
    for (e, vu, vd) in results {
        exc.push(e);
        vxc_up.push(vu);
        vxc_down.push(vd);
    }
    (exc, vxc_up, vxc_down)
}
```

(rayon's `unzip` only handles pairs; for triples, collect then split.)

### Step 3: Parallelize spin-channel diagonalization with `rayon::join`

In `src/scf/mod.rs`, replace the sequential spin-up then spin-down blocks (lines 430-440):

```rust
let (kpoint_results_up, kpoint_results_down) = rayon::join(
    || {
        ctx.kpoints.par_iter().enumerate().map(|(ik, kp)| {
            let mut h = build_hamiltonian_with_v_eff(ctx.basis, &kp.k, &v_eff_up, ctx.grid.dims);
            ctx.vnl_cache[ik].add_to_hamiltonian(&mut h, ctx.crystal, ctx.basis, &kp.k);
            dense::diagonalize_lowest(&h, ctx.params.n_bands)
        }).collect::<Vec<_>>()
    },
    || {
        ctx.kpoints.par_iter().enumerate().map(|(ik, kp)| {
            let mut h = build_hamiltonian_with_v_eff(ctx.basis, &kp.k, &v_eff_down, ctx.grid.dims);
            ctx.vnl_cache[ik].add_to_hamiltonian(&mut h, ctx.crystal, ctx.basis, &kp.k);
            dense::diagonalize_lowest(&h, ctx.params.n_bands)
        }).collect::<Vec<_>>()
    },
);
```

## Verification

```bash
cargo test                          # correctness unchanged
cargo test --test spin_polarization # spin-specific correctness
cargo test --test qe_validation     # physics unchanged
cargo bench --bench scf_benchmarks  # measure speedup
```

## Estimated Effort

Under an hour. The XC change is a drop-in rayon replacement. The `rayon::join` requires verifying `ctx` fields are `Sync` (they are — all shared references).

## Status (2026-04-17)

**Steps 1 + 2 done** (PR #19, merged). `lda_xc_grid` and `lda_xc_spin_grid` parallelized via rayon. See Performance Engineer logbook for measurements.

**Step 3 pending:** `rayon::join` for spin-channel diagonalization in `run_scf_spin`. Independent work — can be picked up at any time.
