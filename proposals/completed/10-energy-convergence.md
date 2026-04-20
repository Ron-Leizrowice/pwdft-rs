# Proposal 10: Energy Convergence Criterion

## Problem

The SCF loop (`src/scf/mod.rs`, line 289) uses a single convergence criterion: the RMS density difference `delta_rho`. This has two problems:

1. **Not what users care about.** Users want converged energies and forces. A density threshold of 1e-6 e/Å^3 doesn't directly map to a known energy precision — the relationship depends on the system size, potential stiffness, and functional.

2. **Not what other codes do.** VASP uses energy difference (`EDIFF`) as its sole criterion. QE's `conv_thr` is an estimated energy error derived from the density residual, not the raw RMS. GPAW checks three criteria simultaneously (energy, density, eigenstates).

Currently, total energy is only computed after convergence is declared (line 306). This means the energy is never monitored during the SCF loop, so oscillating or slowly-converging energies go undetected.

## References

- VASP wiki: [EDIFF](https://www.vasp.at/wiki/index.php/EDIFF) — default 1e-4 eV
- QE docs: `conv_thr` in [pw.x input](https://www.quantum-espresso.org/Doc/INPUT_PW.html) — estimated energy, default 1e-6 Ry
- GPAW docs: [convergence criteria](https://gpaw.readthedocs.io/documentation/convergence.html) — energy + density + eigenstates

## Implementation

### Step 1: Compute energy every iteration

Move the energy computation into the SCF loop. Currently `compute_total_energy` is called once after convergence (line 306). The Ewald energy is constant across iterations, so cache it outside the loop.

In `src/scf/mod.rs`, before the SCF loop:

```rust
// Cache Ewald energy (constant across iterations)
let e_ewald = crate::ewald::ewald_energy(crystal, pseudopotentials);
```

Inside the loop, after occupations are computed (after line 272):

```rust
let e_total = compute_total_energy_fast(
    &eigenvalues_all, &occupations, kpoints,
    &rho_r, &rho_g, &g_squared, &_exc_r, &vxc_r,
    omega, e_ewald,
);
```

The `_fast` variant skips recomputing XC (already have `_exc_r`, `vxc_r` from step 2) and Ewald (cached):

```rust
fn compute_total_energy_fast(
    eigenvalues: &[Vec<f64>],
    occupations: &[Vec<f64>],
    kpoints: &[KPoint],
    rho_r: &[f64],
    rho_g: &[Complex64],
    g_squared: &[f64],
    exc_r: &[f64],
    vxc_r: &[f64],
    omega: f64,
    e_ewald: f64,
) -> f64 {
    let n_grid = rho_g.len();
    let dvol = omega / n_grid as f64;

    let e_band: f64 = /* same as current */;
    let e_hartree: f64 = /* same as current */;
    let e_xc = xc::lda_xc_energy(rho_r, exc_r, omega);
    let e_vxc: f64 = rho_r.iter().zip(vxc_r.iter())
        .map(|(&rho, &vxc)| rho * vxc * dvol).sum();

    e_band - e_hartree + e_xc - e_vxc + e_ewald
}
```

### Step 2: Track energy change

```rust
let mut e_prev: Option<f64> = None;

for iter in 0..params.max_iter {
    // ... existing steps 1-6 ...

    let e_total = compute_total_energy_fast(/* ... */);
    let de = e_prev.map(|ep| (e_total - ep).abs());
    e_prev = Some(e_total);

    // Convergence check: BOTH criteria must be satisfied
    let rho_converged = delta < params.conv_threshold;
    let energy_converged = de.map_or(false, |de| de < params.energy_threshold);

    info!(
        "SCF iter {}: E={:.6} eV  dE={:.2e}  Δρ={:.2e}",
        iter + 1, e_total,
        de.unwrap_or(f64::NAN), delta
    );

    if rho_converged && energy_converged {
        info!("SCF converged after {} iterations", iter + 1);
        // ...
    }
}
```

### Step 3: Input configuration

Add to `ScfConfig` in `src/input.rs`:

```rust
/// Energy convergence threshold in eV. SCF converges when |dE| < threshold.
#[serde(default = "default_energy_thr")]
pub energy_threshold: f64,
```

```rust
fn default_energy_thr() -> f64 { 1e-5 }  // 1e-5 eV ~ 7e-7 Ry
```

```toml
[scf]
energy_threshold = 1e-5    # eV (default)
conv_threshold = 1e-6      # e/Å³ RMS density (existing)
```

Both criteria must be satisfied for convergence. This is stricter than either alone but prevents false convergence where density is stable but energy is still drifting (or vice versa).

### Step 4: Improved logging

Replace the current single-line log (line 287) with a tabular format:

```text
SCF iter  1: E= -215.483721 eV  dE=       N/A  Δρ= 3.2e-02
SCF iter  2: E= -216.012455 eV  dE= 5.29e-01  Δρ= 8.4e-03
SCF iter  3: E= -216.148332 eV  dE= 1.36e-01  Δρ= 2.1e-03
...
SCF iter 12: E= -216.183245 eV  dE= 4.7e-07   Δρ= 8.1e-07  <- converged
```

## Acceptance Criteria

1. **Total energy is computed and logged at every SCF iteration**, not just at convergence.
2. **Energy change `|dE|` is reported** alongside density change at each iteration.
3. **Convergence requires both** `|dE| < energy_threshold` AND `delta_rho < conv_threshold`.
4. **Existing tests pass unchanged** — Si SCF converges to the same energy.
5. **No performance regression:** the extra energy computation per iteration adds < 5% to total SCF time (XC is already computed; Ewald is cached).
6. **New test:** construct a case where density converges but energy has not (or vice versa) and verify that the dual criterion correctly prevents premature convergence.
