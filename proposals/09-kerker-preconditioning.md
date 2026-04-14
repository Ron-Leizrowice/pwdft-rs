# Proposal 09: Kerker Preconditioning

## Problem

Without preconditioning, Anderson mixing treats all Fourier components of the density residual equally. But the Hartree potential amplifies long-wavelength (small-G) residuals since V_H(G) = 4pi*e^2*rho(G)/G^2 — a small density error at G~0 creates a huge potential error, which feeds back into the next density, causing oscillation. This is **charge sloshing** and it makes metallic systems, large cells, and slab geometries fail to converge or require beta < 0.05 with hundreds of iterations.

Kerker preconditioning damps the dangerous long-wavelength components proportionally:

```
P(G) = |G|^2 / (|G|^2 + q_TF^2)
```

where q_TF is the Thomas-Fermi screening wavevector. At G=0, P->0 (full suppression). At large G, P->1 (no change). This keeps the condition number of the SCF map near 1 regardless of system size.

Every production plane-wave code implements this: QE (`mixing_mode='TF'`), VASP (`MIXPRE`), ABINIT (`diemac`/`dielng`), JDFTx (`qKerker`).

## References

- Kerker, G.P., Phys. Rev. B 23, 3082 (1981) — original paper
- QE docs: `mixing_mode` parameter in [pw.x input](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
- VASP wiki: [MIXPRE](https://www.vasp.at/wiki/index.php/MIXPRE)
- Raczkowski et al., Phys. Rev. B 64, 121101 (2001) — local Thomas-Fermi for inhomogeneous systems
- [arXiv:1707.00848](https://arxiv.org/abs/1707.00848) — Kerker for inhomogeneous systems

## Implementation

### Step 1: Preconditioner function

```rust
// src/scf/mixing.rs

/// Kerker preconditioner: damps long-wavelength density residuals to prevent
/// charge sloshing. Returns P(G) = |G|^2 / (|G|^2 + q_tf^2).
fn kerker_weight(g_squared: f64, q_tf_squared: f64) -> f64 {
    if g_squared < 1e-20 {
        0.0  // Suppress G=0 completely
    } else {
        g_squared / (g_squared + q_tf_squared)
    }
}
```

### Step 2: Move mixing to reciprocal space

The mixer currently operates on `rho_r: &[f64]` (real-space density). Kerker preconditioning is defined in reciprocal space, so the mixer needs access to the G-space residual.

Change `AndersonMixer` to store history in G-space and apply preconditioning to residuals:

```rust
pub struct AndersonMixer {
    beta: f64,
    max_history: usize,
    history_in: Vec<Vec<Complex64>>,   // was Vec<f64>
    history_res: Vec<Vec<Complex64>>,  // was Vec<f64>
    q_tf_squared: f64,                 // NEW
    g_squared: Vec<f64>,               // NEW: |G|^2 for each FFT grid point
}

impl AndersonMixer {
    pub fn new(beta: f64, max_history: usize, g_squared: Vec<f64>, q_tf: f64) -> Self {
        Self {
            beta,
            max_history,
            history_in: Vec::new(),
            history_res: Vec::new(),
            q_tf_squared: q_tf * q_tf,
            g_squared,
        }
    }

    /// Mix in reciprocal space with Kerker preconditioning.
    pub fn mix(&mut self, rho_in_g: &[Complex64], rho_out_g: &[Complex64]) -> Vec<Complex64> {
        let n = rho_in_g.len();
        // Preconditioned residual
        let residual: Vec<Complex64> = (0..n)
            .map(|i| {
                let r = rho_out_g[i] - rho_in_g[i];
                r * kerker_weight(self.g_squared[i], self.q_tf_squared)
            })
            .collect();

        // ... rest of Anderson algorithm unchanged, but on Complex64 ...
    }
}
```

### Step 3: Update SCF loop

In `src/scf/mod.rs`, the mixer call (line 321) changes from:

```rust
rho_r = mixer.mix(&rho_r, &rho_r_new);
density_r_to_g(&grid.fft, &rho_r, &mut rho_g);
```

to:

```rust
// Mix in G-space with Kerker preconditioning
let mut rho_out_g = vec![Complex64::new(0.0, 0.0); n_grid];
density_r_to_g(&grid.fft, &rho_r_new, &mut rho_out_g);
rho_g = mixer.mix(&rho_g, &rho_out_g);
// Transform back to real space
rho_r = g_to_real_space(&rho_g, &grid.fft);
```

### Step 4: Input configuration

Add to `ScfConfig` in `src/input.rs`:

```toml
[scf]
mixing_mode = "kerker"     # "plain" | "kerker" | "local-tf"
kerker_q = 1.0             # Thomas-Fermi wavevector in Bohr^-1 (default: auto from avg density)
```

Auto-estimation of q_TF from the average density:
```rust
let rho_avg = n_electrons / omega;  // e/Å^3
let rho_bohr = rho_avg * BOHR3;    // e/Bohr^3
let q_tf = (4.0 * (3.0 * PI * PI * rho_bohr).powf(1.0/3.0) / PI).sqrt();
```

### Future: local Thomas-Fermi (Phase 2)

For slabs and surfaces, the screening varies in space. QE's `local-TF` uses a position-dependent q_TF(r) derived from the local density. This requires computing the preconditioner in real space, which is more expensive but critical for metal/vacuum interfaces.

## Acceptance Criteria

1. **Metallic convergence test:** Run SCF on a metallic system (e.g., Al FCC with smearing) with `mixing_mode="plain"` and `mixing_mode="kerker"`. Kerker should converge in fewer iterations (typically 2-3x fewer) or converge where plain mixing diverges.
2. **Insulator equivalence:** For Si (current test case), Kerker should converge to the same total energy (within 1e-6 eV) as plain mixing.
3. **q_TF sensitivity:** Verify that results are insensitive to q_TF within a factor of 2 of the auto-estimated value.
4. **Plain mixing regression:** `mixing_mode="plain"` reproduces the current behavior exactly (no change to existing tests).
5. **QE comparison:** For a system computed with QE `mixing_mode='TF'`, verify that iteration counts are comparable.
