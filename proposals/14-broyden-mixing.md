# Proposal 14: Broyden Mixing and Adaptive Beta

> **Note:** Line numbers reference the pre-ScfContext codebase. Verify locations before implementing.

## Problem

Anderson/Pulay mixing is effective but not always optimal:

1. **Fixed beta:** The mixing parameter beta = 0.3 is hardcoded for the full calculation. Metals need beta ~ 0.02-0.05 to avoid divergence, while insulators converge fastest at beta ~ 0.5-0.7. A fixed value is either too conservative (slow) or too aggressive (unstable).

2. **No quasi-Newton acceleration:** Anderson mixing finds the optimal linear combination of past residuals, but doesn't build an approximate Jacobian of the SCF map. Broyden's method does, and often converges in fewer iterations for difficult systems.

VASP defaults to modified Broyden (IMIX=4, Johnson 1988 as corrected by Eyert/Kresse). ABINIT offers both Anderson (iscf=7) and Broyden (iscf=6). QE uses Anderson but with Kerker preconditioning which compensates for the lack of Broyden.

## References

- Broyden, C.G., Math. Comp. 19, 577 (1965) — original Broyden's method
- Johnson, D.D., Phys. Rev. B 38, 12807 (1988) — modified Broyden for DFT
- Eyert, V., J. Comp. Phys. 124, 271 (1996) — corrected Johnson algorithm
- Marks, L.D. & Luke, D.R., Phys. Rev. B 78, 075114 (2008) — robust mixing review
- Banerjee et al., J. Chem. Theory Comput. 12, 3053 (2016) — periodic Pulay
- VASP wiki: [IMIX](https://www.vasp.at/wiki/index.php/IMIX)

## Implementation

### Step 1: Mixing trait

Abstract the mixer interface to support multiple algorithms:

```rust
// src/scf/mixing.rs

/// Trait for density mixing algorithms.
pub trait DensityMixer {
    /// Mix input and output densities to produce the next input.
    /// Operates in reciprocal space (Complex64 for Kerker compatibility).
    fn mix(&mut self, rho_in: &[Complex64], rho_out: &[Complex64]) -> Vec<Complex64>;

    /// Reset mixing history (e.g., when restarting).
    fn reset(&mut self);
}
```

### Step 2: Modified Broyden (Johnson/Eyert)

The modified Broyden method avoids storing the full inverse Jacobian by working with a rank-m update built from the iteration history:

```rust
pub struct BroydenMixer {
    beta: f64,
    max_history: usize,
    // History of input vectors and residuals
    delta_rho: Vec<Vec<Complex64>>,   // delta_rho_i = rho_in_{i+1} - rho_in_i
    delta_res: Vec<Vec<Complex64>>,   // delta_R_i = R_{i+1} - R_i
    prev_rho_in: Option<Vec<Complex64>>,
    prev_residual: Option<Vec<Complex64>>,
    // Kerker weights (if enabled)
    kerker_weights: Option<Vec<f64>>,
}

impl DensityMixer for BroydenMixer {
    fn mix(&mut self, rho_in: &[Complex64], rho_out: &[Complex64]) -> Vec<Complex64> {
        let n = rho_in.len();
        let residual: Vec<Complex64> = (0..n)
            .map(|i| rho_out[i] - rho_in[i])
            .collect();

        if let (Some(prev_in), Some(prev_res)) = (&self.prev_rho_in, &self.prev_residual) {
            // Compute differences
            let dr: Vec<Complex64> = (0..n).map(|i| rho_in[i] - prev_in[i]).collect();
            let df: Vec<Complex64> = (0..n).map(|i| residual[i] - prev_res[i]).collect();
            self.delta_rho.push(dr);
            self.delta_res.push(df);

            // Trim history
            if self.delta_rho.len() > self.max_history {
                self.delta_rho.remove(0);
                self.delta_res.remove(0);
            }
        }

        // Build the Broyden update:
        // H_inv * R = beta * R + sum_i gamma_i * (delta_rho_i - beta * delta_res_i)
        // where gamma is determined by solving a small linear system

        let m = self.delta_rho.len();
        let mut result = precondition(&residual, self.beta, &self.kerker_weights);

        if m > 0 {
            // Build overlap matrix: A_{ij} = <delta_res_i | delta_res_j>
            let mut a_mat = vec![0.0; m * m];
            let mut b_vec = vec![Complex64::new(0.0, 0.0); m];

            for i in 0..m {
                b_vec[i] = cdot(&self.delta_res[i], &residual);
                for j in 0..m {
                    a_mat[i * m + j] = cdot(&self.delta_res[i], &self.delta_res[j]).re;
                }
            }

            // Solve for gamma
            let gamma = solve_real_system(&a_mat, &b_vec.iter().map(|c| c.re).collect::<Vec<_>>(), m);

            // Apply Broyden correction
            for i in 0..m {
                for k in 0..n {
                    result[k] += gamma[i] * (self.delta_rho[i][k]
                        - precondition_single(self.delta_res[i][k], self.beta, k, &self.kerker_weights));
                }
            }
        }

        // Update: rho_new = rho_in + H_inv * R
        let rho_new: Vec<Complex64> = (0..n).map(|i| rho_in[i] + result[i]).collect();

        self.prev_rho_in = Some(rho_in.to_vec());
        self.prev_residual = Some(residual);

        rho_new
    }

    fn reset(&mut self) {
        self.delta_rho.clear();
        self.delta_res.clear();
        self.prev_rho_in = None;
        self.prev_residual = None;
    }
}
```

### Step 3: Adaptive beta

Monitor convergence quality and adjust beta:

```rust
pub struct AdaptiveMixer<M: DensityMixer> {
    inner: M,
    beta_min: f64,
    beta_max: f64,
    prev_residual_norm: Option<f64>,
}

impl<M: DensityMixer> DensityMixer for AdaptiveMixer<M> {
    fn mix(&mut self, rho_in: &[Complex64], rho_out: &[Complex64]) -> Vec<Complex64> {
        let res_norm: f64 = rho_out.iter().zip(rho_in.iter())
            .map(|(o, i)| (o - i).norm_sqr()).sum::<f64>().sqrt();

        if let Some(prev_norm) = self.prev_residual_norm {
            let ratio = res_norm / prev_norm;
            if ratio > 1.0 {
                // Residual growing: reduce beta
                self.inner.set_beta((self.inner.beta() * 0.5).max(self.beta_min));
                log::info!("Adaptive mixing: reducing beta to {:.3}", self.inner.beta());
            } else if ratio < 0.5 {
                // Converging well: increase beta cautiously
                self.inner.set_beta((self.inner.beta() * 1.2).min(self.beta_max));
            }
        }
        self.prev_residual_norm = Some(res_norm);

        self.inner.mix(rho_in, rho_out)
    }
}
```

### Step 4: Periodic Pulay restart

An alternative to full adaptive mixing: perform Anderson extrapolation every k-th iteration, with simple linear mixing in between. This is more robust than continuous Pulay for metals:

```rust
pub struct PeriodicPulayMixer {
    anderson: AndersonMixer,
    period: usize,
    iteration: usize,
    beta_linear: f64,
}

impl DensityMixer for PeriodicPulayMixer {
    fn mix(&mut self, rho_in: &[Complex64], rho_out: &[Complex64]) -> Vec<Complex64> {
        self.iteration += 1;
        if self.iteration % self.period == 0 {
            // Pulay extrapolation step
            self.anderson.mix(rho_in, rho_out)
        } else {
            // Simple linear mixing
            let n = rho_in.len();
            (0..n).map(|i| {
                rho_in[i] + self.beta_linear * (rho_out[i] - rho_in[i])
            }).collect()
        }
    }
}
```

### Step 5: Input configuration

```toml
[scf]
mixing_method = "anderson"   # "anderson" | "broyden" | "periodic_pulay"
mixing_beta = 0.3            # initial beta (adaptive adjusts from here)
adaptive_beta = true         # enable adaptive beta adjustment
```

## Acceptance Criteria

1. **Broyden converges for metals:** On a metallic test case (Al FCC), Broyden converges where Anderson with the same beta diverges, or converges in fewer iterations.
2. **Anderson equivalence:** For insulators (Si), Broyden and Anderson produce the same converged energy (within 1e-6 eV) in similar iteration counts.
3. **Adaptive beta works:** On a metallic system, adaptive mixing automatically reduces beta when residuals grow and recovers when convergence stabilizes. Manual beta tuning should not be necessary.
4. **Periodic Pulay robustness:** For a difficult convergence case, periodic Pulay (period=5) converges where continuous Anderson diverges.
5. **Mixing trait abstraction:** All three mixers implement the `DensityMixer` trait. Switching between them requires no changes outside `mixing.rs` and the input parser.
6. **History reset:** `reset()` correctly clears state so that warm-restart from a checkpoint doesn't carry stale Jacobian information.
7. **Backward compatibility:** Default configuration (Anderson, beta=0.3, no adaptive) reproduces current behavior exactly.
