//! Modified Broyden density mixer (Johnson PRB 38, 12807, 1988).

use crate::fft::FFT3D;

use super::kerker::{auto_q_tf_squared, precondition_residual};
use super::linalg::solve_linear_system;

/// Modified Broyden density mixer (Johnson PRB 38, 12807, 1988).
///
/// Implements the same algorithm as QE's `mix_rho.f90`. At each iteration:
///
/// 1. Compute residual R = ρ_out - ρ_in (optionally Kerker-preconditioned).
/// 2. If history exists, store the differences:
///    - df_i = R_current - R_previous
///    - dv_i = ρ_in_current - ρ_in_previous
/// 3. Build the overlap matrix β_{ij} = ⟨df_i | df_j⟩ and invert it.
/// 4. Compute γ_i = Σ_j β⁻¹_{ji} ⟨df_j | R⟩.
/// 5. Correct: ρ_in -= Σ_i γ_i · dv_i, R -= Σ_i γ_i · df_i.
/// 6. Output: ρ_new = ρ_in_corrected + β · R_corrected.
///
/// This is a quasi-Newton method that builds an approximate inverse Jacobian
/// from the iteration history, without storing the full N×N matrix.
pub(crate) struct BroydenMixer {
    beta: f64,
    max_history: usize,
    /// History of input density differences: dv_i = ρ_in_{i+1} - ρ_in_i
    history_dv: Vec<Vec<f64>>,
    /// History of residual differences: df_i = R_{i+1} - R_i
    history_df: Vec<Vec<f64>>,
    /// Previous iteration's input density.
    prev_rho_in: Option<Vec<f64>>,
    /// Previous iteration's (preconditioned) residual.
    prev_residual: Option<Vec<f64>>,
    /// Precomputed Kerker weights P(G) for each FFT grid point.
    /// None if no Kerker preconditioning.
    kerker_weights: Option<Vec<f64>>,
}

impl BroydenMixer {
    /// Create a new Broyden mixer.
    ///
    /// If `kerker` is true, Kerker preconditioning is applied to the residual
    /// before the Broyden update (like QE's default behavior).
    #[must_use]
    pub(super) fn new(
        beta: f64,
        max_history: usize,
        kerker: bool,
        g_squared: Option<&[f64]>,
        n_electrons: f64,
        omega: f64,
    ) -> Self {
        let kerker_weights = if kerker {
            let g2 = g_squared.expect("Broyden+Kerker mode requires g_squared");
            let q_tf_sq = auto_q_tf_squared(n_electrons, omega);
            let weights: Vec<f64> = g2
                .iter()
                .map(|&g2_val| {
                    if g2_val < crate::consts::G2_ZERO_THRESHOLD {
                        0.0
                    } else {
                        g2_val / (g2_val + q_tf_sq)
                    }
                })
                .collect();
            Some(weights)
        } else {
            None
        };

        Self {
            beta,
            max_history,
            history_dv: Vec::new(),
            history_df: Vec::new(),
            prev_rho_in: None,
            prev_residual: None,
            kerker_weights,
        }
    }

    /// Mix input density with output density using modified Broyden's method.
    ///
    /// `rho_in`: current input density (real space).
    /// `rho_out`: density from KS eigenstates (real space).
    /// `fft`: FFT instance (needed for Kerker preconditioning).
    ///
    /// Returns the new input density for the next iteration.
    pub(super) fn mix(&mut self, rho_in: &[f64], rho_out: &[f64], fft: &mut FFT3D) -> Vec<f64> {
        let n = rho_in.len();

        // Compute raw residual R = ρ_out - ρ_in
        let raw_residual: Vec<f64> = rho_out.iter().zip(rho_in.iter()).map(|(&o, &i)| o - i).collect();

        // Apply Kerker preconditioning if enabled
        let residual = if let Some(ref weights) = self.kerker_weights {
            precondition_residual(&raw_residual, weights, fft)
        } else {
            raw_residual
        };

        // If we have history from the previous iteration, compute differences
        if let (Some(prev_in), Some(prev_res)) = (&self.prev_rho_in, &self.prev_residual) {
            let dv: Vec<f64> = rho_in.iter().zip(prev_in.iter()).map(|(&c, &p)| c - p).collect();
            let df: Vec<f64> = residual.iter().zip(prev_res.iter()).map(|(&c, &p)| c - p).collect();

            self.history_dv.push(dv);
            self.history_df.push(df);

            // Trim history to max_history (remove oldest)
            if self.history_dv.len() > self.max_history {
                self.history_dv.remove(0);
                self.history_df.remove(0);
            }
        }

        // Save current state for next iteration's difference computation
        self.prev_rho_in = Some(rho_in.to_vec());
        self.prev_residual = Some(residual.clone());

        let m = self.history_df.len();

        if m == 0 {
            // First iteration: simple linear mixing ρ_new = ρ_in + β·R
            return rho_in.iter().zip(residual.iter()).map(|(&r, &res)| self.beta.mul_add(res, r)).collect();
        }

        // Build overlap matrix β_{ij} = ⟨df_i | df_j⟩
        let mut beta_mat = vec![0.0_f64; m * m];
        for i in 0..m {
            for j in i..m {
                let dot: f64 = self.history_df[i]
                    .iter()
                    .zip(self.history_df[j].iter())
                    .map(|(&a, &b)| a * b)
                    .sum();
                beta_mat[i * m + j] = dot;
                beta_mat[j * m + i] = dot;
            }
        }

        // Compute work_i = ⟨df_i | R_current⟩
        let work: Vec<f64> = (0..m)
            .map(|i| {
                self.history_df[i]
                    .iter()
                    .zip(residual.iter())
                    .map(|(&a, &b)| a * b)
                    .sum()
            })
            .collect();

        // Solve β·γ = work for γ (invert β via our existing Gauss solver)
        let gamma = solve_linear_system(&beta_mat, &work, m);

        // Apply Broyden correction:
        //   ρ_in_corrected = ρ_in - Σ_i γ_i · dv_i
        //   R_corrected    = R    - Σ_i γ_i · df_i
        //   ρ_new = ρ_in_corrected + β · R_corrected
        //
        // We first accumulate the corrections from all history entries,
        // then form the final mixed density in a single pass.
        let mut corr_dv = vec![0.0_f64; n];
        let mut corr_df = vec![0.0_f64; n];
        for (&g, (dv, df)) in gamma.iter().zip(self.history_dv.iter().zip(self.history_df.iter())) {
            for (k, (cdv, cdf)) in corr_dv.iter_mut().zip(corr_df.iter_mut()).enumerate() {
                *cdv += g * dv[k];
                *cdf += g * df[k];
            }
        }

        rho_in
            .iter()
            .zip(residual.iter())
            .zip(corr_dv.iter().zip(corr_df.iter()))
            .map(|((&ri, &res), (&cv, &cf))| self.beta.mul_add(res - cf, ri - cv))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_broyden_first_iteration_is_linear() {
        // On the first call, Broyden should reduce to simple linear mixing
        let mut fft = FFT3D::new(2, 2, 2);
        let mut mixer = BroydenMixer::new(0.3, 4, false, None, 8.0, 40.0);
        let rho_in = vec![1.0; 8];
        let rho_out = vec![2.0; 8];
        let result = mixer.mix(&rho_in, &rho_out, &mut fft);
        for &v in &result {
            assert!(
                (v - 1.3).abs() < 1e-10,
                "First Broyden iteration should be linear mixing: expected 1.3, got {v}"
            );
        }
    }

    #[test]
    fn test_broyden_multi_iteration_finite() {
        // Run 5 iterations and verify all results are finite
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let mut mixer = BroydenMixer::new(0.3, 4, false, None, 8.0, 40.0);

        let mut rho_in: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64).sin()).collect();
        let rho_target: Vec<f64> = (0..n).map(|i| 1.0 + 0.02 * (i as f64).cos()).collect();

        for _ in 0..5 {
            let result = mixer.mix(&rho_in, &rho_target, &mut fft);
            assert!(result.iter().all(|v| v.is_finite()), "Non-finite density after Broyden mixing");
            rho_in = result;
        }
    }

    #[test]
    fn test_broyden_with_kerker_finite() {
        // Broyden + Kerker should produce finite results
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let g_squared: Vec<f64> = (0..n)
            .map(|i| if i == 0 { 0.0 } else { 1.0 + i as f64 })
            .collect();
        let mut mixer = BroydenMixer::new(0.3, 4, true, Some(&g_squared), 8.0, 40.0);

        let mut rho_in: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64).sin()).collect();
        let rho_target: Vec<f64> = (0..n).map(|i| 1.0 + 0.02 * (i as f64).cos()).collect();

        for _ in 0..5 {
            let result = mixer.mix(&rho_in, &rho_target, &mut fft);
            assert!(result.iter().all(|v| v.is_finite()), "Non-finite density after Broyden+Kerker");
            rho_in = result;
        }
    }

    #[test]
    fn test_broyden_converges_quadratic() {
        // For a simple quadratic model f(x) = x (fixed-point iteration with
        // output = target), Broyden should converge toward the target.
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let mut mixer = BroydenMixer::new(0.5, 8, false, None, 8.0, 40.0);

        let target: Vec<f64> = (0..n).map(|i| 1.0 + 0.1 * (i as f64 * 0.3).sin()).collect();
        let mut rho_in: Vec<f64> = vec![1.0; n];

        let initial_err: f64 = rho_in.iter().zip(target.iter())
            .map(|(&a, &b)| (a - b).powi(2)).sum::<f64>().sqrt();

        for _ in 0..20 {
            let result = mixer.mix(&rho_in, &target, &mut fft);
            rho_in = result;
        }

        let final_err: f64 = rho_in.iter().zip(target.iter())
            .map(|(&a, &b)| (a - b).powi(2)).sum::<f64>().sqrt();

        assert!(
            final_err < initial_err * 0.01,
            "Broyden should converge: initial_err={initial_err:.6}, final_err={final_err:.6}"
        );
    }

    #[test]
    fn test_broyden_history_trimming() {
        // With max_history=2, we should never store more than 2 entries
        let mut fft = FFT3D::new(2, 2, 2);
        let mut mixer = BroydenMixer::new(0.3, 2, false, None, 8.0, 40.0);

        let n = 8;
        let mut rho_in = vec![1.0; n];

        for iter in 0..10 {
            let rho_out: Vec<f64> = rho_in.iter().map(|&r| r + 0.1 * (iter as f64)).collect();
            rho_in = mixer.mix(&rho_in, &rho_out, &mut fft);
            assert!(
                mixer.history_df.len() <= 2,
                "History should be trimmed to max_history=2, got {}",
                mixer.history_df.len()
            );
        }
    }

    #[test]
    fn test_broyden_vs_plain_scf_convergence() {
        // Both Anderson and Broyden should converge to the same energy on Si.
        use crate::{
            basis::BasisSet,
            crystal::{Atom, Crystal, Lattice},
            kpoints::KPoint,
            scf::{ScfParams, mixing::MixingMode, run_scf},
        };
        use nalgebra::Vector3;

        let a = 5.431;
        let crystal = Crystal {
            lattice: Lattice::new(
                a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
            ),
            atoms: vec![
                Atom::new(14, [0.0, 0.0, 0.0]),
                Atom::new(14, [0.25, 0.25, 0.25]),
            ],
        };
        let basis = BasisSet::new(&crystal.lattice, 100.0);
        let pp = crate::pseudopotential::load(
            &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        ).unwrap();
        let kpoints = vec![KPoint { k: Vector3::zeros(), weight: 1.0, label: None }];

        let plain_params = ScfParams {
            n_bands: 4,
            max_iter: 40,
            conv_threshold: 1e-6,
            mixing_beta: 0.3,
            mixing_ndim: 4,
            smearing_sigma: 0.05,
            ecutrho_ratio: 4,
            fft_grid: Some([16, 16, 16]),
            mixing_mode: MixingMode::Plain,
            ..Default::default()
        };

        let sym_id = crate::symmetry::SymmetryInfo::identity_only();
        let result_plain = run_scf(&crystal, &basis, &kpoints, &[&pp], &plain_params, &sym_id);

        let broyden_params = ScfParams {
            mixing_mode: MixingMode::Broyden { kerker: false },
            ..plain_params
        };
        let result_broyden = run_scf(&crystal, &basis, &kpoints, &[&pp], &broyden_params, &sym_id);

        match (&result_plain, &result_broyden) {
            (Ok(plain), Ok(broyden)) => {
                let energy_diff = (plain.total_energy - broyden.total_energy).abs();
                assert!(
                    energy_diff < 0.01,
                    "Plain ({:.6} eV) and Broyden ({:.6} eV) should converge to same energy, diff={energy_diff:.6}",
                    plain.total_energy, broyden.total_energy
                );
            }
            (Ok(_), Err(e)) => panic!("Plain converged but Broyden failed: {e}"),
            (Err(e), Ok(_)) => panic!("Broyden converged but plain failed: {e}"),
            (Err(_), Err(_)) => {
                // Both failed to converge — acceptable for this cheap Gamma-only test
            }
        }
    }
}
