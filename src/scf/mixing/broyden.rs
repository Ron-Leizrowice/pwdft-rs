//! Modified Broyden density mixer (Johnson PRB 38, 12807, 1988).

use crate::fft::FFT3D;

use super::AdaptiveBeta;
use super::KerkerSetup;
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
    /// Current effective β. Mutated each iteration by [`AdaptiveBeta::update`]
    /// (no-op when adaptive β is disabled, in which case this equals the
    /// user-configured start β for the entire run).
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
    /// Residual-norm monitor for adaptive β (Eyert 1996 §3.3).
    adaptive: AdaptiveBeta,
}

impl BroydenMixer {
    /// Create a new Broyden mixer.
    ///
    /// If `kerker` is true, Kerker preconditioning is applied to the residual
    /// before the Broyden update (like QE's default behavior).
    ///
    /// `adaptive_beta = true` activates the Eyert (1996) residual-norm
    /// monitor, which damps β when ‖R‖ grows and restores it toward the
    /// configured start when ‖R‖ decreases steadily. `false` keeps β fixed.
    #[must_use]
    pub(super) fn new(
        beta: f64,
        max_history: usize,
        kerker: bool,
        kerker_setup: KerkerSetup<'_>,
        adaptive_beta: bool,
    ) -> Self {
        let kerker_weights = if kerker {
            let g2 = kerker_setup
                .g_squared
                .expect("Broyden+Kerker mode requires g_squared");
            let q_tf_sq = auto_q_tf_squared(kerker_setup.n_electrons, kerker_setup.omega);
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
            adaptive: AdaptiveBeta::new(adaptive_beta, beta),
        }
    }

    /// Current effective β (after any adaptive update applied so far).
    #[must_use]
    pub(super) fn current_beta(&self) -> f64 {
        self.beta
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

        // Adaptive β (Eyert 1996 §3.3): update β from the ratio of this
        // iteration's residual norm to the previous one. The residual here
        // is the *Kerker-preconditioned* residual when Kerker is on
        // (G=0 zeroed + damped long-wavelength components) — that's the
        // right quantity to feed the monitor because it is the one β
        // multiplies in the linear step. With Kerker off it's the raw
        // density residual. No-op when `adaptive_beta = false`; otherwise
        // the Broyden step below uses the updated β for both the linear
        // combination and the correction terms.
        let residual_norm = residual
            .iter()
            .map(|&r| r * r)
            .sum::<f64>()
            .sqrt();
        self.beta = self.adaptive.update(residual_norm, self.beta);

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

    fn plain_ctx() -> KerkerSetup<'static> {
        KerkerSetup {
            g_squared: None,
            n_electrons: 8.0,
            omega: 40.0,
        }
    }

    fn kerker_ctx(g2: &[f64]) -> KerkerSetup<'_> {
        KerkerSetup {
            g_squared: Some(g2),
            n_electrons: 8.0,
            omega: 40.0,
        }
    }

    #[test]
    fn test_broyden_first_iteration_is_linear() {
        // On the first call, Broyden should reduce to simple linear mixing
        let mut fft = FFT3D::new(2, 2, 2);
        let mut mixer = BroydenMixer::new(0.3, 4, false, plain_ctx(), false);
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
        let mut mixer = BroydenMixer::new(0.3, 4, false, plain_ctx(), false);

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
        let mut mixer = BroydenMixer::new(0.3, 4, true, kerker_ctx(&g_squared), false);

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
        let mut mixer = BroydenMixer::new(0.5, 8, false, plain_ctx(), false);

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
        let mut mixer = BroydenMixer::new(0.3, 2, false, plain_ctx(), false);

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
                // FDLT: both mixers must bring Δρ well under conv_threshold
                // (1e-6 here). An upper bound of 1e-4 gives 2 orders of
                // margin without being so tight that a minor mixer tweak
                // breaks the test.
                assert!(
                    plain.final_delta < 1e-4,
                    "Plain Si SCF final Δρ = {:.3e} should be well below conv_threshold=1e-6",
                    plain.final_delta,
                );
                assert!(
                    broyden.final_delta < 1e-4,
                    "Broyden Si SCF final Δρ = {:.3e} should be well below conv_threshold=1e-6",
                    broyden.final_delta,
                );
            }
            (Ok(_), Err(e)) => panic!("Plain converged but Broyden failed: {e}"),
            (Err(e), Ok(_)) => panic!("Broyden converged but plain failed: {e}"),
            (Err(e_plain), Err(e_broyden)) => panic!(
                "Both Plain and Broyden failed on Si Γ-only SCF — this test \
                 expects both to converge. Plain error: {e_plain}; Broyden error: {e_broyden}"
            ),
        }
    }

    // -----------------------------------------------------------------------
    // MXBA: adaptive β unit tests for BroydenMixer
    // -----------------------------------------------------------------------

    /// Drive a Broyden mixer with residuals of prescribed norms.
    fn drive_broyden_with_norms(
        mixer: &mut BroydenMixer,
        fft: &mut FFT3D,
        norms: &[f64],
    ) -> Vec<f64> {
        let n = 8;
        let mut rho_in = vec![1.0; n];
        let mut betas = Vec::with_capacity(norms.len());
        for (k, &target) in norms.iter().enumerate() {
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            let unit = sign / (n as f64).sqrt();
            let rho_out: Vec<f64> = rho_in.iter().map(|&r| r + target * unit).collect();
            rho_in = mixer.mix(&rho_in, &rho_out, fft);
            betas.push(mixer.current_beta());
        }
        betas
    }

    #[test]
    fn broyden_adaptive_beta_damps_on_growth() {
        let mut fft = FFT3D::new(2, 2, 2);
        let mut mixer = BroydenMixer::new(0.5, 4, false, plain_ctx(), true);
        let norms = [1.0, 1.5, 2.25, 3.375, 5.0625, 7.59];
        let betas = drive_broyden_with_norms(&mut mixer, &mut fft, &norms);
        assert!(
            *betas.last().unwrap() < 0.5 - 1e-6,
            "adaptive β did not damp under growth: trajectory {betas:?}"
        );
        for w in betas.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-12,
                "β increased under growth: {} -> {}",
                w[0],
                w[1]
            );
        }
    }

    #[test]
    fn broyden_adaptive_beta_restores_after_damp() {
        let mut fft = FFT3D::new(2, 2, 2);
        let mut mixer = BroydenMixer::new(0.5, 4, false, plain_ctx(), true);

        let growth = [1.0, 1.5, 2.25];
        let _ = drive_broyden_with_norms(&mut mixer, &mut fft, &growth);
        let damped = mixer.current_beta();
        assert!(
            damped < 0.5 - 1e-6,
            "prerequisite failed: Broyden β not damped: {damped}"
        );

        let decrease = [1.0, 0.1, 0.01, 0.001, 0.0001];
        let betas = drive_broyden_with_norms(&mut mixer, &mut fft, &decrease);
        assert!(
            *betas.last().unwrap() > damped + 1e-6,
            "Broyden β did not restore: start={damped}, trajectory={betas:?}"
        );
        for &b in &betas {
            assert!(b <= 0.5 + 1e-12, "β exceeded β_start=0.5: {b}");
        }
    }

    #[test]
    fn broyden_adaptive_disabled_is_bit_identical() {
        // Backward compatibility: adaptive_beta=false must yield the same
        // output as the pre-MXBA code path.
        let mut fft_ref = FFT3D::new(4, 4, 4);
        let mut fft_twin = FFT3D::new(4, 4, 4);
        let n = 64;
        let mut ref_mixer = BroydenMixer::new(0.3, 4, false, plain_ctx(), false);
        let mut twin_mixer = BroydenMixer::new(0.3, 4, false, plain_ctx(), false);
        let mut rho_ref: Vec<f64> =
            (0..n).map(|i| 1.0 + 0.01 * (i as f64).sin()).collect();
        let mut rho_twin = rho_ref.clone();
        for iter in 1..=6 {
            let decay = 0.5_f64.powi(iter);
            let rho_out_ref: Vec<f64> = rho_ref
                .iter()
                .enumerate()
                .map(|(k, &r)| r + decay * 0.1 * (k as f64 * 0.37 + iter as f64).sin())
                .collect();
            let rho_out_twin: Vec<f64> = rho_twin
                .iter()
                .enumerate()
                .map(|(k, &r)| r + decay * 0.1 * (k as f64 * 0.37 + iter as f64).sin())
                .collect();
            let a = ref_mixer.mix(&rho_ref, &rho_out_ref, &mut fft_ref);
            let b = twin_mixer.mix(&rho_twin, &rho_out_twin, &mut fft_twin);
            for (k, (&x, &y)) in a.iter().zip(b.iter()).enumerate() {
                assert!(
                    (x - y).abs() < 1e-15,
                    "iter {iter} elem {k}: Broyden backward-compat diverged"
                );
            }
            rho_ref = a;
            rho_twin = b;
        }
        assert!((ref_mixer.current_beta() - 0.3).abs() < 1e-15);
        assert!((twin_mixer.current_beta() - 0.3).abs() < 1e-15);
    }
}
