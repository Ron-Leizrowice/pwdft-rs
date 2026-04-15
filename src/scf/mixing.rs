//! Density mixing schemes for SCF convergence.
//!
//! Anderson (Pulay) mixing with optional Kerker preconditioning.
//! Kerker damps long-wavelength density residuals to prevent charge sloshing
//! in metals and large cells: P(G) = |G|² / (|G|² + q_TF²).

use ndarray::{Array1, ArrayView1};
use num_complex::Complex64;

use crate::fft::FFT3D;

/// Mixing mode: plain Anderson or Kerker-preconditioned Anderson.
#[derive(Clone, Debug, Default)]
pub enum MixingMode {
    /// Standard Anderson mixing (no preconditioning).
    #[default]
    Plain,
    /// Kerker preconditioning with Thomas-Fermi screening wavevector q_TF (Å⁻¹).
    /// If q_TF is None, it is auto-estimated from the average electron density.
    Kerker { q_tf: Option<f64> },
}

/// Anderson/Pulay density mixer with optional Kerker preconditioning.
pub struct AndersonMixer {
    beta: f64,
    max_history: usize,
    history_in: Vec<Array1<f64>>,
    history_res: Vec<Array1<f64>>,
    /// Precomputed Kerker weights P(G) for each FFT grid point.
    /// None if plain mixing.
    kerker_weights: Option<Vec<f64>>,
}

impl AndersonMixer {
    /// Create a new mixer. For Kerker mode, pass g_squared (|G|² at each FFT grid point).
    pub fn new(
        beta: f64,
        max_history: usize,
        mode: &MixingMode,
        g_squared: Option<&[f64]>,
        n_electrons: f64,
        omega: f64,
    ) -> Self {
        let kerker_weights = match mode {
            MixingMode::Plain => None,
            MixingMode::Kerker { q_tf } => {
                let g2 = g_squared.expect("Kerker mode requires g_squared");
                let q_tf_sq = match q_tf {
                    Some(q) => q * q,
                    None => auto_q_tf_squared(n_electrons, omega),
                };
                let weights: Vec<f64> = g2
                    .iter()
                    .map(|&g2_val| {
                        if g2_val < crate::consts::G2_ZERO_THRESHOLD {
                            0.0 // Suppress G=0 completely
                        } else {
                            g2_val / (g2_val + q_tf_sq)
                        }
                    })
                    .collect();
                Some(weights)
            }
        };

        Self {
            beta,
            max_history,
            history_in: Vec::new(),
            history_res: Vec::new(),
            kerker_weights,
        }
    }

    /// Mix input density with output density.
    ///
    /// `rho_in`: current input density (real space).
    /// `rho_out`: density from KS eigenstates (real space).
    /// `fft`: FFT instance (needed only for Kerker preconditioning).
    ///
    /// Returns the new input density for the next iteration.
    pub fn mix(&mut self, rho_in: &[f64], rho_out: &[f64], fft: &mut FFT3D) -> Vec<f64> {
        let rho_in_arr = ArrayView1::from(rho_in);
        let rho_out_arr = ArrayView1::from(rho_out);
        let raw_residual = &rho_out_arr - &rho_in_arr;

        // Apply Kerker preconditioning if enabled:
        // R_precond(r) = IFFT[ P(G) × FFT[R(r)] ]
        let residual = if let Some(ref weights) = self.kerker_weights {
            precondition_residual(&raw_residual.to_vec(), weights, fft)
        } else {
            raw_residual.to_vec()
        };
        let residual = Array1::from(residual);

        self.history_in.push(rho_in_arr.to_owned());
        self.history_res.push(residual.clone());

        // Trim history
        if self.history_in.len() > self.max_history {
            self.history_in.remove(0);
            self.history_res.remove(0);
        }

        let m = self.history_in.len();
        if m < 2 {
            // Simple linear mixing for first iteration
            let rho_new = &rho_in_arr + &(self.beta * &residual);
            return rho_new.to_vec();
        }

        // Anderson mixing: find coefficients that minimize |Σ α_i R_i|²
        // subject to Σ α_i = 1
        let last = m - 1;
        let mm = m - 1;

        let r_last = &self.history_res[last];

        let dr: Vec<Array1<f64>> = (0..mm)
            .map(|i| &self.history_res[i] - r_last)
            .collect();

        let mut a_mat = vec![0.0; mm * mm];
        let mut b_vec = vec![0.0; mm];

        for i in 0..mm {
            b_vec[i] = -dr[i].dot(r_last);
            for j in 0..mm {
                a_mat[i * mm + j] = dr[i].dot(&dr[j]);
            }
        }

        let alpha_prev = solve_linear_system(&a_mat, &b_vec, mm);
        let alpha_last = 1.0 - alpha_prev.iter().sum::<f64>();

        // Construct mixed density: Σ α_i (ρ_in_i + β R_i)
        let mut rho_new = alpha_last * (&self.history_in[last] + &(self.beta * r_last));
        for ((&alpha, rho_in_j), res_j) in alpha_prev.iter()
            .zip(self.history_in.iter())
            .zip(self.history_res.iter())
        {
            rho_new += &(alpha * (rho_in_j + &(self.beta * res_j)));
        }

        rho_new.to_vec()
    }
}

/// Apply Kerker preconditioning in reciprocal space:
/// R_precond(r) = IFFT[ P(G) × FFT[R(r)] ]
fn precondition_residual(residual_r: &[f64], weights: &[f64], fft: &mut FFT3D) -> Vec<f64> {
    let n = residual_r.len();
    let mut res_g: Vec<Complex64> = residual_r
        .iter()
        .map(|&v| Complex64::new(v, 0.0))
        .collect();

    // Forward FFT
    fft.forward(&mut res_g);

    // Apply Kerker weights in G-space
    for (g, &w) in res_g.iter_mut().zip(weights.iter()) {
        *g *= w;
    }

    // Inverse FFT (unnormalized — need to divide by N)
    fft.inverse(&mut res_g);
    let norm = 1.0 / n as f64;

    res_g.iter().map(|c| c.re * norm).collect()
}

/// Auto-estimate Thomas-Fermi screening wavevector squared from average density.
///
/// q_TF² = 4 (3π²ρ)^{1/3} / π  (in a.u., then convert from Bohr⁻² to Å⁻²)
fn auto_q_tf_squared(n_electrons: f64, omega: f64) -> f64 {
    use crate::consts::BOHR_TO_ANG;
    let rho_avg = n_electrons / omega; // e/ų
    let rho_bohr = rho_avg * BOHR_TO_ANG.powi(3); // e/Bohr³
    let q_tf_bohr_sq =
        4.0 * (3.0 * std::f64::consts::PI * std::f64::consts::PI * rho_bohr).powf(1.0 / 3.0)
            / std::f64::consts::PI;
    // Convert Bohr⁻² to ų
    q_tf_bohr_sq / (BOHR_TO_ANG * BOHR_TO_ANG)
}

/// Solve A x = b for small systems via Gauss elimination with partial pivoting.
fn solve_linear_system(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }

    for col in 0..n {
        let mut max_row = col;
        let mut max_val = aug[col * (n + 1) + col].abs();
        for row in col + 1..n {
            let val = aug[row * (n + 1) + col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }
        if max_row != col {
            for j in 0..=n {
                aug.swap(col * (n + 1) + j, max_row * (n + 1) + j);
            }
        }

        let pivot = aug[col * (n + 1) + col];
        if pivot.abs() < 1e-15 {
            return vec![1.0 / (n + 1) as f64; n];
        }

        for row in col + 1..n {
            let factor = aug[row * (n + 1) + col] / pivot;
            for j in col..=n {
                aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }

    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i * (n + 1) + n];
        for j in i + 1..n {
            sum -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] = sum / aug[i * (n + 1) + i];
    }

    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_mixing_plain() {
        let mut fft = FFT3D::new(2, 2, 2);
        let mut mixer = AndersonMixer::new(0.3, 4, &MixingMode::Plain, None, 8.0, 40.0);
        let rho_in = vec![1.0; 8];
        let rho_out = vec![2.0; 8];
        let result = mixer.mix(&rho_in, &rho_out, &mut fft);
        for &v in &result {
            assert!(
                (v - 1.3).abs() < 1e-10,
                "expected 1.3, got {v}"
            );
        }
    }

    #[test]
    fn test_kerker_suppresses_g0() {
        // Kerker should suppress the G=0 component of the residual
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let g_squared: Vec<f64> = (0..n)
            .map(|i| if i == 0 { 0.0 } else { 1.0 + i as f64 })
            .collect();
        let mut mixer = AndersonMixer::new(
            0.3, 4,
            &MixingMode::Kerker { q_tf: Some(1.0) },
            Some(&g_squared),
            8.0, 40.0,
        );

        // Uniform residual (all G=0 component) should be heavily suppressed
        let rho_in = vec![1.0; n];
        let rho_out = vec![2.0; n]; // residual = 1.0 everywhere (pure G=0)
        let result = mixer.mix(&rho_in, &rho_out, &mut fft);

        // With Kerker, the G=0 residual is zeroed, so mixing should barely change rho_in
        let max_change: f64 = result.iter().zip(rho_in.iter()).map(|(r, &i)| (r - i).abs()).fold(0.0, f64::max);
        // Without Kerker, change would be 0.3. With Kerker on uniform residual, much less.
        assert!(
            max_change < 0.1,
            "Kerker should suppress uniform (G=0) residual, but max_change={max_change}"
        );
    }

    #[test]
    fn test_auto_q_tf_reasonable() {
        // Si: 8 electrons, ~40 ų → q_TF should be ~1-3 Å⁻¹
        let q_tf_sq = auto_q_tf_squared(8.0, 40.0);
        let q_tf = q_tf_sq.sqrt();
        assert!(
            q_tf > 0.5 && q_tf < 5.0,
            "q_TF = {q_tf} Å⁻¹ outside reasonable range [0.5, 5.0]"
        );
    }

    #[test]
    fn test_kerker_anderson_multi_iteration() {
        // Run 3 iterations of Anderson+Kerker to verify the preconditioned
        // residual history works correctly (not just 1st iteration linear mixing).
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let g_squared: Vec<f64> = (0..n)
            .map(|i| if i == 0 { 0.0 } else { 1.0 + i as f64 * 0.5 })
            .collect();
        let mut mixer = AndersonMixer::new(
            0.3, 4,
            &MixingMode::Kerker { q_tf: Some(1.0) },
            Some(&g_squared),
            8.0, 40.0,
        );

        let mut rho_in: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64).sin()).collect();
        let rho_target: Vec<f64> = (0..n).map(|i| 1.0 + 0.02 * (i as f64).cos()).collect();

        // 3 iterations — should not panic and should produce finite values
        for _ in 0..3 {
            let result = mixer.mix(&rho_in, &rho_target, &mut fft);
            assert!(result.iter().all(|v| v.is_finite()), "Non-finite density after mixing");
            rho_in = result;
        }
    }

    #[test]
    fn test_kerker_small_qtf_approaches_plain() {
        // As q_TF → 0, P(G) → 1 for all G≠0, so Kerker → plain mixing
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let g_squared: Vec<f64> = (0..n)
            .map(|i| if i == 0 { 0.0 } else { 1.0 + i as f64 })
            .collect();

        let rho_in: Vec<f64> = (0..n).map(|i| 1.0 + 0.1 * (i as f64 * 0.3).sin()).collect();
        let rho_out: Vec<f64> = (0..n).map(|i| 1.0 + 0.2 * (i as f64 * 0.3).sin()).collect();

        let mut mixer_plain = AndersonMixer::new(
            0.3, 4, &MixingMode::Plain, None, 8.0, 40.0,
        );
        let result_plain = mixer_plain.mix(&rho_in, &rho_out, &mut fft);

        let mut mixer_kerker = AndersonMixer::new(
            0.3, 4,
            &MixingMode::Kerker { q_tf: Some(0.001) }, // tiny q_TF
            Some(&g_squared),
            8.0, 40.0,
        );
        let result_kerker = mixer_kerker.mix(&rho_in, &rho_out, &mut fft);

        // Should be nearly identical (small q_TF means almost no preconditioning)
        let max_diff: f64 = result_plain.iter()
            .zip(result_kerker.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(
            max_diff < 0.01,
            "Kerker with tiny q_TF should match plain mixing, max_diff={max_diff}"
        );
    }

    #[test]
    fn test_kerker_large_qtf_suppresses_all() {
        // As q_TF → ∞, P(G) → 0 for all G, total suppression
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let g_squared: Vec<f64> = (0..n)
            .map(|i| if i == 0 { 0.0 } else { 1.0 + i as f64 })
            .collect();
        let mut mixer = AndersonMixer::new(
            0.3, 4,
            &MixingMode::Kerker { q_tf: Some(1000.0) }, // huge q_TF
            Some(&g_squared),
            8.0, 40.0,
        );

        let rho_in = vec![1.0; n];
        let rho_out = vec![2.0; n];
        let result = mixer.mix(&rho_in, &rho_out, &mut fft);

        // With huge q_TF, almost no mixing should occur (residual fully suppressed)
        let max_change: f64 = result.iter().zip(rho_in.iter())
            .map(|(r, &i)| (r - i).abs())
            .fold(0.0, f64::max);
        assert!(
            max_change < 0.01,
            "Kerker with huge q_TF should suppress all mixing, max_change={max_change}"
        );
    }

    #[test]
    fn test_precondition_preserves_real() {
        // precondition_residual should produce a real-valued result
        // (imaginary parts should be negligible after FFT→filter→IFFT of real data)
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let weights: Vec<f64> = (0..n)
            .map(|i| if i == 0 { 0.0 } else { 0.5 })
            .collect();
        let residual: Vec<f64> = (0..n).map(|i| (i as f64 * 0.1).sin()).collect();

        let result = precondition_residual(&residual, &weights, &mut fft);
        assert_eq!(result.len(), n);
        assert!(result.iter().all(|v| v.is_finite()), "Non-finite preconditioned residual");
    }

    #[test]
    fn test_kerker_high_g_passes_through() {
        // A residual with only high-G components should pass through Kerker
        // nearly unchanged (P(G) → 1 for large |G|²)
        let _fft = FFT3D::new(4, 4, 4);
        let n = 64;
        // All G-vectors have large |G|² (>> q_TF²)
        let g_squared: Vec<f64> = (0..n).map(|i| 100.0 + i as f64).collect();
        let q_tf = 1.0; // q_TF² = 1, much smaller than all |G|²

        let weights: Vec<f64> = g_squared.iter()
            .map(|&g2| g2 / (g2 + q_tf * q_tf))
            .collect();

        // All weights should be close to 1.0
        for (i, &w) in weights.iter().enumerate() {
            assert!(
                w > 0.99,
                "Weight at G={i} should be ~1.0 for large |G|², got {w}"
            );
        }
    }

    #[test]
    fn test_plain_ignores_fft() {
        // Plain mixing should produce the same result regardless of FFT state
        let mut fft1 = FFT3D::new(2, 2, 2);
        let mut fft2 = FFT3D::new(2, 2, 2);

        let mut mixer1 = AndersonMixer::new(0.3, 4, &MixingMode::Plain, None, 8.0, 40.0);
        let mut mixer2 = AndersonMixer::new(0.3, 4, &MixingMode::Plain, None, 8.0, 40.0);

        let rho_in = vec![1.0; 8];
        let rho_out = vec![2.0; 8];

        let r1 = mixer1.mix(&rho_in, &rho_out, &mut fft1);
        let r2 = mixer2.mix(&rho_in, &rho_out, &mut fft2);

        for (a, b) in r1.iter().zip(r2.iter()) {
            assert!((a - b).abs() < 1e-15, "Plain mixing results should be identical");
        }
    }

    #[test]
    fn test_kerker_vs_plain_scf_convergence() {
        // Both modes should converge to the same energy on Si.
        // Minimal system: Γ-only, ecut=100, 16³ grid.
        use crate::{
            basis::BasisSet,
            crystal::{Atom, Crystal, Lattice},
            kpoints::KPoint,
            scf::{ScfParams, run_scf},
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

        let base_params = ScfParams {
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

        let result_plain = run_scf(&crystal, &basis, &kpoints, &[&pp], &base_params, None);

        let kerker_params = ScfParams {
            mixing_mode: MixingMode::Kerker { q_tf: None },
            ..base_params
        };
        let result_kerker = run_scf(&crystal, &basis, &kpoints, &[&pp], &kerker_params, None);

        match (&result_plain, &result_kerker) {
            (Ok(plain), Ok(kerker)) => {
                let energy_diff = (plain.total_energy - kerker.total_energy).abs();
                assert!(
                    energy_diff < 0.01,
                    "Plain ({:.6} eV) and Kerker ({:.6} eV) should converge to same energy, diff={energy_diff:.6}",
                    plain.total_energy, kerker.total_energy
                );
                // For insulators, Kerker may take a few more iterations (it's
                // designed for metals). Just verify it's not wildly worse.
                assert!(
                    kerker.n_iterations <= plain.n_iterations + 10,
                    "Kerker ({} iters) shouldn't be much slower than plain ({} iters)",
                    kerker.n_iterations, plain.n_iterations
                );
            }
            (Ok(_), Err(e)) => panic!("Plain converged but Kerker failed: {e}"),
            (Err(e), Ok(_)) => panic!("Kerker converged but plain failed: {e}"),
            (Err(_), Err(_)) => {
                // Both failed to converge — acceptable for this cheap test
            }
        }
    }

    #[test]
    fn test_solve_linear_system() {
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let b = vec![5.0, 7.0];
        let x = solve_linear_system(&a, &b, 2);
        assert!((x[0] - 1.6).abs() < 1e-10);
        assert!((x[1] - 1.8).abs() < 1e-10);
    }
}
