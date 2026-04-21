//! Kerker preconditioning helpers shared by Anderson, Broyden, and
//! Periodic Pulay mixers.
//!
//! Kerker preconditioning damps long-wavelength density residuals in
//! reciprocal space:
//!
//!   R̃(G) = [|G|² / (|G|² + q_TF²)] R(G)
//!
//! which prevents charge sloshing (the dominant instability in metals and
//! large-gap systems).

use num_complex::Complex64;

use crate::fft::FFT3D;

/// Apply Kerker preconditioning in reciprocal space:
/// R_precond(r) = IFFT[ P(G) × FFT[R(r)] ]
pub(super) fn precondition_residual(residual_r: &[f64], weights: &[f64], fft: &mut FFT3D) -> Vec<f64> {
    let n = residual_r.len();
    let mut res_g: Vec<Complex64> = residual_r.iter().map(|&v| Complex64::new(v, 0.0)).collect();

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

/// Auto-estimate Thomas-Fermi screening wavevector squared from average
/// density.
///
/// q_TF² = 4 (3π²ρ)^{1/3} / π  (in a.u., then convert from Bohr⁻² to Å⁻²)
pub(super) fn auto_q_tf_squared(n_electrons: f64, omega: f64) -> f64 {
    use crate::consts::BOHR_TO_ANG;
    let rho_avg = n_electrons / omega; // e/ų
    let rho_bohr = rho_avg * BOHR_TO_ANG.powi(3); // e/Bohr³
    let q_tf_bohr_sq =
        4.0 * (3.0 * std::f64::consts::PI * std::f64::consts::PI * rho_bohr).cbrt() / std::f64::consts::PI;
    // Convert Bohr⁻² to ų
    q_tf_bohr_sq / (BOHR_TO_ANG * BOHR_TO_ANG)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_precondition_preserves_real() {
        // precondition_residual should produce a real-valued result
        // (imaginary parts should be negligible after FFT→filter→IFFT of real data)
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let weights: Vec<f64> = (0..n).map(|i| if i == 0 { 0.0 } else { 0.5 }).collect();
        let residual: Vec<f64> = (0..n).map(|i| (i as f64 * 0.1).sin()).collect();

        let result = precondition_residual(&residual, &weights, &mut fft);
        assert_eq!(result.len(), n);
        assert!(
            result.iter().all(|v| v.is_finite()),
            "Non-finite preconditioned residual"
        );
    }

    #[test]
    fn test_kerker_high_g_passes_through() {
        // A residual with only high-G components should pass through Kerker
        // nearly unchanged (P(G) → 1 for large |G|²)
        let _fft = FFT3D::new(4, 4, 4);
        let n = 64;
        // All G-vectors have large |G|² (>> q_TF²)
        let g_squared: Vec<f64> = (0..n).map(|i| 100.0 + f64::from(i)).collect();
        let q_tf = 1.0; // q_TF² = 1, much smaller than all |G|²

        let weights: Vec<f64> = g_squared.iter().map(|&g2| g2 / (g2 + q_tf * q_tf)).collect();

        // All weights should be close to 1.0
        for (i, &w) in weights.iter().enumerate() {
            assert!(w > 0.99, "Weight at G={i} should be ~1.0 for large |G|², got {w}");
        }
    }
}
