//! Density mixing schemes for SCF convergence.
//!
//! Four mixing algorithms are available:
//!
//! - **Anderson (Pulay/DIIS):** finds the optimal linear combination of past
//!   residuals to minimize the residual norm. Default method.
//!
//! - **Modified Broyden (Johnson PRB 38, 12807):** builds an approximate inverse
//!   Jacobian from the history of density residuals. Same algorithm as QE's
//!   `mix_rho.f90` and VASP's IMIX=4. Often converges faster for difficult
//!   systems (metals, large cells, charge sloshing).
//!
//! - **Periodic Pulay (Banerjee, Suryanarayana, Pask, JCTC 12, 3053 (2016)):**
//!   plain linear mixing on every iteration *except* every k-th, where
//!   Anderson/DIIS extrapolation is performed using the accumulated history.
//!   Avoids divergence on early iterations (when history is too short for DIIS
//!   to be reliable) while still getting the acceleration on later iterations.
//!   Paper reports 30–50% iteration-count reduction on transition-metal-oxide
//!   cases versus continuous Anderson.
//!
//! - **Kerker preconditioning:** can be combined with Anderson, Broyden, or
//!   Periodic Pulay. Damps long-wavelength density residuals to prevent charge
//!   sloshing: P(G) = |G|² / (|G|² + q_TF²).
//!
//! Use [`Mixer`] as the unified interface — it dispatches to the right algorithm
//! based on [`MixingMode`].
//!
//! ## Adaptive mixing β
//!
//! All three mixer algorithms support an optional **residual-norm monitor**
//! (Eyert, J. Comp. Phys. 124, 271 (1996), §3.3) that adjusts the mixing
//! parameter β between iterations based on the ratio of successive residual
//! norms:
//!
//! - If `‖R_i‖ / ‖R_{i-1}‖ > growth_threshold` (residual grew): damp
//!   `β ← max(β · damp_factor, β_min)`.
//! - If the ratio is `< restore_threshold` for `restore_window` consecutive
//!   iterations (steady convergence): restore
//!   `β ← min(β / damp_factor, β_start)` toward the user-configured start.
//! - Otherwise: β unchanged (hysteresis band).
//!
//! Defaults: `growth_threshold = 1.2`, `damp_factor = 0.7`,
//! `restore_threshold = 0.5`, `restore_window = 3`,
//! `β_min = max(0.05 · β_start, 0.01)`.
//!
//! The feature is gated on the `adaptive_beta` flag threaded through
//! `Mixer::new`. When `false` (the default for backwards compatibility) the
//! mixers behave bit-for-bit identically to the pre-MXBA implementation.

mod anderson;
mod broyden;
mod kerker;
mod linalg;

use crate::fft::FFT3D;

use anderson::{AndersonMixer, PeriodicPulayMixer};
use broyden::BroydenMixer;

/// Parameters needed to build Kerker-preconditioning weights.
///
/// Bundled so that mixer constructors don't grow their argument lists when
/// unrelated features (e.g. adaptive β) are added. These three quantities
/// always travel together — Kerker's `q_TF` is either user-provided or
/// auto-estimated from `(n_electrons, omega)`, and the weights are evaluated
/// on the FFT grid's `|G|²` array.
#[derive(Copy, Clone, Debug)]
pub(super) struct KerkerSetup<'a> {
    /// |G|² on the FFT grid; required when any Kerker variant is active.
    pub g_squared: Option<&'a [f64]>,
    /// Total valence electron count (for auto-q_TF estimation).
    pub n_electrons: f64,
    /// Unit-cell volume in Å³ (for auto-q_TF estimation).
    pub omega: f64,
}

/// Residual-norm monitor for adaptive mixing β.
///
/// Tracks the L2 norm of successive residuals and proposes a β update each
/// iteration per the Eyert (1996, §3.3) rule. The struct is held by both
/// [`AndersonMixer`] and [`BroydenMixer`]; when `enabled == false` the
/// monitor is a no-op (`update` returns the current β unchanged and touches
/// no internal state), which restores the pre-MXBA fixed-β behavior
/// bit-for-bit.
///
/// References:
/// - Eyert, *J. Comp. Phys.* **124**, 271 (1996), §3.3.
/// - Banerjee, Suryanarayana, Pask, *JCTC* **12**, 3053 (2016), §3.
#[derive(Clone, Debug)]
pub(super) struct AdaptiveBeta {
    /// If false, `update` returns β unchanged — no state is mutated.
    enabled: bool,
    /// User-configured starting β (upper clamp when restoring).
    beta_start: f64,
    /// Lower clamp — `max(0.05 · β_start, 0.01)`.
    beta_min: f64,
    /// Damp β when `‖R_i‖ / ‖R_{i-1}‖` exceeds this (>1 → residual grew).
    growth_threshold: f64,
    /// Restore β when `‖R_i‖ / ‖R_{i-1}‖` stays below this for
    /// `restore_window` iterations (strong monotone decrease).
    restore_threshold: f64,
    /// Multiplicative factor applied on growth (also 1/x when restoring).
    damp_factor: f64,
    /// Consecutive "good" iterations before restoring β.
    restore_window: usize,
    /// Previous residual norm; `None` on the first iteration.
    prev_norm: Option<f64>,
    /// Counter of consecutive iterations with strong residual decrease.
    restore_streak: usize,
}

impl AdaptiveBeta {
    /// Construct a monitor. Use `enabled = false` to disable (β never changes).
    #[must_use]
    pub(super) fn new(enabled: bool, beta_start: f64) -> Self {
        Self {
            enabled,
            beta_start,
            beta_min: (beta_start * 0.05).max(0.01),
            growth_threshold: 1.2,
            restore_threshold: 0.5,
            damp_factor: 0.7,
            restore_window: 3,
            prev_norm: None,
            restore_streak: 0,
        }
    }

    /// Update state with the current residual norm and return the new β.
    ///
    /// Called once per mixer iteration immediately after the (optionally
    /// preconditioned) residual has been computed. When `enabled == false`,
    /// returns `current_beta` unchanged and does not touch internal state.
    pub(super) fn update(&mut self, residual_norm: f64, current_beta: f64) -> f64 {
        if !self.enabled {
            return current_beta;
        }
        let prev = self.prev_norm.replace(residual_norm);
        let Some(prev) = prev else {
            // First iteration — no ratio yet, β unchanged.
            return current_beta;
        };
        // Degenerate case: non-positive prev norm or non-finite current norm.
        // Keep β, don't divide by zero.
        if prev <= 0.0 || !residual_norm.is_finite() {
            return current_beta;
        }
        let ratio = residual_norm / prev;

        if ratio > self.growth_threshold {
            // Residual grew — damp β toward β_min.
            self.restore_streak = 0;
            (current_beta * self.damp_factor).max(self.beta_min)
        } else if ratio < self.restore_threshold {
            self.restore_streak += 1;
            if self.restore_streak >= self.restore_window {
                // Sustained strong decrease — restore β toward β_start.
                self.restore_streak = 0;
                (current_beta / self.damp_factor).min(self.beta_start)
            } else {
                current_beta
            }
        } else {
            // Ratio in the hysteresis band [restore, growth]: do nothing.
            self.restore_streak = 0;
            current_beta
        }
    }
}

/// Mixing mode: plain Anderson, Kerker-preconditioned Anderson, Broyden, or
/// Periodic Pulay.
#[derive(Clone, Debug, Default)]
pub enum MixingMode {
    /// Standard Anderson mixing (no preconditioning).
    #[default]
    Plain,
    /// Kerker preconditioning with Thomas-Fermi screening wavevector q_TF (Å⁻¹).
    /// If q_TF is None, it is auto-estimated from the average electron density.
    Kerker { q_tf: Option<f64> },
    /// Modified Broyden mixing (Johnson PRB 38, 12807, 1988).
    ///
    /// Builds an approximate inverse Jacobian from the history of density
    /// residuals. Often converges faster and more robustly than Anderson for
    /// difficult systems (metals, large cells). This is the default in VASP
    /// (IMIX=4) and the algorithm used by QE's `mix_rho.f90`.
    ///
    /// Optionally combined with Kerker preconditioning.
    Broyden { kerker: bool },
    /// Periodic Pulay mixing (Banerjee et al., JCTC 12, 3053 (2016)).
    ///
    /// Plain linear mixing with `β` on every iteration except every
    /// `period`-th, where Anderson/DIIS extrapolation is performed using the
    /// accumulated history. Optional Kerker preconditioning is applied to the
    /// residual on both linear and DIIS steps.
    ///
    /// Period 3–5 is recommended for LDA/GGA on insulators and semiconductors;
    /// 5–8 for metals (per the paper). Default: 3.
    PeriodicPulay {
        /// k in the paper — do a DIIS step every k-th iteration.
        period: usize,
        /// Whether to apply Kerker preconditioning to the residual.
        kerker: bool,
    },
}

/// Unified mixer that dispatches to Anderson, Broyden, or Periodic Pulay
/// based on `MixingMode`.
///
/// This avoids the need for a trait object or generic parameter in the SCF loop.
pub(crate) enum Mixer {
    Anderson(AndersonMixer),
    Broyden(BroydenMixer),
    PeriodicPulay(PeriodicPulayMixer),
}

impl Mixer {
    /// Create a mixer from the given parameters and mixing mode.
    ///
    /// `adaptive_beta` enables the Eyert residual-norm monitor (see module
    /// docs). When `false`, β stays fixed at the supplied value for the whole
    /// run, which reproduces the pre-MXBA behavior bit-for-bit.
    #[must_use]
    pub(crate) fn new(
        beta: f64,
        max_history: usize,
        mode: &MixingMode,
        g_squared: Option<&[f64]>,
        n_electrons: f64,
        omega: f64,
        adaptive_beta: bool,
    ) -> Self {
        let kerker = KerkerSetup {
            g_squared,
            n_electrons,
            omega,
        };
        match mode {
            MixingMode::Plain | MixingMode::Kerker { .. } => Mixer::Anderson(AndersonMixer::new(
                beta,
                max_history,
                mode,
                kerker,
                adaptive_beta,
            )),
            MixingMode::Broyden { kerker: use_kerker } => Mixer::Broyden(BroydenMixer::new(
                beta,
                max_history,
                *use_kerker,
                kerker,
                adaptive_beta,
            )),
            MixingMode::PeriodicPulay {
                period,
                kerker: use_kerker,
            } => Mixer::PeriodicPulay(PeriodicPulayMixer::new(
                beta,
                max_history,
                *period,
                *use_kerker,
                kerker,
                adaptive_beta,
            )),
        }
    }

    /// Mix input and output densities to produce the next input density.
    pub(crate) fn mix(&mut self, rho_in: &[f64], rho_out: &[f64], fft: &mut FFT3D) -> Vec<f64> {
        match self {
            Mixer::Anderson(m) => m.mix(rho_in, rho_out, fft),
            Mixer::Broyden(m) => m.mix(rho_in, rho_out, fft),
            Mixer::PeriodicPulay(m) => m.mix(rho_in, rho_out, fft),
        }
    }

    /// Return the mixer's current effective β.
    ///
    /// With adaptive β disabled this is constant (equal to the value passed
    /// to [`Mixer::new`]); with adaptive β enabled it reflects the most
    /// recent update from the residual-norm monitor. Useful for logging
    /// at each SCF iteration.
    #[must_use]
    pub(crate) fn current_beta(&self) -> f64 {
        match self {
            Mixer::Anderson(m) => m.current_beta(),
            Mixer::Broyden(m) => m.current_beta(),
            Mixer::PeriodicPulay(m) => m.current_beta(),
        }
    }
}

#[cfg(test)]
mod adaptive_beta_tests {
    use super::*;

    #[test]
    fn disabled_monitor_never_changes_beta() {
        // With `enabled = false`, update() is a pure pass-through.
        let mut ab = AdaptiveBeta::new(false, 0.3);
        // Feed a wildly growing residual sequence that would normally trigger
        // heavy damping; with the monitor off β must stay at 0.3.
        let norms = [1.0, 10.0, 100.0, 1000.0];
        let mut beta = 0.3;
        for &n in &norms {
            beta = ab.update(n, beta);
            assert!(
                (beta - 0.3).abs() < 1e-15,
                "disabled monitor must not mutate β; got {beta}"
            );
        }
        assert!(ab.prev_norm.is_none(), "disabled monitor must not touch prev_norm");
    }

    #[test]
    fn growing_residual_damps_beta_monotonically() {
        // ‖R_i‖ / ‖R_{i-1}‖ = 1.5 > growth_threshold (1.2) every step.
        // β should be multiplied by damp_factor (0.7) on every update
        // until clamped to β_min = max(0.05·β_start, 0.01) = 0.05·0.5 = 0.025.
        let beta_start = 0.5_f64;
        let mut ab = AdaptiveBeta::new(true, beta_start);
        let beta_min = (beta_start * 0.05).max(0.01);
        let mut beta = beta_start;

        // Seed prev_norm = 1.0 on the first call (no damping yet).
        beta = ab.update(1.0, beta);
        assert!((beta - beta_start).abs() < 1e-15, "first iter should be no-op");

        let mut residual = 1.0_f64;
        let mut history = vec![beta];
        for _ in 0..15 {
            residual *= 1.5;
            beta = ab.update(residual, beta);
            history.push(beta);
        }
        // Every β should be <= its predecessor (monotone non-increasing).
        for w in history.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-15,
                "β not monotonically decreasing under growth: {} -> {}",
                w[0],
                w[1]
            );
        }
        // Eventually clamped at β_min.
        assert!(
            (*history.last().unwrap() - beta_min).abs() < 1e-12,
            "β did not clamp to β_min={beta_min}; final={}",
            history.last().unwrap()
        );
    }

    #[test]
    fn shrinking_residual_restores_beta_toward_start() {
        // Start from a damped β and feed strongly decreasing residuals.
        // After `restore_window = 3` consecutive "good" iterations, β should
        // be multiplied by 1/damp_factor (≈1.4286).
        let beta_start = 0.5_f64;
        let mut ab = AdaptiveBeta::new(true, beta_start);
        let mut beta = 0.05_f64;
        beta = ab.update(1.0, beta); // seed (no ratio yet)
        assert!((beta - 0.05).abs() < 1e-15, "first update must be no-op");
        beta = ab.update(0.1, beta); // ratio 0.10 — good, streak=1
        assert!((beta - 0.05).abs() < 1e-15, "no restore yet (streak=1)");
        beta = ab.update(0.01, beta); // ratio 0.10 — good, streak=2
        assert!((beta - 0.05).abs() < 1e-15, "no restore yet (streak=2)");
        beta = ab.update(0.001, beta); // ratio 0.10 — good, streak=3 → fire
        let expected = (0.05_f64 / 0.7).min(beta_start);
        assert!(
            (beta - expected).abs() < 1e-12,
            "restore should have fired at streak=3: expected {expected}, got {beta}"
        );
    }

    #[test]
    fn hysteresis_band_holds_beta_steady() {
        // Ratio in [restore_threshold, growth_threshold] = [0.5, 1.2]: β unchanged.
        let mut ab = AdaptiveBeta::new(true, 0.5);
        let mut beta = 0.3_f64;
        beta = ab.update(1.0, beta); // seed
        // Ratios: 0.9, 1.0, 0.8 (all inside the band).
        for next in [0.9, 0.9, 0.72] {
            beta = ab.update(next, beta);
            assert!(
                (beta - 0.3).abs() < 1e-15,
                "β should be unchanged in hysteresis band, got {beta}"
            );
        }
    }

    #[test]
    fn restore_cannot_exceed_beta_start() {
        // Repeatedly feed great ratios — β should clamp at β_start, not blow up.
        let beta_start = 0.5_f64;
        let mut ab = AdaptiveBeta::new(true, beta_start);
        let mut beta = beta_start;
        beta = ab.update(1.0, beta); // seed
        for k in 0..30 {
            let norm = 0.1_f64.powi(k + 1);
            beta = ab.update(norm, beta);
            assert!(
                beta <= beta_start + 1e-15,
                "β exceeded β_start ceiling: {beta} > {beta_start}"
            );
        }
    }

    #[test]
    fn growth_then_restore_cycle() {
        // Simulate: growth, growth, decrease, decrease, decrease → β should
        // drop twice, then restore once after 3 consecutive good iters.
        let beta_start = 0.4_f64;
        let damp = 0.7_f64;
        let mut ab = AdaptiveBeta::new(true, beta_start);
        let mut beta = beta_start;

        beta = ab.update(1.0, beta);     // seed
        beta = ab.update(2.0, beta);     // ratio 2.0 > 1.2 → damp
        assert!((beta - (beta_start * damp)).abs() < 1e-12);
        beta = ab.update(4.0, beta);     // ratio 2.0 > 1.2 → damp again
        let after_two_damps = beta_start * damp * damp;
        assert!((beta - after_two_damps).abs() < 1e-12);

        // Now three strong decreases.
        beta = ab.update(0.4, beta);     // ratio 0.1, streak=1
        assert!((beta - after_two_damps).abs() < 1e-12);
        beta = ab.update(0.04, beta);    // ratio 0.1, streak=2
        assert!((beta - after_two_damps).abs() < 1e-12);
        beta = ab.update(0.004, beta);   // ratio 0.1, streak=3 → restore
        let expected = (after_two_damps / damp).min(beta_start);
        assert!(
            (beta - expected).abs() < 1e-12,
            "expected restore to {expected}, got {beta}"
        );
    }

    #[test]
    fn mixer_current_beta_matches_construction_without_adaptive() {
        // Smoke test: Mixer::current_beta returns the user-configured β when
        // adaptive_beta is off, for every mixer variant.
        let mut fft = FFT3D::new(2, 2, 2);
        let modes = [
            MixingMode::Plain,
            MixingMode::Broyden { kerker: false },
            MixingMode::PeriodicPulay {
                period: 3,
                kerker: false,
            },
        ];
        for mode in &modes {
            let mut mixer = Mixer::new(0.42, 4, mode, None, 8.0, 40.0, false);
            assert!(
                (mixer.current_beta() - 0.42).abs() < 1e-15,
                "pre-mix current_beta mismatch for mode {mode:?}"
            );
            let _ = mixer.mix(&[1.0; 8], &[2.0; 8], &mut fft);
            assert!(
                (mixer.current_beta() - 0.42).abs() < 1e-15,
                "post-mix β changed with adaptive off for mode {mode:?}"
            );
        }
    }

    #[test]
    fn anderson_adaptive_backward_compat_bit_identical() {
        // With adaptive_beta = false the mixer's output must be bit-identical
        // to the pre-MXBA code path (which used a fixed self.beta). This is
        // the regression guard for backward compatibility.
        let mut fft_off = FFT3D::new(4, 4, 4);
        let mut fft_on_then_off = FFT3D::new(4, 4, 4);
        let n = 64;

        // mixer_ref: adaptive off, driven through several iterations.
        let mut ref_mixer = Mixer::new(0.3, 4, &MixingMode::Plain, None, 8.0, 40.0, false);
        let mut twin_mixer = Mixer::new(0.3, 4, &MixingMode::Plain, None, 8.0, 40.0, false);
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
            let new_ref = ref_mixer.mix(&rho_ref, &rho_out_ref, &mut fft_off);
            let new_twin = twin_mixer.mix(&rho_twin, &rho_out_twin, &mut fft_on_then_off);
            for (k, (&a, &b)) in new_ref.iter().zip(new_twin.iter()).enumerate() {
                assert!(
                    (a - b).abs() < 1e-15,
                    "iter {iter} elem {k}: backward-compat diverged"
                );
            }
            rho_ref = new_ref;
            rho_twin = new_twin;
        }
        // Both mixers should report the same constant β.
        assert!((ref_mixer.current_beta() - 0.3).abs() < 1e-15);
        assert!((twin_mixer.current_beta() - 0.3).abs() < 1e-15);
    }
}
