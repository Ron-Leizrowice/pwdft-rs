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
//! Use `Mixer` as the unified interface — it dispatches to the right algorithm
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
/// no internal state).
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
            let new_beta = (current_beta * self.damp_factor).max(self.beta_min);
            if (new_beta - current_beta).abs() > f64::EPSILON * current_beta.max(1.0) {
                // debug! (not info!) because under adverse conditions this
                // fires every iteration; the per-iteration β line already
                // surfaces the running value.
                log::debug!(
                    "AdaptiveBeta: damp β {current_beta:.4} → {new_beta:.4} \
                     (residual ratio {ratio:.3} > growth {growth:.3})",
                    growth = self.growth_threshold,
                );
            }
            new_beta
        } else if ratio < self.restore_threshold {
            self.restore_streak += 1;
            if self.restore_streak >= self.restore_window {
                // Sustained strong decrease — restore β toward β_start.
                let window = self.restore_streak;
                self.restore_streak = 0;
                let new_beta = (current_beta / self.damp_factor).min(self.beta_start);
                if (new_beta - current_beta).abs() > f64::EPSILON * current_beta.max(1.0) {
                    log::debug!(
                        "AdaptiveBeta: restore β {current_beta:.4} → {new_beta:.4} \
                         (after {window} consecutive ratios < restore {restore:.3})",
                        restore = self.restore_threshold,
                    );
                }
                new_beta
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
    /// Standard Anderson / Pulay DIIS mixing (no preconditioning).
    ///
    /// **Known limitation:** plain Anderson stalls on wide-gap insulators
    /// (C diamond, polar oxides) with residual plateaus at
    /// Δρ ≈ 10⁻⁵ – 10⁻⁶ that do not decay within a normal `max_iter`
    /// budget. Combine with [`MixingMode::Kerker`] preconditioning or
    /// switch to [`MixingMode::Broyden`] / [`MixingMode::PeriodicPulay`]
    /// for those systems. Pinned by a mixer-robustness regression test
    /// on C diamond (see `tests/mixer_robustness.rs`).
    #[default]
    Plain,
    /// Kerker preconditioning with Thomas-Fermi screening wavevector q_TF (Å⁻¹).
    /// If q_TF is None, it is auto-estimated from the average electron density.
    Kerker { q_tf: Option<f64> },
    /// Modified Broyden mixing (Johnson PRB 38, 12807, 1988).
    ///
    /// Builds an approximate inverse Jacobian from the history of density
    /// residuals. Often converges faster and more robustly than Anderson for
    /// difficult systems (metals, large cells).
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
    /// docs). When `false`, β stays fixed at the supplied value for the
    /// whole run.
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
    pub(crate) fn current_beta(&self) -> f64 {
        match self {
            Mixer::Anderson(m) => m.current_beta(),
            Mixer::Broyden(m) => m.current_beta(),
            Mixer::PeriodicPulay(m) => m.current_beta(),
        }
    }

    /// The chosen Thomas-Fermi q_TF and how it was picked, when Kerker
    /// preconditioning is active. Plain mixers return `None`.
    fn kerker_q_tf(&self) -> Option<anderson::KerkerQtf> {
        match self {
            Mixer::Anderson(m) => m.kerker_q_tf,
            Mixer::Broyden(m) => m.kerker_q_tf,
            Mixer::PeriodicPulay(m) => m.kerker_q_tf(),
        }
    }

    /// Emit a single `info!` line describing this mixer's active
    /// configuration: algorithm, β, history depth, adaptive-β flag, and
    /// Kerker status (with q_TF). Call once, at SCF start, after the mixer
    /// is constructed. Callers may pass a short `tag` (e.g. `"charge"` /
    /// `"magnetization"`) to distinguish the two mixers of the spin driver
    /// in the log; empty string for a single-channel run.
    pub(crate) fn log_init(
        &self,
        mode: &MixingMode,
        max_history: usize,
        adaptive: bool,
        tag: &str,
    ) {
        let prefix = if tag.is_empty() {
            "Mixer".to_string()
        } else {
            format!("Mixer[{tag}]")
        };
        log::info!(
            "{prefix}: {algo}  β={beta:.3}  history={hist}  adaptive_β={adaptive_flag}  \
             kerker={kerker}",
            algo = mode_label(mode),
            beta = self.current_beta(),
            hist = max_history,
            adaptive_flag = if adaptive { "on" } else { "off" },
            kerker = kerker_summary(self.kerker_q_tf(), mode_wants_kerker(mode)),
        );
    }
}

/// One-line `info!` label for a mixing mode. Owned String so PeriodicPulay
/// can carry its `period`; short &'static str round-trip via `Cow` would
/// complicate the format call site for no real win.
fn mode_label(mode: &MixingMode) -> String {
    match mode {
        MixingMode::Plain => "Anderson".to_string(),
        MixingMode::Kerker { .. } => "Anderson + Kerker".to_string(),
        MixingMode::Broyden { kerker: false } => "Broyden".to_string(),
        MixingMode::Broyden { kerker: true } => "Broyden + Kerker".to_string(),
        MixingMode::PeriodicPulay { period, kerker: false } => {
            format!("PeriodicPulay(k={period})")
        }
        MixingMode::PeriodicPulay { period, kerker: true } => {
            format!("PeriodicPulay(k={period}) + Kerker")
        }
    }
}

/// Does this mode request Kerker preconditioning?
fn mode_wants_kerker(mode: &MixingMode) -> bool {
    matches!(
        mode,
        MixingMode::Kerker { .. }
            | MixingMode::Broyden { kerker: true }
            | MixingMode::PeriodicPulay { kerker: true, .. }
    )
}

/// Short summary of the Kerker status for [`Mixer::log_init`]. When Kerker
/// was requested by the mode but the mixer reports no q_TF, we still log
/// `on (q_TF=?)` instead of silently falling back to "off" — that would
/// hide a configuration bug.
fn kerker_summary(q_tf: Option<anderson::KerkerQtf>, wants: bool) -> String {
    match (q_tf, wants) {
        (None, false) => "off".to_string(),
        (None, true) => "on (q_TF=?)".to_string(),
        (Some(k), _) => {
            let origin = if k.user_supplied { "user" } else { "auto" };
            format!("on (q_TF={origin}, {q:.3} Å⁻¹)", q = k.q_tf)
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
    fn adaptive_beta_fe_failure_trajectory_floors_to_beta_min() {
        // Cheap deterministic companion to `tests/mxba_adaptive_beta_fe.rs`.
        //
        // Approximates the Fe BCC CCMX residual trajectory that drove MXBA's
        // `adaptive_beta = false` default. See `tests/mxba_adaptive_beta_fe.rs`
        // for the real 80-iter SCF run and
        // `proposals/completed/MXBA-adaptive-mixing-beta.md` for the paper
        // reference (Eyert, J. Comp. Phys. 124, 271 (1996), §3.3).
        //
        // Synthetic sequence — 80 iterations with MXBA defaults (β_start=0.3,
        // β_min=0.015, growth=1.2, restore=0.5, damp=0.7, window=3):
        //
        //   iters 1–10: flat plateau at 0.34 (ratio ≈ 1.0, inside the
        //   hysteresis band [0.5, 1.2] → no β change). Mirrors the observed
        //   "monitor silent" phase of iters 1–9 in the SCF log.
        //
        //   iters 11–80: residual oscillates between 0.34·1.3 and 0.34·0.9,
        //   producing alternating ratios of ~1.444 (damp fires: β ← 0.7·β,
        //   clamped at β_min) and ~0.692 (inside band, streak resets to 0).
        //
        // Under this pattern the 3-iter streak of ratios < 0.5 that would
        // restore β is *unreachable*, so β monotonically ratchets down and
        // floors at β_min long before iter 80.
        //
        // When MXB2 lands a fix to the Eyert tuning (e.g. "require a 2-iter
        // growth streak before damping"), this test must update with it —
        // either (a) assert β recovers above β_min under the same sequence,
        // or (b) gain an `#[ignore]` marker matching
        // `tests/mxba_adaptive_beta_fe.rs` if the failure mode persists.
        let beta_start = 0.3_f64;
        let beta_min = (beta_start * 0.05).max(0.01); // 0.015 with these defaults
        let mut ab = AdaptiveBeta::new(true, beta_start);
        let mut beta = beta_start;

        // Sanity-check the monitor's configured band before driving it.
        assert!((ab.beta_min - beta_min).abs() < 1e-15);
        assert!((ab.growth_threshold - 1.2).abs() < 1e-15);
        assert!((ab.restore_threshold - 0.5).abs() < 1e-15);
        assert!((ab.damp_factor - 0.7).abs() < 1e-15);
        assert_eq!(ab.restore_window, 3);

        let base = 0.34_f64;
        let plateau_len = 10;
        let oscillation_multipliers = [1.3_f64, 0.9_f64];
        let total_iters: usize = 80;

        // Track the maximum run of consecutive sub-restore-threshold ratios so
        // we can pin the "restore cannot fire" property below.
        let mut max_good_streak: usize = 0;
        let mut cur_good_streak: usize = 0;
        let mut prev_residual: Option<f64> = None;

        for iter in 1..=total_iters {
            let residual = if iter <= plateau_len {
                base
            } else {
                let osc_idx = (iter - plateau_len - 1) % oscillation_multipliers.len();
                base * oscillation_multipliers[osc_idx]
            };
            if let Some(prev) = prev_residual {
                let ratio = residual / prev;
                if ratio < ab.restore_threshold {
                    cur_good_streak += 1;
                    max_good_streak = max_good_streak.max(cur_good_streak);
                } else {
                    cur_good_streak = 0;
                }
            }
            prev_residual = Some(residual);
            beta = ab.update(residual, beta);
        }

        // Core assertion: β has floored at β_min well before iter 80.
        assert!(
            (beta - beta_min).abs() <= 0.01 * beta_min,
            "β did not reach β_min within 80 iters of Fe-trajectory: β={beta}, β_min={beta_min}"
        );

        // Defense-in-depth: the documented failure mode requires that the
        // `restore_window = 3` streak of ratios < 0.5 is never reached, so β
        // has no way to recover during the SCF. If a future tuning change
        // relaxes `restore_threshold`, this assertion will flag that the
        // synthetic sequence also needs updating.
        assert!(
            max_good_streak < ab.restore_window,
            "synthetic sequence accidentally triggers a restore streak ({max_good_streak} >= \
             restore_window {}); rework multipliers to preserve the failure mode",
            ab.restore_window
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
            (0..n).map(|i| 1.0 + 0.01 * f64::from(i).sin()).collect();
        let mut rho_twin = rho_ref.clone();

        for iter in 1..=6 {
            let decay = 0.5_f64.powi(iter);
            let rho_out_ref: Vec<f64> = rho_ref
                .iter()
                .enumerate()
                .map(|(k, &r)| r + decay * 0.1 * (k as f64 * 0.37 + f64::from(iter)).sin())
                .collect();
            let rho_out_twin: Vec<f64> = rho_twin
                .iter()
                .enumerate()
                .map(|(k, &r)| r + decay * 0.1 * (k as f64 * 0.37 + f64::from(iter)).sin())
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
