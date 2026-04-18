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

mod anderson;
mod broyden;
mod kerker;
mod linalg;

use crate::fft::FFT3D;

use anderson::{AndersonMixer, PeriodicPulayMixer};
use broyden::BroydenMixer;

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
    #[must_use]
    pub(crate) fn new(
        beta: f64,
        max_history: usize,
        mode: &MixingMode,
        g_squared: Option<&[f64]>,
        n_electrons: f64,
        omega: f64,
    ) -> Self {
        match mode {
            MixingMode::Plain | MixingMode::Kerker { .. } => {
                Mixer::Anderson(AndersonMixer::new(
                    beta, max_history, mode, g_squared, n_electrons, omega,
                ))
            }
            MixingMode::Broyden { kerker } => {
                Mixer::Broyden(BroydenMixer::new(
                    beta, max_history, *kerker, g_squared, n_electrons, omega,
                ))
            }
            MixingMode::PeriodicPulay { period, kerker } => {
                Mixer::PeriodicPulay(PeriodicPulayMixer::new(
                    beta,
                    max_history,
                    *period,
                    *kerker,
                    g_squared,
                    n_electrons,
                    omega,
                ))
            }
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
}
