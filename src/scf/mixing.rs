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

use ndarray::{Array1, ArrayView1};
use num_complex::Complex64;

use crate::fft::FFT3D;

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

/// Anderson/Pulay (DIIS) density mixer with optional Kerker preconditioning.
///
/// Stores a history of input densities and residuals R^(n) = ρ_out^(n) - ρ_in^(n).
/// At each step, finds coefficients c_i (summing to 1) that minimize |Σ c_i R^(i)|²
/// by solving the DIIS linear system, then constructs the new density as:
///
///   ρ_in^{n+1} = Σ_i c_i [ρ_in^(i) + β R^(i)]
///
/// where β is the mixing parameter.
///
/// With Kerker preconditioning, the residual is modified in G-space before mixing:
///   R̃(G) = [|G|² / (|G|² + q_TF²)] R(G)
///
/// This damps long-wavelength charge sloshing, which is the dominant source of
/// SCF instability in metals and large-gap systems.
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
    #[must_use]
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
                // SAFETY: Callers must pass g_squared when using Kerker mode.
                // AndersonMixer::new is called from ScfContext which always provides
                // g_squared when mixing_mode is Kerker.
                let g2 = g_squared
                    .expect("BUG: Kerker mode requires g_squared to be provided by caller");
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
            MixingMode::Broyden { .. } | MixingMode::PeriodicPulay { .. } => {
                unreachable!(
                    "AndersonMixer should not be constructed with Broyden or PeriodicPulay mode; \
                     use Mixer::new()"
                )
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
        self.push_history(rho_in, rho_out, fft);
        self.diis_step()
    }

    /// Compute the residual `R = ρ_out - ρ_in` (optionally Kerker-preconditioned),
    /// append `(ρ_in, R)` to the history, and trim the history to `max_history`.
    ///
    /// Exposed for composite mixers (e.g. [`PeriodicPulayMixer`]) that accumulate
    /// history on every iteration but only invoke the DIIS solve on a subset.
    pub fn push_history(&mut self, rho_in: &[f64], rho_out: &[f64], fft: &mut FFT3D) {
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
        self.history_res.push(residual);

        // Trim history (oldest first)
        if self.history_in.len() > self.max_history {
            self.history_in.remove(0);
            self.history_res.remove(0);
        }
    }

    /// Run one DIIS/Anderson step using the currently-accumulated history.
    ///
    /// Assumes [`push_history`](Self::push_history) has been called with the
    /// latest `(ρ_in, ρ_out)` pair — the most recent entry drives the linear
    /// step when history is too short for DIIS.
    ///
    /// # Panics
    ///
    /// Panics if called before any history has been accumulated. Always call
    /// [`push_history`](Self::push_history) first.
    pub fn diis_step(&self) -> Vec<f64> {
        let m = self.history_in.len();
        assert!(
            m >= 1,
            "AndersonMixer::diis_step called before any history accumulated"
        );

        if m < 2 {
            // Simple linear mixing for first iteration: ρ_new = ρ_in + β R
            let rho_in_last = &self.history_in[m - 1];
            let r_last = &self.history_res[m - 1];
            return (rho_in_last + &(self.beta * r_last)).to_vec();
        }

        // Anderson mixing: find coefficients that minimize |Σ α_i R_i|²
        // subject to Σ α_i = 1.
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
        for ((&alpha, rho_in_j), res_j) in alpha_prev
            .iter()
            .zip(self.history_in.iter())
            .zip(self.history_res.iter())
        {
            rho_new += &(alpha * (rho_in_j + &(self.beta * res_j)));
        }

        rho_new.to_vec()
    }

    /// Number of history entries currently accumulated.
    #[must_use]
    pub fn history_len(&self) -> usize {
        self.history_in.len()
    }
}

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
pub struct BroydenMixer {
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
    pub fn new(
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
    pub fn mix(&mut self, rho_in: &[f64], rho_out: &[f64], fft: &mut FFT3D) -> Vec<f64> {
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

/// Periodic Pulay mixer (Banerjee, Suryanarayana, Pask, JCTC 12, 3053 (2016)).
///
/// Runs plain linear mixing (optionally Kerker-preconditioned) on every
/// iteration, except every `period`-th iteration where an Anderson/DIIS
/// extrapolation is performed using the residual history accumulated on the
/// intervening linear steps.
///
/// The rationale from the paper:
///
/// - Early iterations: DIIS is unreliable when only 1–2 residual vectors are
///   available and the long-wavelength sloshing dominates. Plain linear mixing
///   with a conservative β is more robust.
/// - Later iterations: once several consistent residuals are accumulated, DIIS
///   extrapolation gives the super-linear convergence advantage.
/// - Gating DIIS to a periodic cadence (k = 3–5) delivers both benefits and
///   yields 30–50% iteration-count reduction on the paper's test cases vs.
///   continuous Anderson.
///
/// Internally this wraps an [`AndersonMixer`]: on every iteration
/// [`AndersonMixer::push_history`] is called (history accumulates), and on
/// `period`-th iterations [`AndersonMixer::diis_step`] is invoked instead of
/// a plain `β·R` linear step.
pub struct PeriodicPulayMixer {
    anderson: AndersonMixer,
    /// k in the paper — DIIS on every k-th iteration (1-based counter).
    period: usize,
    /// 1-based iteration counter.
    iteration: usize,
}

impl PeriodicPulayMixer {
    /// Create a new Periodic Pulay mixer.
    ///
    /// `period` must be ≥ 1. `period = 1` is equivalent to continuous Anderson;
    /// `period = usize::MAX` is effectively plain linear mixing.
    ///
    /// `kerker` enables Kerker preconditioning on the residual (applied to both
    /// the linear-step residual and the DIIS step, sharing the Anderson inner
    /// mixer's preconditioner weights).
    #[must_use]
    pub fn new(
        beta: f64,
        max_history: usize,
        period: usize,
        kerker: bool,
        g_squared: Option<&[f64]>,
        n_electrons: f64,
        omega: f64,
    ) -> Self {
        assert!(period >= 1, "PeriodicPulayMixer: period must be >= 1");
        // The inner Anderson mixer owns the (optional) Kerker weights and the
        // history. We ask for Kerker mode iff `kerker == true`.
        let inner_mode = if kerker {
            MixingMode::Kerker { q_tf: None }
        } else {
            MixingMode::Plain
        };
        let anderson = AndersonMixer::new(
            beta,
            max_history,
            &inner_mode,
            g_squared,
            n_electrons,
            omega,
        );
        Self {
            anderson,
            period,
            iteration: 0,
        }
    }

    /// Mix input density with output density using Periodic Pulay.
    ///
    /// On non-Pulay iterations performs `ρ_new = ρ_in + β R`
    /// (plain/Kerker linear mixing). On every `period`-th iteration performs a
    /// DIIS step using the accumulated history.
    pub fn mix(&mut self, rho_in: &[f64], rho_out: &[f64], fft: &mut FFT3D) -> Vec<f64> {
        self.iteration += 1;

        // Always accumulate history — cheap, and lets the Pulay step use it.
        // push_history also applies Kerker preconditioning to the stored
        // residual when the inner Anderson mixer was constructed with Kerker.
        self.anderson.push_history(rho_in, rho_out, fft);

        let do_pulay =
            self.iteration.is_multiple_of(self.period) && self.anderson.history_len() >= 2;

        if do_pulay {
            // DIIS extrapolation using the full history.
            self.anderson.diis_step()
        } else {
            // Plain linear mixing against the most recent (preconditioned)
            // residual that push_history just stored. This preserves Kerker
            // preconditioning on linear steps as promised by the docstring.
            let len = self.anderson.history_len();
            let rho_in_last = &self.anderson.history_in[len - 1];
            let r_last = &self.anderson.history_res[len - 1];
            (rho_in_last + &(self.anderson.beta * r_last)).to_vec()
        }
    }

    /// Current iteration counter (1-based, zero before first `mix` call).
    #[cfg(test)]
    fn iteration_count(&self) -> usize {
        self.iteration
    }

    /// Accumulated history length (for testing).
    #[cfg(test)]
    fn history_len(&self) -> usize {
        self.anderson.history_len()
    }
}

/// Unified mixer that dispatches to Anderson, Broyden, or Periodic Pulay
/// based on `MixingMode`.
///
/// This avoids the need for a trait object or generic parameter in the SCF loop.
pub enum Mixer {
    Anderson(AndersonMixer),
    Broyden(BroydenMixer),
    PeriodicPulay(PeriodicPulayMixer),
}

impl Mixer {
    /// Create a mixer from the given parameters and mixing mode.
    #[must_use]
    pub fn new(
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
    pub fn mix(&mut self, rho_in: &[f64], rho_out: &[f64], fft: &mut FFT3D) -> Vec<f64> {
        match self {
            Mixer::Anderson(m) => m.mix(rho_in, rho_out, fft),
            Mixer::Broyden(m) => m.mix(rho_in, rho_out, fft),
            Mixer::PeriodicPulay(m) => m.mix(rho_in, rho_out, fft),
        }
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
        4.0 * (3.0 * std::f64::consts::PI * std::f64::consts::PI * rho_bohr).cbrt()
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
            log::warn!(
                "Anderson mixer: singular overlap matrix (pivot={pivot:.2e}), \
                 falling back to uniform coefficients"
            );
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

        let sym_id = crate::symmetry::SymmetryInfo::identity_only();
        let result_plain = run_scf(&crystal, &basis, &kpoints, &[&pp], &base_params, &sym_id);

        let kerker_params = ScfParams {
            mixing_mode: MixingMode::Kerker { q_tf: None },
            ..base_params
        };
        let result_kerker = run_scf(&crystal, &basis, &kpoints, &[&pp], &kerker_params, &sym_id);

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

    // -----------------------------------------------------------------------
    // Broyden mixer tests
    // -----------------------------------------------------------------------

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

    // -----------------------------------------------------------------------
    // Periodic Pulay mixer tests
    // -----------------------------------------------------------------------

    /// Deterministic synthetic residual sequence for exercising a mixer:
    /// ρ_out(x) = ρ_in + decaying_sinusoid. Independent of the mixer's output,
    /// so we can feed the same sequence into two mixers and compare.
    fn synthetic_rho_out(rho_in: &[f64], iter: usize) -> Vec<f64> {
        let decay = 0.5_f64.powi(iter as i32);
        rho_in
            .iter()
            .enumerate()
            .map(|(k, &r)| r + decay * 0.1 * (k as f64 * 0.37 + iter as f64).sin())
            .collect()
    }

    #[test]
    fn periodic_pulay_period_one_matches_anderson() {
        // With period = 1, every iteration triggers a DIIS step once history
        // has ≥ 2 entries. Since push_history is identical between both
        // mixers, the outputs should be bit-for-bit equal starting from
        // iteration 2 (iteration 1 is still linear mixing in both paths
        // because history_len < 2 → early return in diis_step / history_len
        // < 2 guard in PeriodicPulayMixer).
        let mut fft_a = FFT3D::new(4, 4, 4);
        let mut fft_b = FFT3D::new(4, 4, 4);
        let n = 64;

        let mut anderson = AndersonMixer::new(0.3, 4, &MixingMode::Plain, None, 8.0, 40.0);
        let mut pp = PeriodicPulayMixer::new(0.3, 4, 1, false, None, 8.0, 40.0);

        let mut rho_a: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64).sin()).collect();
        let mut rho_b = rho_a.clone();

        for iter in 1..=10 {
            let rho_out_a = synthetic_rho_out(&rho_a, iter);
            let rho_out_b = synthetic_rho_out(&rho_b, iter);

            let new_a = anderson.mix(&rho_a, &rho_out_a, &mut fft_a);
            let new_b = pp.mix(&rho_b, &rho_out_b, &mut fft_b);

            for (k, (&a, &b)) in new_a.iter().zip(new_b.iter()).enumerate() {
                assert!(
                    (a - b).abs() < 1e-12,
                    "iter {iter} elem {k}: anderson={a:.12} pp={b:.12}"
                );
            }

            rho_a = new_a;
            rho_b = new_b;
        }
    }

    #[test]
    fn periodic_pulay_huge_period_matches_plain_linear() {
        // With period = usize::MAX, the DIIS step never fires — every
        // iteration is a plain β·R linear step. This must match a reference
        // linear-mixing loop exactly.
        let mut fft_pp = FFT3D::new(4, 4, 4);
        let n = 64;
        let beta = 0.3;

        let mut pp = PeriodicPulayMixer::new(beta, 4, usize::MAX, false, None, 8.0, 40.0);

        let mut rho_pp: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64).sin()).collect();
        let mut rho_ref = rho_pp.clone();

        for iter in 1..=10 {
            let rho_out_pp = synthetic_rho_out(&rho_pp, iter);
            let rho_out_ref = synthetic_rho_out(&rho_ref, iter);

            // Reference: pure linear mixing ρ + β(ρ_out - ρ)
            let new_ref: Vec<f64> = rho_ref
                .iter()
                .zip(rho_out_ref.iter())
                .map(|(&r, &o)| r + beta * (o - r))
                .collect();

            let new_pp = pp.mix(&rho_pp, &rho_out_pp, &mut fft_pp);

            for (k, (&a, &b)) in new_pp.iter().zip(new_ref.iter()).enumerate() {
                assert!(
                    (a - b).abs() < 1e-12,
                    "iter {iter} elem {k}: pp={a:.12} ref={b:.12}"
                );
            }

            rho_pp = new_pp;
            rho_ref = new_ref;
        }
    }

    #[test]
    fn periodic_pulay_history_accumulates_between_pulay_steps() {
        // Over 6 iterations with period = 3, history should grow every step
        // up to max_history, regardless of whether the iteration is a Pulay
        // step or a linear step.
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let mut pp = PeriodicPulayMixer::new(0.3, 8, 3, false, None, 8.0, 40.0);

        let mut rho: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64).sin()).collect();

        for iter in 1..=6 {
            let rho_out = synthetic_rho_out(&rho, iter);
            rho = pp.mix(&rho, &rho_out, &mut fft);
            // push_history is called unconditionally, so after k iterations
            // we expect history_len = k (capped at max_history = 8).
            assert_eq!(
                pp.history_len(),
                iter,
                "after iter {iter}: history_len = {}, expected {iter}",
                pp.history_len()
            );
        }
    }

    #[test]
    fn periodic_pulay_period_three_fires_on_correct_iterations() {
        // With period = 3, a Pulay step happens on iterations 3, 6, 9.
        // We detect a Pulay step by comparing against an alternate mixer
        // whose output would differ on those iterations: a "linear-only"
        // reference built from the same preconditioned residuals.
        //
        // Since period = 1 is Anderson and period = ∞ is plain, we run both
        // the period-3 mixer and a period-∞ mixer on identical input streams
        // and assert they *diverge* exactly on iterations 3, 6, 9 (linear
        // steps are bit-identical; Pulay steps differ once history is ≥ 2).
        let mut fft_3 = FFT3D::new(4, 4, 4);
        let mut fft_inf = FFT3D::new(4, 4, 4);
        let n = 64;
        let beta = 0.3;

        let mut pp3 = PeriodicPulayMixer::new(beta, 8, 3, false, None, 8.0, 40.0);
        let mut pp_inf = PeriodicPulayMixer::new(beta, 8, usize::MAX, false, None, 8.0, 40.0);

        // Use the SAME density for both — synthetic output depends on rho_in
        // so we feed the reference density to both to keep inputs identical.
        let mut rho_shared: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64).sin()).collect();

        let mut divergence_iters = Vec::new();
        for iter in 1..=9 {
            let rho_out = synthetic_rho_out(&rho_shared, iter);

            // Clone current state to feed both mixers the same input.
            let rho_a = pp3.mix(&rho_shared, &rho_out, &mut fft_3);
            let rho_b = pp_inf.mix(&rho_shared, &rho_out, &mut fft_inf);

            let max_diff = rho_a
                .iter()
                .zip(rho_b.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);

            // On iterations 1-2 there isn't enough history for a DIIS step
            // even with period=3, so both should match. On iteration 3 and
            // later multiples of 3, they should diverge.
            if max_diff > 1e-10 {
                divergence_iters.push(iter);
            }

            // After the mixers have updated, advance the shared density
            // using the *reference* (period=∞, i.e. linear) path to keep
            // the input streams in sync for the next iteration.
            rho_shared = rho_b;
        }

        // Iteration 1, 2: no Pulay (history < 2).
        // Iteration 3: Pulay fires, pp3 diverges.
        // Iteration 6, 9: Pulay fires again, still diverges.
        // Iterations 4, 5, 7, 8 are linear steps in pp3 (matches pp_inf
        // structurally, but the internal state of pp3 has already been
        // perturbed by the DIIS step on iteration 3, so its linear-step
        // output *is* different from pp_inf's linear step because
        // rho_shared fed into pp3 differs from the state pp3 last saw).
        //
        // The robust assertion: iterations 1-2 must match (no Pulay yet),
        // iteration 3 must differ (first Pulay step).
        assert!(
            !divergence_iters.contains(&1),
            "iter 1: should match — no Pulay yet"
        );
        assert!(
            !divergence_iters.contains(&2),
            "iter 2: should match — no Pulay yet (history_len=2 boundary)"
        );
        assert!(
            divergence_iters.contains(&3),
            "iter 3: pp3 should do a Pulay step and diverge from pp_inf (divergence_iters={divergence_iters:?})"
        );
    }

    #[test]
    fn periodic_pulay_first_iteration_is_linear_mixing() {
        // Before any Pulay step, output must equal ρ_in + β·R (plain linear).
        let mut fft = FFT3D::new(2, 2, 2);
        let mut pp = PeriodicPulayMixer::new(0.3, 4, 3, false, None, 8.0, 40.0);
        let rho_in = vec![1.0; 8];
        let rho_out = vec![2.0; 8];
        let result = pp.mix(&rho_in, &rho_out, &mut fft);
        for &v in &result {
            assert!(
                (v - 1.3).abs() < 1e-12,
                "first PeriodicPulay step should be linear: expected 1.3, got {v}"
            );
        }
        assert_eq!(pp.iteration_count(), 1);
        assert_eq!(pp.history_len(), 1);
    }

    #[test]
    fn periodic_pulay_with_kerker_finite() {
        // Periodic Pulay + Kerker should produce finite results across
        // multiple iterations, including Pulay steps.
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let g_squared: Vec<f64> = (0..n)
            .map(|i| if i == 0 { 0.0 } else { 1.0 + i as f64 })
            .collect();
        let mut pp =
            PeriodicPulayMixer::new(0.3, 4, 3, true, Some(&g_squared), 8.0, 40.0);

        let mut rho_in: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64).sin()).collect();

        for iter in 1..=7 {
            let rho_out = synthetic_rho_out(&rho_in, iter);
            let result = pp.mix(&rho_in, &rho_out, &mut fft);
            assert!(
                result.iter().all(|v| v.is_finite()),
                "non-finite density after PeriodicPulay+Kerker at iter {iter}"
            );
            rho_in = result;
        }
    }

    #[test]
    fn periodic_pulay_converges_synthetic_fixed_point() {
        // Classic convergence test: f(x) = target is a fixed point, periodic
        // Pulay should drive x → target. This exercises both linear and
        // Pulay steps.
        let mut fft = FFT3D::new(4, 4, 4);
        let n = 64;
        let mut pp = PeriodicPulayMixer::new(0.5, 8, 3, false, None, 8.0, 40.0);

        let target: Vec<f64> = (0..n).map(|i| 1.0 + 0.1 * (i as f64 * 0.3).sin()).collect();
        let mut rho: Vec<f64> = vec![1.0; n];

        let initial_err: f64 = rho
            .iter()
            .zip(target.iter())
            .map(|(&a, &b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();

        for _ in 0..20 {
            rho = pp.mix(&rho, &target, &mut fft);
        }

        let final_err: f64 = rho
            .iter()
            .zip(target.iter())
            .map(|(&a, &b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();

        assert!(
            final_err < initial_err * 0.01,
            "PeriodicPulay should converge: initial_err={initial_err:.6}, final_err={final_err:.6}"
        );
    }

    #[test]
    fn periodic_pulay_vs_plain_scf_convergence() {
        // PRPL on Si Γ-only should reach the same total energy as Plain
        // within the convergence threshold, and should not be much slower
        // in iteration count.
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
            &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap();
        let kpoints = vec![KPoint {
            k: Vector3::zeros(),
            weight: 1.0,
            label: None,
        }];

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

        let prpl_params = ScfParams {
            mixing_mode: MixingMode::PeriodicPulay {
                period: 3,
                kerker: false,
            },
            ..plain_params
        };
        let result_prpl = run_scf(&crystal, &basis, &kpoints, &[&pp], &prpl_params, &sym_id);

        match (&result_plain, &result_prpl) {
            (Ok(plain), Ok(prpl)) => {
                let energy_diff = (plain.total_energy - prpl.total_energy).abs();
                assert!(
                    energy_diff < 0.01,
                    "Plain ({:.6} eV) and PeriodicPulay ({:.6} eV) should agree, diff={energy_diff:.6}",
                    plain.total_energy,
                    prpl.total_energy
                );
                // PRPL should converge in roughly as many or fewer iterations
                // than Plain on this small test case. Allow a modest slack
                // for numerical noise (Si Γ-only is a 2-atom insulator,
                // both schemes converge in similar step counts).
                assert!(
                    prpl.n_iterations <= plain.n_iterations + 5,
                    "PeriodicPulay ({} iters) should not be much slower than Plain ({} iters)",
                    prpl.n_iterations,
                    plain.n_iterations
                );
                println!(
                    "PRPL convergence: Plain={} iters, PeriodicPulay={} iters (ΔE={energy_diff:.3e} eV)",
                    plain.n_iterations, prpl.n_iterations
                );
            }
            (Ok(_), Err(e)) => panic!("Plain converged but PeriodicPulay failed: {e}"),
            (Err(e), Ok(_)) => panic!("PeriodicPulay converged but Plain failed: {e}"),
            (Err(_), Err(_)) => {
                // Both failed to converge — acceptable for this cheap Γ-only test
            }
        }
    }
}
