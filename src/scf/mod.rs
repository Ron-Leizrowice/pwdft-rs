//! Self-consistent field (SCF) driver and supporting machinery.
//!
//! [`run_scf`] is the public entry point: it validates its
//! [`ScfParams`], dispatches to the non-spin or spin-polarized driver,
//! and returns a converged [`ScfResult`] with eigenvalues, density, and
//! [`EnergyComponents`]. The dispatch is thin — the real work lives in
//! the private `driver` / `driver_spin` modules.
//!
//! Supporting submodules:
//! - [`mixing`] — density mixing schemes (Anderson, Broyden, Kerker).
//! - [`smearing`] — finite-temperature occupation functions.
//! - [`initial_density`] — superposition-of-atomic-densities (SAD) start.
//! - [`density`] — density reconstruction from occupied orbitals.
//!
//! Internal helpers (context, energy, grid, potentials, report) are
//! crate-private and shared between the two drivers.

pub(crate) mod context;
pub mod density;
mod driver;
mod driver_spin;
pub(crate) mod energy;
pub(crate) mod grid;
pub mod initial_density;
pub mod mixing;
pub(crate) mod potentials;
mod report;
pub mod smearing;

use num_complex::Complex64;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    error::{PwdftError, Result},
    kpoints::KPoint,
    pseudopotential::PseudopotentialData,
};

pub use self::energy::EnergyComponents;
use self::grid::FftGrid;

/// Parameters controlling the self-consistent field iteration.
///
/// The SCF loop solves the Kohn-Sham equations iteratively:
/// 1. Construct `V_eff = V_local + V_Hartree[ρ_val] + V_xc[ρ_val + ρ_core]`
///    (the core charge ρ_core enters only V_xc, via the NLCC path; it is
///    *not* added to the Hartree source and *not* counted as valence).
/// 2. Diagonalize H = T + V_eff + V_NL at each k-point
/// 3. Compute occupations from eigenvalues (Fermi-Dirac or other smearing)
/// 4. Reconstruct density ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²
/// 5. Mix input and output densities (Anderson/Pulay) and repeat
///
/// Convergence requires both density (Δρ < conv_threshold) and
/// energy (ΔE < energy_threshold) criteria to be met.
///
/// Nonlinear core correction (NLCC; Louie, Froyen, Cohen, *Phys. Rev. B*
/// **26**, 1738 (1982)) is enabled automatically whenever the UPF file
/// has `core_correction="T"`. The invariants (ρ_core in XC only, not
/// Hartree; not counted as valence; split evenly between spin channels
/// in LSDA) are pinned by unit tests in `src/pseudopotential/upf/` and
/// by an integration regression guard for BCC Fe.
#[derive(Clone)]
pub struct ScfParams {
    pub n_bands: usize,
    pub max_iter: usize,
    /// Density convergence threshold (RMS, e/ų).
    pub conv_threshold: f64,
    /// Energy convergence threshold (eV). Both density AND energy must converge.
    pub energy_threshold: f64,
    pub mixing_beta: f64,
    pub mixing_ndim: usize,
    pub smearing_sigma: f64,
    /// Smearing scheme for occupation numbers.
    pub smearing_scheme: smearing::SmearingScheme,
    /// Charge density cutoff as multiple of wavefunction cutoff.
    pub ecutrho_ratio: u32,
    /// Explicit FFT grid dimensions. If set, overrides ecutrho_ratio.
    pub fft_grid: Option<[usize; 3]>,
    /// Mixing mode: plain Anderson or Kerker-preconditioned.
    pub mixing_mode: mixing::MixingMode,
    /// Enable adaptive mixing β (Eyert 1996, §3.3 residual-norm monitor).
    ///
    /// When `true`, β is damped when the residual norm grows and restored
    /// toward `mixing_beta` when it decreases steadily. When `false`
    /// (default), β stays fixed at `mixing_beta` for the entire run.
    /// See `src/scf/mixing/mod.rs` module docs for the full rule and
    /// defaults.
    pub adaptive_beta: bool,
    /// Number of spin channels: 1 (unpolarized) or 2 (collinear spin-polarized).
    pub nspin: usize,
    /// Starting magnetization per atom type (fractional, -1 to 1).
    /// Maps from element symbol to magnetization. Empty = non-magnetic.
    pub starting_magnetization: std::collections::HashMap<String, f64>,
    /// Fixed total magnetization (n_up - n_down) in electrons.
    /// If None, magnetization is determined self-consistently.
    pub tot_magnetization: Option<f64>,
    /// Eigensolver backend for per-k-point diagonalization.
    ///
    /// `Dense` (default): full `faer::SelfAdjointEigen` O(n³).
    /// `Iterative`: faer's implicitly-restarted Arnoldi partial solver,
    /// computing only the lowest `n_bands` eigenpairs. Falls back to dense
    /// when the problem is too small for Arnoldi to be profitable or when
    /// iteration fails to converge in the restart budget.
    pub eigensolver: crate::eigensolver::EigensolverKind,
    /// Enable subspace-diagonalization warm-start for the dense eigensolver.
    ///
    /// When `true` and `eigensolver == Dense`, the driver caches the
    /// previous iteration's eigenvectors per k-point and uses them as a
    /// Rayleigh–Ritz subspace for the next SCF iteration's diagonalization
    /// ([`crate::eigensolver::dense::diagonalize_subspace`]). The subspace
    /// path automatically falls back to the full dense diagonalization
    /// whenever the per-eigenpair residual gate
    /// ([`crate::eigensolver::dense::WFRX_RESIDUAL_TOL`]) is tripped, so
    /// correctness is never sacrificed. When `eigensolver == Iterative`,
    /// this flag has no effect (the iterative path has its own warm-start
    /// via `v0`).
    ///
    /// Default `false` (opt-in).
    pub wfrx_subspace: bool,
    /// Exchange-correlation functional selector.
    ///
    /// Only [`crate::settings::XcFunctional::Pz`] (Perdew-Zunger LDA) is
    /// actually implemented today. Any other variant causes `run_scf` to
    /// fail fast with [`PwdftError::NotImplemented`] so that a YAML typo
    /// cannot silently produce LDA results under a GGA or hybrid label.
    pub xc_functional: crate::settings::XcFunctional,
    /// Width (Å) of the Gaussian model atomic charge used by the initial
    /// Superposition-of-Atomic-Densities (SAD) guess when a pseudopotential
    /// lacks `PP_RHOATOM`.
    ///
    /// Default is 1.0 Å (see `scf::initial_density::DEFAULT_GAUSSIAN_SIGMA`
    /// for the canonical constant). The SCF refines the initial guess
    /// away in the first few iterations, so this value does not affect
    /// the converged density; it is a sensitivity-study knob — a wider
    /// sigma damps high-G content in the starting density and can slightly
    /// change iteration count for ill-conditioned systems.
    pub gaussian_sigma: f64,
}

impl ScfParams {
    /// Validate parameters before starting an SCF calculation.
    ///
    /// # Errors
    ///
    /// Returns [`PwdftError::InvalidParam`] when any of the following is
    /// true:
    /// - `n_bands == 0`
    /// - `conv_threshold <= 0`
    /// - `mixing_beta` outside the half-open interval `(0, 1]`
    /// - `smearing_sigma < 0`
    /// - `ecutrho_ratio < 1`
    /// - `nspin` is neither 1 nor 2
    /// - `mixing_mode` is `PeriodicPulay` with `period == 0`
    /// - `gaussian_sigma` is non-finite or non-positive
    pub fn validate(&self) -> Result<()> {
        if self.n_bands == 0 {
            return Err(PwdftError::InvalidParam {
                name: "n_bands",
                reason: "must be > 0".into(),
            });
        }
        if self.conv_threshold <= 0.0 {
            return Err(PwdftError::InvalidParam {
                name: "conv_threshold",
                reason: "must be positive".into(),
            });
        }
        if self.mixing_beta <= 0.0 || self.mixing_beta > 1.0 {
            return Err(PwdftError::InvalidParam {
                name: "mixing_beta",
                reason: format!("must be in (0, 1], got {}", self.mixing_beta),
            });
        }
        if self.smearing_sigma < 0.0 {
            return Err(PwdftError::InvalidParam {
                name: "smearing_sigma",
                reason: "must be non-negative".into(),
            });
        }
        if self.ecutrho_ratio < 1 {
            return Err(PwdftError::InvalidParam {
                name: "ecutrho_ratio",
                reason: format!("must be >= 1, got {}", self.ecutrho_ratio),
            });
        }
        if self.nspin != 1 && self.nspin != 2 {
            return Err(PwdftError::InvalidParam {
                name: "nspin",
                reason: format!("must be 1 or 2, got {}", self.nspin),
            });
        }
        // Periodic Pulay: period k must be >= 1, else PeriodicPulayMixer::new
        // would panic (`period >= 1` assertion). Surface this as InvalidParam
        // so YAML parsing / programmatic callers get a structured error.
        if let mixing::MixingMode::PeriodicPulay { period, .. } = self.mixing_mode
            && period == 0
        {
            return Err(PwdftError::InvalidParam {
                name: "pulay_period",
                reason: "must be >= 1 for PeriodicPulay mixing".into(),
            });
        }
        // Gaussian sigma for the SAD initial density (CFGN Phase 1) must be
        // positive and finite. A NaN or non-positive value would propagate
        // into `exp(-|G|^2 * sigma^2 / 2)` and poison the entire starting
        // density.
        if !self.gaussian_sigma.is_finite() || self.gaussian_sigma <= 0.0 {
            return Err(PwdftError::InvalidParam {
                name: "gaussian_sigma",
                reason: format!(
                    "must be positive and finite, got {}",
                    self.gaussian_sigma
                ),
            });
        }
        Ok(())
    }
}

impl Default for ScfParams {
    fn default() -> Self {
        Self {
            n_bands: 8,
            max_iter: 100,
            conv_threshold: 1e-6,
            energy_threshold: 1e-5,
            mixing_beta: 0.3,
            mixing_ndim: 8,
            smearing_sigma: 0.01,
            smearing_scheme: smearing::SmearingScheme::FermiDirac,
            ecutrho_ratio: 4,
            fft_grid: None,
            mixing_mode: mixing::MixingMode::Plain,
            adaptive_beta: false,
            nspin: 1,
            starting_magnetization: std::collections::HashMap::new(),
            tot_magnetization: None,
            eigensolver: crate::eigensolver::EigensolverKind::default(),
            wfrx_subspace: false,
            xc_functional: crate::settings::XcFunctional::default(),
            gaussian_sigma: initial_density::DEFAULT_GAUSSIAN_SIGMA,
        }
    }
}

/// Output of a converged SCF calculation.
///
/// All energies are in eV. The four energy quantities are:
/// - `total_energy`: F = E_band - E_H + E_xc - E_vxc + E_ewald + V_local(G=0)·N_el − TS
///   (the Mermin free-energy functional — matches QE's `! total energy`
///   line from `pw.x` output, which is also F = E − TS).
/// - `free_energy`: same quantity as `total_energy`; retained as an
///   explicit alias for backward-compatible callers.
/// - `harris_foulkes_energy`: HF estimator of the same Mermin F.
/// - `energy_sigma0`: E₀ = (E_internal + F)/2, a σ → 0 extrapolation.
///
/// Note: prior to TSEN (2026-04-19) `total_energy` omitted the −TS
/// contribution. The pre-TSEN internal energy `E = F + TS` is recoverable
/// as `total_energy + entropy_ts`.
#[derive(Debug)]
pub struct ScfResult {
    /// Mermin free energy `F = E − TS` (includes −TS; matches QE's
    /// `! total energy` line). Pre-TSEN internal energy is
    /// `total_energy + entropy_ts`.
    pub total_energy: f64,
    /// Harris-Foulkes estimator of the Mermin free energy.
    ///
    /// Uses the input density for Hartree/XC corrections but output eigenvalues.
    /// Stationary at self-consistency: `|E_HF − F_KS| → 0` quadratically.
    /// Serves as a convergence quality indicator.
    pub harris_foulkes_energy: f64,
    /// Alias for `total_energy`. Both fields carry the Mermin free energy
    /// `F = E − TS`; the alias is retained so callers migrating from the
    /// pre-TSEN API (where `total_energy` was the internal `E` and
    /// `free_energy` was `E − TS`) can continue to read `free_energy`
    /// and receive the physically correct quantity.
    pub free_energy: f64,
    /// σ → 0 extrapolated energy `E₀ = (E_internal + F)/2 = F + TS/2`.
    /// Best estimate of the T = 0 total energy at a fixed k-grid.
    pub energy_sigma0: f64,
    /// Positive entropy contribution `T · S` in eV. The `−TS` term
    /// folded into `total_energy` has the opposite sign.
    pub entropy_ts: f64,
    /// Eigenvalues indexed as `[spin_k_index][band]`.
    /// For nspin=1: length = n_kpoints. For nspin=2: length = 2 * n_kpoints
    /// (spin-up k-points followed by spin-down k-points).
    pub eigenvalues: Vec<Vec<f64>>,
    pub fermi_energy: f64,
    pub n_iterations: usize,
    /// The last Δρ the SCF saw.
    ///
    /// On successful convergence, this is the value that fell below
    /// `ScfParams::conv_threshold`. For nspin=2, it is the per-channel
    /// max(‖Δρ↑‖, ‖Δρ↓‖) per-spin convergence criterion (see
    /// `scf::driver_spin`).
    ///
    /// On `max_iter` exhaustion the driver returns
    /// `PwdftError::ConvergenceFailure { delta, .. }` instead of an
    /// `ScfResult`, so this field only ever carries a converged value.
    /// It is exposed primarily for pathology-specific regression guards
    /// — e.g. a Fe limit cycle pinned Δρ at ~0.254, a failure mode the
    /// `n_iterations < max_iter` guard alone cannot catch if `max_iter`
    /// is relaxed.
    pub final_delta: f64,
    pub rho_g: Vec<Complex64>,
    /// Total magnetization M = ∫(ρ_up - ρ_down)dr in μB (Bohr magnetons).
    /// Zero for nspin=1.
    pub magnetization: f64,
    /// Number of spin channels (1 or 2).
    pub nspin: usize,
    /// Per-term energy breakdown (kinetic, local, non-local, Hartree,
    /// XC, Ewald) plus the Harris-Foulkes stationary estimator.
    pub components: EnergyComponents,
}

/// Run the self-consistent field loop.
///
/// Validates parameters and dispatches to the non-spin (`nspin == 1`) or
/// spin-polarized (`nspin == 2`) driver. The two drivers live in
/// `scf::driver` and `scf::driver_spin`; this front door centralizes
/// argument validation.
///
/// # Errors
///
/// - `PwdftError::InvalidParam` if [`ScfParams::validate`] rejects the
///   params (see that method's own `# Errors` section for the per-field
///   variant map).
/// - `PwdftError::InvalidCrystal` if `crystal.atoms` is empty, if
///   `kpoints` is empty, or if the lattice volume is effectively zero
///   (< 1e-10 ų).
/// - `PwdftError::NotImplemented` if `params.xc_functional` is an
///   unsupported variant (anything other than Perdew-Zunger LDA today);
///   surfaced by `XcEvaluator::from_settings` so a YAML typo fails fast
///   before compute work starts.
/// - `PwdftError::ConvergenceFailure` from the selected driver when the
///   SCF loop exhausts `max_iter` without satisfying both density and
///   energy thresholds.
/// - Any error propagated from the driver (eigensolver failure, NaN
///   density, etc.) — see `scf::driver::run_scf_unpolarized` and
///   `scf::driver_spin::run_scf_spin`.
pub fn run_scf(
    crystal: &Crystal,
    basis: &BasisSet,
    kpoints: &[KPoint],
    pseudopotentials: &[&PseudopotentialData],
    params: &ScfParams,
    symmetry: &crate::symmetry::SymmetryInfo,
) -> Result<ScfResult> {
    params.validate()?;
    // GGAP Phase A dispatch: construct the XC evaluator up-front so that
    // unsupported functional labels (pbe0, hse06) fail fast at SCF entry,
    // before any compute work — the same XCNI guarantee. `Pbe` constructs
    // successfully here but its evaluation returns NotImplemented from
    // inside the driver; Phase B replaces that branch with real PBE.
    //
    // The evaluator is a *data* enum. This shape is load-bearing for HYBR
    // (PBE0/HSE06) because hybrid functionals need (ρ, ψ) access during
    // Hamiltonian construction, which a closure-shaped dispatch could not
    // reach. See `proposals/HYBR-hybrid-functional-support.md` §3.
    let xc_evaluator =
        crate::potential::xc::XcEvaluator::from_settings(params.xc_functional)?;
    if crystal.atoms.is_empty() {
        return Err(PwdftError::InvalidCrystal {
            reason: "at least one atom is required",
        });
    }
    if kpoints.is_empty() {
        return Err(PwdftError::InvalidCrystal {
            reason: "at least one k-point is required",
        });
    }
    let omega = crystal.lattice.volume();
    if omega < 1e-10 {
        return Err(PwdftError::InvalidCrystal {
            reason: "lattice has zero or near-zero volume",
        });
    }

    if params.nspin == 2 {
        driver_spin::run_scf_spin(
            crystal, basis, kpoints, pseudopotentials, params, symmetry, xc_evaluator,
        )
    } else {
        driver::run_scf_unpolarized(
            crystal, basis, kpoints, pseudopotentials, params, symmetry, xc_evaluator,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_rejects_zero_pulay_period() {
        // Nit from PR #39 review: ScfParams::validate() must reject
        // pulay_period == 0 with a structured error, rather than letting
        // the PeriodicPulayMixer::new `period >= 1` assertion panic at
        // SCF setup time. ERR2 P1.c migrated this path from the catch-all
        // `InvalidInput` to the structured `InvalidParam` variant.
        let params = ScfParams {
            mixing_mode: mixing::MixingMode::PeriodicPulay {
                period: 0,
                kerker: false,
            },
            ..Default::default()
        };
        let err = params.validate().expect_err("expected InvalidParam error");
        match err {
            PwdftError::InvalidParam { name, reason } => {
                assert_eq!(name, "pulay_period");
                assert!(
                    reason.contains("PeriodicPulay"),
                    "reason should mention PeriodicPulay, got: {reason}"
                );
            }
            other => panic!("expected InvalidParam, got: {other:?}"),
        }

        // And the Kerker variant too.
        let params_k = ScfParams {
            mixing_mode: mixing::MixingMode::PeriodicPulay {
                period: 0,
                kerker: true,
            },
            ..Default::default()
        };
        assert!(matches!(
            params_k.validate(),
            Err(PwdftError::InvalidParam { name, .. }) if name == "pulay_period"
        ));

        // Sanity: period = 1 must pass.
        let ok = ScfParams {
            mixing_mode: mixing::MixingMode::PeriodicPulay {
                period: 1,
                kerker: false,
            },
            ..Default::default()
        };
        assert!(ok.validate().is_ok());
    }

    /// ERR2 P1.c: pin the `Display` output of the new `InvalidParam`
    /// variant end-to-end through `ScfParams::validate` for the
    /// `mixing_beta` path (which interpolates the offending value).
    #[test]
    fn validate_mixing_beta_display_matches_invalid_param() {
        let params = ScfParams {
            mixing_beta: -0.1,
            ..Default::default()
        };
        let err = params
            .validate()
            .expect_err("expected InvalidParam error for mixing_beta = -0.1");
        assert_eq!(
            format!("{err}"),
            "invalid parameter mixing_beta: must be in (0, 1], got -0.1"
        );
    }

    /// ERR2 P1.c: every PARAM-cluster site in `ScfParams::validate`
    /// must return `InvalidParam` with its expected `name`. Pins all
    /// eight scf/mod.rs sites so a regression to `InvalidInput` fails
    /// here rather than silently downgrading.
    #[test]
    fn validate_param_cluster_names() {
        // n_bands
        let p = ScfParams {
            n_bands: 0,
            ..Default::default()
        };
        assert!(matches!(
            p.validate(),
            Err(PwdftError::InvalidParam { name, .. }) if name == "n_bands"
        ));

        // conv_threshold
        let p = ScfParams {
            conv_threshold: 0.0,
            ..Default::default()
        };
        assert!(matches!(
            p.validate(),
            Err(PwdftError::InvalidParam { name, .. }) if name == "conv_threshold"
        ));

        // mixing_beta (too large)
        let p = ScfParams {
            mixing_beta: 1.5,
            ..Default::default()
        };
        assert!(matches!(
            p.validate(),
            Err(PwdftError::InvalidParam { name, .. }) if name == "mixing_beta"
        ));

        // smearing_sigma
        let p = ScfParams {
            smearing_sigma: -0.1,
            ..Default::default()
        };
        assert!(matches!(
            p.validate(),
            Err(PwdftError::InvalidParam { name, .. }) if name == "smearing_sigma"
        ));

        // ecutrho_ratio
        let p = ScfParams {
            ecutrho_ratio: 0,
            ..Default::default()
        };
        assert!(matches!(
            p.validate(),
            Err(PwdftError::InvalidParam { name, .. }) if name == "ecutrho_ratio"
        ));

        // nspin
        let p = ScfParams {
            nspin: 3,
            ..Default::default()
        };
        assert!(matches!(
            p.validate(),
            Err(PwdftError::InvalidParam { name, .. }) if name == "nspin"
        ));

        // Sanity: a default-constructed params passes validation.
        assert!(ScfParams::default().validate().is_ok());
    }

    /// ERR2 P1.b: the three CRYSTAL-cluster preconditions in `run_scf`
    /// (empty atom list, empty k-point list, degenerate lattice volume)
    /// must surface as `PwdftError::InvalidCrystal`, not the catch-all
    /// `InvalidInput`. Exercises all three branches and pins both the
    /// discriminant and the `Display` output.
    #[test]
    fn run_scf_rejects_crystal_shape_errors() {
        use crate::crystal::{Atom, Lattice};
        use crate::symmetry::SymmetryInfo;
        use nalgebra::Vector3;

        let si_a = 5.431_f64;
        let lattice = Lattice::new(
            si_a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            si_a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            si_a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        );
        let basis = BasisSet::new(&lattice, 5.0);
        let params = ScfParams::default();
        let sym = SymmetryInfo::identity_only();

        // Branch 1: empty atoms.
        let crystal_no_atoms = Crystal {
            lattice: lattice.clone(),
            atoms: vec![],
        };
        let kpts = vec![KPoint {
            k: Vector3::new(0.0, 0.0, 0.0),
            weight: 1.0,
            label: None,
        }];
        let err = run_scf(&crystal_no_atoms, &basis, &kpts, &[], &params, &sym)
            .expect_err("empty atoms must fail");
        match &err {
            PwdftError::InvalidCrystal { reason } => {
                assert!(
                    reason.contains("atom"),
                    "InvalidCrystal reason should mention atom, got: {reason}"
                );
            }
            other => panic!("expected InvalidCrystal, got: {other:?}"),
        }
        assert_eq!(
            format!("{err}"),
            "invalid crystal input: at least one atom is required"
        );

        // Branch 2: empty k-point list.
        let crystal_ok = Crystal {
            lattice,
            atoms: vec![Atom::new(14, [0.0, 0.0, 0.0])],
        };
        let err = run_scf(&crystal_ok, &basis, &[], &[], &params, &sym)
            .expect_err("empty kpoints must fail");
        assert!(matches!(
            err,
            PwdftError::InvalidCrystal {
                reason: "at least one k-point is required"
            }
        ));

        // Branch 3: degenerate lattice (zero volume from co-planar vectors).
        let bad_lattice = Lattice::new(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        );
        let crystal_degenerate = Crystal {
            lattice: bad_lattice,
            atoms: vec![Atom::new(14, [0.0, 0.0, 0.0])],
        };
        let err = run_scf(&crystal_degenerate, &basis, &kpts, &[], &params, &sym)
            .expect_err("zero volume must fail");
        assert!(matches!(
            err,
            PwdftError::InvalidCrystal {
                reason: "lattice has zero or near-zero volume"
            }
        ));
    }
}
