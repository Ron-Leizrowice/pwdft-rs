//! Comprehensive settings/configuration for plane-wave DFT calculations.
//!
//! Parsed from YAML files via `serde_yaml_ng`. This is the primary (and only)
//! input configuration path, wired into `main.rs`.
//!
//! Defaults follow Quantum ESPRESSO conventions where applicable.

use std::collections::HashMap;
use std::path::Path;

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::{
    crystal::{Atom, Crystal, Lattice},
    error::{PwdftError, Result},
    kpoints::{HighSymPoint, KGridShift},
    scf::{smearing::SmearingScheme, ScfParams},
};

// ---------------------------------------------------------------------------
// Top-level Settings
// ---------------------------------------------------------------------------

/// Complete calculation settings parsed from a YAML file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Settings {
    /// Crystal structure: lattice vectors and atomic positions.
    pub system: SystemSettings,

    /// Plane-wave basis parameters.
    #[serde(default)]
    pub basis: BasisSettings,

    /// k-point sampling.
    pub kpoints: KPointSettings,

    /// Self-consistent field iteration control.
    #[serde(default)]
    pub scf: ScfSettings,

    /// Electronic structure parameters (mixing, smearing, occupations).
    #[serde(default)]
    pub electrons: ElectronSettings,

    /// Exchange-correlation functional.
    #[serde(default)]
    pub xc: XcSettings,

    /// Crystal symmetry handling.
    #[serde(default)]
    pub symmetry: SymmetrySettings,

    /// Pseudopotential file paths keyed by element symbol.
    #[serde(default)]
    pub pseudopotentials: PseudopotentialSettings,

    /// Output verbosity and write flags.
    #[serde(default)]
    pub output: OutputSettings,

    /// Initial-density (SAD) parameters.
    ///
    /// Exposes the `gaussian_sigma` knob for the starting guess. Omitting
    /// this block from YAML uses the default (see `InitialDensitySettings`).
    #[serde(default)]
    pub initial_density: InitialDensitySettings,
}

// ---------------------------------------------------------------------------
// Sub-structs
// ---------------------------------------------------------------------------

/// Crystal structure definition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemSettings {
    /// Lattice vectors in Angstroms: `[[ax,ay,az],[bx,by,bz],[cx,cy,cz]]`.
    pub lattice: [[f64; 3]; 3],

    /// Atomic positions in fractional (crystal) coordinates.
    #[serde(default)]
    pub atoms: Vec<AtomSetting>,
}

/// A single atom: element symbol and fractional position.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtomSetting {
    pub symbol: String,
    pub position: [f64; 3],
}

/// Plane-wave basis set parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct BasisSettings {
    /// Wavefunction kinetic-energy cutoff in eV.
    ///
    /// `None` means "resolve from the per-element recommended cutoff
    /// table"; a concrete value bypasses the lookup and is used
    /// verbatim. Callers must route through
    /// [`Settings::resolve_ecutwfc`] to turn `None` into a concrete
    /// eV value — `None` never leaks into the SCF layer.
    pub ecutwfc: Option<f64>,
    /// Charge-density cutoff ratio: ecutrho = ecutrho_ratio * ecutwfc.
    /// Default `4` is appropriate for norm-conserving pseudopotentials.
    pub ecutrho_ratio: u32,
    /// Explicit FFT grid dimensions [n1, n2, n3]. If set, overrides ecutrho_ratio.
    pub fft_grid: Option<[usize; 3]>,
}

impl Default for BasisSettings {
    fn default() -> Self {
        Self {
            ecutwfc: None,
            ecutrho_ratio: 4,
            fft_grid: None,
        }
    }
}

/// k-point sampling configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum KPointSettings {
    /// Monkhorst-Pack uniform grid.
    #[serde(rename = "monkhorst_pack")]
    MonkhorstPack {
        /// Grid dimensions [n1, n2, n3].
        grid: [u32; 3],
        /// Shift convention (Γ-centered vs MP-1976 shifted).
        ///
        /// Defaults to `gamma_centered` — the uniform mesh includes the
        /// Γ point (fractional coords `i/N` for `i ∈ 0..N`). Use `mp1976`
        /// for the original Monkhorst & Pack 1976 half-shift `(i + ½)/N`,
        /// or the free-form `[k1, k2, k3]` integers (each 0 or 1) for a
        /// per-axis half-shift.
        #[serde(default)]
        shift: KGridShift,
    },

    /// High-symmetry band path for band-structure calculations.
    #[serde(rename = "band_path")]
    BandPath {
        /// Ordered list of high-symmetry points defining the path.
        path: Vec<PathPointSetting>,
        /// Number of k-points per segment (default 50).
        #[serde(default = "KPointSettings::default_band_npoints")]
        npoints: usize,
    },
}

/// A high-symmetry point on a band path.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathPointSetting {
    /// Label (e.g. "G", "X", "L").
    pub label: String,
    /// Fractional reciprocal-space coordinates.
    pub frac: [f64; 3],
}

/// SCF iteration parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ScfSettings {
    /// Maximum number of SCF iterations.
    pub max_iter: usize,
    /// Convergence threshold: RMS density change (e/Å³).
    pub conv_threshold: f64,
    /// Energy convergence threshold (eV). Both density AND energy must converge.
    pub energy_threshold: f64,
    /// Number of Kohn-Sham bands. `None` = automatic from n_electrons/2 + padding.
    pub n_bands: Option<usize>,
    /// Eigensolver backend for per-k-point diagonalization.
    pub eigensolver: EigensolverType,
    /// Enable subspace-diagonalization warm-start for the dense
    /// eigensolver: project H into the previous iteration's eigenvector
    /// subspace before the full diagonalization (opt-in; ignored when
    /// `eigensolver == iterative`). Default `false`.
    pub wfrx_subspace: bool,
}

impl Default for ScfSettings {
    fn default() -> Self {
        Self {
            max_iter: 100,
            conv_threshold: 1e-6,
            energy_threshold: 1e-5,
            n_bands: None,
            eigensolver: EigensolverType::default(),
            wfrx_subspace: false,
        }
    }
}

/// Eigensolver backend selection (YAML-friendly adapter for
/// [`crate::eigensolver::EigensolverKind`]).
///
/// The default is `Dense`; flip to `Iterative` to activate the partial
/// Arnoldi / Krylov-Schur path that solves only for the lowest `n_bands`
/// eigenpairs (typically 3-10× faster at `n_pw ≥ 200`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum EigensolverType {
    /// Full dense Hermitian eigendecomposition (current default).
    #[default]
    Dense,
    /// ITEV: iterative partial Krylov-Schur eigensolver.
    Iterative,
}

impl From<EigensolverType> for crate::eigensolver::EigensolverKind {
    fn from(kind: EigensolverType) -> Self {
        match kind {
            EigensolverType::Dense => Self::Dense,
            EigensolverType::Iterative => Self::Iterative,
        }
    }
}

/// Electronic-structure parameters: mixing, smearing, occupations.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ElectronSettings {
    /// Density mixing parameter (`0 < beta <= 1`).
    pub mixing_beta: f64,
    /// Number of past densities kept for Anderson/Pulay mixing.
    pub mixing_ndim: usize,
    /// Smearing scheme for partial occupations.
    pub smearing: SmearingScheme,
    /// Smearing width in eV.
    pub smearing_width: f64,
    /// Occupation scheme.
    pub occupations: OccupationType,
    /// Mixing preconditioning mode: plain Anderson or Kerker-preconditioned.
    pub mixing_mode: MixingModeType,
    /// Periodic Pulay period k (Banerjee et al., JCTC 12, 3053 (2016)).
    ///
    /// Only consulted when `mixing_mode == periodic_pulay` or
    /// `periodic_pulay_kerker`; ignored otherwise. Must be ≥ 1. Default: 3.
    pub pulay_period: usize,
    /// Enable adaptive mixing β (Eyert 1996, §3.3 residual-norm monitor).
    ///
    /// When `true`, β is damped when the residual norm grows and restored
    /// toward `mixing_beta` when it decreases steadily for three
    /// consecutive iterations. Default `false` uses the fixed `mixing_beta`.
    /// See `src/scf/mixing/mod.rs` module docs for the rule and thresholds.
    pub adaptive_beta: bool,
    /// Number of spin channels: 1 (unpolarized) or 2 (collinear spin-polarized).
    pub nspin: usize,
    /// Starting magnetization per atom type (fractional, -1 to 1).
    /// Maps from element symbol to magnetization. Empty = non-magnetic.
    pub starting_magnetization: HashMap<String, f64>,
    /// Fixed total magnetization (n_up - n_down) in electrons.
    /// If None, magnetization is determined self-consistently.
    pub tot_magnetization: Option<f64>,
}

impl Default for ElectronSettings {
    fn default() -> Self {
        Self {
            mixing_beta: 0.3,
            mixing_ndim: 8,
            smearing: SmearingScheme::default(),
            smearing_width: 0.05,
            occupations: OccupationType::default(),
            mixing_mode: MixingModeType::default(),
            pulay_period: 3,
            adaptive_beta: false,
            nspin: 1,
            starting_magnetization: HashMap::new(),
            tot_magnetization: None,
        }
    }
}

/// How occupation numbers are determined.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum OccupationType {
    /// Use smearing to determine occupations (metals / finite-T).
    #[default]
    Smearing,
    /// Fixed integer occupations (insulators at T=0).
    Fixed,
}

/// Mixing preconditioning mode for SCF density mixing.
///
/// This is the serde-friendly adapter for [`crate::scf::mixing::MixingMode`].
/// The SCF-internal variants carry run-time parameters (auto-estimated q_tf,
/// configured Pulay period) that are injected during conversion; we keep this
/// flat enum for YAML parsing and convert via
/// [`MixingModeType::to_scf_mode`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum MixingModeType {
    /// Standard Anderson mixing (no preconditioning).
    #[default]
    Plain,
    /// Kerker preconditioning.
    Kerker,
    /// Modified Broyden mixing (Johnson PRB 38, 12807, 1988).
    Broyden,
    /// Modified Broyden mixing with Kerker preconditioning.
    BroydenKerker,
    /// Periodic Pulay mixing (Banerjee et al., JCTC 12, 3053 (2016)).
    ///
    /// Plain linear mixing except every `pulay_period`-th iteration, where a
    /// DIIS extrapolation is performed over the accumulated history.
    PeriodicPulay,
    /// Periodic Pulay mixing with Kerker preconditioning on the residual.
    PeriodicPulayKerker,
}

impl MixingModeType {
    /// Convert to the runtime `MixingMode`, supplying `pulay_period` for the
    /// Periodic Pulay variants (ignored for others).
    #[must_use]
    pub fn to_scf_mode(self, pulay_period: usize) -> crate::scf::mixing::MixingMode {
        use crate::scf::mixing::MixingMode;
        match self {
            MixingModeType::Plain => MixingMode::Plain,
            MixingModeType::Kerker => MixingMode::Kerker { q_tf: None },
            MixingModeType::Broyden => MixingMode::Broyden { kerker: false },
            MixingModeType::BroydenKerker => MixingMode::Broyden { kerker: true },
            MixingModeType::PeriodicPulay => MixingMode::PeriodicPulay {
                period: pulay_period,
                kerker: false,
            },
            MixingModeType::PeriodicPulayKerker => MixingMode::PeriodicPulay {
                period: pulay_period,
                kerker: true,
            },
        }
    }
}

// Back-compat: preserve the simple `.into()` for modes that don't consume
// `pulay_period`. For Periodic Pulay, callers must use `to_scf_mode(period)`
// explicitly so the period is threaded through from settings.
impl From<MixingModeType> for crate::scf::mixing::MixingMode {
    fn from(mode: MixingModeType) -> Self {
        // Default period of 3 matches the `ElectronSettings` default and the
        // paper's recommendation; callers that want a non-default period use
        // `to_scf_mode` instead.
        mode.to_scf_mode(3)
    }
}

/// Exchange-correlation functional specification.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct XcSettings {
    /// Functional name. Currently implemented: `"pz"` (Perdew-Zunger LDA).
    /// The YAML parser also accepts `"pbe"`, `"pbe0"`, `"hse06"`, but these
    /// are rejected with [`crate::error::PwdftError::NotImplemented`] at
    /// SCF entry — a YAML typo cannot silently produce LDA numbers under
    /// a GGA or hybrid label.
    pub functional: XcFunctional,
}

/// Supported exchange-correlation functionals.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum XcFunctional {
    /// Perdew-Zunger LDA (Ceperley-Alder).
    #[default]
    Pz,
    /// Perdew-Burke-Ernzerhof GGA.
    ///
    /// **Not yet implemented.** Returns
    /// [`crate::error::PwdftError::NotImplemented`] at SCF entry so a
    /// YAML typo cannot silently produce LDA numbers under a PBE label.
    Pbe,
    /// PBE0 hybrid functional.
    ///
    /// **Not yet implemented.** Returns
    /// [`crate::error::PwdftError::NotImplemented`] at SCF entry.
    Pbe0,
    /// Heyd-Scuseria-Ernzerhof screened hybrid.
    ///
    /// **Not yet implemented.** Returns
    /// [`crate::error::PwdftError::NotImplemented`] at SCF entry.
    Hse06,
}

/// Crystal symmetry settings.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct SymmetrySettings {
    /// Whether to detect and exploit crystal symmetry.
    pub enabled: bool,
    /// Whether to apply time-reversal symmetry (`k → -k`).
    pub time_reversal: bool,
    /// Tolerance for symmetry detection in fractional coordinates.
    pub tolerance: f64,
}

impl Default for SymmetrySettings {
    fn default() -> Self {
        Self {
            enabled: true,
            time_reversal: true,
            tolerance: 1e-5,
        }
    }
}

/// Pseudopotential file paths keyed by element symbol.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(transparent)]
pub struct PseudopotentialSettings {
    /// Map from element symbol (e.g. "Si") to file path.
    pub files: HashMap<String, String>,
}

/// Verbosity level for output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Verbosity {
    Low,
    #[default]
    Normal,
    High,
}

/// Output and I/O settings.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct OutputSettings {
    /// Verbosity level.
    pub verbosity: Verbosity,
    /// Whether to write the converged charge density to a file.
    pub write_density: bool,
    /// Whether to write band-structure eigenvalues.
    pub write_bands: bool,
}

impl Default for OutputSettings {
    fn default() -> Self {
        Self {
            verbosity: Verbosity::Normal,
            write_density: false,
            write_bands: true,
        }
    }
}

/// Initial-density (Superposition of Atomic Densities) parameters.
///
/// The SAD initial guess uses a pseudopotential's `PP_RHOATOM` when
/// available; otherwise it falls back to a Gaussian model of width
/// `gaussian_sigma` (Å) per atom. This block exposes that width as a
/// configurable YAML knob.
///
/// Default is the canonical constant
/// `scf::initial_density::DEFAULT_GAUSSIAN_SIGMA` (1.0 Å).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default)]
pub struct InitialDensitySettings {
    /// Gaussian model-charge width in Å.
    ///
    /// Must be positive and finite. Typical values: 0.5–2.0 Å. Wider
    /// sigma smooths the initial high-|G| content; narrower sigma gives
    /// a sharper but numerically stiffer starting density.
    pub gaussian_sigma: f64,
}

impl Default for InitialDensitySettings {
    fn default() -> Self {
        Self {
            gaussian_sigma: crate::scf::initial_density::DEFAULT_GAUSSIAN_SIGMA,
        }
    }
}

// ---------------------------------------------------------------------------
// Default value functions for serde (only where #[serde(default)] on struct
// level doesn't work — e.g. enum variants with non-Default values)
// ---------------------------------------------------------------------------

impl KPointSettings {
    fn default_band_npoints() -> usize { 50 }
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

impl Settings {
    /// Parse settings from a YAML string.
    ///
    /// # Errors
    ///
    /// Returns [`PwdftError::Parse`] wrapping the underlying
    /// `serde_yaml_ng` error when the YAML is syntactically invalid or
    /// fails to deserialize into [`Settings`] (unknown field, wrong type,
    /// missing required key, etc.).
    pub fn from_yaml_str(s: &str) -> Result<Self> {
        serde_yaml_ng::from_str(s).map_err(|e| PwdftError::Parse(format!("YAML parse error: {e}")))
    }

    /// Parse settings from a YAML file on disk.
    ///
    /// # Errors
    ///
    /// - `PwdftError::Io` (from the `?` on `std::fs::read_to_string`) if
    ///   the file cannot be opened or read.
    /// - Any `PwdftError::Parse` forwarded from [`Self::from_yaml_str`].
    pub fn from_yaml_file(path: &Path) -> Result<Self> {
        let contents = std::fs::read_to_string(path)?;
        Self::from_yaml_str(&contents)
    }

    // -----------------------------------------------------------------------
    // Conversion helpers — bridge to existing crate types
    // -----------------------------------------------------------------------

    /// Build a `Crystal` from the system settings.
    ///
    /// # Errors
    /// Returns `PwdftError::InvalidInput` if any atom has an unrecognized element symbol.
    pub fn to_crystal(&self) -> Result<Crystal> {
        let [a, b, c] = self.system.lattice;
        let lattice = Lattice::new(
            Vector3::new(a[0], a[1], a[2]),
            Vector3::new(b[0], b[1], b[2]),
            Vector3::new(c[0], c[1], c[2]),
        );

        let atoms = self
            .system
            .atoms
            .iter()
            .map(|ai| {
                let elem = crate::atoms::Element::iter()
                    .find(|e| e.symbol() == ai.symbol)
                    .ok_or_else(|| {
                        PwdftError::InvalidInput(format!("unknown element: {}", ai.symbol))
                    })?;
                Ok(Atom::new(elem.atomic_number(), ai.position))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Crystal { atoms, lattice })
    }

    /// Build `ScfParams` from the SCF + electron + basis settings.
    ///
    /// `n_bands_fallback` is used when `scf.n_bands` is `None` (auto mode).
    #[must_use]
    pub fn to_scf_params(&self, n_bands_fallback: usize) -> ScfParams {
        ScfParams {
            n_bands: self.scf.n_bands.unwrap_or(n_bands_fallback),
            max_iter: self.scf.max_iter,
            conv_threshold: self.scf.conv_threshold,
            energy_threshold: self.scf.energy_threshold,
            mixing_beta: self.electrons.mixing_beta,
            mixing_ndim: self.electrons.mixing_ndim,
            smearing_sigma: self.electrons.smearing_width,
            smearing_scheme: self.electrons.smearing,
            ecutrho_ratio: self.basis.ecutrho_ratio,
            fft_grid: self.basis.fft_grid,
            mixing_mode: self
                .electrons
                .mixing_mode
                .to_scf_mode(self.electrons.pulay_period),
            adaptive_beta: self.electrons.adaptive_beta,
            nspin: self.electrons.nspin,
            starting_magnetization: self.electrons.starting_magnetization.clone(),
            tot_magnetization: self.electrons.tot_magnetization,
            eigensolver: self.scf.eigensolver.into(),
            wfrx_subspace: self.scf.wfrx_subspace,
            xc_functional: self.xc.functional,
            gaussian_sigma: self.initial_density.gaussian_sigma,
        }
    }

    /// Resolve the wavefunction energy cutoff in eV for this crystal.
    ///
    /// - If `basis.ecutwfc` is set in YAML, it is returned verbatim and
    ///   no warning is emitted.
    /// - Otherwise the per-element recommended cutoff table is consulted
    ///   and the maximum across all species in the crystal is used. A
    ///   [`log::warn!`] is emitted naming the driving species so users
    ///   cannot silently run with a sub-convergence cutoff.
    /// - If the crystal contains only elements outside the table, the
    ///   caller receives an [`PwdftError::InvalidParam`] asking them to
    ///   set `basis.ecutwfc` explicitly.
    ///
    /// # Errors
    ///
    /// Returns [`PwdftError::InvalidParam`] when no species in
    /// `crystal` has a tabulated recommended cutoff and `basis.ecutwfc`
    /// is unset — the calculation cannot proceed with an unknown cutoff.
    pub fn resolve_ecutwfc(&self, crystal: &Crystal) -> Result<f64> {
        use crate::pseudopotential::recommended_ecut::{
            recommended_ecut_for_crystal, EcutVariant,
        };

        if let Some(user) = self.basis.ecutwfc {
            return Ok(user);
        }

        match recommended_ecut_for_crystal(crystal, EcutVariant::Standard) {
            Some(rec) => {
                let symbol = crate::atoms::Element::iter()
                    .find(|e| e.atomic_number() == rec.z)
                    .map_or_else(|| format!("Z={}", rec.z), |e| e.symbol().to_string());
                log::warn!(
                    "basis.ecutwfc not set in YAML; using {:.2} eV recommended for {} (PseudoDojo .standard, max over species). Set `basis.ecutwfc` in YAML to override.",
                    rec.ev,
                    symbol,
                );
                Ok(rec.ev)
            }
            None => Err(PwdftError::InvalidParam {
                name: "basis.ecutwfc",
                reason: "not set and no species in the crystal has a tabulated recommended cutoff; set `basis.ecutwfc` explicitly in YAML".to_string(),
            }),
        }
    }

    /// Extract the Monkhorst-Pack grid dimensions, if configured.
    #[must_use]
    pub fn mp_grid(&self) -> Option<[u32; 3]> {
        match &self.kpoints {
            KPointSettings::MonkhorstPack { grid, .. } => Some(*grid),
            _ => None,
        }
    }

    /// Extract the Monkhorst-Pack shift convention, if configured.
    #[must_use]
    pub fn mp_shift(&self) -> Option<KGridShift> {
        match &self.kpoints {
            KPointSettings::MonkhorstPack { shift, .. } => Some(*shift),
            _ => None,
        }
    }

    /// Extract the high-symmetry path for band-structure calculations.
    #[must_use]
    pub fn to_high_sym_path(&self) -> Option<Vec<HighSymPoint>> {
        match &self.kpoints {
            KPointSettings::BandPath { path, .. } => Some(
                path.iter()
                    .map(|p| HighSymPoint {
                        label: p.label.clone(),
                        frac: p.frac,
                    })
                    .collect(),
            ),
            _ => None,
        }
    }

    /// Number of k-points per band-path segment (if band_path mode).
    #[must_use]
    pub fn band_path_npoints(&self) -> Option<usize> {
        match &self.kpoints {
            KPointSettings::BandPath { npoints, .. } => Some(*npoints),
            _ => None,
        }
    }

    /// Build `SymmetryInfo` from these settings and a crystal.
    ///
    /// When `symmetry.enabled = false`, this returns the trivial group
    /// (identity-only, no time reversal) via
    /// [`crate::symmetry::SymmetryInfo::identity_only`] rather than `None`.
    /// Downstream code treats that as "no symmetrization to apply",
    /// bit-identically to the legacy `Option::None` path.
    #[must_use]
    pub fn to_symmetry_info(
        &self,
        crystal: &Crystal,
    ) -> crate::symmetry::SymmetryInfo {
        if !self.symmetry.enabled {
            return crate::symmetry::SymmetryInfo::identity_only();
        }
        let mut info =
            crate::symmetry::SymmetryInfo::from_crystal(crystal, self.symmetry.tolerance);
        info.has_time_reversal = self.symmetry.time_reversal;
        info
    }

    /// Return the pseudopotential file path for a given element symbol.
    #[must_use]
    pub fn pseudopotential_path(&self, symbol: &str) -> Option<&str> {
        self.pseudopotentials.files.get(symbol).map(|s| s.as_str())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal valid YAML with only required fields.
    const MINIMAL_YAML: &str = r#"
system:
  lattice:
    - [0.0, 2.7155, 2.7155]
    - [2.7155, 0.0, 2.7155]
    - [2.7155, 2.7155, 0.0]
  atoms:
    - symbol: Si
      position: [0.0, 0.0, 0.0]
    - symbol: Si
      position: [0.25, 0.25, 0.25]

basis:
  ecutwfc: 204.09

kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]
"#;

    /// Full YAML exercising every field.
    const FULL_YAML: &str = r#"
system:
  lattice:
    - [0.0, 2.7155, 2.7155]
    - [2.7155, 0.0, 2.7155]
    - [2.7155, 2.7155, 0.0]
  atoms:
    - symbol: Si
      position: [0.0, 0.0, 0.0]
    - symbol: Si
      position: [0.25, 0.25, 0.25]

basis:
  ecutwfc: 204.09
  ecutrho_ratio: 4
  fft_grid: [24, 24, 24]

kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]

scf:
  max_iter: 100
  conv_threshold: 1.0e-6
  energy_threshold: 1.0e-5
  n_bands: 8

electrons:
  mixing_beta: 0.3
  mixing_ndim: 8
  smearing: fermi_dirac
  smearing_width: 0.05
  occupations: smearing
  mixing_mode: plain
  nspin: 1
  starting_magnetization: {}
  tot_magnetization: null

xc:
  functional: pz

symmetry:
  enabled: true
  time_reversal: true
  tolerance: 1.0e-5

pseudopotentials:
  Si: "../pseudopotentials/nc/lda/Si.upf"

output:
  verbosity: normal
  write_density: false
  write_bands: true
"#;

    const BAND_PATH_YAML: &str = r#"
system:
  lattice:
    - [0.0, 2.7155, 2.7155]
    - [2.7155, 0.0, 2.7155]
    - [2.7155, 2.7155, 0.0]
  atoms:
    - symbol: Si
      position: [0.0, 0.0, 0.0]
    - symbol: Si
      position: [0.25, 0.25, 0.25]

basis:
  ecutwfc: 200.0

kpoints:
  type: band_path
  path:
    - { label: "G", frac: [0.0, 0.0, 0.0] }
    - { label: "X", frac: [0.5, 0.0, 0.5] }
    - { label: "G", frac: [0.0, 0.0, 0.0] }
    - { label: "L", frac: [0.5, 0.5, 0.5] }
  npoints: 50
"#;

    #[test]
    fn parse_minimal_yaml() {
        let s = Settings::from_yaml_str(MINIMAL_YAML).unwrap();
        assert_eq!(s.system.atoms.len(), 2);
        assert_eq!(s.system.atoms[0].symbol, "Si");
        assert!((s.basis.ecutwfc.unwrap() - 204.09).abs() < f64::EPSILON);
        assert!(s.mp_grid().is_some());
        assert_eq!(s.mp_grid().unwrap(), [4, 4, 4]);
        // MPSH default: Γ-centered (matches QE's `automatic / … 0 0 0`).
        assert_eq!(s.mp_shift(), Some(KGridShift::GammaCentered));
    }

    #[test]
    fn parse_monkhorst_pack_shift_gamma_centered_explicit() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]
  shift: gamma_centered
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.mp_shift(), Some(KGridShift::GammaCentered));
    }

    #[test]
    fn parse_monkhorst_pack_shift_mp1976() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]
  shift: mp1976
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.mp_shift(), Some(KGridShift::MP1976));
    }

    #[test]
    fn parse_monkhorst_pack_shift_custom() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]
  shift:
    custom: [1, 0, 1]
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.mp_shift(), Some(KGridShift::Custom([1, 0, 1])));
    }

    #[test]
    fn parse_monkhorst_pack_shift_missing_defaults_to_gamma() {
        // MPSH: missing `shift:` in YAML must default to Γ-centered
        // (QE convention). This is the default-behavior-drift guard.
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.mp_shift(), Some(KGridShift::GammaCentered));
    }

    #[test]
    fn parse_full_yaml() {
        let s = Settings::from_yaml_str(FULL_YAML).unwrap();
        assert_eq!(s.scf.max_iter, 100);
        assert_eq!(s.scf.n_bands, Some(8));
        assert!((s.scf.energy_threshold - 1e-5).abs() < 1e-15);
        assert_eq!(s.basis.fft_grid, Some([24, 24, 24]));
        assert!((s.electrons.mixing_beta - 0.3).abs() < f64::EPSILON);
        assert_eq!(s.electrons.mixing_ndim, 8);
        assert_eq!(s.electrons.smearing, SmearingScheme::FermiDirac);
        assert!((s.electrons.smearing_width - 0.05).abs() < 1e-15);
        assert_eq!(s.electrons.occupations, OccupationType::Smearing);
        assert_eq!(s.electrons.mixing_mode, MixingModeType::Plain);
        assert_eq!(s.electrons.nspin, 1);
        assert!(s.electrons.starting_magnetization.is_empty());
        assert!(s.electrons.tot_magnetization.is_none());
        assert_eq!(s.xc.functional, XcFunctional::Pz);
        assert!(s.symmetry.enabled);
        assert!(s.symmetry.time_reversal);
        assert!((s.symmetry.tolerance - 1e-5).abs() < 1e-15);
        assert_eq!(
            s.pseudopotentials.files.get("Si").unwrap(),
            "../pseudopotentials/nc/lda/Si.upf"
        );
        assert_eq!(s.output.verbosity, Verbosity::Normal);
        assert!(!s.output.write_density);
        assert!(s.output.write_bands);
    }

    #[test]
    fn parse_band_path() {
        let s = Settings::from_yaml_str(BAND_PATH_YAML).unwrap();
        let path = s.to_high_sym_path().unwrap();
        assert_eq!(path.len(), 4);
        assert_eq!(path[0].label, "G");
        assert_eq!(path[1].label, "X");
        assert_eq!(s.band_path_npoints(), Some(50));
    }

    #[test]
    fn defaults_match_qe_conventions() {
        let s = Settings::from_yaml_str(MINIMAL_YAML).unwrap();

        assert_eq!(s.basis.ecutrho_ratio, 4);
        assert!(s.basis.fft_grid.is_none());
        assert_eq!(s.scf.max_iter, 100);
        assert!((s.scf.conv_threshold - 1e-6).abs() < 1e-15);
        assert!((s.scf.energy_threshold - 1e-5).abs() < 1e-15);
        assert!(s.scf.n_bands.is_none());
        assert!((s.electrons.mixing_beta - 0.3).abs() < 1e-15);
        assert_eq!(s.electrons.mixing_ndim, 8);
        assert_eq!(s.electrons.smearing, SmearingScheme::FermiDirac);
        assert!((s.electrons.smearing_width - 0.05).abs() < 1e-15);
        assert_eq!(s.electrons.occupations, OccupationType::Smearing);
        assert_eq!(s.electrons.mixing_mode, MixingModeType::Plain);
        assert_eq!(s.electrons.nspin, 1);
        assert!(s.electrons.starting_magnetization.is_empty());
        assert!(s.electrons.tot_magnetization.is_none());
        assert_eq!(s.xc.functional, XcFunctional::Pz);
        assert!(s.symmetry.enabled);
        assert!(s.symmetry.time_reversal);
        assert!((s.symmetry.tolerance - 1e-5).abs() < 1e-15);
        assert_eq!(s.output.verbosity, Verbosity::Normal);
        assert!(!s.output.write_density);
        assert!(s.output.write_bands);
    }

    #[test]
    fn to_crystal_produces_correct_structure() {
        let s = Settings::from_yaml_str(FULL_YAML).unwrap();
        let crystal = s.to_crystal().unwrap();
        assert_eq!(crystal.atoms.len(), 2);
        assert_eq!(crystal.atoms[0].z, 14);
        assert_eq!(crystal.atoms[1].z, 14);
        for (&got, &want) in crystal.atoms[1].position.iter().zip([0.25, 0.25, 0.25].iter()) {
            assert!((got - want).abs() < f64::EPSILON);
        }
        let lat = &crystal.lattice;
        assert!((lat.a[0] - 0.0).abs() < 1e-10);
        assert!((lat.a[1] - 2.7155).abs() < 1e-10);
    }

    #[test]
    fn to_scf_params_merges_sections() {
        let s = Settings::from_yaml_str(FULL_YAML).unwrap();
        let params = s.to_scf_params(4);
        assert_eq!(params.n_bands, 8);
        assert_eq!(params.max_iter, 100);
        assert!((params.mixing_beta - 0.3).abs() < 1e-15);
        assert_eq!(params.ecutrho_ratio, 4);
        assert!((params.energy_threshold - 1e-5).abs() < 1e-15);
        assert_eq!(params.fft_grid, Some([24, 24, 24]));
        assert_eq!(
            params.smearing_scheme,
            crate::scf::smearing::SmearingScheme::FermiDirac
        );
        assert!(matches!(
            params.mixing_mode,
            crate::scf::mixing::MixingMode::Plain
        ));
        assert_eq!(params.nspin, 1);
        assert!(params.starting_magnetization.is_empty());
        assert!(params.tot_magnetization.is_none());
    }

    #[test]
    fn to_scf_params_threads_xc_functional() {
        // XCNI: the XC functional must flow from YAML → ScfParams so the
        // SCF entry dispatch can fail-fast on non-LDA options.
        let s_lda = Settings::from_yaml_str(FULL_YAML).unwrap();
        let p_lda = s_lda.to_scf_params(4);
        assert_eq!(p_lda.xc_functional, XcFunctional::Pz);

        // `xc_functional: pbe` must propagate through, so the SCF layer
        // can trigger NotImplemented at entry.
        let yaml_pbe = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
xc:
  functional: pbe
"#;
        let s_pbe = Settings::from_yaml_str(yaml_pbe).unwrap();
        let p_pbe = s_pbe.to_scf_params(4);
        assert_eq!(p_pbe.xc_functional, XcFunctional::Pbe);
    }

    #[test]
    fn to_scf_params_uses_fallback_n_bands() {
        let s = Settings::from_yaml_str(MINIMAL_YAML).unwrap();
        let params = s.to_scf_params(12);
        assert_eq!(params.n_bands, 12);
    }

    #[test]
    fn pseudopotential_path_lookup() {
        let s = Settings::from_yaml_str(FULL_YAML).unwrap();
        assert_eq!(
            s.pseudopotential_path("Si"),
            Some("../pseudopotentials/nc/lda/Si.upf")
        );
        assert_eq!(s.pseudopotential_path("Ge"), None);
    }

    #[test]
    fn smearing_types_roundtrip() {
        for variant in [
            SmearingScheme::FermiDirac,
            SmearingScheme::Gaussian,
            SmearingScheme::MethfesselPaxton,
            SmearingScheme::Cold,
            SmearingScheme::Fixed,
        ] {
            let yaml = serde_yaml_ng::to_string(&variant).unwrap();
            let parsed: SmearingScheme = serde_yaml_ng::from_str(&yaml).unwrap();
            assert_eq!(parsed, variant);
        }
    }

    #[test]
    fn xc_functional_roundtrip() {
        for variant in [XcFunctional::Pz, XcFunctional::Pbe, XcFunctional::Pbe0, XcFunctional::Hse06] {
            let yaml = serde_yaml_ng::to_string(&variant).unwrap();
            let parsed: XcFunctional = serde_yaml_ng::from_str(&yaml).unwrap();
            assert_eq!(parsed, variant);
        }
    }

    #[test]
    fn verbosity_roundtrip() {
        for variant in [Verbosity::Low, Verbosity::Normal, Verbosity::High] {
            let yaml = serde_yaml_ng::to_string(&variant).unwrap();
            let parsed: Verbosity = serde_yaml_ng::from_str(&yaml).unwrap();
            assert_eq!(parsed, variant);
        }
    }

    #[test]
    fn occupation_type_roundtrip() {
        for variant in [OccupationType::Smearing, OccupationType::Fixed] {
            let yaml = serde_yaml_ng::to_string(&variant).unwrap();
            let parsed: OccupationType = serde_yaml_ng::from_str(&yaml).unwrap();
            assert_eq!(parsed, variant);
        }
    }

    #[test]
    fn mixing_mode_type_roundtrip() {
        for variant in [
            MixingModeType::Plain,
            MixingModeType::Kerker,
            MixingModeType::Broyden,
            MixingModeType::BroydenKerker,
            MixingModeType::PeriodicPulay,
            MixingModeType::PeriodicPulayKerker,
        ] {
            let yaml = serde_yaml_ng::to_string(&variant).unwrap();
            let parsed: MixingModeType = serde_yaml_ng::from_str(&yaml).unwrap();
            assert_eq!(parsed, variant);
        }
    }

    #[test]
    fn mixing_mode_type_converts_to_scf_mixing_mode() {
        use crate::scf::mixing::MixingMode;

        let plain: MixingMode = MixingModeType::Plain.into();
        assert!(matches!(plain, MixingMode::Plain));

        let kerker: MixingMode = MixingModeType::Kerker.into();
        // Kerker conversion should default q_tf to None (auto-estimated at runtime)
        assert!(matches!(kerker, MixingMode::Kerker { q_tf: None }));

        let broyden: MixingMode = MixingModeType::Broyden.into();
        assert!(matches!(broyden, MixingMode::Broyden { kerker: false }));

        let broyden_kerker: MixingMode = MixingModeType::BroydenKerker.into();
        assert!(matches!(broyden_kerker, MixingMode::Broyden { kerker: true }));

        // Periodic Pulay picks up the supplied period.
        let pp = MixingModeType::PeriodicPulay.to_scf_mode(5);
        assert!(matches!(
            pp,
            MixingMode::PeriodicPulay {
                period: 5,
                kerker: false
            }
        ));

        let pp_k = MixingModeType::PeriodicPulayKerker.to_scf_mode(7);
        assert!(matches!(
            pp_k,
            MixingMode::PeriodicPulay {
                period: 7,
                kerker: true
            }
        ));

        // Default period via `Into` is 3 (matches `ElectronSettings` default).
        let pp_default: MixingMode = MixingModeType::PeriodicPulay.into();
        assert!(matches!(
            pp_default,
            MixingMode::PeriodicPulay {
                period: 3,
                kerker: false
            }
        ));
    }

    #[test]
    fn periodic_pulay_mixing_mode_parses() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
electrons:
  mixing_mode: periodic_pulay
  pulay_period: 4
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.electrons.mixing_mode, MixingModeType::PeriodicPulay);
        assert_eq!(s.electrons.pulay_period, 4);
    }

    #[test]
    fn periodic_pulay_kerker_mixing_mode_parses() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
electrons:
  mixing_mode: periodic_pulay_kerker
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(
            s.electrons.mixing_mode,
            MixingModeType::PeriodicPulayKerker
        );
        // Default period
        assert_eq!(s.electrons.pulay_period, 3);
    }

    #[test]
    fn adaptive_beta_defaults_off_and_parses_from_yaml() {
        // Default: absent from YAML → `false` (preserves pre-MXBA behaviour).
        let yaml_no_key = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
electrons:
  mixing_beta: 0.25
"#;
        let s = Settings::from_yaml_str(yaml_no_key).unwrap();
        assert!(!s.electrons.adaptive_beta, "default must be false");
        // And the ScfParams round-trip preserves that.
        let params = s.to_scf_params(8);
        assert!(!params.adaptive_beta);

        // Explicit opt-in in YAML.
        let yaml_on = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
electrons:
  adaptive_beta: true
"#;
        let s_on = Settings::from_yaml_str(yaml_on).unwrap();
        assert!(s_on.electrons.adaptive_beta);
        assert!(s_on.to_scf_params(8).adaptive_beta);
    }

    // -----------------------------------------------------------------------
    // CFGN Phase 1 — initial_density.gaussian_sigma plumbing tests
    // -----------------------------------------------------------------------

    #[test]
    fn initial_density_gaussian_sigma_default_matches_pre_cfgn_constant() {
        // Defaulting path: no `initial_density` block in YAML. The field
        // must land at the pre-CFGN hardcoded 1.0 Å, i.e.
        // `DEFAULT_GAUSSIAN_SIGMA`, preserving bit-identical starting
        // densities for every existing YAML input.
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        let expected = crate::scf::initial_density::DEFAULT_GAUSSIAN_SIGMA;
        assert!(
            (s.initial_density.gaussian_sigma - expected).abs() < 1e-15,
            "default gaussian_sigma {} != DEFAULT_GAUSSIAN_SIGMA {}",
            s.initial_density.gaussian_sigma,
            expected,
        );
        // And the value must survive the Settings -> ScfParams plumbing.
        let p = s.to_scf_params(8);
        assert!(
            (p.gaussian_sigma - expected).abs() < 1e-15,
            "ScfParams.gaussian_sigma {} != default {}",
            p.gaussian_sigma,
            expected,
        );
    }

    #[test]
    fn initial_density_gaussian_sigma_override_plumbs_through_to_scf_params() {
        // Override path: a non-default value in YAML must flow unchanged
        // into `ScfParams`. This is the CFGN Phase 1 plumbing test.
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
initial_density:
  gaussian_sigma: 0.75
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert!(
            (s.initial_density.gaussian_sigma - 0.75).abs() < 1e-15,
            "YAML override not picked up: {}",
            s.initial_density.gaussian_sigma,
        );
        let p = s.to_scf_params(8);
        assert!(
            (p.gaussian_sigma - 0.75).abs() < 1e-15,
            "override not threaded into ScfParams: {}",
            p.gaussian_sigma,
        );
    }

    #[test]
    fn initial_density_settings_roundtrip_via_yaml() {
        // Serialize -> parse should preserve the value exactly.
        let mut s = Settings::from_yaml_str(FULL_YAML).unwrap();
        s.initial_density.gaussian_sigma = 1.25;
        let yaml = serde_yaml_ng::to_string(&s).unwrap();
        let restored = Settings::from_yaml_str(&yaml).unwrap();
        assert!(
            (restored.initial_density.gaussian_sigma - 1.25).abs() < 1e-15
        );
    }

    #[test]
    fn scf_params_validate_rejects_non_positive_gaussian_sigma() {
        // `ScfParams::validate()` must reject sigma <= 0 or non-finite
        // sigma with a structured error — these values would poison
        // the SAD Gaussian exp(-|G|^2 sigma^2 / 2) term. ERR2 P1.c
        // migrated this path from the catch-all `InvalidInput` to the
        // structured `InvalidParam` variant.
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let params = ScfParams {
                gaussian_sigma: bad,
                ..Default::default()
            };
            match params.validate() {
                Err(PwdftError::InvalidParam { name, reason }) => {
                    assert_eq!(name, "gaussian_sigma");
                    assert!(
                        reason.contains("positive"),
                        "reason should mention positivity, got: {reason}"
                    );
                }
                other => panic!(
                    "expected InvalidParam for gaussian_sigma = {bad}, got: {other:?}"
                ),
            }
        }

        // Sanity: a positive finite value passes.
        let ok = ScfParams {
            gaussian_sigma: 1.0,
            ..Default::default()
        };
        assert!(ok.validate().is_ok());
    }

    #[test]
    fn full_settings_roundtrip() {
        let original = Settings::from_yaml_str(FULL_YAML).unwrap();
        let serialized = serde_yaml_ng::to_string(&original).unwrap();
        let restored = Settings::from_yaml_str(&serialized).unwrap();

        assert_eq!(original.system.atoms.len(), restored.system.atoms.len());
        assert_eq!(original.basis.ecutwfc, restored.basis.ecutwfc);
        assert_eq!(original.scf.max_iter, restored.scf.max_iter);
        assert_eq!(original.scf.n_bands, restored.scf.n_bands);
        assert!((original.electrons.mixing_beta - restored.electrons.mixing_beta).abs() < f64::EPSILON);
        assert_eq!(original.xc.functional, restored.xc.functional);
        assert_eq!(original.symmetry.enabled, restored.symmetry.enabled);
        assert_eq!(original.output.verbosity, restored.output.verbosity);
    }

    #[test]
    fn invalid_yaml_returns_parse_error() {
        let bad = "this is not: [valid: yaml: {{{}}}";
        let err = Settings::from_yaml_str(bad).unwrap_err();
        match err {
            PwdftError::Parse(msg) => assert!(msg.contains("YAML")),
            other => panic!("expected Parse error, got: {other:?}"),
        }
    }

    #[test]
    fn missing_required_field_returns_error() {
        let bad = r#"
system:
  atoms: []
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
"#;
        assert!(Settings::from_yaml_str(bad).is_err());
    }

    #[test]
    fn partial_scf_uses_defaults() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
scf:
  max_iter: 50
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.scf.max_iter, 50);
        assert!((s.scf.conv_threshold - 1e-6).abs() < 1e-15);
        assert!(s.scf.n_bands.is_none());
    }

    #[test]
    fn partial_electrons_uses_defaults() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
electrons:
  mixing_beta: 0.7
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert!((s.electrons.mixing_beta - 0.7).abs() < 1e-15);
        assert_eq!(s.electrons.mixing_ndim, 8);
        assert_eq!(s.electrons.smearing, SmearingScheme::FermiDirac);
    }

    #[test]
    fn broyden_mixing_mode_parses() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
electrons:
  mixing_mode: broyden
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.electrons.mixing_mode, MixingModeType::Broyden);
    }

    #[test]
    fn broyden_kerker_mixing_mode_parses() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
electrons:
  mixing_mode: broyden_kerker
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.electrons.mixing_mode, MixingModeType::BroydenKerker);
    }

    #[test]
    fn gaussian_smearing_parses() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
electrons:
  smearing: gaussian
  smearing_width: 0.1
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.electrons.smearing, SmearingScheme::Gaussian);
        assert!((s.electrons.smearing_width - 0.1).abs() < 1e-15);
    }

    #[test]
    fn pbe_functional_parses() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
xc:
  functional: pbe
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.xc.functional, XcFunctional::Pbe);
    }

    #[test]
    fn symmetry_disabled_parses() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
symmetry:
  enabled: false
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert!(!s.symmetry.enabled);
        assert!(s.symmetry.time_reversal);
        assert!((s.symmetry.tolerance - 1e-5).abs() < 1e-15);
    }

    #[test]
    fn fixed_occupations_with_fixed_smearing() {
        let yaml = r#"
system:
  lattice: [[1,0,0],[0,1,0],[0,0,1]]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
electrons:
  smearing: fixed
  occupations: fixed
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert_eq!(s.electrons.smearing, SmearingScheme::Fixed);
        assert_eq!(s.electrons.occupations, OccupationType::Fixed);
    }

    // -----------------------------------------------------------------------
    // ECUT — per-species recommended-ecutwfc defaulting
    // -----------------------------------------------------------------------

    #[test]
    fn resolve_ecutwfc_returns_user_value_when_set() {
        // When YAML sets `basis.ecutwfc`, the resolver returns it verbatim
        // and does not consult the per-element table.
        let s = Settings::from_yaml_str(FULL_YAML).unwrap();
        let crystal = s.to_crystal().unwrap();
        let resolved = s.resolve_ecutwfc(&crystal).unwrap();
        assert!((resolved - 204.09).abs() < f64::EPSILON);
    }

    #[test]
    fn resolve_ecutwfc_falls_back_to_table_for_silicon() {
        // Omit `basis.ecutwfc` from YAML. Si (Z=14) is in the table at
        // 12 Ha ≈ 326.5 eV. Resolver must produce that value.
        use crate::consts::HA_TO_EV;
        let yaml = r#"
system:
  lattice:
    - [0.0, 2.7155, 2.7155]
    - [2.7155, 0.0, 2.7155]
    - [2.7155, 2.7155, 0.0]
  atoms:
    - symbol: Si
      position: [0.0, 0.0, 0.0]
    - symbol: Si
      position: [0.25, 0.25, 0.25]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        assert!(s.basis.ecutwfc.is_none(), "ecutwfc must be None when omitted");
        let crystal = s.to_crystal().unwrap();
        let resolved = s.resolve_ecutwfc(&crystal).unwrap();
        let expected = 12.0 * HA_TO_EV;
        assert!(
            (resolved - expected).abs() < 1e-9,
            "Si default ecutwfc {resolved} != expected {expected}"
        );
    }

    #[test]
    fn resolve_ecutwfc_takes_max_across_species_iron_wins_over_oxygen() {
        // Mixed Fe + O crystal: Fe is 30 Ha (~816 eV), O is 24 Ha
        // (~653 eV). The resolver must pick Fe's higher cutoff.
        use crate::consts::HA_TO_EV;
        let yaml = r#"
system:
  lattice:
    - [4.3, 0.0, 0.0]
    - [0.0, 4.3, 0.0]
    - [0.0, 0.0, 4.3]
  atoms:
    - symbol: Fe
      position: [0.0, 0.0, 0.0]
    - symbol: O
      position: [0.5, 0.5, 0.5]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        let crystal = s.to_crystal().unwrap();
        let resolved = s.resolve_ecutwfc(&crystal).unwrap();
        let fe_ev = 30.0 * HA_TO_EV;
        assert!(
            (resolved - fe_ev).abs() < 1e-9,
            "Fe/O default ecutwfc {resolved} must match Fe {fe_ev}, not O {}",
            24.0 * HA_TO_EV,
        );
    }

    #[test]
    fn resolve_ecutwfc_returns_invalid_param_when_table_has_no_match() {
        // ERR2 P1.c: when `basis.ecutwfc` is unset and every species in
        // the crystal is outside the recommended-ecut table, the caller
        // must receive `PwdftError::InvalidParam { name: "basis.ecutwfc",
        // .. }` rather than the catch-all `InvalidInput`. The table
        // lookup returns `None` for synthetic Z=200, so we build a
        // Crystal directly rather than round-tripping through YAML
        // (the YAML element parser rejects unknown symbols upstream).
        use crate::crystal::{Atom, Crystal, Lattice};
        use nalgebra::Vector3;
        let lat = Lattice::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(0.0, 5.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        );
        let crystal = Crystal {
            atoms: vec![Atom::new(200, [0.0, 0.0, 0.0])],
            lattice: lat,
        };
        // A minimal Settings with no `basis.ecutwfc` set. We construct
        // via a YAML hop that omits `basis`.
        let yaml = r#"
system:
  lattice:
    - [5.0, 0.0, 0.0]
    - [0.0, 5.0, 0.0]
    - [0.0, 0.0, 5.0]
  atoms:
    - symbol: Si
      position: [0.0, 0.0, 0.0]
kpoints:
  type: monkhorst_pack
  grid: [1, 1, 1]
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        // Hand-replace the crystal with the synthetic Z=200 one so the
        // table lookup falls through.
        let err = s.resolve_ecutwfc(&crystal).expect_err(
            "resolve_ecutwfc must fail when no species has a recommended cutoff",
        );
        match err {
            PwdftError::InvalidParam { name, reason } => {
                assert_eq!(name, "basis.ecutwfc");
                assert!(
                    reason.contains("tabulated"),
                    "reason should mention the tabulated lookup, got: {reason}"
                );
            }
            other => panic!("expected InvalidParam, got: {other:?}"),
        }
    }

    #[test]
    fn resolve_ecutwfc_defaults_to_carbon_for_diamond() {
        // Diamond (C, Z=6) is harder than Si. .standard = 18 Ha.
        use crate::consts::HA_TO_EV;
        let yaml = r#"
system:
  lattice:
    - [0.0, 1.784, 1.784]
    - [1.784, 0.0, 1.784]
    - [1.784, 1.784, 0.0]
  atoms:
    - symbol: C
      position: [0.0, 0.0, 0.0]
    - symbol: C
      position: [0.25, 0.25, 0.25]
kpoints:
  type: monkhorst_pack
  grid: [2, 2, 2]
"#;
        let s = Settings::from_yaml_str(yaml).unwrap();
        let crystal = s.to_crystal().unwrap();
        let resolved = s.resolve_ecutwfc(&crystal).unwrap();
        let expected = 18.0 * HA_TO_EV;
        assert!(
            (resolved - expected).abs() < 1e-9,
            "C default ecutwfc {resolved} != expected {expected}"
        );
    }
}
