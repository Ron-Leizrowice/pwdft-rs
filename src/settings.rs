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
    kpoints::HighSymPoint,
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
    pub ecutwfc: f64,
    /// Charge-density cutoff ratio: ecutrho = ecutrho_ratio * ecutwfc.
    /// QE default for norm-conserving PPs is 4.
    pub ecutrho_ratio: u32,
    /// Explicit FFT grid dimensions [n1, n2, n3]. If set, overrides ecutrho_ratio.
    pub fft_grid: Option<[usize; 3]>,
}

impl Default for BasisSettings {
    fn default() -> Self {
        Self {
            ecutwfc: 204.09,
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
    /// Maximum number of SCF iterations (QE: electron_maxstep).
    pub max_iter: usize,
    /// Convergence threshold: RMS density change (e/Å³).
    pub conv_threshold: f64,
    /// Energy convergence threshold (eV). Both density AND energy must converge.
    pub energy_threshold: f64,
    /// Number of Kohn-Sham bands. `None` = automatic from n_electrons/2 + padding.
    pub n_bands: Option<usize>,
    /// Eigensolver backend for per-k-point diagonalization.
    pub eigensolver: EigensolverType,
}

impl Default for ScfSettings {
    fn default() -> Self {
        Self {
            max_iter: 100,
            conv_threshold: 1e-6,
            energy_threshold: 1e-5,
            n_bands: None,
            eigensolver: EigensolverType::default(),
        }
    }
}

/// Eigensolver backend selection (YAML-friendly adapter for
/// [`crate::eigensolver::EigensolverKind`]).
///
/// The default is `Dense` for compatibility; flip to `Iterative` to
/// activate the ITEV partial Arnoldi path (typically 3-10× faster at
/// n_pw ≥ 200).
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
    /// Density mixing parameter (0 < beta <= 1). QE: mixing_beta.
    pub mixing_beta: f64,
    /// Number of past densities kept for Anderson/Pulay mixing. QE: mixing_ndim.
    pub mixing_ndim: usize,
    /// Smearing scheme for partial occupations.
    pub smearing: SmearingScheme,
    /// Smearing width in eV (QE: degauss, but QE uses Ry internally).
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
    /// consecutive iterations. Default `false` preserves the pre-MXBA
    /// fixed-β behaviour for existing inputs. See
    /// `src/scf/mixing/mod.rs` module docs for the rule and thresholds.
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
    /// Functional name. Currently supported: "pz" (LDA).
    /// Future: "pbe" (GGA), "pbe0", "hse06".
    pub functional: XcFunctional,
}

/// Supported exchange-correlation functionals.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum XcFunctional {
    /// Perdew-Zunger LDA (Ceperley-Alder). QE: input_dft = 'PZ'.
    #[default]
    Pz,
    /// Perdew-Burke-Ernzerhof GGA. QE: input_dft = 'PBE'.
    Pbe,
    /// PBE0 hybrid functional.
    Pbe0,
    /// Heyd-Scuseria-Ernzerhof screened hybrid.
    Hse06,
}

/// Crystal symmetry settings.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct SymmetrySettings {
    /// Whether to detect and exploit crystal symmetry. QE: nosym = .false.
    pub enabled: bool,
    /// Whether to apply time-reversal symmetry (k -> -k). QE: noinv = .false.
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
    pub fn from_yaml_str(s: &str) -> Result<Self> {
        serde_yaml_ng::from_str(s).map_err(|e| PwdftError::Parse(format!("YAML parse error: {e}")))
    }

    /// Parse settings from a YAML file on disk.
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
                let elem = crate::atoms::from_symbol(&ai.symbol).ok_or_else(|| {
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
        }
    }

    /// Extract the wavefunction energy cutoff in eV.
    #[must_use]
    pub fn ecutwfc(&self) -> f64 {
        self.basis.ecutwfc
    }

    /// Extract the Monkhorst-Pack grid dimensions, if configured.
    #[must_use]
    pub fn mp_grid(&self) -> Option<[u32; 3]> {
        match &self.kpoints {
            KPointSettings::MonkhorstPack { grid } => Some(*grid),
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
        assert!((s.basis.ecutwfc - 204.09).abs() < f64::EPSILON);
        assert!(s.mp_grid().is_some());
        assert_eq!(s.mp_grid().unwrap(), [4, 4, 4]);
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

    #[test]
    fn full_settings_roundtrip() {
        let original = Settings::from_yaml_str(FULL_YAML).unwrap();
        let serialized = serde_yaml_ng::to_string(&original).unwrap();
        let restored = Settings::from_yaml_str(&serialized).unwrap();

        assert_eq!(original.system.atoms.len(), restored.system.atoms.len());
        assert!((original.basis.ecutwfc - restored.basis.ecutwfc).abs() < f64::EPSILON);
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
}
