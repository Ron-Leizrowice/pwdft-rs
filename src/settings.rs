//! Comprehensive settings/configuration for plane-wave DFT calculations.
//!
//! Parsed from YAML files via `serde_yaml_ng`. This is a parallel configuration path
//! alongside the existing TOML-based `input.rs` -- both remain fully functional.
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
    scf::ScfParams,
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
    /// Lattice vectors in Angstroms: [[ax,ay,az],[bx,by,bz],[cx,cy,cz]].
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
}

impl Default for BasisSettings {
    fn default() -> Self {
        Self {
            ecutwfc: 204.09,
            ecutrho_ratio: 4,
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
        /// Number of k-points per segment.
        #[serde(default = "default_band_npoints")]
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
    /// Number of Kohn-Sham bands. `None` = automatic from n_electrons/2 + padding.
    pub n_bands: Option<usize>,
}

impl Default for ScfSettings {
    fn default() -> Self {
        Self {
            max_iter: 100,
            conv_threshold: 1e-6,
            n_bands: None,
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
    pub smearing: SmearingType,
    /// Smearing width in eV (QE: degauss, but QE uses Ry internally).
    pub smearing_width: f64,
    /// Occupation scheme.
    pub occupations: OccupationType,
}

impl Default for ElectronSettings {
    fn default() -> Self {
        Self {
            mixing_beta: 0.3,
            mixing_ndim: 8,
            smearing: SmearingType::default(),
            smearing_width: 0.05,
            occupations: OccupationType::default(),
        }
    }
}

/// Smearing function for partial occupation numbers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum SmearingType {
    /// Fermi-Dirac (finite-temperature) smearing. QE: 'fermi-dirac' / 'f-d'.
    #[default]
    FermiDirac,
    /// Gaussian smearing. QE: 'gaussian' / 'gauss'.
    Gaussian,
    /// Methfessel-Paxton first-order smearing. QE: 'methfessel-paxton' / 'm-p'.
    MethfesselPaxton,
    /// Marzari-Vanderbilt-DeVita-Payne cold smearing. QE: 'cold' / 'm-v'.
    Cold,
    /// Fixed occupations (no smearing, insulator mode).
    Fixed,
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

fn default_band_npoints() -> usize {
    50
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
    pub fn to_crystal(&self) -> Crystal {
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
                let elem = crate::atoms::Element::from_symbol(&ai.symbol)
                    .unwrap_or_else(|| panic!("unknown element: {}", ai.symbol));
                Atom::new(elem.atomic_number(), ai.position)
            })
            .collect();

        Crystal { atoms, lattice }
    }

    /// Build `ScfParams` from the SCF + electron + basis settings.
    ///
    /// `n_bands_fallback` is used when `scf.n_bands` is `None` (auto mode).
    pub fn to_scf_params(&self, n_bands_fallback: usize) -> ScfParams {
        ScfParams {
            n_bands: self.scf.n_bands.unwrap_or(n_bands_fallback),
            max_iter: self.scf.max_iter,
            conv_threshold: self.scf.conv_threshold,
            mixing_beta: self.electrons.mixing_beta,
            mixing_ndim: self.electrons.mixing_ndim,
            smearing_sigma: self.electrons.smearing_width,
            ecutrho_ratio: self.basis.ecutrho_ratio,
            fft_grid: None,
            ..Default::default()
        }
    }

    /// Extract the wavefunction energy cutoff in eV.
    pub fn ecutwfc(&self) -> f64 {
        self.basis.ecutwfc
    }

    /// Extract the Monkhorst-Pack grid dimensions, if configured.
    pub fn mp_grid(&self) -> Option<[u32; 3]> {
        match &self.kpoints {
            KPointSettings::MonkhorstPack { grid } => Some(*grid),
            _ => None,
        }
    }

    /// Extract the high-symmetry path for band-structure calculations.
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
    pub fn band_path_npoints(&self) -> Option<usize> {
        match &self.kpoints {
            KPointSettings::BandPath { npoints, .. } => Some(*npoints),
            _ => None,
        }
    }

    /// Build `SymmetryInfo` from these settings and a crystal.
    pub fn to_symmetry_info(
        &self,
        crystal: &Crystal,
    ) -> Option<crate::symmetry::SymmetryInfo> {
        if !self.symmetry.enabled {
            return None;
        }
        let mut info =
            crate::symmetry::SymmetryInfo::from_crystal(crystal, self.symmetry.tolerance);
        info.has_time_reversal = self.symmetry.time_reversal;
        Some(info)
    }

    /// Return the pseudopotential file path for a given element symbol.
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

kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]

scf:
  max_iter: 100
  conv_threshold: 1.0e-6
  n_bands: 8

electrons:
  mixing_beta: 0.3
  mixing_ndim: 8
  smearing: fermi_dirac
  smearing_width: 0.05
  occupations: smearing

xc:
  functional: pz

symmetry:
  enabled: true
  time_reversal: true
  tolerance: 1.0e-5

pseudopotentials:
  Si: "../pseudopotentials/Si.UPF"

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
        assert_eq!(s.basis.ecutwfc, 204.09);
        assert!(s.mp_grid().is_some());
        assert_eq!(s.mp_grid().unwrap(), [4, 4, 4]);
    }

    #[test]
    fn parse_full_yaml() {
        let s = Settings::from_yaml_str(FULL_YAML).unwrap();
        assert_eq!(s.scf.max_iter, 100);
        assert_eq!(s.scf.n_bands, Some(8));
        assert_eq!(s.electrons.mixing_beta, 0.3);
        assert_eq!(s.electrons.mixing_ndim, 8);
        assert_eq!(s.electrons.smearing, SmearingType::FermiDirac);
        assert!((s.electrons.smearing_width - 0.05).abs() < 1e-15);
        assert_eq!(s.electrons.occupations, OccupationType::Smearing);
        assert_eq!(s.xc.functional, XcFunctional::Pz);
        assert!(s.symmetry.enabled);
        assert!(s.symmetry.time_reversal);
        assert!((s.symmetry.tolerance - 1e-5).abs() < 1e-15);
        assert_eq!(
            s.pseudopotentials.files.get("Si").unwrap(),
            "../pseudopotentials/Si.UPF"
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
        assert_eq!(s.scf.max_iter, 100);
        assert!((s.scf.conv_threshold - 1e-6).abs() < 1e-15);
        assert!(s.scf.n_bands.is_none());
        assert!((s.electrons.mixing_beta - 0.3).abs() < 1e-15);
        assert_eq!(s.electrons.mixing_ndim, 8);
        assert_eq!(s.electrons.smearing, SmearingType::FermiDirac);
        assert!((s.electrons.smearing_width - 0.05).abs() < 1e-15);
        assert_eq!(s.electrons.occupations, OccupationType::Smearing);
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
        let crystal = s.to_crystal();
        assert_eq!(crystal.atoms.len(), 2);
        assert_eq!(crystal.atoms[0].z, 14);
        assert_eq!(crystal.atoms[1].z, 14);
        assert_eq!(crystal.atoms[1].position, [0.25, 0.25, 0.25]);
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
            Some("../pseudopotentials/Si.UPF")
        );
        assert_eq!(s.pseudopotential_path("Ge"), None);
    }

    #[test]
    fn smearing_types_roundtrip() {
        for variant in [
            SmearingType::FermiDirac,
            SmearingType::Gaussian,
            SmearingType::MethfesselPaxton,
            SmearingType::Cold,
            SmearingType::Fixed,
        ] {
            let yaml = serde_yaml_ng::to_string(&variant).unwrap();
            let parsed: SmearingType = serde_yaml_ng::from_str(&yaml).unwrap();
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
    fn full_settings_roundtrip() {
        let original = Settings::from_yaml_str(FULL_YAML).unwrap();
        let serialized = serde_yaml_ng::to_string(&original).unwrap();
        let restored = Settings::from_yaml_str(&serialized).unwrap();

        assert_eq!(original.system.atoms.len(), restored.system.atoms.len());
        assert_eq!(original.basis.ecutwfc, restored.basis.ecutwfc);
        assert_eq!(original.scf.max_iter, restored.scf.max_iter);
        assert_eq!(original.scf.n_bands, restored.scf.n_bands);
        assert_eq!(original.electrons.mixing_beta, restored.electrons.mixing_beta);
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
        assert_eq!(s.electrons.smearing, SmearingType::FermiDirac);
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
        assert_eq!(s.electrons.smearing, SmearingType::Gaussian);
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
        assert_eq!(s.electrons.smearing, SmearingType::Fixed);
        assert_eq!(s.electrons.occupations, OccupationType::Fixed);
    }
}
