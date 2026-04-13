use nalgebra::Vector3;
use serde::Deserialize;

use crate::{
    crystal::{Atom, Crystal, Lattice},
    error::{PwdftError, Result},
    kpoints::HighSymPoint,
};

/// Top-level input file structure.
#[derive(Debug, Deserialize)]
pub struct InputFile {
    pub system: SystemConfig,
    pub kpoints: KPointsConfig,
}

#[derive(Debug, Deserialize)]
pub struct SystemConfig {
    /// Lattice vectors as [[ax,ay,az],[bx,by,bz],[cx,cy,cz]].
    pub lattice: [[f64; 3]; 3],
    /// Atoms: each has symbol and fractional position.
    #[serde(default)]
    pub atoms: Vec<AtomInput>,
    /// Plane-wave energy cutoff in eV.
    pub ecut: f64,
    /// Number of bands to compute (default: auto from n_electrons/2 + padding).
    pub n_bands: Option<usize>,
}

#[derive(Debug, Deserialize)]
pub struct AtomInput {
    pub symbol: String,
    pub position: [f64; 3],
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type")]
pub enum KPointsConfig {
    #[serde(rename = "monkhorst_pack")]
    MonkhorstPack { grid: [u32; 3] },
    #[serde(rename = "band_path")]
    BandPath {
        /// Path segments as list of {label, frac} entries.
        points: Vec<PathPoint>,
        /// Points per segment.
        npoints: usize,
    },
}

#[derive(Debug, Deserialize)]
pub struct PathPoint {
    pub label: String,
    pub frac: [f64; 3],
}

impl InputFile {
    pub fn from_str(s: &str) -> Result<Self> {
        toml::from_str(s).map_err(|e| PwdftError::Parse(e.to_string()))
    }

    pub fn from_file(path: &std::path::Path) -> Result<Self> {
        let contents = std::fs::read_to_string(path)?;
        Self::from_str(&contents)
    }

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

    pub fn to_high_sym_path(&self) -> Option<Vec<HighSymPoint>> {
        match &self.kpoints {
            KPointsConfig::BandPath { points, .. } => Some(
                points
                    .iter()
                    .map(|p| HighSymPoint {
                        label: p.label.clone(),
                        frac: p.frac,
                    })
                    .collect(),
            ),
            _ => None,
        }
    }

    pub fn band_path_npoints(&self) -> Option<usize> {
        match &self.kpoints {
            KPointsConfig::BandPath { npoints, .. } => Some(*npoints),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_band_path_input() {
        let input = r#"
[system]
lattice = [
    [0.0, 2.7155, 2.7155],
    [2.7155, 0.0, 2.7155],
    [2.7155, 2.7155, 0.0],
]
ecut = 200.0
n_bands = 8

[[system.atoms]]
symbol = "Si"
position = [0.0, 0.0, 0.0]

[[system.atoms]]
symbol = "Si"
position = [0.25, 0.25, 0.25]

[kpoints]
type = "band_path"
npoints = 50

[[kpoints.points]]
label = "Γ"
frac = [0.0, 0.0, 0.0]

[[kpoints.points]]
label = "X"
frac = [0.5, 0.0, 0.5]

[[kpoints.points]]
label = "Γ"
frac = [0.0, 0.0, 0.0]

[[kpoints.points]]
label = "L"
frac = [0.5, 0.5, 0.5]
"#;
        let config = InputFile::from_str(input).unwrap();
        assert_eq!(config.system.ecut, 200.0);
        assert_eq!(config.system.atoms.len(), 2);
        assert_eq!(config.system.n_bands, Some(8));

        let crystal = config.to_crystal();
        assert_eq!(crystal.atoms.len(), 2);

        let path = config.to_high_sym_path().unwrap();
        assert_eq!(path.len(), 4);
        assert_eq!(path[0].label, "Γ");
        assert_eq!(path[1].label, "X");
    }
}
