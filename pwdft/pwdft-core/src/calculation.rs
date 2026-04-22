use std::collections::HashMap;

use elements_rs::Element;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    error::PwdftError,
    kpoints::KPointSet,
    pseudopotential::UpfPseudoPotential,
    scf::ScfParams,
    settings::{InputSettings, PredefinedXcFunctionals},
    symmetry::SymmetryInfo,
};

pub enum Calculation {
    Scf(ScfCalculation),
    BandStructure(BandStructureCalculation),
}

pub struct ScfCalculation {
    pub crystal: Crystal,
    pub kpoints: KPointSet,
    pub symmetry: SymmetryInfo,
    pub xc_functional: PredefinedXcFunctionals,
    pub pseudopotentials: HashMap<Element, UpfPseudoPotential>,
    pub basis: BasisSet,
    /// IBZ-reduced Monkhorst-Pack grid. Provenance
    /// (`SamplingKind::Irreducible`) carries the parent full-grid
    /// dimensions + shift for density symmetrization.
    pub scf_settings: ScfParams,
}

pub struct BandStructureCalculation {
    pub crystal: Crystal,
    pub basis: BasisSet,
    /// Band path. Provenance is `SamplingKind::BandPath { distances }`;
    /// pull the distances via `kpoints.distances().unwrap()` when
    /// writing the TSV.
    pub kpoints: KPointSet,
    pub n_bands: usize,
}

impl TryFrom<InputSettings> for ScfCalculation {
    type Error = PwdftError;

    fn try_from(input: InputSettings) -> Result<Self, Self::Error> {
        let crystal = Crystal::from(input.system);
        let basis_set = BasisSet::new(&crystal.lattice, input.basis.ecutwfc);
    }
}

impl ScfCalculation {
    fn load_pseudopotentials(&self) -> Result<(), PwdftError> {}
}
