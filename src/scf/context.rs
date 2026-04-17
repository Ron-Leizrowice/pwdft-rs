//! SCF calculation context — immutable state shared across iterations.
//!
//! `ScfContext` holds everything computed once at setup time: crystal,
//! basis, grid, pseudopotentials, cached V_NL, Ewald energy, etc.
//! The SCF loop methods take `&self` plus mutable density state.

use log::info;
use num_complex::Complex64;
use rayon::prelude::*;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    error::{PwdftError, Result},
    kpoints::KPoint,
    potential::nonlocal::NonlocalPotential,
    pseudopotential::PseudopotentialData,
};

use super::grid::FftGrid;
use super::potentials;
use super::ScfParams;

/// Immutable context for an SCF calculation.
///
/// Created once before the SCF loop. Methods compute potentials,
/// build Hamiltonians, and evaluate energies without owning
/// the mutable density state.
pub(crate) struct ScfContext<'a> {
    pub crystal: &'a Crystal,
    pub basis: &'a BasisSet,
    pub kpoints: &'a [KPoint],
    pub pseudopotentials: &'a [&'a PseudopotentialData],
    pub params: &'a ScfParams,
    pub symmetry: &'a crate::symmetry::SymmetryInfo,

    // Precomputed (immutable across iterations)
    pub grid: FftGrid,
    pub g_to_fft: Vec<usize>,
    pub g_squared: Vec<f64>,
    pub v_local_fft: Vec<Complex64>,
    pub v_local_g0: f64,
    pub vnl_cache: Vec<NonlocalPotential>,
    pub rho_core_r: Vec<f64>,
    pub e_ewald: f64,
    pub omega: f64,
    pub n_electrons: f64,
    pub n_grid: usize,
    pub kpt_weights: Vec<f64>,
    pub spin_factor: f64,
}

impl<'a> ScfContext<'a> {
    /// Build the SCF context: precompute everything that doesn't change
    /// between iterations.
    ///
    /// # Errors
    /// Returns `PwdftError::MissingPseudopotential` if any atom lacks a loaded PP.
    pub fn new(
        crystal: &'a Crystal,
        basis: &'a BasisSet,
        kpoints: &'a [KPoint],
        pseudopotentials: &'a [&'a PseudopotentialData],
        params: &'a ScfParams,
        symmetry: &'a crate::symmetry::SymmetryInfo,
    ) -> Result<Self> {
        let omega = crystal.lattice.volume();
        let n_electrons: f64 = crystal
            .atoms
            .iter()
            .map(|a| {
                crate::pseudopotential::find_for_atom(a.z, pseudopotentials)
                    .map(|pp| pp.z_valence)
                    .ok_or_else(|| PwdftError::MissingPseudopotential(
                        format!("Z={} not found in loaded pseudopotentials", a.z)
                    ))
            })
            .collect::<Result<Vec<f64>>>()?
            .into_iter()
            .sum();

        info!("SCF: {n_electrons} electrons, {omega:.3} ų cell volume, nspin={}", params.nspin);

        let mut grid = FftGrid::new(basis, &crystal.lattice, params.ecutrho_ratio, params.fft_grid);
        let n_grid = grid.total_size();
        let [nx, ny, nz] = grid.dims;
        info!("FFT grid: {nx}×{ny}×{nz} = {n_grid} points");

        let g_to_fft = grid.basis_to_fft(basis);

        // V_local with G=0 excluded
        let mut v_local_fft = potentials::compute_v_local(crystal, &grid, pseudopotentials, omega)?;
        let v_local_g0 = v_local_fft[0].re;
        v_local_fft[0] = Complex64::new(0.0, 0.0);
        info!("V_local(G=0) = {v_local_g0:.6} eV (excluded from Hamiltonian)");

        // Precompute |G|²
        let dims = grid.dims;
        let recip = grid.recip.clone();
        let g_squared: Vec<f64> = (0..n_grid)
            .into_par_iter()
            .map(|idx| super::grid::g_vector_at_dims(idx, dims, &recip).norm_squared())
            .collect();

        // NLCC core density
        let rho_core_r = potentials::compute_core_density(crystal, &mut grid, pseudopotentials);
        if !rho_core_r.is_empty() {
            let core_min = rho_core_r.iter().copied().fold(f64::INFINITY, f64::min);
            let core_max = rho_core_r.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            info!("NLCC core density: min={core_min:.4e} max={core_max:.4e}");
        }

        // Cache V_NL per k-point
        let vnl_cache: Result<Vec<NonlocalPotential>> = kpoints
            .par_iter()
            .map(|kp| NonlocalPotential::new(crystal, basis, &kp.k, pseudopotentials))
            .collect();
        let vnl_cache = vnl_cache?;

        let e_ewald = crate::ewald::ewald_energy(crystal, pseudopotentials);
        let kpt_weights: Vec<f64> = kpoints.iter().map(|kp| kp.weight).collect();
        let spin_factor = 2.0 / params.nspin as f64;

        Ok(Self {
            crystal,
            basis,
            kpoints,
            pseudopotentials,
            params,
            symmetry,
            grid,
            g_to_fft,
            g_squared,
            v_local_fft,
            v_local_g0,
            vnl_cache,
            rho_core_r,
            e_ewald,
            omega,
            n_electrons,
            n_grid,
            kpt_weights,
            spin_factor,
        })
    }
}
