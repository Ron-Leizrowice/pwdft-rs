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

/// Per-calculation SCF context: immutable physics inputs plus
/// preallocated scratch that is reused across SCF iterations.
///
/// Created once before the SCF loop. Most fields are genuinely
/// immutable (crystal, basis, precomputed V_local, VNL cache, Ewald
/// energy, etc.). [`Self::h_scratch`] is the one exception: it is a
/// caller-owned per-(spin, k-point) `faer::Mat` scratch buffer that
/// the driver mutably borrows each iteration to assemble the
/// kinetic + V_eff + V_NL Hamiltonian in place. Callers
/// hold `&mut ScfContext` through the SCF loop, which lets them
/// `par_iter_mut()` over `h_scratch` while immutably borrowing the
/// rest of the context.
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

    /// Per-(spin, k-point) scratch Hamiltonian matrices, length
    /// `params.nspin * kpoints.len()`. Indexed as
    /// `ispin * n_k + ik` (row-major `[spin][k]`). Each entry is an
    /// `n_pw × n_pw` `faer::Mat<Complex64>` that the driver
    /// **fully overwrites** (via
    /// [`super::potentials::fill_hamiltonian_with_v_eff`]) before
    /// the non-local KB term accumulates on top, so no zero-fill
    /// is needed between iterations. Resident footprint:
    /// `nspin · n_k · n_pw² · 16` bytes (84 MB for nspin=1,
    /// n_k=10, n_pw=725 — the production Si 4×4×4 @ ecut=400 point).
    /// Replaces the per-iteration `Mat::zeros(n, n)` that previously
    /// burned the same 84 MB as transient allocations every SCF step.
    pub h_scratch: Vec<faer::Mat<Complex64>>,
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

        // Allocate per-(spin, k-point) Hamiltonian scratch (ALOC F-5).
        // One n_pw × n_pw Complex64 Mat per spin channel per k-point,
        // zero-initialised; `fill_hamiltonian_with_v_eff` fully
        // overwrites every entry before use, so the zero-init is just
        // the default for new faer matrices (no per-iteration clear).
        let n_pw = basis.len();
        let nspin = params.nspin;
        let n_k = kpoints.len();
        let h_scratch: Vec<faer::Mat<Complex64>> = (0..nspin * n_k)
            .map(|_| faer::Mat::<Complex64>::zeros(n_pw, n_pw))
            .collect();
        let mb_resident = (nspin * n_k * n_pw * n_pw * 16) as f64 / (1024.0 * 1024.0);
        info!(
            "Hamiltonian scratch: {} × {n_pw}² Mat<Complex64> ({mb_resident:.1} MB resident)",
            nspin * n_k,
        );

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
            h_scratch,
        })
    }
}
