//! SCF calculation context — immutable state shared across iterations.
//!
//! `ScfContext` holds everything computed once at setup time: crystal,
//! basis, grid, pseudopotentials, cached V_NL, Ewald energy, etc.
//! The SCF loop methods take `&self` plus mutable density state.

use std::collections::HashMap;

use elements_rs::Element;
use log::info;
use num_complex::Complex64;
use rayon::prelude::*;

use super::{ScfParams, grid::FftGrid, potentials};
use crate::{
    basis::BasisSet,
    crystal::Crystal,
    error::Result,
    fft,
    kpoints::KPoint,
    potential::{nonlocal::NonlocalPotential, xc::XcEvaluator},
    pseudopotential::UpfPseudoPotential,
    symmetry::SymmetryInfo,
};

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
    pub pseudopotentials: &'a HashMap<Element, UpfPseudoPotential>,
    pub params: &'a ScfParams,
    pub symmetry: &'a SymmetryInfo,

    // Precomputed (immutable across iterations)
    pub grid: FftGrid,
    pub g_to_fft: Vec<usize>,
    pub g_squared: Vec<f64>,
    /// Per-FFT-grid-point reciprocal-space vectors `G` in Å⁻¹, in the
    /// FFT-aligned ordering of `scf::grid::g_vector_at_dims`. Cached so
    /// GGA density-gradient evaluators (GGAP Phase A.1) don't rebuild
    /// the Vec<[f64;3]> every SCF iteration.
    pub g_vectors: Vec<[f64; 3]>,
    pub v_local_fft: Vec<Complex64>,
    /// Diagnostic copy of `V_local(G=0)` in eV — the uniform-background
    /// DC offset carried by every Kohn-Sham eigenvalue under the
    /// QE-compatible gauge. Logged at `ScfContext::new` time and kept
    /// on the struct for downstream diagnostic callers; the SCF loop
    /// itself does not read it.
    #[allow(
        dead_code,
        reason = "diagnostic-only field; consumed via info! log at construction site"
    )]
    pub v_local_g0: f64,
    pub vnl_cache: Vec<NonlocalPotential>,
    pub rho_core_r: Vec<f64>,
    /// Pre-FFT'd ∇ρ_core on the FFT grid (length-3 arrays). `Some`
    /// whenever `rho_core_r` is non-empty *and* the active XC
    /// functional needs a gradient (GGA); `None` otherwise. Cached at
    /// construction time because the core density is fixed per geometry.
    pub rho_core_grad_r: Option<Vec<[f64; 3]>>,
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
    /// Returns `PwdftError::MissingPseudopotential` if any atom lacks a loaded
    /// PP.
    pub fn new(
        crystal: &'a Crystal,
        basis: &'a BasisSet,
        kpoints: &'a [KPoint],
        pseudopotentials: &'a HashMap<Element, UpfPseudoPotential>,
        params: &'a ScfParams,
        symmetry: &'a SymmetryInfo,
    ) -> Result<Self> {
        let omega = crystal.lattice.volume();
        let n_electrons: f64 = crystal
            .atoms
            .iter()
            .map(|a| pseudopotentials.get(&a.symbol).expect("BUG: Missing PP").z_valence)
            .sum();

        info!(
            "SCF: {n_electrons} electrons, {omega:.3} ų cell volume, nspin={}",
            params.nspin
        );

        let mut grid = FftGrid::new(basis, &crystal.lattice, params.ecutrho_ratio, params.fft_grid);
        let n_grid = grid.total_size();
        let [nx, ny, nz] = grid.dims;
        info!("FFT grid: {nx}×{ny}×{nz} = {n_grid} points");

        let g_to_fft = grid.basis_to_fft(basis);

        // V_local in G-space. The G=0 component is kept on the Hamiltonian
        // diagonal (via `v_local_fft[0]`) so every Kohn-Sham eigenvalue
        // carries the uniform-background DC offset, matching QE's
        // convention (`qe-7.5/PW/src/setlocal.f90:91-96`: `v_of_0 =
        // DBLE(aux(1))` is recorded for diagnostics but `aux(1)` stays on
        // `vltot(r)`). The `v_local_g0` stash is retained for logging and
        // diagnostic scripts; it is no longer subtracted out of the
        // Hamiltonian and re-added to the total energy.
        let v_local_fft = potentials::compute_v_local(crystal, &grid, pseudopotentials, omega)?;
        let v_local_g0 = v_local_fft[0].re;
        info!("V_local(G=0) = {v_local_g0:.6} eV (on Hamiltonian diagonal)");

        // Precompute |G|² and the full G-vector cache. `g_vectors` is
        // reused by GGA gradient / divergence FFTs (GGAP Phase A.1);
        // building it once avoids rebuilding every SCF iteration.
        let dims = grid.dims;
        let recip = grid.recip.clone();
        let g_vectors: Vec<[f64; 3]> = (0..n_grid)
            .into_par_iter()
            .map(|idx| {
                let g = super::grid::g_vector_at_dims(idx, dims, &recip);
                [g.x, g.y, g.z]
            })
            .collect();
        let g_squared: Vec<f64> = g_vectors
            .par_iter()
            .map(|g| g[0] * g[0] + g[1] * g[1] + g[2] * g[2])
            .collect();

        // NLCC core density
        let rho_core_r = potentials::compute_core_density(crystal, &mut grid, pseudopotentials);
        if !rho_core_r.is_empty() {
            let core_min = rho_core_r.iter().copied().fold(f64::INFINITY, f64::min);
            let core_max = rho_core_r.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            info!("NLCC core density: min={core_min:.4e} max={core_max:.4e}");
        }

        // Precompute ∇ρ_core once per calculation (it is time-independent
        // across the SCF loop). Only populate when the active functional
        // needs gradients (GGA) and NLCC is active. For LDA this stays
        // `None` and the LDA path pays zero FFT cost.
        let xc_needs_gradient = XcEvaluator::from_settings(params.xc_functional)
            .map(|xc| xc.needs_gradient())
            .unwrap_or(false);
        let rho_core_grad_r: Option<Vec<[f64; 3]>> = if xc_needs_gradient && !rho_core_r.is_empty() {
            Some(fft::compute_density_gradient(&rho_core_r, &mut grid.fft, &g_vectors))
        } else {
            None
        };

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
            g_vectors,
            v_local_fft,
            v_local_g0,
            vnl_cache,
            rho_core_r,
            rho_core_grad_r,
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
