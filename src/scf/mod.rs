pub(crate) mod context;
pub mod density;
pub(crate) mod energy;
pub(crate) mod grid;
pub mod initial_density;
pub mod mixing;
pub(crate) mod potentials;
pub mod smearing;

use log::info;
use num_complex::Complex64;
use rayon::prelude::*;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    eigensolver::dense,
    error::{PwdftError, Result},
    kpoints::KPoint,
    potential::{nonlocal::NonlocalPotential, xc},
    pseudopotential::PseudopotentialData,
};

use self::energy::{
    add_core_density, assemble_v_eff, band_energy, density_diff,
    density_r_to_g, hartree_energy, hartree_on_fft_grid, real_to_g_space,
    total_energy, xc_energy_corrected,
};
use self::grid::{FftGrid, g_vector_at_dims};
use self::potentials::{build_hamiltonian_with_v_eff, compute_core_density, compute_v_local};

/// Parameters for an SCF calculation.
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
    /// Number of spin channels: 1 (unpolarized) or 2 (collinear spin-polarized).
    pub nspin: usize,
    /// Starting magnetization per atom type (fractional, -1 to 1).
    /// Maps from element symbol to magnetization. Empty = non-magnetic.
    pub starting_magnetization: std::collections::HashMap<String, f64>,
    /// Fixed total magnetization (n_up - n_down) in electrons.
    /// If None, magnetization is determined self-consistently.
    pub tot_magnetization: Option<f64>,
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
            nspin: 1,
            starting_magnetization: std::collections::HashMap::new(),
            tot_magnetization: None,
        }
    }
}

/// Result of an SCF calculation.
pub struct ScfResult {
    /// Kohn-Sham total energy (no entropy).
    pub total_energy: f64,
    /// Free energy F = E - TS (Mermin functional, variational quantity).
    pub free_energy: f64,
    /// Sigma→0 extrapolated energy E₀ = (E + F) / 2.
    pub energy_sigma0: f64,
    /// Entropy contribution T*S in eV.
    pub entropy_ts: f64,
    /// Eigenvalues indexed as [spin_k_index][band].
    /// For nspin=1: length = n_kpoints. For nspin=2: length = 2 * n_kpoints
    /// (spin-up k-points followed by spin-down k-points).
    pub eigenvalues: Vec<Vec<f64>>,
    pub fermi_energy: f64,
    pub n_iterations: usize,
    pub rho_g: Vec<Complex64>,
    /// Total magnetization M = ∫(ρ_up - ρ_down)dr in μB (Bohr magnetons).
    /// Zero for nspin=1.
    pub magnetization: f64,
    /// Number of spin channels (1 or 2).
    pub nspin: usize,
}

/// Run the self-consistent field loop.
///
/// Dispatches to `run_scf_spin` for nspin=2.
pub fn run_scf(
    crystal: &Crystal,
    basis: &BasisSet,
    kpoints: &[KPoint],
    pseudopotentials: &[&PseudopotentialData],
    params: &ScfParams,
    symmetry: Option<&crate::symmetry::SymmetryInfo>,
) -> Result<ScfResult> {
    if params.nspin == 2 {
        return run_scf_spin(crystal, basis, kpoints, pseudopotentials, params, symmetry);
    }

    let omega = crystal.lattice.volume().abs();
    let n_electrons: f64 = crystal
        .atoms
        .iter()
        .map(|a| crate::pseudopotential::find_for_atom(a.z, pseudopotentials).z_valence)
        .sum();

    info!("SCF: {n_electrons} electrons, {omega:.3} ų cell volume");

    let mut grid = FftGrid::new(basis, &crystal.lattice, params.ecutrho_ratio, params.fft_grid);
    let n_grid = grid.total_size();
    let [nx, ny, nz] = grid.dims;
    info!("FFT grid: {nx}×{ny}×{nz} = {n_grid} points (ecutrho_ratio={})", params.ecutrho_ratio);

    let g_to_fft = grid.basis_to_fft(basis);

    // Precompute local pseudopotential on the FULL FFT grid
    let mut v_local_fft = compute_v_local(crystal, &grid, pseudopotentials, omega);

    // Store V_local(G=0) separately and zero it in the FFT grid.
    // V_local(G=0) is an arbitrary constant (depends on PP construction) that shifts
    // all eigenvalues equally. QE excludes it from the Hamiltonian and adds it to the
    // total energy as v_of_0 * n_electrons. We follow the same convention.
    let v_local_g0 = v_local_fft[0].re;
    v_local_fft[0] = Complex64::new(0.0, 0.0);
    info!("V_local(G=0) = {v_local_g0:.6} eV (excluded from Hamiltonian)");

    // Try to initialize GPU if compiled with gpu feature
    #[cfg(feature = "gpu")]
    let mut gpu = crate::gpu::GpuAccelerator::try_new();
    #[cfg(feature = "gpu")]
    if gpu.is_some() {
        info!("GPU acceleration enabled for grid operations");
    }
    // Precompute |G|² for each FFT grid point (used by Hartree, both CPU and GPU)
    let dims = grid.dims;
    let recip = grid.recip.clone();
    let g_squared: Vec<f64> = (0..n_grid)
        .into_par_iter()
        .map(|idx| g_vector_at_dims(idx, dims, &recip).norm_squared())
        .collect();

    // Prepare persistent GPU buffers if GPU is available
    #[cfg(feature = "gpu")]
    if let Some(ref mut g) = gpu {
        g.prepare_buffers(n_grid, &g_squared);
    }

    // Initial density: superposition of atomic densities (SAD)
    let init_config = initial_density::InitialDensityConfig::non_magnetic(crystal.atoms.len());
    let mut rho_r = initial_density::generate_initial_density(
        crystal, &mut grid, pseudopotentials, n_electrons, &init_config,
    );
    let mut rho_g = vec![Complex64::new(0.0, 0.0); n_grid];
    density_r_to_g(&mut grid.fft, &rho_r, &mut rho_g);
    info!("Initial density: superposition of atomic densities (Gaussian model)");

    let mut mixer = mixing::AndersonMixer::new(
        params.mixing_beta,
        params.mixing_ndim,
        &params.mixing_mode,
        Some(&g_squared),
        n_electrons,
        omega,
    );
    let mut eigenvalues_all: Vec<Vec<f64>>;
    let mut fermi_energy;

    // Cache Ewald energy (constant across iterations)
    let e_ewald = crate::ewald::ewald_energy(crystal, pseudopotentials);

    // Precompute non-local projectors per k-point (depends only on k+G, not density)
    let vnl_cache: Vec<NonlocalPotential> = kpoints
        .par_iter()
        .map(|kp| NonlocalPotential::new(crystal, basis, &kp.k, pseudopotentials))
        .collect();
    let mut e_prev: Option<f64> = None;

    // Precompute NLCC core density on real-space grid (constant across iterations).
    // ρ_core(r) is added to ρ_valence(r) before XC evaluation.
    let rho_core_r = compute_core_density(crystal, &mut grid, pseudopotentials);
    if !rho_core_r.is_empty() {
        let core_min = rho_core_r.iter().copied().fold(f64::INFINITY, f64::min);
        let core_max = rho_core_r.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let core_sum: f64 = rho_core_r.iter().sum::<f64>() * omega / n_grid as f64;
        let _has_nan = rho_core_r.iter().any(|v| v.is_nan());
        info!("NLCC core density: min={core_min:.4e} max={core_max:.4e} integral={core_sum:.4}");
    }

    for iter in 0..params.max_iter {
        // Steps 1-3: Hartree, XC, V_eff assembly.
        // GPU path uses f32 for Hartree, XC, and V_eff; CPU path uses f64 + rayon.
        // XC FFT normalization is shared between both paths.

        // 1. Hartree potential
        #[cfg(feature = "gpu")]
        let v_h_fft = if let Some(ref gpu) = gpu {
            let fourpi_e2 = 4.0 * std::f64::consts::PI * hartree::E2;
            gpu.hartree_potential(&rho_g, &g_squared, fourpi_e2)
        } else {
            hartree_on_fft_grid(&rho_g, &g_squared)
        };
        #[cfg(not(feature = "gpu"))]
        let v_h_fft = hartree_on_fft_grid(&rho_g, &g_squared);

        // 2. XC potential: compute in real space, FFT to G-space
        // NLCC: add core density to valence density for XC evaluation
        let rho_for_xc = add_core_density(&rho_r, &rho_core_r);
        #[cfg(feature = "gpu")]
        let (_exc_r, vxc_r) = if let Some(ref gpu) = gpu {
            gpu.lda_xc(&rho_for_xc)
        } else {
            xc::lda_xc_grid(&rho_for_xc)
        };
        #[cfg(not(feature = "gpu"))]
        let (_exc_r, vxc_r) = xc::lda_xc_grid(&rho_for_xc);

        let vxc_g = real_to_g_space(&vxc_r, &mut grid.fft);

        // 3. V_eff = V_local + V_H + V_xc
        #[cfg(feature = "gpu")]
        let v_eff_fft = if let Some(ref gpu) = gpu {
            gpu.v_eff_assembly(&v_local_fft, &v_h_fft, &vxc_g)
        } else {
            assemble_v_eff(&v_local_fft, &v_h_fft, &vxc_g)
        };
        #[cfg(not(feature = "gpu"))]
        let v_eff_fft = assemble_v_eff(&v_local_fft, &v_h_fft, &vxc_g);

        // 4. Solve eigenvalue problem at each k-point (parallel over k-points)
        let kpoint_results: Vec<_> = kpoints
            .par_iter()
            .enumerate()
            .map(|(ik, kp)| {
                let mut h = build_hamiltonian_with_v_eff(basis, &kp.k, &v_eff_fft, grid.dims);
                vnl_cache[ik].add_to_hamiltonian(&mut h, crystal, basis, &kp.k);
                dense::diagonalize_lowest(&h, params.n_bands)
            })
            .collect();

        eigenvalues_all = kpoint_results.iter().map(|r| r.eigenvalues.clone()).collect();
        let all_kpoint_wavefns: Vec<_> = kpoint_results.into_iter().map(|r| r.eigenvectors).collect();

        // 5. Occupations (configurable smearing scheme)
        let spin_factor = 2.0; // nspin=1: each state holds 2 electrons
        let kpt_weights: Vec<f64> = kpoints.iter().map(|kp| kp.weight).collect();
        fermi_energy = smearing::find_fermi_energy(
            &eigenvalues_all,
            &kpt_weights,
            n_electrons,
            params.smearing_sigma,
            params.smearing_scheme,
            spin_factor,
        );
        let occupations: Vec<Vec<f64>> = eigenvalues_all
            .iter()
            .map(|evs| {
                evs.iter()
                    .map(|&e| smearing::occupation(params.smearing_scheme, e, fermi_energy, params.smearing_sigma, spin_factor))
                    .collect()
            })
            .collect();

        // 6. New density
        let mut rho_r_new = density::compute_density(
            basis, kpoints, &all_kpoint_wavefns, &occupations, &g_to_fft, &mut grid.fft,
            n_electrons, omega,
        );

        // 6b. Symmetrize density if symmetry info is available
        if let Some(symm) = symmetry {
            crate::symmetry::density::symmetrize_density(&mut rho_r_new, grid.dims, symm);
        }

        // 7. Convergence check (dual criterion: density AND energy)
        let delta = density_diff(&rho_r, &rho_r_new, omega, n_grid);

        // Compute energy every iteration for convergence monitoring
        let mut rho_g_new = vec![Complex64::new(0.0, 0.0); n_grid];
        for (i, &r) in rho_r_new.iter().enumerate() {
            rho_g_new[i] = Complex64::new(r, 0.0);
        }
        grid.fft.forward(&mut rho_g_new);
        let fft_norm = 1.0 / n_grid as f64;
        for v in &mut rho_g_new {
            *v *= fft_norm;
        }

        let rho_new_for_xc = add_core_density(&rho_r_new, &rho_core_r);
        let (exc_r, vxc_r_energy) = xc::lda_xc_grid(&rho_new_for_xc);

        let e_total = total_energy(
            band_energy(&eigenvalues_all, &occupations, &kpt_weights),
            hartree_energy(&rho_g_new, &g_squared, omega),
            xc_energy_corrected(&rho_new_for_xc, &rho_r_new, &exc_r, &vxc_r_energy, omega),
            e_ewald,
        ) + v_local_g0 * n_electrons;

        let de = e_prev.map(|ep| (e_total - ep).abs());
        e_prev = Some(e_total);

        let rho_converged = delta < params.conv_threshold;
        let energy_converged = de.is_some_and(|de| de < params.energy_threshold);

        info!(
            "SCF iter {:>3}: E={:.6} eV  dE={:>10}  Δρ={:.2e}",
            iter + 1, e_total,
            de.map_or("N/A".to_string(), |de| format!("{de:.2e}")),
            delta
        );

        if rho_converged && energy_converged {
            info!("SCF converged after {} iterations", iter + 1);
            rho_g = rho_g_new;
            let rho_g_basis: Vec<Complex64> = g_to_fft.iter().map(|&idx| rho_g[idx]).collect();

            // Entropy and free energy
            let ts = smearing::entropy_ts(
                &eigenvalues_all, &kpt_weights, fermi_energy,
                params.smearing_sigma, params.smearing_scheme, spin_factor,
            );
            let free_energy = e_total - ts;
            let energy_sigma0 = (e_total + free_energy) / 2.0;

            info!("Energy (E):      {e_total:.6} eV");
            info!("Free energy (F): {free_energy:.6} eV");
            info!("E sigma→0 (E₀):  {energy_sigma0:.6} eV");
            if ts.abs() > 1e-8 {
                info!(
                    "Entropy (-TS):   {:.6} eV ({:.3} meV/atom)",
                    -ts, -ts * 1000.0 / crystal.atoms.len() as f64
                );
            }

            return Ok(ScfResult {
                total_energy: e_total,
                free_energy,
                energy_sigma0,
                entropy_ts: ts,
                eigenvalues: eigenvalues_all,
                fermi_energy,
                n_iterations: iter + 1,
                rho_g: rho_g_basis,
                magnetization: 0.0,
                nspin: 1,
            });
        }

        // 8. Mix (Kerker preconditioning applied inside if enabled)
        rho_r = mixer.mix(&rho_r, &rho_r_new, &mut grid.fft);
        density_r_to_g(&mut grid.fft, &rho_r, &mut rho_g);
    }

    Err(PwdftError::ConvergenceFailure {
        iterations: params.max_iter,
        delta: density_diff(&rho_r, &rho_r, omega, n_grid),
    })
}

/// Compute local pseudopotential V_local(G) on the FULL FFT grid.
/// Spin-polarized SCF loop (nspin=2).
///
/// Two spin channels with independent densities, XC potentials, and Hamiltonians.
/// Hartree and V_local are spin-independent. V_xc is spin-dependent via LSDA.
fn run_scf_spin(
    crystal: &Crystal,
    basis: &BasisSet,
    kpoints: &[KPoint],
    pseudopotentials: &[&PseudopotentialData],
    params: &ScfParams,
    symmetry: Option<&crate::symmetry::SymmetryInfo>,
) -> Result<ScfResult> {
    let omega = crystal.lattice.volume().abs();
    let n_electrons: f64 = crystal
        .atoms
        .iter()
        .map(|a| crate::pseudopotential::find_for_atom(a.z, pseudopotentials).z_valence)
        .sum();

    info!("Spin-polarized SCF: {n_electrons} electrons, nspin=2");

    let mut grid = FftGrid::new(basis, &crystal.lattice, params.ecutrho_ratio, params.fft_grid);
    let n_grid = grid.total_size();
    let [nx, ny, nz] = grid.dims;
    info!("FFT grid: {nx}×{ny}×{nz} = {n_grid} points");

    let g_to_fft = grid.basis_to_fft(basis);
    let mut v_local_fft = compute_v_local(crystal, &grid, pseudopotentials, omega);
    let v_local_g0 = v_local_fft[0].re;
    v_local_fft[0] = Complex64::new(0.0, 0.0);
    info!("V_local(G=0) = {v_local_g0:.6} eV (excluded from Hamiltonian)");

    let dims = grid.dims;
    let recip = grid.recip.clone();
    let g_squared: Vec<f64> = (0..n_grid)
        .into_par_iter()
        .map(|idx| g_vector_at_dims(idx, dims, &recip).norm_squared())
        .collect();

    let rho_core_r = compute_core_density(crystal, &mut grid, pseudopotentials);

    // Determine initial spin split from starting_magnetization
    let per_atom_mag: Vec<f64> = crystal.atoms.iter().map(|a| {
        let sym = crate::atoms::from_z(a.z).map(|e| e.symbol().to_string()).unwrap_or_default();
        *params.starting_magnetization.get(&sym).unwrap_or(&0.0)
    }).collect();

    // Determine n_up, n_down
    let (n_up, n_down) = if let Some(tot_mag) = params.tot_magnetization {
        let n_up = (n_electrons + tot_mag) / 2.0;
        let n_down = (n_electrons - tot_mag) / 2.0;
        info!("Fixed magnetization: n_up={n_up:.2}, n_down={n_down:.2}");
        (n_up, n_down)
    } else {
        (n_electrons / 2.0, n_electrons / 2.0) // initial guess; self-consistent
    };

    // Initial spin density via SAD with magnetic moments
    let init_config = initial_density::InitialDensityConfig {
        magnetic_moments: per_atom_mag.clone(),
        gaussian_sigma: None,
    };
    // Generate total density then split
    let rho_total = initial_density::generate_initial_density(
        crystal, &mut grid, pseudopotentials, n_electrons, &init_config,
    );
    // Split: rho_up = (1+m)/2 * rho, rho_down = (1-m)/2 * rho
    // For a uniform initial split, m=0 → equal channels
    let avg_mag: f64 = per_atom_mag.iter().sum::<f64>() / per_atom_mag.len().max(1) as f64;
    let mut rho_up_r: Vec<f64> = rho_total.iter().map(|&r| r * (1.0 + avg_mag) / 2.0).collect();
    let mut rho_down_r: Vec<f64> = rho_total.iter().map(|&r| r * (1.0 - avg_mag) / 2.0).collect();

    let spin_factor = 1.0; // nspin=2: each state holds 1 electron
    let kpt_weights: Vec<f64> = kpoints.iter().map(|kp| kp.weight).collect();
    let _n_kpts = kpoints.len();

    // Precompute non-local projectors per k-point
    let vnl_cache: Vec<NonlocalPotential> = kpoints
        .par_iter()
        .map(|kp| NonlocalPotential::new(crystal, basis, &kp.k, pseudopotentials))
        .collect();

    // Two mixers (one per spin channel)
    let mut mixer_up = mixing::AndersonMixer::new(
        params.mixing_beta, params.mixing_ndim, &params.mixing_mode,
        Some(&g_squared), n_up, omega,
    );
    let mut mixer_down = mixing::AndersonMixer::new(
        params.mixing_beta, params.mixing_ndim, &params.mixing_mode,
        Some(&g_squared), n_down, omega,
    );

    let e_ewald = crate::ewald::ewald_energy(crystal, pseudopotentials);
    let mut e_prev: Option<f64> = None;

    for iter in 0..params.max_iter {
        // Total density for Hartree (spin-independent)
        let rho_total_r: Vec<f64> = rho_up_r.iter().zip(rho_down_r.iter())
            .map(|(&u, &d)| u + d).collect();
        let mut rho_total_g = vec![Complex64::new(0.0, 0.0); n_grid];
        density_r_to_g(&mut grid.fft, &rho_total_r, &mut rho_total_g);

        // 1. Hartree from total density
        let v_h_fft = hartree_on_fft_grid(&rho_total_g, &g_squared);

        // 2. Spin-dependent XC
        let rho_up_xc = add_core_density(&rho_up_r, &rho_core_r.iter().map(|&c| c / 2.0).collect::<Vec<_>>());
        let rho_down_xc = add_core_density(&rho_down_r, &rho_core_r.iter().map(|&c| c / 2.0).collect::<Vec<_>>());
        let (exc_r, vxc_up_r, vxc_down_r) = xc::lda_xc_spin_grid(&rho_up_xc, &rho_down_xc);

        let vxc_up_g = real_to_g_space(&vxc_up_r, &mut grid.fft);
        let vxc_down_g = real_to_g_space(&vxc_down_r, &mut grid.fft);

        // 3. Two V_eff: V_local + V_H + V_xc_σ
        let v_eff_up: Vec<Complex64> = v_local_fft.iter().zip(v_h_fft.iter()).zip(vxc_up_g.iter())
            .map(|((&vl, &vh), &vxc)| vl + vh + vxc).collect();
        let v_eff_down: Vec<Complex64> = v_local_fft.iter().zip(v_h_fft.iter()).zip(vxc_down_g.iter())
            .map(|((&vl, &vh), &vxc)| vl + vh + vxc).collect();

        // 4. Diagonalize both spins at each k-point
        let grid_dims = grid.dims;
        let kpoint_results_up: Vec<_> = kpoints.par_iter().enumerate().map(|(ik, kp)| {
            let mut h = build_hamiltonian_with_v_eff(basis, &kp.k, &v_eff_up, grid_dims);
            vnl_cache[ik].add_to_hamiltonian(&mut h, crystal, basis, &kp.k);
            dense::diagonalize_lowest(&h, params.n_bands)
        }).collect();

        let kpoint_results_down: Vec<_> = kpoints.par_iter().enumerate().map(|(ik, kp)| {
            let mut h = build_hamiltonian_with_v_eff(basis, &kp.k, &v_eff_down, grid_dims);
            vnl_cache[ik].add_to_hamiltonian(&mut h, crystal, basis, &kp.k);
            dense::diagonalize_lowest(&h, params.n_bands)
        }).collect();

        // Flatten eigenvalues: [up_k0, up_k1, ..., down_k0, down_k1, ...]
        let eig_up: Vec<Vec<f64>> = kpoint_results_up.iter().map(|r| r.eigenvalues.clone()).collect();
        let eig_down: Vec<Vec<f64>> = kpoint_results_down.iter().map(|r| r.eigenvalues.clone()).collect();
        let eigenvalues_all: Vec<Vec<f64>> = eig_up.iter().chain(eig_down.iter()).cloned().collect();
        let weights_all: Vec<f64> = kpt_weights.iter().chain(kpt_weights.iter()).copied().collect();

        let wfn_up: Vec<_> = kpoint_results_up.into_iter().map(|r| r.eigenvectors).collect();
        let wfn_down: Vec<_> = kpoint_results_down.into_iter().map(|r| r.eigenvectors).collect();

        // 5. Fermi energy and occupations
        let (fermi_energy, occ_up, occ_down) = if let Some(tot_mag) = params.tot_magnetization {
            // Fixed magnetization: separate Fermi energies per spin
            let n_up_target = (n_electrons + tot_mag) / 2.0;
            let n_down_target = (n_electrons - tot_mag) / 2.0;

            let ef_up = smearing::find_fermi_energy(
                &eig_up, &kpt_weights, n_up_target,
                params.smearing_sigma, params.smearing_scheme, spin_factor,
            );
            let ef_down = smearing::find_fermi_energy(
                &eig_down, &kpt_weights, n_down_target,
                params.smearing_sigma, params.smearing_scheme, spin_factor,
            );

            let occ_up: Vec<Vec<f64>> = eig_up.iter().map(|evs| {
                evs.iter().map(|&e| smearing::occupation(params.smearing_scheme, e, ef_up, params.smearing_sigma, spin_factor)).collect()
            }).collect();
            let occ_down: Vec<Vec<f64>> = eig_down.iter().map(|evs| {
                evs.iter().map(|&e| smearing::occupation(params.smearing_scheme, e, ef_down, params.smearing_sigma, spin_factor)).collect()
            }).collect();

            // Report average Fermi energy
            ((ef_up + ef_down) / 2.0, occ_up, occ_down)
        } else {
            // Free magnetization: single Fermi energy for both spins
            let fermi_energy = smearing::find_fermi_energy(
                &eigenvalues_all, &weights_all, n_electrons,
                params.smearing_sigma, params.smearing_scheme, spin_factor,
            );

            let occ_up: Vec<Vec<f64>> = eig_up.iter().map(|evs| {
                evs.iter().map(|&e| smearing::occupation(params.smearing_scheme, e, fermi_energy, params.smearing_sigma, spin_factor)).collect()
            }).collect();
            let occ_down: Vec<Vec<f64>> = eig_down.iter().map(|evs| {
                evs.iter().map(|&e| smearing::occupation(params.smearing_scheme, e, fermi_energy, params.smearing_sigma, spin_factor)).collect()
            }).collect();

            (fermi_energy, occ_up, occ_down)
        };

        // 6. Reconstruct spin densities
        let n_el_up: f64 = occ_up.iter().zip(kpt_weights.iter())
            .flat_map(|(occs, &w)| occs.iter().map(move |&f| f * w))
            .sum();
        let n_el_down: f64 = occ_down.iter().zip(kpt_weights.iter())
            .flat_map(|(occs, &w)| occs.iter().map(move |&f| f * w))
            .sum();

        let rho_up_new = density::compute_density(
            basis, kpoints, &wfn_up, &occ_up, &g_to_fft, &mut grid.fft, n_el_up, omega,
        );
        let rho_down_new = density::compute_density(
            basis, kpoints, &wfn_down, &occ_down, &g_to_fft, &mut grid.fft, n_el_down, omega,
        );

        // Symmetrize each channel
        let mut rho_up_sym = rho_up_new;
        let mut rho_down_sym = rho_down_new;
        if let Some(symm) = symmetry {
            crate::symmetry::density::symmetrize_density(&mut rho_up_sym, grid.dims, symm);
            crate::symmetry::density::symmetrize_density(&mut rho_down_sym, grid.dims, symm);
        }

        // 7. Convergence check
        let rho_total_new: Vec<f64> = rho_up_sym.iter().zip(rho_down_sym.iter())
            .map(|(&u, &d)| u + d).collect();
        let delta = density_diff(&rho_total_r, &rho_total_new, omega, n_grid);

        // Energy from new density
        let mut rho_total_new_g = vec![Complex64::new(0.0, 0.0); n_grid];
        density_r_to_g(&mut grid.fft, &rho_total_new, &mut rho_total_new_g);

        let rho_xc_total: Vec<f64> = add_core_density(&rho_total_new, &rho_core_r);
        let occ_all: Vec<Vec<f64>> = occ_up.iter().chain(occ_down.iter()).cloned().collect();

        // Spin XC double-counting: E_vxc = ∫(V_xc_up ρ_up + V_xc_down ρ_down) dr
        let dvol = omega / n_grid as f64;
        let e_vxc_spin: f64 = rho_up_sym.iter().zip(vxc_up_r.iter())
            .zip(rho_down_sym.iter().zip(vxc_down_r.iter()))
            .map(|((&ru, &vu), (&rd, &vd))| (ru * vu + rd * vd) * dvol)
            .sum();
        let e_xc = xc::lda_xc_energy(&rho_xc_total, &exc_r, omega);
        let e_xc_corrected = e_xc - e_vxc_spin;

        let e_total = total_energy(
            band_energy(&eigenvalues_all, &occ_all, &weights_all),
            hartree_energy(&rho_total_new_g, &g_squared, omega),
            e_xc_corrected,
            e_ewald,
        ) + v_local_g0 * n_electrons;

        let de = e_prev.map(|ep| (e_total - ep).abs());
        e_prev = Some(e_total);

        let rho_converged = delta < params.conv_threshold;
        let energy_converged = de.is_some_and(|de| de < params.energy_threshold);

        let mag = (n_el_up - n_el_down).abs();
        info!(
            "SCF iter {:>3}: E={:.6} eV  dE={:>10}  Δρ={:.2e}  M={:.3} μB",
            iter + 1, e_total,
            de.map_or("N/A".to_string(), |de| format!("{de:.2e}")),
            delta, mag
        );

        if rho_converged && energy_converged {
            info!("SCF converged after {} iterations", iter + 1);
            info!("Magnetization: {mag:.4} μB ({n_el_up:.4} up, {n_el_down:.4} down)");

            let rho_g_basis: Vec<Complex64> = g_to_fft.iter().map(|&idx| rho_total_new_g[idx]).collect();

            let ts = smearing::entropy_ts(
                &eigenvalues_all, &weights_all, fermi_energy,
                params.smearing_sigma, params.smearing_scheme, spin_factor,
            );
            let free_energy = e_total - ts;
            let energy_sigma0 = (e_total + free_energy) / 2.0;

            return Ok(ScfResult {
                total_energy: e_total,
                free_energy,
                energy_sigma0,
                entropy_ts: ts,
                eigenvalues: eigenvalues_all,
                fermi_energy,
                n_iterations: iter + 1,
                rho_g: rho_g_basis,
                magnetization: mag,
                nspin: 2,
            });
        }

        // 8. Mix each spin channel independently
        rho_up_r = mixer_up.mix(&rho_up_r, &rho_up_sym, &mut grid.fft);
        rho_down_r = mixer_down.mix(&rho_down_r, &rho_down_sym, &mut grid.fft);
    }

    Err(PwdftError::ConvergenceFailure {
        iterations: params.max_iter,
        delta: 0.0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    
    use approx::relative_eq;
    use crate::fft::FFT3D;

    #[test]
    fn test_real_to_g_space_dc_component() {
        // A constant real-space function f(r) = C should give
        // F(G=0) = C and F(G≠0) = 0.
        let mut fft = FFT3D::new(8, 8, 8);
        let n = fft.total_size();
        let c = 3.5;
        let data_r = vec![c; n];
        let data_g = real_to_g_space(&data_r, &mut fft);

        // G=0 component (index 0) should be C
        assert!(
            relative_eq!(data_g[0].re, c, epsilon = 1e-10),
            "DC component: expected {c}, got {}", data_g[0].re
        );
        assert!(data_g[0].im.abs() < 1e-10);

        // All other G-components should be ~0
        for (i, &v) in data_g.iter().enumerate().skip(1) {
            assert!(
                v.norm() < 1e-10,
                "G≠0 component at {i}: expected ~0, got {v}"
            );
        }
    }

    #[test]
    fn test_real_to_g_space_roundtrip() {
        let mut fft = FFT3D::new(8, 8, 8);
        let n = fft.total_size();
        let data_r: Vec<f64> = (0..n).map(|i| (i as f64 * 0.1).sin()).collect();
        let data_g = real_to_g_space(&data_r, &mut fft);

        // Inverse FFT should recover original (unnormalized → need N factor)
        let mut data_back = data_g;
        fft.inverse(&mut data_back);
        // real_to_g_space divides by N, inverse multiplies by N → should recover original
        for (i, (&orig, &back)) in data_r.iter().zip(data_back.iter()).enumerate() {
            assert!(
                relative_eq!(orig, back.re, epsilon = 1e-10),
                "Roundtrip failed at {i}: original={orig}, recovered={}", back.re
            );
            assert!(back.im.abs() < 1e-10, "Imaginary part at {i}: {}", back.im);
        }
    }

    #[test]
    fn test_assemble_v_eff_adds_correctly() {
        let n = 100;
        let v1: Vec<Complex64> = (0..n).map(|i| Complex64::new(i as f64, 0.0)).collect();
        let v2: Vec<Complex64> = (0..n).map(|i| Complex64::new(0.0, i as f64 * 0.1)).collect();
        let v3: Vec<Complex64> = (0..n).map(|i| Complex64::new(-(i as f64) * 0.5, 0.0)).collect();

        let result = assemble_v_eff(&v1, &v2, &v3);

        for i in 0..n {
            let expected = v1[i] + v2[i] + v3[i];
            assert!(
                (result[i] - expected).norm() < 1e-14,
                "V_eff mismatch at {i}: expected {expected}, got {}", result[i]
            );
        }
    }

    #[test]
    fn test_hartree_on_fft_grid_g0_zero() {
        // V_H(G=0) should be zero (no divergence)
        let rho_g = vec![Complex64::new(1.0, 0.0); 10];
        let g_squared = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let v_h = hartree_on_fft_grid(&rho_g, &g_squared);
        assert!(v_h[0].norm() < 1e-15, "V_H(G=0) should be zero, got {}", v_h[0]);
        // V_H(G≠0) should be finite and positive real for positive ρ
        for &v in &v_h[1..] {
            assert!(v.re > 0.0, "V_H should be positive for positive ρ: {v}");
        }
    }

    #[test]
    fn test_density_diff_identical() {
        let rho = vec![1.0; 100];
        let diff = density_diff(&rho, &rho, 40.0, 100);
        assert!(diff < 1e-15, "Identical densities should give zero diff: {diff}");
    }

    #[test]
    fn test_density_diff_known() {
        let omega = 40.0;
        let n = 100;
        let rho_a = vec![1.0; n];
        let rho_b = vec![2.0; n];
        // diff = sqrt(Σ(1.0)² × dvol / omega) = sqrt(n × dvol / omega) = sqrt(dvol × n / omega)
        // dvol = omega / n = 0.4
        // diff = sqrt(0.4 * 100 / 40) = sqrt(1.0) = 1.0
        let diff = density_diff(&rho_a, &rho_b, omega, n);
        assert!(
            relative_eq!(diff, 1.0, epsilon = 1e-10),
            "Expected diff=1.0, got {diff}"
        );
    }
}
