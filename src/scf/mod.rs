pub mod density;
pub mod initial_density;
pub mod mixing;
pub mod smearing;

use log::info;
use nalgebra::Vector3;
use num_complex::Complex64;
use rayon::prelude::*;

use crate::{
    basis::BasisSet,
    consts::HBAR2_OVER_2M,
    crystal::Crystal,
    eigensolver::dense,
    error::{PwdftError, Result},
    fft::{fft_grid_size, FFT3D},
    kpoints::KPoint,
    potential::{hartree, nonlocal::NonlocalPotential, xc},
    pseudopotential::PseudopotentialData,
};

/// Parameters for an SCF calculation.
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
    pub eigenvalues: Vec<Vec<f64>>,
    pub fermi_energy: f64,
    pub n_iterations: usize,
    pub rho_g: Vec<Complex64>,
}

/// Describes the FFT grid and provides index mapping.
struct FftGrid {
    dims: [usize; 3],
    fft: FFT3D,
    /// Reciprocal lattice vectors for computing G from Miller indices.
    recip: crate::crystal::Lattice,
}

impl FftGrid {
    /// Create an FFT grid sized for the charge density.
    ///
    /// `ecutrho_ratio`: multiplier for the wavefunction cutoff (default 4, matching
    /// QE's ecutrho = 4*ecutwfc for NC PPs). The grid must accommodate G-vectors
    /// up to sqrt(ecutrho_ratio) × G_max in each direction. Higher ratio = more
    /// accurate V_xc (nonlinear function of ρ) but larger grid.
    fn new(
        basis: &BasisSet,
        lattice: &crate::crystal::Lattice,
        ecutrho_ratio: u32,
        explicit_dims: Option<[usize; 3]>,
    ) -> Self {
        let dims = if let Some(d) = explicit_dims {
            d
        } else {
            let miller = basis.miller_indices();
            let n_max: Vec<i32> = (0..3)
                .map(|dim| miller.iter().map(|m| m[dim].abs()).max().unwrap_or(0))
                .collect();
            // Grid must accommodate G-vectors up to sqrt(ratio) × G_max.
            // For ratio=4 (QE default), this is 2× G_max in each direction.
            let scale = (ecutrho_ratio as f64).sqrt().ceil() as i32;
            [
                fft_grid_size(scale * n_max[0]),
                fft_grid_size(scale * n_max[1]),
                fft_grid_size(scale * n_max[2]),
            ]
        };
        let fft = FFT3D::new(dims[0], dims[1], dims[2]);
        let recip = lattice.reciprocal();
        Self { dims, fft, recip }
    }

    fn total_size(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    fn miller_to_idx(&self, n1: i32, n2: i32, n3: i32) -> usize {
        miller_to_idx(self.dims, n1, n2, n3)
    }

    /// Compute the Cartesian G-vector for FFT grid index.
    fn g_vector_at(&self, idx: usize) -> Vector3<f64> {
        g_vector_at_dims(idx, self.dims, &self.recip)
    }

    /// Mapping from basis G-vectors to FFT grid indices.
    fn basis_to_fft(&self, basis: &BasisSet) -> Vec<usize> {
        basis
            .miller_indices()
            .iter()
            .map(|&[n1, n2, n3]| self.miller_to_idx(n1, n2, n3))
            .collect()
    }
}

/// Run the self-consistent field loop.
///
/// If `symmetry` is provided, the charge density is symmetrized after each
/// SCF step to enforce crystal symmetry. K-point reduction should be done
/// by the caller before passing `kpoints`.
/// Compute G-vector from FFT grid index (standalone, safe to call from parallel contexts).
fn g_vector_at_dims(idx: usize, dims: [usize; 3], recip: &crate::crystal::Lattice) -> Vector3<f64> {
    let [nx, ny, nz] = dims;
    let i1 = idx / (ny * nz);
    let i2 = (idx / nz) % ny;
    let i3 = idx % nz;
    let n1 = if i1 > nx / 2 { i1 as i32 - nx as i32 } else { i1 as i32 };
    let n2 = if i2 > ny / 2 { i2 as i32 - ny as i32 } else { i2 as i32 };
    let n3 = if i3 > nz / 2 { i3 as i32 - nz as i32 } else { i3 as i32 };
    n1 as f64 * recip.a + n2 as f64 * recip.b + n3 as f64 * recip.c
}

pub fn run_scf(
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

    info!("SCF: {n_electrons} electrons, {omega:.3} ų cell volume");

    let mut grid = FftGrid::new(basis, &crystal.lattice, params.ecutrho_ratio, params.fft_grid);
    let n_grid = grid.total_size();
    let [nx, ny, nz] = grid.dims;
    info!("FFT grid: {nx}×{ny}×{nz} = {n_grid} points (ecutrho_ratio={})", params.ecutrho_ratio);

    let g_to_fft = grid.basis_to_fft(basis);

    // Precompute local pseudopotential on the FULL FFT grid
    let mut v_local_fft = compute_v_local_on_fft_grid(crystal, &grid, pseudopotentials, omega);

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
    let mut e_prev: Option<f64> = None;

    // Precompute NLCC core density on real-space grid (constant across iterations).
    // ρ_core(r) is added to ρ_valence(r) before XC evaluation.
    let rho_core_r = compute_core_density(crystal, &mut grid, pseudopotentials);

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
            .map(|kp| {
                let mut h = build_hamiltonian_with_v_eff(basis, &kp.k, &v_eff_fft, grid.dims);

                // Add non-local pseudopotential
                let vnl = NonlocalPotential::new(crystal, basis, &kp.k, pseudopotentials);
                vnl.add_to_hamiltonian(&mut h, crystal, basis, &kp.k);

                dense::diagonalize_lowest(&h, params.n_bands)
            })
            .collect();

        eigenvalues_all = kpoint_results.iter().map(|r| r.eigenvalues.clone()).collect();
        let all_kpoint_wavefns: Vec<_> = kpoint_results.into_iter().map(|r| r.eigenvectors).collect();

        // 5. Occupations (configurable smearing scheme)
        fermi_energy = smearing::find_fermi_energy(
            &eigenvalues_all,
            kpoints,
            n_electrons,
            params.smearing_sigma,
            params.smearing_scheme,
        );
        let occupations: Vec<Vec<f64>> = eigenvalues_all
            .iter()
            .map(|evs| {
                evs.iter()
                    .map(|&e| smearing::occupation(params.smearing_scheme, e, fermi_energy, params.smearing_sigma))
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
        // NLCC: E_xc uses total density (val+core), E_vxc uses valence only
        let e_total = compute_total_energy_from_components(
            &eigenvalues_all, &occupations, kpoints,
            &rho_new_for_xc, // total density for E_xc = ∫ ε_xc × (ρ_val+ρ_core) dr
            &rho_r_new,       // valence density for E_vxc = ∫ V_xc × ρ_val dr
            &rho_g_new, &g_squared, &exc_r, &vxc_r_energy, omega, e_ewald,
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
                &eigenvalues_all, kpoints, fermi_energy,
                params.smearing_sigma, params.smearing_scheme,
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
fn compute_v_local_on_fft_grid(
    crystal: &Crystal,
    grid: &FftGrid,
    pseudopotentials: &[&PseudopotentialData],
    omega: f64,
) -> Vec<Complex64> {
    let n_grid = grid.total_size();

    // Precompute per-atom data to avoid repeated lookups in the hot loop
    let atom_data: Vec<(Vector3<f64>, &PseudopotentialData)> = crystal
        .atoms
        .iter()
        .map(|atom| {
            let pp = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials);
            (atom.cart_position(&crystal.lattice), pp)
        })
        .collect();

    let dims = grid.dims;
    let recip = grid.recip.clone();
    (0..n_grid)
        .into_par_iter()
        .map(|idx| {
            let g = g_vector_at_dims(idx, dims, &recip);
            let g_norm = g.norm();
            let mut v = Complex64::new(0.0, 0.0);

            for &(ref tau, pp) in &atom_data {
                let phase = -g.dot(tau);
                let sf = Complex64::new(phase.cos(), phase.sin());
                let v_form = pp.v_local_of_g(g_norm, omega);
                v += sf * v_form;
            }
            v
        })
        .collect()
}

/// Compute Hartree potential on the FULL FFT grid.
fn hartree_on_fft_grid(rho_g: &[Complex64], g_squared: &[f64]) -> Vec<Complex64> {
    let fourpi_e2 = 4.0 * std::f64::consts::PI * hartree::E2;

    rho_g
        .par_iter()
        .zip(g_squared.par_iter())
        .map(|(&rho, &g2)| {
            if g2 > 1e-20 {
                rho * fourpi_e2 / g2
            } else {
                Complex64::new(0.0, 0.0)
            }
        })
        .collect()
}

/// Convert real-space array to G-space with FFT normalization.
fn real_to_g_space(data_r: &[f64], fft: &mut FFT3D) -> Vec<Complex64> {
    let n = data_r.len();
    let mut data_g: Vec<Complex64> = data_r.iter().map(|&v| Complex64::new(v, 0.0)).collect();
    fft.forward(&mut data_g);
    let norm = 1.0 / n as f64;
    for v in &mut data_g {
        *v *= norm;
    }
    data_g
}

/// Assemble V_eff = V_local + V_H + V_xc (CPU path).
fn assemble_v_eff(
    v_local: &[Complex64],
    v_h: &[Complex64],
    v_xc: &[Complex64],
) -> Vec<Complex64> {
    v_local
        .par_iter()
        .zip(v_h.par_iter())
        .zip(v_xc.par_iter())
        .map(|((&vl, &vh), &vxc)| vl + vh + vxc)
        .collect()
}

/// Map Miller indices to a flat FFT grid index (standalone, for use in parallel contexts).
fn miller_to_idx(dims: [usize; 3], n1: i32, n2: i32, n3: i32) -> usize {
    let i1 = ((n1 % dims[0] as i32) + dims[0] as i32) as usize % dims[0];
    let i2 = ((n2 % dims[1] as i32) + dims[1] as i32) as usize % dims[1];
    let i3 = ((n3 % dims[2] as i32) + dims[2] as i32) as usize % dims[2];
    i1 * dims[1] * dims[2] + i2 * dims[2] + i3
}

/// Build Hamiltonian: kinetic + V_eff(G-G') looked up from FFT grid.
fn build_hamiltonian_with_v_eff(
    basis: &BasisSet,
    k: &Vector3<f64>,
    v_eff_fft: &[Complex64],
    grid_dims: [usize; 3],
) -> faer::Mat<Complex64> {
    let n = basis.len();
    let mut h = faer::Mat::<Complex64>::zeros(n, n);

    // Kinetic (diagonal)
    for (i, g) in basis.g_vectors().iter().enumerate() {
        let ke = HBAR2_OVER_2M * (k + g).norm_squared();
        h[(i, i)] = Complex64::new(ke, 0.0);
    }

    // Potential: V_{G,G'} = V_eff(G-G') from FFT grid
    let miller_idx = basis.miller_indices();
    for i in 0..n {
        for j in 0..n {
            let dn1 = miller_idx[i][0] - miller_idx[j][0];
            let dn2 = miller_idx[i][1] - miller_idx[j][1];
            let dn3 = miller_idx[i][2] - miller_idx[j][2];
            let fft_idx = miller_to_idx(grid_dims, dn1, dn2, dn3);
            h[(i, j)] += v_eff_fft[fft_idx];
        }
    }

    h
}

/// FFT density from real space to G-space (normalized).
fn density_r_to_g(fft: &mut FFT3D, rho_r: &[f64], rho_g: &mut [Complex64]) {
    for (i, &r) in rho_r.iter().enumerate() {
        rho_g[i] = Complex64::new(r, 0.0);
    }
    fft.forward(rho_g);
    let norm = 1.0 / fft.total_size() as f64;
    for v in rho_g.iter_mut() {
        *v *= norm;
    }
}

fn density_diff(rho_old: &[f64], rho_new: &[f64], omega: f64, n_grid: usize) -> f64 {
    let dvol = omega / n_grid as f64;
    let sum_sq: f64 = rho_old
        .iter()
        .zip(rho_new.iter())
        .map(|(&a, &b)| (a - b).powi(2) * dvol)
        .sum();
    (sum_sq / omega).sqrt()
}

/// Add NLCC core density to valence density for XC evaluation.
/// Returns rho_val if no core density is present (empty vec).
fn add_core_density(rho_val: &[f64], rho_core: &[f64]) -> Vec<f64> {
    if rho_core.is_empty() {
        rho_val.to_vec()
    } else {
        rho_val
            .iter()
            .zip(rho_core.iter())
            .map(|(&v, &c)| v + c)
            .collect()
    }
}

/// Compute NLCC core density on the real-space FFT grid.
///
/// For each atom with NLCC, computes ρ_core(G) via spherical Bessel
/// transform of PP_NLCC, accumulates with structure factors, then
/// inverse FFTs to real space. Returns empty vec if no PP has NLCC.
fn compute_core_density(
    crystal: &Crystal,
    grid: &mut FftGrid,
    pseudopotentials: &[&PseudopotentialData],
) -> Vec<f64> {
    // Check if any PP has NLCC
    let any_nlcc = pseudopotentials.iter().any(|pp| pp.has_nlcc);
    if !any_nlcc {
        return vec![];
    }

    let n_grid = grid.total_size();
    let omega = crystal.lattice.volume().abs();
    let mut rho_core_g = vec![Complex64::new(0.0, 0.0); n_grid];

    for atom in &crystal.atoms {
        let pp = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials);
        if !pp.has_nlcc || pp.core_charge.is_empty() {
            continue;
        }

        let tau = atom.cart_position(&crystal.lattice);

        // Bessel transform of core charge at each G-vector
        for (idx, rho_g_val) in rho_core_g.iter_mut().enumerate() {
            let g = grid.g_vector_at(idx);
            let g_norm = g.norm();

            // ∫ [4πr²ρ_core(r)] j₀(|G|r) dr
            let mut integral = 0.0;
            for ((&rho_c, &r), &dr) in pp
                .core_charge
                .iter()
                .zip(pp.r_grid.iter())
                .zip(pp.rab.iter())
            {
                let gr = g_norm * r;
                let j0 = if gr < 1e-10 {
                    1.0 - gr * gr / 6.0
                } else {
                    gr.sin() / gr
                };
                integral += rho_c * j0 * dr;
            }

            // Structure factor and normalization
            let phase = -g.dot(&tau);
            let sf = Complex64::new(phase.cos(), phase.sin());
            *rho_g_val += sf * (integral / omega);
        }
    }

    // Inverse FFT to real space
    grid.fft.inverse(&mut rho_core_g);

    // Extract real part (imaginary should be negligible)
    rho_core_g.iter().map(|c| c.re).collect()
}

/// Compute total energy reusing pre-computed XC and cached Ewald.
/// Called every SCF iteration for energy convergence monitoring.
///
/// `rho_xc_r`: density for XC energy (ρ_val + ρ_core if NLCC, else ρ_val).
/// `rho_val_r`: valence density only (for E_vxc double-counting correction).
#[allow(clippy::too_many_arguments)]
fn compute_total_energy_from_components(
    eigenvalues: &[Vec<f64>],
    occupations: &[Vec<f64>],
    kpoints: &[KPoint],
    rho_xc_r: &[f64],
    rho_val_r: &[f64],
    rho_g: &[Complex64],
    g_squared: &[f64],
    exc_r: &[f64],
    vxc_r: &[f64],
    omega: f64,
    e_ewald: f64,
) -> f64 {
    let n_grid = rho_g.len();

    let e_band: f64 = eigenvalues
        .iter()
        .zip(occupations.iter())
        .zip(kpoints.iter())
        .map(|((evs, occs), kp)| {
            evs.iter().zip(occs.iter()).map(|(&e, &f)| f * kp.weight * e).sum::<f64>()
        })
        .sum();

    let fourpi_e2 = 4.0 * std::f64::consts::PI * hartree::E2;
    let e_hartree: f64 = rho_g
        .iter()
        .zip(g_squared.iter())
        .map(|(rho, &g2)| {
            if g2 > 1e-20 { rho.norm_sqr() * fourpi_e2 / g2 } else { 0.0 }
        })
        .sum::<f64>()
        * 0.5
        * omega;

    // E_xc = ∫ ε_xc(ρ_total) × ρ_total dr  (uses val+core for NLCC)
    let e_xc = xc::lda_xc_energy(rho_xc_r, exc_r, omega);
    let dvol = omega / n_grid as f64;
    // E_vxc = ∫ V_xc(ρ_total) × ρ_val dr  (only valence for double-counting)
    let e_vxc: f64 = rho_val_r
        .iter()
        .zip(vxc_r.iter())
        .map(|(&rho, &vxc)| rho * vxc * dvol)
        .sum();

    e_band - e_hartree + e_xc - e_vxc + e_ewald
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

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
