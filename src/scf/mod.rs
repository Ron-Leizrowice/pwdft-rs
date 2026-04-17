pub(crate) mod context;
pub mod density;
pub(crate) mod energy;
pub(crate) mod grid;
pub mod initial_density;
pub mod mixing;
pub(crate) mod potentials;
pub mod smearing;

use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use num_complex::Complex64;
use rayon::prelude::*;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    eigensolver::dense,
    error::{PwdftError, Result},
    kpoints::KPoint,
    potential::xc,
    pseudopotential::PseudopotentialData,
};

use self::energy::{
    add_core_density, assemble_v_eff, band_energy, density_diff,
    density_r_to_g, harris_foulkes_energy, hartree_energy, hartree_on_fft_grid,
    kinetic_expectation, local_pp_energy_grid, nonlocal_expectation,
    real_to_g_space, total_energy, xc_energy_bare, xc_energy_corrected,
};
use self::grid::FftGrid;
use self::potentials::build_hamiltonian_with_v_eff;

/// Parameters controlling the self-consistent field iteration.
///
/// The SCF loop solves the Kohn-Sham equations iteratively:
/// 1. Construct V_eff = V_local + V_Hartree[ρ] + V_xc[ρ]
/// 2. Diagonalize H = T + V_eff + V_NL at each k-point
/// 3. Compute occupations from eigenvalues (Fermi-Dirac or other smearing)
/// 4. Reconstruct density ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²
/// 5. Mix input and output densities (Anderson/Pulay) and repeat
///
/// Convergence requires both density (Δρ < conv_threshold) and
/// energy (ΔE < energy_threshold) criteria to be met.
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

impl ScfParams {
    /// Validate parameters before starting an SCF calculation.
    pub fn validate(&self) -> Result<()> {
        if self.n_bands == 0 {
            return Err(PwdftError::InvalidInput("n_bands must be > 0".into()));
        }
        if self.conv_threshold <= 0.0 {
            return Err(PwdftError::InvalidInput("conv_threshold must be positive".into()));
        }
        if self.mixing_beta <= 0.0 || self.mixing_beta > 1.0 {
            return Err(PwdftError::InvalidInput(
                format!("mixing_beta must be in (0, 1], got {}", self.mixing_beta),
            ));
        }
        if self.smearing_sigma < 0.0 {
            return Err(PwdftError::InvalidInput("smearing_sigma must be non-negative".into()));
        }
        if self.ecutrho_ratio < 1 {
            return Err(PwdftError::InvalidInput(
                format!("ecutrho_ratio must be >= 1, got {}", self.ecutrho_ratio),
            ));
        }
        if self.nspin != 1 && self.nspin != 2 {
            return Err(PwdftError::InvalidInput(
                format!("nspin must be 1 or 2, got {}", self.nspin),
            ));
        }
        Ok(())
    }
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

/// Per-component energy decomposition of a converged SCF total energy.
///
/// All values in eV. Identity (at convergence):
/// ```text
/// E_total = e_kinetic
///         + e_local
///         + e_local_g0_shift   (= V_local(G=0) · N_el)
///         + e_nonlocal
///         + e_hartree
///         + e_xc
///         + e_ewald
/// ```
/// and by the Kohn-Sham double-counting identity:
/// ```text
/// e_band = e_kinetic + e_local + e_nonlocal + 2·e_hartree + e_vxc
/// ```
/// where `e_vxc = ∫ρ(r)·V_xc(r)dr`. The `V_local(G=0)·N_el` background shift
/// is the compensating term for zeroing the G=0 component of the local
/// pseudopotential in the Hamiltonian; see `scf::context::ScfContext::new`.
///
/// Mirrors QE's `pw.x` standard-output decomposition:
/// ```text
///   one-electron contribution = e_kinetic + e_local + e_nonlocal + e_local_g0_shift
///   hartree    contribution = e_hartree
///   xc         contribution = e_xc
///   ewald      contribution = e_ewald
/// ```
/// Intended for validation (see proposal VGC5) rather than routine SCF use.
/// Computed on the final iteration by one extra pass over wavefunctions,
/// V_local on the FFT grid, and the V_NL operator.
#[derive(Debug, Clone)]
pub struct EnergyComponents {
    /// Band energy: Σ_{n,k} f_{n,k} w_k ε_{n,k}.
    pub e_band: f64,
    /// Kinetic: Σ_{n,k} f·w·⟨ψ|T|ψ⟩ = Σ_{n,k} f·w·Σ_G |c_G|² · ℏ²/(2m)·|k+G|².
    pub e_kinetic: f64,
    /// Local PP (G ≠ 0): ∫ρ(r)·V_local(r)dr on the FFT grid (G=0 excluded).
    pub e_local: f64,
    /// Local PP G=0 compensating shift: V_local(G=0)·N_el.
    /// Constant background subtracted from `v_local_fft` at setup to keep the
    /// Hamiltonian diagonal finite.
    pub e_local_g0_shift: f64,
    /// Non-local (KB separable): Σ_{n,k} f·w·⟨ψ|V_NL|ψ⟩.
    pub e_nonlocal: f64,
    /// Hartree: (Ω/2) Σ_G |ρ(G)|² · 4πe²/|G|² (from OUTPUT density).
    pub e_hartree: f64,
    /// XC energy: ∫ρ(r)·ε_xc(r)dr (from OUTPUT density; same sign as QE's
    /// "xc contribution"). NLCC: ρ here is ρ_val + ρ_core.
    pub e_xc: f64,
    /// Ewald ion-ion energy (spin- and density-independent).
    pub e_ewald: f64,
}

/// Output of a converged SCF calculation.
///
/// All energies are in eV. The three energy quantities are:
/// - `total_energy`: E = E_band - E_H + E_xc - E_vxc + E_ewald + V_local(G=0)·N_el
/// - `free_energy`: F = E - TS (Mermin functional, variational at finite σ)
/// - `energy_sigma0`: E₀ = (E + F)/2 (best estimate of T=0 energy)
pub struct ScfResult {
    /// Kohn-Sham total energy (no entropy).
    pub total_energy: f64,
    /// Harris-Foulkes energy (double-counting from input density).
    ///
    /// Uses the input density for Hartree/XC corrections but output eigenvalues.
    /// Stationary at self-consistency: |E_HF - E_KS| -> 0 quadratically.
    /// Serves as a convergence quality indicator.
    pub harris_foulkes_energy: f64,
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
    /// Per-term energy breakdown (VGC5 diagnostic).
    pub components: EnergyComponents,
}

/// Compute occupation numbers for all k-points from eigenvalues and Fermi energy.
fn compute_occupations(
    eigenvalues: &[Vec<f64>],
    scheme: smearing::SmearingScheme,
    fermi_energy: f64,
    sigma: f64,
    spin_factor: f64,
) -> Vec<Vec<f64>> {
    eigenvalues
        .iter()
        .map(|evs| {
            evs.iter()
                .map(|&e| smearing::occupation(scheme, e, fermi_energy, sigma, spin_factor))
                .collect()
        })
        .collect()
}

fn scf_progress_bar(max_iter: usize) -> ProgressBar {
    let pb = ProgressBar::new(max_iter as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "SCF [{bar:30}] {pos}/{len}  {msg}  [{elapsed_precise} elapsed]",
        )
        // SAFETY: This is a static, valid template string -- with_template cannot fail.
        .expect("BUG: invalid progress bar template")
        .progress_chars("##-"),
    );
    pb
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
    symmetry: &crate::symmetry::SymmetryInfo,
) -> Result<ScfResult> {
    params.validate()?;
    if crystal.atoms.is_empty() {
        return Err(PwdftError::InvalidInput("at least one atom is required".into()));
    }
    if kpoints.is_empty() {
        return Err(PwdftError::InvalidInput("at least one k-point is required".into()));
    }
    let omega = crystal.lattice.volume();
    if omega < 1e-10 {
        return Err(PwdftError::InvalidInput("lattice has zero or near-zero volume".into()));
    }

    if params.nspin == 2 {
        return run_scf_spin(crystal, basis, kpoints, pseudopotentials, params, symmetry);
    }

    let mut ctx = context::ScfContext::new(
        crystal, basis, kpoints, pseudopotentials, params, symmetry,
    )?;

    // Try to initialize GPU if compiled with gpu feature
    #[cfg(feature = "gpu")]
    let mut gpu = crate::gpu::GpuAccelerator::try_new();
    #[cfg(feature = "gpu")]
    if gpu.is_some() {
        info!("GPU acceleration enabled for grid operations");
    }
    // Prepare persistent GPU buffers if GPU is available
    #[cfg(feature = "gpu")]
    if let Some(ref mut g) = gpu {
        g.prepare_buffers(ctx.n_grid, &ctx.g_squared);
    }

    // Initial density: superposition of atomic densities (SAD)
    let init_config = initial_density::InitialDensityConfig::non_magnetic(ctx.crystal.atoms.len());
    let mut rho_r = initial_density::generate_initial_density(
        ctx.crystal, &mut ctx.grid, ctx.pseudopotentials, ctx.n_electrons, &init_config,
    );
    let mut rho_g = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
    density_r_to_g(&mut ctx.grid.fft, &rho_r, &mut rho_g);
    info!("Initial density: superposition of atomic densities (Gaussian model)");

    let mut mixer = mixing::Mixer::new(
        ctx.params.mixing_beta,
        ctx.params.mixing_ndim,
        &ctx.params.mixing_mode,
        Some(&ctx.g_squared),
        ctx.n_electrons,
        ctx.omega,
    );
    let mut eigenvalues_all: Vec<Vec<f64>>;
    let mut fermi_energy;
    let mut e_prev: Option<f64> = None;
    let mut last_delta = f64::INFINITY;
    let pb = scf_progress_bar(ctx.params.max_iter);

    for iter in 0..ctx.params.max_iter {
        // Steps 1-3: Hartree, XC, V_eff assembly.
        // GPU path uses f32 for Hartree, XC, and V_eff; CPU path uses f64 + rayon.
        // XC FFT normalization is shared between both paths.

        // 1. Hartree potential
        #[cfg(feature = "gpu")]
        let v_h_fft = if let Some(ref gpu) = gpu {
            let fourpi_e2 = 4.0 * std::f64::consts::PI * crate::consts::E2_COULOMB;
            gpu.hartree_potential(&rho_g, &ctx.g_squared, fourpi_e2)
        } else {
            hartree_on_fft_grid(&rho_g, &ctx.g_squared)
        };
        #[cfg(not(feature = "gpu"))]
        let v_h_fft = hartree_on_fft_grid(&rho_g, &ctx.g_squared);

        // 2. XC potential: compute in real space, FFT to G-space
        // NLCC: add core density to valence density for XC evaluation
        let rho_for_xc = add_core_density(&rho_r, &ctx.rho_core_r);
        #[cfg(feature = "gpu")]
        let (exc_r_in, vxc_r) = if let Some(ref gpu) = gpu {
            gpu.lda_xc(&rho_for_xc)
        } else {
            xc::lda_xc_grid(&rho_for_xc)
        };
        #[cfg(not(feature = "gpu"))]
        let (exc_r_in, vxc_r) = xc::lda_xc_grid(&rho_for_xc);

        let vxc_g = real_to_g_space(&vxc_r, &mut ctx.grid.fft);

        // 3. V_eff = V_local + V_H + V_xc
        #[cfg(feature = "gpu")]
        let v_eff_fft = if let Some(ref gpu) = gpu {
            gpu.v_eff_assembly(&ctx.v_local_fft, &v_h_fft, &vxc_g)
        } else {
            assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_g)
        };
        #[cfg(not(feature = "gpu"))]
        let v_eff_fft = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_g);

        // 4. Solve eigenvalue problem at each k-point (parallel over k-points)
        let kpoint_results: Result<Vec<_>> = ctx.kpoints
            .par_iter()
            .enumerate()
            .map(|(ik, kp)| {
                let mut h = build_hamiltonian_with_v_eff(ctx.basis, &kp.k, &v_eff_fft, ctx.grid.dims);
                ctx.vnl_cache[ik].add_to_hamiltonian(&mut h, ctx.crystal, ctx.basis, &kp.k);
                dense::diagonalize_lowest(&h, ctx.params.n_bands)
            })
            .collect();
        let kpoint_results = kpoint_results?;

        eigenvalues_all = kpoint_results.iter().map(|r| r.eigenvalues.clone()).collect();
        let all_kpoint_wavefns: Vec<_> = kpoint_results.into_iter().map(|r| r.eigenvectors).collect();

        // 5. Occupations (configurable smearing scheme)
        fermi_energy = smearing::find_fermi_energy(
            &eigenvalues_all,
            &ctx.kpt_weights,
            ctx.n_electrons,
            ctx.params.smearing_sigma,
            ctx.params.smearing_scheme,
            ctx.spin_factor,
        );
        let occupations = compute_occupations(
            &eigenvalues_all, ctx.params.smearing_scheme, fermi_energy,
            ctx.params.smearing_sigma, ctx.spin_factor,
        );

        // 6. New density
        let mut rho_r_new = density::compute_density(
            &mut density::DensityGrid {
                basis: ctx.basis, g_to_fft: &ctx.g_to_fft, fft: &mut ctx.grid.fft,
                n_electrons: ctx.n_electrons, omega: ctx.omega,
            },
            ctx.kpoints, &all_kpoint_wavefns, &occupations,
        );

        // 6b. Symmetrize density. `symmetrize_density` short-circuits for a
        //     trivial group (identity-only), so this call is a no-op when the
        //     user disabled symmetry — bit-identical to the legacy skip.
        crate::symmetry::density::symmetrize_density(&mut rho_r_new, ctx.grid.dims, ctx.symmetry);

        // 7. Convergence check (dual criterion: density AND energy)
        let delta = density_diff(&rho_r, &rho_r_new, ctx.omega, ctx.n_grid);
        last_delta = delta;

        // Compute energy every iteration for convergence monitoring
        let mut rho_g_new = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
        density_r_to_g(&mut ctx.grid.fft, &rho_r_new, &mut rho_g_new);

        let rho_new_for_xc = add_core_density(&rho_r_new, &ctx.rho_core_r);
        let (exc_r, vxc_r_energy) = xc::lda_xc_grid(&rho_new_for_xc);

        let e_band = band_energy(&eigenvalues_all, &occupations, &ctx.kpt_weights);

        // Kohn-Sham energy: double-counting from OUTPUT density
        let e_total = total_energy(
            e_band,
            hartree_energy(&rho_g_new, &ctx.g_squared, ctx.omega),
            xc_energy_corrected(&rho_new_for_xc, &rho_r_new, &exc_r, &vxc_r_energy, ctx.omega),
            ctx.e_ewald,
        ) + ctx.v_local_g0 * ctx.n_electrons;

        // Harris-Foulkes energy: double-counting from INPUT density
        // rho_g, rho_for_xc, exc_r_in, vxc_r are all from the input density
        let e_harris = harris_foulkes_energy(
            e_band,
            hartree_energy(&rho_g, &ctx.g_squared, ctx.omega),
            xc_energy_corrected(&rho_for_xc, &rho_r, &exc_r_in, &vxc_r, ctx.omega),
            ctx.e_ewald,
        ) + ctx.v_local_g0 * ctx.n_electrons;

        let hf_diff = (e_harris - e_total).abs();

        let de = e_prev.map(|ep| (e_total - ep).abs());
        e_prev = Some(e_total);

        let rho_converged = delta < ctx.params.conv_threshold;
        let energy_converged = de.is_some_and(|de| de < ctx.params.energy_threshold);

        pb.set_position((iter + 1) as u64);
        pb.set_message(format!(
            "E={e_total:.4} eV  Δρ={delta:.1e}"
        ));
        info!(
            "SCF iter {:>3}: E_KS={:.6} eV  E_HF={:.6} eV  |HF-KS|={:.2e}  dE={:>10}  Δρ={:.2e}",
            iter + 1, e_total, e_harris, hf_diff,
            de.map_or("N/A".to_string(), |de| format!("{de:.2e}")),
            delta
        );

        // Warn if density converged but Harris-Foulkes difference is large
        if rho_converged && energy_converged && hf_diff > 0.01 {
            log::warn!(
                "Density converged but |E_HF - E_KS| = {hf_diff:.2e} eV — \
                 energy may not be reliable. Consider tightening conv_threshold."
            );
        }

        if rho_converged && energy_converged {
            pb.finish_and_clear();
            info!("SCF converged after {} iterations", iter + 1);
            rho_g = rho_g_new;
            let rho_g_basis: Vec<Complex64> = ctx.g_to_fft.iter().map(|&idx| rho_g[idx]).collect();

            // Entropy and free energy
            let ts = smearing::entropy_ts(
                &eigenvalues_all, &ctx.kpt_weights, fermi_energy,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );
            let free_energy = e_total - ts;
            let energy_sigma0 = f64::midpoint(e_total, free_energy);

            // -------------------------------------------------------------
            // VGC5 per-component decomposition (diagnostic).
            // Computed once at convergence; final wavefunctions are alive here.
            // -------------------------------------------------------------
            let k_vecs: Vec<nalgebra::Vector3<f64>> =
                ctx.kpoints.iter().map(|kp| kp.k).collect();

            let e_kinetic = kinetic_expectation(
                ctx.basis, &k_vecs, &ctx.kpt_weights,
                &all_kpoint_wavefns, &occupations,
            );

            // V_local in real space (G=0 already zeroed in v_local_fft).
            let mut v_local_cplx = ctx.v_local_fft.clone();
            ctx.grid.fft.inverse(&mut v_local_cplx);
            let v_local_r: Vec<f64> = v_local_cplx.iter().map(|c| c.re).collect();
            let e_local = local_pp_energy_grid(&rho_r_new, &v_local_r, ctx.omega, ctx.n_grid);
            let e_local_g0_shift = ctx.v_local_g0 * ctx.n_electrons;

            let e_nonlocal = nonlocal_expectation(
                ctx.basis, ctx.crystal, &k_vecs, &ctx.kpt_weights,
                &all_kpoint_wavefns, &occupations, &ctx.vnl_cache,
            );

            // NB: `rho_g` was assigned `rho_g_new` in the convergence branch
            // above; use `rho_g` here (the final-iteration density in G-space).
            let e_hartree_term = hartree_energy(&rho_g, &ctx.g_squared, ctx.omega);
            let e_xc_term = xc_energy_bare(&rho_new_for_xc, &exc_r, ctx.omega);
            let e_ewald_term = ctx.e_ewald;

            let components = EnergyComponents {
                e_band,
                e_kinetic,
                e_local,
                e_local_g0_shift,
                e_nonlocal,
                e_hartree: e_hartree_term,
                e_xc: e_xc_term,
                e_ewald: e_ewald_term,
            };

            let e_sum = e_kinetic + e_local + e_local_g0_shift + e_nonlocal
                + e_hartree_term + e_xc_term + e_ewald_term;
            info!("Energy (E_KS):   {e_total:.6} eV");
            info!("Harris-Foulkes:  {e_harris:.6} eV  (|HF-KS|={hf_diff:.2e})");
            info!("Free energy (F): {free_energy:.6} eV");
            info!("E sigma→0 (E₀):  {energy_sigma0:.6} eV");
            if ts.abs() > 1e-8 {
                info!(
                    "Entropy (-TS):   {:.6} eV ({:.3} meV/atom)",
                    -ts, -ts * 1000.0 / ctx.crystal.atoms.len() as f64
                );
            }
            info!("--- Per-component energies (eV) ---");
            info!("  E_band       = {e_band:.6}");
            info!("  E_kinetic    = {e_kinetic:.6}");
            info!("  E_local      = {e_local:.6}");
            info!("  E_local(G=0) = {e_local_g0_shift:.6}  (= V_loc(G=0)·N_el, N_el={:.3})", ctx.n_electrons);
            info!("  E_nonlocal   = {e_nonlocal:.6}");
            info!("  E_hartree    = {e_hartree_term:.6}");
            info!("  E_xc         = {e_xc_term:.6}");
            info!("  E_ewald      = {e_ewald_term:.6}");
            info!("  E_sum(comp)  = {e_sum:.6}   (vs E_KS {e_total:.6}, Δ={:.2e})",
                  e_sum - e_total);

            return Ok(ScfResult {
                total_energy: e_total,
                harris_foulkes_energy: e_harris,
                free_energy,
                energy_sigma0,
                entropy_ts: ts,
                eigenvalues: eigenvalues_all,
                fermi_energy,
                n_iterations: iter + 1,
                rho_g: rho_g_basis,
                magnetization: 0.0,
                nspin: 1,
                components,
            });
        }

        // 8. Mix (Kerker preconditioning applied inside if enabled)
        rho_r = mixer.mix(&rho_r, &rho_r_new, &mut ctx.grid.fft);
        density_r_to_g(&mut ctx.grid.fft, &rho_r, &mut rho_g);
    }

    pb.abandon_with_message("did not converge");
    Err(PwdftError::ConvergenceFailure {
        iterations: ctx.params.max_iter,
        delta: last_delta,
    })
}

/// Spin-polarized SCF loop (nspin=2).
///
/// Two spin channels with independent densities, XC potentials, and
/// Hamiltonians. Hartree and V_local are spin-independent (computed from
/// total density ρ↑ + ρ↓). V_xc is spin-dependent via LSDA.
/// NLCC core charge is split equally between channels: ρ_core/2 per spin.
fn run_scf_spin(
    crystal: &Crystal,
    basis: &BasisSet,
    kpoints: &[KPoint],
    pseudopotentials: &[&PseudopotentialData],
    params: &ScfParams,
    symmetry: &crate::symmetry::SymmetryInfo,
) -> Result<ScfResult> {
    let mut ctx = context::ScfContext::new(crystal, basis, kpoints, pseudopotentials, params, symmetry)?;

    // Determine initial spin split from starting_magnetization
    let per_atom_mag: Vec<f64> = ctx.crystal.atoms.iter().map(|a| {
        let sym = crate::atoms::from_z(a.z).map(|e| e.symbol().to_string()).unwrap_or_default();
        *ctx.params.starting_magnetization.get(&sym).unwrap_or(&0.0)
    }).collect();

    // Determine n_up, n_down
    let (n_up, n_down) = if let Some(tot_mag) = ctx.params.tot_magnetization {
        let n_up = f64::midpoint(ctx.n_electrons, tot_mag);
        let n_down = (ctx.n_electrons - tot_mag) / 2.0;
        info!("Fixed magnetization: n_up={n_up:.2}, n_down={n_down:.2}");
        (n_up, n_down)
    } else {
        (ctx.n_electrons / 2.0, ctx.n_electrons / 2.0) // initial guess; self-consistent
    };

    // Initial spin density via SAD with magnetic moments
    let init_config = initial_density::InitialDensityConfig {
        magnetic_moments: per_atom_mag.clone(),
        gaussian_sigma: None,
    };
    // Generate total density then split
    let rho_total = initial_density::generate_initial_density(
        ctx.crystal, &mut ctx.grid, ctx.pseudopotentials, ctx.n_electrons, &init_config,
    );
    // Split: rho_up = (1+m)/2 * rho, rho_down = (1-m)/2 * rho
    // For a uniform initial split, m=0 -> equal channels
    let avg_mag: f64 = per_atom_mag.iter().sum::<f64>() / per_atom_mag.len().max(1) as f64;
    let mut rho_up_r: Vec<f64> = rho_total.iter().map(|&r| r * (1.0 + avg_mag) / 2.0).collect();
    let mut rho_down_r: Vec<f64> = rho_total.iter().map(|&r| r * (1.0 - avg_mag) / 2.0).collect();

    // Two mixers (one per spin channel)
    let mut mixer_up = mixing::Mixer::new(
        ctx.params.mixing_beta, ctx.params.mixing_ndim, &ctx.params.mixing_mode,
        Some(&ctx.g_squared), n_up, ctx.omega,
    );
    let mut mixer_down = mixing::Mixer::new(
        ctx.params.mixing_beta, ctx.params.mixing_ndim, &ctx.params.mixing_mode,
        Some(&ctx.g_squared), n_down, ctx.omega,
    );

    let mut e_prev: Option<f64> = None;
    let mut last_delta = f64::INFINITY;
    let pb = scf_progress_bar(ctx.params.max_iter);

    // Precompute half core density (constant across iterations)
    let rho_core_half: Vec<f64> = ctx.rho_core_r.iter().map(|&c| c / 2.0).collect();

    for iter in 0..ctx.params.max_iter {
        // Total density for Hartree (spin-independent)
        let rho_total_r: Vec<f64> = rho_up_r.iter().zip(rho_down_r.iter())
            .map(|(&u, &d)| u + d).collect();
        let mut rho_total_g = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
        density_r_to_g(&mut ctx.grid.fft, &rho_total_r, &mut rho_total_g);

        // 1. Hartree from total density
        let v_h_fft = hartree_on_fft_grid(&rho_total_g, &ctx.g_squared);

        // 2. Spin-dependent XC
        let rho_up_xc = add_core_density(&rho_up_r, &rho_core_half);
        let rho_down_xc = add_core_density(&rho_down_r, &rho_core_half);
        let (exc_r, vxc_up_r, vxc_down_r) = xc::lda_xc_spin_grid(&rho_up_xc, &rho_down_xc);

        let vxc_up_g = real_to_g_space(&vxc_up_r, &mut ctx.grid.fft);
        let vxc_down_g = real_to_g_space(&vxc_down_r, &mut ctx.grid.fft);

        // 3. Two V_eff: V_local + V_H + V_xc_sigma
        let v_eff_up = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_up_g);
        let v_eff_down = assemble_v_eff(&ctx.v_local_fft, &v_h_fft, &vxc_down_g);

        // 4. Diagonalize both spins at each k-point
        let kpoint_results_up: Result<Vec<_>> = ctx.kpoints.par_iter().enumerate().map(|(ik, kp)| {
            let mut h = build_hamiltonian_with_v_eff(ctx.basis, &kp.k, &v_eff_up, ctx.grid.dims);
            ctx.vnl_cache[ik].add_to_hamiltonian(&mut h, ctx.crystal, ctx.basis, &kp.k);
            dense::diagonalize_lowest(&h, ctx.params.n_bands)
        }).collect();
        let kpoint_results_up = kpoint_results_up?;

        let kpoint_results_down: Result<Vec<_>> = ctx.kpoints.par_iter().enumerate().map(|(ik, kp)| {
            let mut h = build_hamiltonian_with_v_eff(ctx.basis, &kp.k, &v_eff_down, ctx.grid.dims);
            ctx.vnl_cache[ik].add_to_hamiltonian(&mut h, ctx.crystal, ctx.basis, &kp.k);
            dense::diagonalize_lowest(&h, ctx.params.n_bands)
        }).collect();
        let kpoint_results_down = kpoint_results_down?;

        // Flatten eigenvalues: [up_k0, up_k1, ..., down_k0, down_k1, ...]
        let eig_up: Vec<Vec<f64>> = kpoint_results_up.iter().map(|r| r.eigenvalues.clone()).collect();
        let eig_down: Vec<Vec<f64>> = kpoint_results_down.iter().map(|r| r.eigenvalues.clone()).collect();
        let eigenvalues_all: Vec<Vec<f64>> = eig_up.iter().chain(eig_down.iter()).cloned().collect();
        let weights_all: Vec<f64> = ctx.kpt_weights.iter().chain(ctx.kpt_weights.iter()).copied().collect();

        let wfn_up: Vec<_> = kpoint_results_up.into_iter().map(|r| r.eigenvectors).collect();
        let wfn_down: Vec<_> = kpoint_results_down.into_iter().map(|r| r.eigenvectors).collect();

        // 5. Fermi energy and occupations
        let (fermi_energy, occ_up, occ_down) = if let Some(tot_mag) = ctx.params.tot_magnetization {
            // Fixed magnetization: separate Fermi energies per spin
            let n_up_target = f64::midpoint(ctx.n_electrons, tot_mag);
            let n_down_target = (ctx.n_electrons - tot_mag) / 2.0;

            let ef_up = smearing::find_fermi_energy(
                &eig_up, &ctx.kpt_weights, n_up_target,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );
            let ef_down = smearing::find_fermi_energy(
                &eig_down, &ctx.kpt_weights, n_down_target,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );

            let occ_up = compute_occupations(
                &eig_up, ctx.params.smearing_scheme, ef_up,
                ctx.params.smearing_sigma, ctx.spin_factor,
            );
            let occ_down = compute_occupations(
                &eig_down, ctx.params.smearing_scheme, ef_down,
                ctx.params.smearing_sigma, ctx.spin_factor,
            );

            // Report average Fermi energy
            (f64::midpoint(ef_up, ef_down), occ_up, occ_down)
        } else {
            // Free magnetization: single Fermi energy for both spins
            let fermi_energy = smearing::find_fermi_energy(
                &eigenvalues_all, &weights_all, ctx.n_electrons,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );

            let occ_up = compute_occupations(
                &eig_up, ctx.params.smearing_scheme, fermi_energy,
                ctx.params.smearing_sigma, ctx.spin_factor,
            );
            let occ_down = compute_occupations(
                &eig_down, ctx.params.smearing_scheme, fermi_energy,
                ctx.params.smearing_sigma, ctx.spin_factor,
            );

            (fermi_energy, occ_up, occ_down)
        };

        // 6. Reconstruct spin densities
        let n_el_up: f64 = occ_up.iter().zip(ctx.kpt_weights.iter())
            .flat_map(|(occs, &w)| occs.iter().map(move |&f| f * w))
            .sum();
        let n_el_down: f64 = occ_down.iter().zip(ctx.kpt_weights.iter())
            .flat_map(|(occs, &w)| occs.iter().map(move |&f| f * w))
            .sum();

        let rho_up_new = density::compute_density(
            &mut density::DensityGrid {
                basis: ctx.basis, g_to_fft: &ctx.g_to_fft, fft: &mut ctx.grid.fft,
                n_electrons: n_el_up, omega: ctx.omega,
            },
            ctx.kpoints, &wfn_up, &occ_up,
        );
        let rho_down_new = density::compute_density(
            &mut density::DensityGrid {
                basis: ctx.basis, g_to_fft: &ctx.g_to_fft, fft: &mut ctx.grid.fft,
                n_electrons: n_el_down, omega: ctx.omega,
            },
            ctx.kpoints, &wfn_down, &occ_down,
        );

        // Symmetrize each channel. Trivial (identity-only) groups
        // short-circuit inside `symmetrize_density`, preserving the legacy
        // "no symmetrization" behavior bit-identically.
        let mut rho_up_sym = rho_up_new;
        let mut rho_down_sym = rho_down_new;
        crate::symmetry::density::symmetrize_density(&mut rho_up_sym, ctx.grid.dims, ctx.symmetry);
        crate::symmetry::density::symmetrize_density(&mut rho_down_sym, ctx.grid.dims, ctx.symmetry);

        // 7. Convergence check
        let rho_total_new: Vec<f64> = rho_up_sym.iter().zip(rho_down_sym.iter())
            .map(|(&u, &d)| u + d).collect();
        // SPNC: use per-spin max, not total-density diff. A spin-flip fluctuation
        // (+ε in rho_up, −ε in rho_down) is invisible to the total but keeps
        // zeta_in ≠ zeta_out, which leaves E_xc[rho, zeta] inconsistent and
        // spoils the O(Δρ²) convergence of |E_HF − E_KS|. Per-spin max is
        // strictly stronger than total and keeps the scalar `conv_threshold`
        // semantics unchanged. See proposals/SPNC-spin-per-density-convergence.md.
        let delta_up = density_diff(&rho_up_r, &rho_up_sym, ctx.omega, ctx.n_grid);
        let delta_down = density_diff(&rho_down_r, &rho_down_sym, ctx.omega, ctx.n_grid);
        let delta = delta_up.max(delta_down);
        last_delta = delta;

        // Energy from new density
        let mut rho_total_new_g = vec![Complex64::new(0.0, 0.0); ctx.n_grid];
        density_r_to_g(&mut ctx.grid.fft, &rho_total_new, &mut rho_total_new_g);

        let rho_xc_total: Vec<f64> = add_core_density(&rho_total_new, &ctx.rho_core_r);
        let occ_all: Vec<Vec<f64>> = occ_up.iter().chain(occ_down.iter()).cloned().collect();

        // SPXC fix: recompute spin-polarized XC from OUTPUT spin densities for E_KS.
        // Previously `exc_r` (INPUT density, depends on zeta_in) was integrated against
        // `rho_xc_total` (OUTPUT total density), which is not physically meaningful and
        // introduced an O(delta_rho) error that spoiled the quadratic convergence of
        // |E_HF - E_KS|. The non-spin run_scf already recomputes XC from the output
        // density (see `let (exc_r, vxc_r_energy) = xc::lda_xc_grid(&rho_new_for_xc)`).
        // INPUT-based quantities (exc_r, vxc_up_r, vxc_down_r) remain for E_HF below.
        let rho_up_xc_out = add_core_density(&rho_up_sym, &rho_core_half);
        let rho_down_xc_out = add_core_density(&rho_down_sym, &rho_core_half);
        let (exc_r_out, vxc_up_r_out, vxc_down_r_out) =
            xc::lda_xc_spin_grid(&rho_up_xc_out, &rho_down_xc_out);

        // Spin XC double-counting (OUTPUT density): E_vxc = integral(V_xc_up rho_up_out + V_xc_down rho_down_out) dr
        // Both V_xc and rho_sigma here come from the OUTPUT (symmetrized) density,
        // matching the non-spin run_scf convention.
        let dvol = ctx.omega / ctx.n_grid as f64;
        let e_vxc_spin_out: f64 = rho_up_sym.iter().zip(vxc_up_r_out.iter())
            .zip(rho_down_sym.iter().zip(vxc_down_r_out.iter()))
            .map(|((&ru, &vu), (&rd, &vd))| (ru * vu + rd * vd) * dvol)
            .sum();
        let e_xc_out = xc::lda_xc_energy(&rho_xc_total, &exc_r_out, ctx.omega);
        let e_xc_corrected_out = e_xc_out - e_vxc_spin_out;

        let e_band = band_energy(&eigenvalues_all, &occ_all, &weights_all);

        // Kohn-Sham energy: double-counting from OUTPUT density
        let e_total = total_energy(
            e_band,
            hartree_energy(&rho_total_new_g, &ctx.g_squared, ctx.omega),
            e_xc_corrected_out,
            ctx.e_ewald,
        ) + ctx.v_local_g0 * ctx.n_electrons;

        // Harris-Foulkes energy: double-counting from INPUT density
        // rho_total_g, rho_up_xc, rho_down_xc, exc_r, vxc_up_r, vxc_down_r
        // are all from the input density
        let rho_xc_total_in: Vec<f64> = add_core_density(&rho_total_r, &ctx.rho_core_r);
        let e_xc_in = xc::lda_xc_energy(&rho_xc_total_in, &exc_r, ctx.omega);
        let e_vxc_spin_in: f64 = rho_up_r.iter().zip(vxc_up_r.iter())
            .zip(rho_down_r.iter().zip(vxc_down_r.iter()))
            .map(|((&ru, &vu), (&rd, &vd))| (ru * vu + rd * vd) * dvol)
            .sum();
        let e_xc_corrected_in = e_xc_in - e_vxc_spin_in;

        let e_harris = harris_foulkes_energy(
            e_band,
            hartree_energy(&rho_total_g, &ctx.g_squared, ctx.omega),
            e_xc_corrected_in,
            ctx.e_ewald,
        ) + ctx.v_local_g0 * ctx.n_electrons;

        let hf_diff = (e_harris - e_total).abs();

        let de = e_prev.map(|ep| (e_total - ep).abs());
        e_prev = Some(e_total);

        let rho_converged = delta < ctx.params.conv_threshold;
        let energy_converged = de.is_some_and(|de| de < ctx.params.energy_threshold);

        let mag = (n_el_up - n_el_down).abs();
        pb.set_position((iter + 1) as u64);
        pb.set_message(format!(
            "E={e_total:.4} eV  Δρ={delta:.1e}  M={mag:.2} μB"
        ));
        info!(
            "SCF iter {:>3}: E_KS={:.6} eV  E_HF={:.6} eV  |HF-KS|={:.2e}  dE={:>10}  Δρ={:.2e} (↑{:.2e} ↓{:.2e})  M={:.3} μB",
            iter + 1, e_total, e_harris, hf_diff,
            de.map_or("N/A".to_string(), |de| format!("{de:.2e}")),
            delta, delta_up, delta_down, mag
        );

        // Warn if density converged but Harris-Foulkes difference is large
        if rho_converged && energy_converged && hf_diff > 0.01 {
            log::warn!(
                "Density converged but |E_HF - E_KS| = {hf_diff:.2e} eV — \
                 energy may not be reliable. Consider tightening conv_threshold."
            );
        }

        if rho_converged && energy_converged {
            pb.finish_and_clear();
            info!("SCF converged after {} iterations", iter + 1);
            info!("Magnetization: {mag:.4} μB ({n_el_up:.4} up, {n_el_down:.4} down)");

            let rho_g_basis: Vec<Complex64> = ctx.g_to_fft.iter().map(|&idx| rho_total_new_g[idx]).collect();

            let ts = smearing::entropy_ts(
                &eigenvalues_all, &weights_all, fermi_energy,
                ctx.params.smearing_sigma, ctx.params.smearing_scheme, ctx.spin_factor,
            );
            let free_energy = e_total - ts;
            let energy_sigma0 = f64::midpoint(e_total, free_energy);

            // -------------------------------------------------------------
            // VGC5 per-component decomposition (diagnostic, spin-polarized).
            // Kinetic / non-local: sum over both spin channels.
            // Local / Hartree / Ewald: built from total density (spin-
            //   independent operators).
            // XC: ∫(ρ_up+ρ_down+ρ_core)·ε_xc(ρ_up,ρ_down)dr (bare, from OUTPUT).
            // -------------------------------------------------------------
            let k_vecs: Vec<nalgebra::Vector3<f64>> =
                ctx.kpoints.iter().map(|kp| kp.k).collect();

            let e_kin_up = kinetic_expectation(
                ctx.basis, &k_vecs, &ctx.kpt_weights, &wfn_up, &occ_up,
            );
            let e_kin_down = kinetic_expectation(
                ctx.basis, &k_vecs, &ctx.kpt_weights, &wfn_down, &occ_down,
            );
            let e_kinetic = e_kin_up + e_kin_down;

            let mut v_local_cplx = ctx.v_local_fft.clone();
            ctx.grid.fft.inverse(&mut v_local_cplx);
            let v_local_r: Vec<f64> = v_local_cplx.iter().map(|c| c.re).collect();
            let e_local = local_pp_energy_grid(&rho_total_new, &v_local_r, ctx.omega, ctx.n_grid);
            let e_local_g0_shift = ctx.v_local_g0 * ctx.n_electrons;

            let e_nl_up = nonlocal_expectation(
                ctx.basis, ctx.crystal, &k_vecs, &ctx.kpt_weights,
                &wfn_up, &occ_up, &ctx.vnl_cache,
            );
            let e_nl_down = nonlocal_expectation(
                ctx.basis, ctx.crystal, &k_vecs, &ctx.kpt_weights,
                &wfn_down, &occ_down, &ctx.vnl_cache,
            );
            let e_nonlocal = e_nl_up + e_nl_down;

            let e_hartree_term = hartree_energy(&rho_total_new_g, &ctx.g_squared, ctx.omega);
            let e_xc_term = xc_energy_bare(&rho_xc_total, &exc_r_out, ctx.omega);
            let e_ewald_term = ctx.e_ewald;

            let components = EnergyComponents {
                e_band,
                e_kinetic,
                e_local,
                e_local_g0_shift,
                e_nonlocal,
                e_hartree: e_hartree_term,
                e_xc: e_xc_term,
                e_ewald: e_ewald_term,
            };

            let e_sum = e_kinetic + e_local + e_local_g0_shift + e_nonlocal
                + e_hartree_term + e_xc_term + e_ewald_term;
            info!("Energy (E_KS):   {e_total:.6} eV");
            info!("Harris-Foulkes:  {e_harris:.6} eV  (|HF-KS|={hf_diff:.2e})");
            info!("Free energy (F): {free_energy:.6} eV");
            info!("E sigma→0 (E₀):  {energy_sigma0:.6} eV");
            info!("--- Per-component energies (eV) ---");
            info!("  E_band       = {e_band:.6}");
            info!("  E_kinetic    = {e_kinetic:.6}");
            info!("  E_local      = {e_local:.6}");
            info!("  E_local(G=0) = {e_local_g0_shift:.6}  (= V_loc(G=0)·N_el, N_el={:.3})", ctx.n_electrons);
            info!("  E_nonlocal   = {e_nonlocal:.6}");
            info!("  E_hartree    = {e_hartree_term:.6}");
            info!("  E_xc         = {e_xc_term:.6}");
            info!("  E_ewald      = {e_ewald_term:.6}");
            info!("  E_sum(comp)  = {e_sum:.6}   (vs E_KS {e_total:.6}, Δ={:.2e})",
                  e_sum - e_total);

            return Ok(ScfResult {
                total_energy: e_total,
                harris_foulkes_energy: e_harris,
                free_energy,
                energy_sigma0,
                entropy_ts: ts,
                eigenvalues: eigenvalues_all,
                fermi_energy,
                n_iterations: iter + 1,
                rho_g: rho_g_basis,
                magnetization: mag,
                nspin: 2,
                components,
            });
        }

        // 8. Mix each spin channel independently
        rho_up_r = mixer_up.mix(&rho_up_r, &rho_up_sym, &mut ctx.grid.fft);
        rho_down_r = mixer_down.mix(&rho_down_r, &rho_down_sym, &mut ctx.grid.fft);
    }

    pb.abandon_with_message("did not converge");
    Err(PwdftError::ConvergenceFailure {
        iterations: ctx.params.max_iter,
        delta: last_delta,
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
