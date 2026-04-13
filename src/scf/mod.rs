pub mod density;
pub mod mixing;
pub mod smearing;

use log::info;
use nalgebra::Vector3;
use num_complex::Complex64;

use crate::{
    basis::BasisSet,
    consts::HBAR2_OVER_2M,
    crystal::Crystal,
    eigensolver::dense,
    error::{PwdftError, Result},
    fft::{fft_grid_size, FFT3D},
    kpoints::KPoint,
    potential::{hartree, local::LocalPotential, nonlocal::NonlocalPotential, xc},
    pseudopotential::PseudopotentialData,
};

/// Parameters for an SCF calculation.
pub struct ScfParams {
    pub n_bands: usize,
    pub max_iter: usize,
    pub conv_threshold: f64,
    pub mixing_beta: f64,
    pub smearing_sigma: f64,
}

impl Default for ScfParams {
    fn default() -> Self {
        Self {
            n_bands: 8,
            max_iter: 100,
            conv_threshold: 1e-6,
            mixing_beta: 0.3,
            smearing_sigma: 0.01, // eV
        }
    }
}

/// Result of an SCF calculation.
pub struct ScfResult {
    /// Converged total energy (eV).
    pub total_energy: f64,
    /// Eigenvalues at each k-point (eV).
    pub eigenvalues: Vec<Vec<f64>>,
    /// Fermi energy (eV).
    pub fermi_energy: f64,
    /// Number of SCF iterations performed.
    pub n_iterations: usize,
    /// Final charge density in G-space.
    pub rho_g: Vec<Complex64>,
}

/// Run the self-consistent field loop.
pub fn run_scf(
    crystal: &Crystal,
    basis: &BasisSet,
    kpoints: &[KPoint],
    pseudopotentials: &[&PseudopotentialData],
    params: &ScfParams,
) -> Result<ScfResult> {
    let omega = crystal.lattice.volume().abs();
    let n_electrons: f64 = crystal
        .atoms
        .iter()
        .map(|a| {
            let pp = pseudopotentials
                .iter()
                .find(|pp| {
                    crate::atoms::Element::from_symbol(&pp.element)
                        .map_or(false, |e| e.atomic_number() == a.z)
                })
                .unwrap();
            pp.z_valence
        })
        .sum();

    info!("SCF: {n_electrons} electrons, {omega:.3} ų cell volume");

    // Determine FFT grid dimensions from basis
    let miller = basis.miller_indices();
    let n_max: Vec<i32> = (0..3)
        .map(|dim| miller.iter().map(|m| m[dim].abs()).max().unwrap_or(0))
        .collect();
    let grid = [
        fft_grid_size(n_max[0]),
        fft_grid_size(n_max[1]),
        fft_grid_size(n_max[2]),
    ];
    let fft = FFT3D::new(grid[0], grid[1], grid[2]);
    let n_grid = fft.total_size();
    info!("FFT grid: {}×{}×{} = {} points", grid[0], grid[1], grid[2], n_grid);

    // Precompute local pseudopotential
    let v_local = LocalPotential::new(crystal, basis, pseudopotentials);

    // Mapping from G-vector Miller indices to FFT grid index
    let g_to_fft: Vec<usize> = miller
        .iter()
        .map(|&[n1, n2, n3]| {
            let i1 = ((n1 % grid[0] as i32) + grid[0] as i32) as usize % grid[0];
            let i2 = ((n2 % grid[1] as i32) + grid[1] as i32) as usize % grid[1];
            let i3 = ((n3 % grid[2] as i32) + grid[2] as i32) as usize % grid[2];
            i1 * grid[1] * grid[2] + i2 * grid[2] + i3
        })
        .collect();

    // Initial density guess: uniform density
    let rho_init = n_electrons / omega;
    let mut rho_r = vec![rho_init; n_grid];
    let mut rho_g = vec![Complex64::new(0.0, 0.0); n_grid];

    // FFT rho_r → rho_g
    density_r_to_g(&fft, &rho_r, &mut rho_g);

    let mut mixer = mixing::AndersonMixer::new(params.mixing_beta, 4, n_grid);
    let mut eigenvalues_all = Vec::new();
    let mut fermi_energy;

    for iter in 0..params.max_iter {
        // Build effective potential on the FFT grid
        // 1. Hartree potential (G-space)
        let v_h_g = hartree_on_basis(&rho_g, basis, &g_to_fft, n_grid);

        // 2. XC potential (real space)
        let (_exc_r, vxc_r) = xc::lda_xc_grid(&rho_r);

        // FFT V_xc to G-space
        let mut vxc_g = vec![Complex64::new(0.0, 0.0); n_grid];
        for (i, &v) in vxc_r.iter().enumerate() {
            vxc_g[i] = Complex64::new(v, 0.0);
        }
        fft.forward(&mut vxc_g);
        // Normalize: FFT convention
        let norm = 1.0 / n_grid as f64;
        for v in &mut vxc_g {
            *v *= norm;
        }

        // Build V_eff(G) on the basis G-vectors: V_local + V_H + V_xc
        let n_pw = basis.len();
        let mut v_eff_basis: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n_pw];
        for ig in 0..n_pw {
            let fft_idx = g_to_fft[ig];
            v_eff_basis[ig] = v_local.v_of_g(ig) + v_h_g[fft_idx] + vxc_g[fft_idx];
        }

        // Solve eigenvalue problem at each k-point
        eigenvalues_all.clear();
        let mut all_kpoint_wavefns = Vec::new();

        for kp in kpoints {
            let mut h = build_hamiltonian_with_potential(basis, &kp.k, &v_eff_basis);

            // Add non-local pseudopotential (k-dependent)
            let vnl = NonlocalPotential::new(crystal, basis, &kp.k, pseudopotentials);
            vnl.add_to_hamiltonian(&mut h, crystal, basis, &kp.k);

            let result = dense::diagonalize_lowest(&h, params.n_bands);
            eigenvalues_all.push(result.eigenvalues);
            all_kpoint_wavefns.push(result.eigenvectors);
        }

        // Compute occupations (Fermi-Dirac smearing)
        fermi_energy = smearing::find_fermi_energy(
            &eigenvalues_all,
            kpoints,
            n_electrons,
            params.smearing_sigma,
        );
        let occupations: Vec<Vec<f64>> = eigenvalues_all
            .iter()
            .map(|evs| {
                evs.iter()
                    .map(|&e| smearing::fermi_dirac(e, fermi_energy, params.smearing_sigma))
                    .collect()
            })
            .collect();

        // Compute new density
        let rho_r_new = density::compute_density(
            basis,
            kpoints,
            &all_kpoint_wavefns,
            &occupations,
            &g_to_fft,
            &fft,
            n_electrons,
            omega,
        );

        // Check convergence
        let delta = density_diff(&rho_r, &rho_r_new, omega, n_grid);

        info!(
            "SCF iter {}: E_fermi = {:.6} eV, delta_rho = {:.2e}",
            iter + 1,
            fermi_energy,
            delta
        );

        if delta < params.conv_threshold {
            info!("SCF converged after {} iterations", iter + 1);
            // Compute final rho_g
            rho_r = rho_r_new;
            density_r_to_g(&fft, &rho_r, &mut rho_g);
            // Extract rho_g on basis
            let rho_g_basis: Vec<Complex64> = g_to_fft.iter().map(|&idx| rho_g[idx]).collect();

            let total_energy = compute_total_energy(
                &eigenvalues_all,
                &occupations,
                kpoints,
                &rho_r,
                &rho_g_basis,
                basis,
                crystal,
                pseudopotentials,
                omega,
                n_grid,
            );

            return Ok(ScfResult {
                total_energy,
                eigenvalues: eigenvalues_all,
                fermi_energy,
                n_iterations: iter + 1,
                rho_g: rho_g_basis,
            });
        }

        // Mix densities
        rho_r = mixer.mix(&rho_r, &rho_r_new);
        density_r_to_g(&fft, &rho_r, &mut rho_g);
    }

    Err(PwdftError::ConvergenceFailure {
        iterations: params.max_iter,
        delta: density_diff(&rho_r, &rho_r, omega, n_grid),
    })
}

/// Build the local Hamiltonian (kinetic + local potential) at k-point k.
/// Non-local potential is added separately via NonlocalPotential::add_to_hamiltonian.
fn build_hamiltonian_with_potential(
    basis: &BasisSet,
    k: &Vector3<f64>,
    v_eff_basis: &[Complex64],
) -> nalgebra::DMatrix<Complex64> {
    let n = basis.len();
    let mut h = nalgebra::DMatrix::zeros(n, n);

    // Kinetic energy (diagonal)
    for (i, g) in basis.g_vectors().iter().enumerate() {
        let kpg = k + g;
        let ke = HBAR2_OVER_2M * kpg.norm_squared();
        h[(i, i)] = Complex64::new(ke, 0.0);
    }

    // Potential: V_{G,G'} = V_eff(G-G') looked up via basis Miller index difference
    let miller = basis.miller_indices();
    for i in 0..n {
        for j in 0..n {
            let dn = [
                miller[i][0] - miller[j][0],
                miller[i][1] - miller[j][1],
                miller[i][2] - miller[j][2],
            ];
            // Map difference to FFT grid index
            // We need grid dims — extract from g_to_fft mapping
            // This is inefficient but correct. Phase 4 will use FFT-based H|ψ⟩.
            if let Some(idx) = basis.index_of(dn[0], dn[1], dn[2]) {
                h[(i, j)] += v_eff_basis[idx];
            }
            // If G-G' is outside the basis, the potential contribution is zero
            // (this is the approximation inherent in the plane-wave cutoff)
        }
    }

    h
}

/// FFT density from real space to G-space.
fn density_r_to_g(fft: &FFT3D, rho_r: &[f64], rho_g: &mut [Complex64]) {
    for (i, &r) in rho_r.iter().enumerate() {
        rho_g[i] = Complex64::new(r, 0.0);
    }
    fft.forward(rho_g);
    let norm = 1.0 / fft.total_size() as f64;
    for v in rho_g.iter_mut() {
        *v *= norm;
    }
}

/// Compute Hartree potential on the FFT grid from rho_g.
fn hartree_on_basis(
    rho_g: &[Complex64],
    basis: &BasisSet,
    g_to_fft: &[usize],
    n_fft: usize,
) -> Vec<Complex64> {
    let fourpi_e2 = 4.0 * std::f64::consts::PI * hartree::E2;
    let mut v_h = vec![Complex64::new(0.0, 0.0); n_fft];

    // Only compute Hartree at G-vectors in the basis (the rest are zero)
    for (ig, g) in basis.g_vectors().iter().enumerate() {
        let g2 = g.norm_squared();
        let fft_idx = g_to_fft[ig];
        if g2 > 1e-20 {
            v_h[fft_idx] = rho_g[fft_idx] * fourpi_e2 / g2;
        }
    }

    v_h
}

/// Compute RMS density difference.
fn density_diff(rho_old: &[f64], rho_new: &[f64], omega: f64, n_grid: usize) -> f64 {
    let dvol = omega / n_grid as f64;
    let sum_sq: f64 = rho_old
        .iter()
        .zip(rho_new.iter())
        .map(|(&a, &b)| (a - b) * (a - b) * dvol)
        .sum();
    (sum_sq / omega).sqrt()
}

/// Compute total energy.
fn compute_total_energy(
    eigenvalues: &[Vec<f64>],
    occupations: &[Vec<f64>],
    kpoints: &[KPoint],
    rho_r: &[f64],
    rho_g_basis: &[Complex64],
    basis: &BasisSet,
    crystal: &Crystal,
    pseudopotentials: &[&PseudopotentialData],
    omega: f64,
    n_grid: usize,
) -> f64 {
    // Band energy: E_band = Σ_{n,k} f_{n,k} w_k ε_{n,k}
    let e_band: f64 = eigenvalues
        .iter()
        .zip(occupations.iter())
        .zip(kpoints.iter())
        .map(|((evs, occs), kp)| {
            evs.iter()
                .zip(occs.iter())
                .map(|(&e, &f)| f * kp.weight * e)
                .sum::<f64>()
        })
        .sum();

    // Hartree energy (double-counting correction)
    let e_hartree = hartree::hartree_energy(rho_g_basis, basis.g_vectors(), omega);

    // XC energy and potential integral
    let (exc_r, vxc_r) = xc::lda_xc_grid(rho_r);
    let e_xc = xc::lda_xc_energy(rho_r, &exc_r, omega);
    let dvol = omega / n_grid as f64;
    let e_vxc: f64 = rho_r
        .iter()
        .zip(vxc_r.iter())
        .map(|(&rho, &vxc)| rho * vxc * dvol)
        .sum();

    // Ewald energy (ion-ion)
    let e_ewald = crate::ewald::ewald_energy(crystal, pseudopotentials);

    // Total: E = E_band - E_H + E_xc - E_vxc + E_ewald
    let e_total = e_band - e_hartree + e_xc - e_vxc + e_ewald;

    info!("Energy components: E_band={e_band:.6}, E_H={e_hartree:.6}, E_xc={e_xc:.6}, E_vxc={e_vxc:.6}, E_ewald={e_ewald:.6}");
    info!("Total energy: {e_total:.6} eV");

    e_total
}
