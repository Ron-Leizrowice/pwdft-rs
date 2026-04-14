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
    potential::{hartree, nonlocal::NonlocalPotential, xc},
    pseudopotential::PseudopotentialData,
};

/// Parameters for an SCF calculation.
pub struct ScfParams {
    pub n_bands: usize,
    pub max_iter: usize,
    pub conv_threshold: f64,
    pub mixing_beta: f64,
    pub mixing_ndim: usize,
    pub smearing_sigma: f64,
    /// Charge density cutoff as multiple of wavefunction cutoff.
    /// Controls FFT grid density. QE default is 4 for NC PPs.
    pub ecutrho_ratio: u32,
    /// Explicit FFT grid dimensions. If set, overrides ecutrho_ratio.
    pub fft_grid: Option<[usize; 3]>,
}

impl Default for ScfParams {
    fn default() -> Self {
        Self {
            n_bands: 8,
            max_iter: 100,
            conv_threshold: 1e-6,
            mixing_beta: 0.3,
            mixing_ndim: 8,
            smearing_sigma: 0.01,
            ecutrho_ratio: 4,
            fft_grid: None,
        }
    }
}

/// Result of an SCF calculation.
pub struct ScfResult {
    pub total_energy: f64,
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

    /// Map Miller indices (n1, n2, n3) to a flat FFT grid index.
    fn miller_to_idx(&self, n1: i32, n2: i32, n3: i32) -> usize {
        let i1 = ((n1 % self.dims[0] as i32) + self.dims[0] as i32) as usize % self.dims[0];
        let i2 = ((n2 % self.dims[1] as i32) + self.dims[1] as i32) as usize % self.dims[1];
        let i3 = ((n3 % self.dims[2] as i32) + self.dims[2] as i32) as usize % self.dims[2];
        i1 * self.dims[1] * self.dims[2] + i2 * self.dims[2] + i3
    }

    /// Compute the Cartesian G-vector for FFT grid index.
    fn g_vector_at(&self, idx: usize) -> Vector3<f64> {
        let [nx, ny, nz] = self.dims;
        let i1 = idx / (ny * nz);
        let i2 = (idx / nz) % ny;
        let i3 = idx % nz;
        // Map FFT index back to Miller index (centered)
        let n1 = if i1 > nx / 2 { i1 as i32 - nx as i32 } else { i1 as i32 };
        let n2 = if i2 > ny / 2 { i2 as i32 - ny as i32 } else { i2 as i32 };
        let n3 = if i3 > nz / 2 { i3 as i32 - nz as i32 } else { i3 as i32 };
        n1 as f64 * self.recip.a + n2 as f64 * self.recip.b + n3 as f64 * self.recip.c
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
        .map(|a| {
            pseudopotentials
                .iter()
                .find(|pp| {
                    crate::atoms::Element::from_symbol(&pp.element)
                        .map_or(false, |e| e.atomic_number() == a.z)
                })
                .unwrap()
                .z_valence
        })
        .sum();

    info!("SCF: {n_electrons} electrons, {omega:.3} ų cell volume");

    let grid = FftGrid::new(basis, &crystal.lattice, params.ecutrho_ratio, params.fft_grid);
    let n_grid = grid.total_size();
    let [nx, ny, nz] = grid.dims;
    info!("FFT grid: {nx}×{ny}×{nz} = {n_grid} points (ecutrho_ratio={})", params.ecutrho_ratio);

    let g_to_fft = grid.basis_to_fft(basis);

    // Precompute local pseudopotential on the FULL FFT grid
    let v_local_fft = compute_v_local_on_fft_grid(crystal, &grid, pseudopotentials, omega);

    // Initial density: uniform
    let rho_init = n_electrons / omega;
    let mut rho_r = vec![rho_init; n_grid];
    let mut rho_g = vec![Complex64::new(0.0, 0.0); n_grid];
    density_r_to_g(&grid.fft, &rho_r, &mut rho_g);

    let mut mixer = mixing::AndersonMixer::new(params.mixing_beta, params.mixing_ndim, n_grid);
    let mut eigenvalues_all = Vec::new();
    let mut fermi_energy;

    for iter in 0..params.max_iter {
        // 1. Hartree potential on FULL FFT grid
        let v_h_fft = hartree_on_fft_grid(&rho_g, &grid);

        // 2. XC potential in real space → FFT to G-space
        let (_exc_r, vxc_r) = xc::lda_xc_grid(&rho_r);
        let mut vxc_g = vec![Complex64::new(0.0, 0.0); n_grid];
        for (i, &v) in vxc_r.iter().enumerate() {
            vxc_g[i] = Complex64::new(v, 0.0);
        }
        grid.fft.forward(&mut vxc_g);
        let norm = 1.0 / n_grid as f64;
        for v in &mut vxc_g {
            *v *= norm;
        }

        // 3. V_eff on the FULL FFT grid: V_local + V_H + V_xc
        let mut v_eff_fft = vec![Complex64::new(0.0, 0.0); n_grid];
        for i in 0..n_grid {
            v_eff_fft[i] = v_local_fft[i] + v_h_fft[i] + vxc_g[i];
        }

        // 4. Solve eigenvalue problem at each k-point
        eigenvalues_all.clear();
        let mut all_kpoint_wavefns = Vec::new();

        for kp in kpoints {
            let mut h = build_hamiltonian_with_v_eff(basis, &kp.k, &v_eff_fft, &grid);

            // Add non-local pseudopotential
            let vnl = NonlocalPotential::new(crystal, basis, &kp.k, pseudopotentials);
            vnl.add_to_hamiltonian(&mut h, crystal, basis, &kp.k);

            let result = dense::diagonalize_lowest(&h, params.n_bands);
            eigenvalues_all.push(result.eigenvalues);
            all_kpoint_wavefns.push(result.eigenvectors);
        }

        // 5. Occupations
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

        // 6. New density
        let mut rho_r_new = density::compute_density(
            basis, kpoints, &all_kpoint_wavefns, &occupations, &g_to_fft, &grid.fft,
            n_electrons, omega,
        );

        // 6b. Symmetrize density if symmetry info is available
        if let Some(symm) = symmetry {
            crate::symmetry::density::symmetrize_density(&mut rho_r_new, grid.dims, symm);
        }

        // 7. Convergence check
        let delta = density_diff(&rho_r, &rho_r_new, omega, n_grid);
        info!("SCF iter {}: E_fermi = {:.6} eV, delta_rho = {:.2e}", iter + 1, fermi_energy, delta);

        if delta < params.conv_threshold {
            info!("SCF converged after {} iterations", iter + 1);
            rho_r = rho_r_new;
            density_r_to_g(&grid.fft, &rho_r, &mut rho_g);
            let rho_g_basis: Vec<Complex64> = g_to_fft.iter().map(|&idx| rho_g[idx]).collect();

            // Diagnostics: print V_eff at key G-vectors for comparison with QE
            for &[n1, n2, n3] in &[[0,0,0], [1,0,0], [1,1,0], [1,1,1], [2,0,0]] {
                let fft_idx = grid.miller_to_idx(n1, n2, n3);
                info!(
                    "V_eff(G=({},{},{})) = {:+.6} {:+.6}i eV  (loc={:+.6} H={:+.6} xc={:+.6})",
                    n1, n2, n3,
                    v_eff_fft[fft_idx].re, v_eff_fft[fft_idx].im,
                    v_local_fft[fft_idx].re, v_h_fft[fft_idx].re, vxc_g[fft_idx].re
                );
            }

            let total_energy = compute_total_energy(
                &eigenvalues_all, &occupations, kpoints, &rho_r, &rho_g,
                &grid, crystal, pseudopotentials, omega,
            );

            return Ok(ScfResult {
                total_energy,
                eigenvalues: eigenvalues_all,
                fermi_energy,
                n_iterations: iter + 1,
                rho_g: rho_g_basis,
            });
        }

        // 8. Mix
        rho_r = mixer.mix(&rho_r, &rho_r_new);
        density_r_to_g(&grid.fft, &rho_r, &mut rho_g);
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
    let mut v_local = vec![Complex64::new(0.0, 0.0); n_grid];

    for idx in 0..n_grid {
        let g = grid.g_vector_at(idx);
        let g_norm = g.norm();

        for atom in &crystal.atoms {
            let pp = pseudopotentials
                .iter()
                .find(|pp| {
                    crate::atoms::Element::from_symbol(&pp.element)
                        .map_or(false, |e| e.atomic_number() == atom.z)
                })
                .unwrap();

            let tau = atom.cart_position(&crystal.lattice);
            let phase = -g.dot(&tau);
            let sf = Complex64::new(phase.cos(), phase.sin());
            let v_form = pp.v_local_of_g(g_norm, omega);

            v_local[idx] += sf * v_form;
        }
    }

    v_local
}

/// Compute Hartree potential on the FULL FFT grid.
fn hartree_on_fft_grid(rho_g: &[Complex64], grid: &FftGrid) -> Vec<Complex64> {
    let fourpi_e2 = 4.0 * std::f64::consts::PI * hartree::E2;
    let n_grid = grid.total_size();
    let mut v_h = vec![Complex64::new(0.0, 0.0); n_grid];

    for idx in 0..n_grid {
        let g = grid.g_vector_at(idx);
        let g2 = g.norm_squared();
        if g2 > 1e-20 {
            v_h[idx] = rho_g[idx] * fourpi_e2 / g2;
        }
    }

    v_h
}

/// Build Hamiltonian: kinetic + V_eff(G-G') looked up from FFT grid.
fn build_hamiltonian_with_v_eff(
    basis: &BasisSet,
    k: &Vector3<f64>,
    v_eff_fft: &[Complex64],
    grid: &FftGrid,
) -> nalgebra::DMatrix<Complex64> {
    let n = basis.len();
    let mut h = nalgebra::DMatrix::zeros(n, n);

    // Kinetic (diagonal)
    for (i, g) in basis.g_vectors().iter().enumerate() {
        let ke = HBAR2_OVER_2M * (k + g).norm_squared();
        h[(i, i)] = Complex64::new(ke, 0.0);
    }

    // Potential: V_{G,G'} = V_eff(G-G') from FFT grid
    let miller = basis.miller_indices();
    for i in 0..n {
        for j in 0..n {
            let dn1 = miller[i][0] - miller[j][0];
            let dn2 = miller[i][1] - miller[j][1];
            let dn3 = miller[i][2] - miller[j][2];
            let fft_idx = grid.miller_to_idx(dn1, dn2, dn3);
            h[(i, j)] += v_eff_fft[fft_idx];
        }
    }

    h
}

/// FFT density from real space to G-space (normalized).
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

fn density_diff(rho_old: &[f64], rho_new: &[f64], omega: f64, n_grid: usize) -> f64 {
    let dvol = omega / n_grid as f64;
    let sum_sq: f64 = rho_old
        .iter()
        .zip(rho_new.iter())
        .map(|(&a, &b)| (a - b).powi(2) * dvol)
        .sum();
    (sum_sq / omega).sqrt()
}

/// Compute total energy with proper double-counting corrections.
fn compute_total_energy(
    eigenvalues: &[Vec<f64>],
    occupations: &[Vec<f64>],
    kpoints: &[KPoint],
    rho_r: &[f64],
    rho_g: &[Complex64],
    grid: &FftGrid,
    crystal: &Crystal,
    pseudopotentials: &[&PseudopotentialData],
    omega: f64,
) -> f64 {
    let n_grid = grid.total_size();

    // Band energy
    let e_band: f64 = eigenvalues
        .iter()
        .zip(occupations.iter())
        .zip(kpoints.iter())
        .map(|((evs, occs), kp)| {
            evs.iter().zip(occs.iter()).map(|(&e, &f)| f * kp.weight * e).sum::<f64>()
        })
        .sum();

    // Hartree energy on full FFT grid
    let fourpi_e2 = 4.0 * std::f64::consts::PI * hartree::E2;
    let e_hartree: f64 = (0..n_grid)
        .map(|idx| {
            let g2 = grid.g_vector_at(idx).norm_squared();
            if g2 > 1e-20 {
                rho_g[idx].norm_sqr() * fourpi_e2 / g2
            } else {
                0.0
            }
        })
        .sum::<f64>()
        * 0.5
        * omega;

    // XC energy and potential integral
    let (exc_r, vxc_r) = xc::lda_xc_grid(rho_r);
    let e_xc = xc::lda_xc_energy(rho_r, &exc_r, omega);
    let dvol = omega / n_grid as f64;
    let e_vxc: f64 = rho_r
        .iter()
        .zip(vxc_r.iter())
        .map(|(&rho, &vxc)| rho * vxc * dvol)
        .sum();

    let e_ewald = crate::ewald::ewald_energy(crystal, pseudopotentials);

    let e_total = e_band - e_hartree + e_xc - e_vxc + e_ewald;

    info!("Energy: band={e_band:.6} H={e_hartree:.6} xc={e_xc:.6} vxc={e_vxc:.6} ewald={e_ewald:.6}");
    info!("Total energy: {e_total:.6} eV");

    e_total
}
