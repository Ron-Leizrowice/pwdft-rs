//! Initial electron density for SCF convergence.
//!
//! Provides Superposition of Atomic Densities (SAD) using either:
//! - Atomic charge density from the pseudopotential file (PP_RHOATOM)
//! - Gaussian model charges when PP_RHOATOM is unavailable
//!
//! The initial density is constructed in G-space on the FFT grid,
//! then inverse-FFT'd to real space.
//!
//! For spin-polarized calculations, the initial spin density is:
//!   ρ_up(r)   = (1 + m_i)/2 × ρ_atom_i(r)
//!   ρ_down(r) = (1 - m_i)/2 × ρ_atom_i(r)
//! where m_i ∈ [-1, 1] is the initial magnetic moment fraction of atom i.

use num_complex::Complex64;

use crate::{
    crystal::Crystal,
    pseudopotential::PseudopotentialData,
};

use super::FftGrid;

/// Default Gaussian width parameter for the SAD model charge (Å).
///
/// Controls the spatial extent of the initial density guess when a
/// pseudopotential file does not ship with `PP_RHOATOM`. A value of
/// ~1.0 Å gives a physically reasonable distribution for most elements;
/// the SCF refines this away in the first few iterations, so the exact
/// number has weak effect on the converged density. For numerical
/// sensitivity studies (a wider sigma smooths the initial high-frequency
/// G content; a narrower sigma makes the initial Delta-rho larger) the
/// value is exposed as the YAML input `initial_density.gaussian_sigma`.
///
/// The `pub(crate)` visibility is deliberate: this is the single source
/// of truth for the default in both `Settings` (YAML deserialization)
/// and `ScfParams` (programmatic construction). Changing the value here
/// changes both entry points at once.
pub(crate) const DEFAULT_GAUSSIAN_SIGMA: f64 = 1.0;

/// Configuration for initial density generation.
pub struct InitialDensityConfig {
    /// Initial magnetic moment fraction per atom, in [-1, 1].
    /// Length must equal crystal.atoms.len(). Default: all zeros (non-magnetic).
    pub magnetic_moments: Vec<f64>,
    /// Gaussian width override (Å). None uses the default.
    pub gaussian_sigma: Option<f64>,
}

impl InitialDensityConfig {
    /// Non-magnetic default: all moments zero.
    pub fn non_magnetic(n_atoms: usize) -> Self {
        Self {
            magnetic_moments: vec![0.0; n_atoms],
            gaussian_sigma: None,
        }
    }
}

/// Public diagnostic wrapper around the private SCF-driver SAD
/// assembler for integration tests that need the output on an
/// explicitly-sized FFT grid (e.g. to compare against a Python
/// reference computed on the same grid).
///
/// Returns `(dims, rho_r)`: the chosen FFT grid dimensions (echoed back
/// so the caller can verify) and the flat row-major real-space density
/// in e/Å³ on that grid, post-clamp and post-renormalization (i.e. the
/// same object the SCF driver would hand off to its hot loop).
/// Integration tests can then shell-average `rho_r` around each atom
/// and compare bin-by-bin against an independent reference.
///
/// The caller supplies the crystal, pseudopotentials, electron count,
/// an `ecutwfc` (eV) used only to seed a dummy [`crate::basis::BasisSet`]
/// (FFT grid dims come from `explicit_dims` when supplied), and an
/// [`InitialDensityConfig`]. Pass `explicit_dims = Some([nx, ny, nz])` to
/// force a specific grid; the wrapper then constructs the same internal
/// grid state the SCF driver would use for that size.
pub fn build_sad_density_for_diagnostic(
    crystal: &Crystal,
    pseudopotentials: &[&PseudopotentialData],
    n_electrons: f64,
    ecutwfc_ev: f64,
    ecutrho_ratio: u32,
    explicit_dims: Option<[usize; 3]>,
    config: &InitialDensityConfig,
) -> ([usize; 3], Vec<f64>) {
    let basis = crate::basis::BasisSet::new(&crystal.lattice, ecutwfc_ev);
    let mut grid = FftGrid::new(&basis, &crystal.lattice, ecutrho_ratio, explicit_dims);
    let dims = grid.dims;
    let rho_r = generate_initial_density(crystal, &mut grid, pseudopotentials, n_electrons, config);
    (dims, rho_r)
}

/// Extra statistics from the SAD assembly — returned by
/// [`build_sad_density_for_diagnostic_verbose`] so integration tests
/// can separately diagnose the Bessel + IFFT, the negative-density
/// clamp, and the renormalization-to-N_el steps.
#[derive(Debug, Clone, Copy)]
pub struct SadDiagnosticStats {
    /// ∫ρ d³r immediately after the inverse FFT, before any clamp or
    /// renormalization is applied. Ideally `= n_electrons` to double
    /// precision if the Bessel transform and FFT conventions are
    /// consistent.
    pub integrated_pre_clamp: f64,
    /// Total *negative* mass: Σ_{j: ρ_j < 0} ρ_j · dV  (≤ 0 by construction).
    /// Large magnitude here means the clamp step removes a physically
    /// meaningful amount of charge and the subsequent renormalization
    /// deforms the density non-uniformly relative to a no-clamp pipeline.
    /// The initial-density clamp is a production-only safety net against
    /// Gibbs ringing at the atom cores; a no-op on all systems that have
    /// strictly non-negative radial atomic densities on the FFT grid.
    pub negative_mass_clamped: f64,
    /// Renormalization scale factor applied after the clamp:
    /// `rho *= n_electrons / integral_post_clamp`. Value `= 1.0` means
    /// clamp removed zero net charge and ρ was already normalized.
    pub renorm_scale: f64,
}

/// Verbose variant of [`build_sad_density_for_diagnostic`] that returns
/// diagnostic statistics alongside the clamped-and-renormalized ρ(r),
/// plus the intermediate ρ(r) BEFORE the clamp+renorm so integration
/// tests can attribute mismatches to the upstream (Bessel+FFT) vs.
/// downstream (clamp+renorm) parts of the pipeline.
///
/// Returns `(dims, rho_pre_clamp, rho_final, stats)`.
///
/// # Panics
///
/// Panics if any atom in `crystal` does not have a matching
/// pseudopotential in `pseudopotentials` — the caller is expected to
/// pre-validate, mirroring the production path through the SCF driver.
pub fn build_sad_density_for_diagnostic_verbose(
    crystal: &Crystal,
    pseudopotentials: &[&PseudopotentialData],
    n_electrons: f64,
    ecutwfc_ev: f64,
    ecutrho_ratio: u32,
    explicit_dims: Option<[usize; 3]>,
    config: &InitialDensityConfig,
) -> ([usize; 3], Vec<f64>, Vec<f64>, SadDiagnosticStats) {
    let basis = crate::basis::BasisSet::new(&crystal.lattice, ecutwfc_ev);
    let mut grid = FftGrid::new(&basis, &crystal.lattice, ecutrho_ratio, explicit_dims);
    let dims = grid.dims;

    // Manually walk the SAD path so we can snapshot the intermediate.
    let omega = crystal.lattice.volume();
    let n_grid = grid.total_size();
    let sigma = config.gaussian_sigma.unwrap_or(DEFAULT_GAUSSIAN_SIGMA);

    let mut rho_g = vec![Complex64::new(0.0, 0.0); n_grid];
    for atom in &crystal.atoms {
        // Mirrors the sibling `generate_initial_density` precondition:
        // callers must pre-validate PPs. The diagnostic entry point is
        // documented as such in the `# Panics` section above.
        #[expect(
            clippy::expect_used,
            reason = "BUG: atom has no matching pseudopotential — pre-validation is the caller's contract"
        )]
        let pp = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials)
            .expect("BUG: atom has no matching pseudopotential (should have been validated at startup)");
        let tau = atom.cart_position(&crystal.lattice);
        let z_val = pp.z_valence;
        if pp.has_rho_atom() {
            add_atomic_density_from_pp(pp, &tau, &grid, omega, &mut rho_g);
        } else {
            add_gaussian_density(&tau, z_val, sigma, &grid, omega, &mut rho_g);
        }
    }
    grid.fft.inverse(&mut rho_g);
    let rho_pre_clamp: Vec<f64> = rho_g.iter().map(|c| c.re).collect();

    let dvol = omega / n_grid as f64;
    let integrated_pre_clamp: f64 = rho_pre_clamp.iter().sum::<f64>() * dvol;
    let negative_mass_clamped: f64 = rho_pre_clamp
        .iter()
        .filter(|v| **v < 0.0)
        .copied()
        .sum::<f64>()
        * dvol;

    // Apply the production clamp + renorm on a fresh copy.
    let mut rho_final: Vec<f64> = rho_pre_clamp.clone();
    for v in &mut rho_final {
        if *v < 0.0 {
            *v = 0.0;
        }
    }
    let post_clamp_integral: f64 = rho_final.iter().sum::<f64>() * dvol;
    let renorm_scale = if post_clamp_integral.abs() > 1e-15 {
        n_electrons / post_clamp_integral
    } else {
        1.0
    };
    for v in &mut rho_final {
        *v *= renorm_scale;
    }

    let stats = SadDiagnosticStats {
        integrated_pre_clamp,
        negative_mass_clamped,
        renorm_scale,
    };
    (dims, rho_pre_clamp, rho_final, stats)
}

/// Generate initial charge density on the FFT real-space grid via
/// Superposition of Atomic Densities (SAD).
///
/// Uses PP_RHOATOM if available; falls back to Gaussian model charges.
/// Returns ρ(r) on the FFT grid, normalized to integrate to n_electrons.
pub(super) fn generate_initial_density(
    crystal: &Crystal,
    grid: &mut FftGrid,
    pseudopotentials: &[&PseudopotentialData],
    n_electrons: f64,
    config: &InitialDensityConfig,
) -> Vec<f64> {
    let omega = crystal.lattice.volume();
    let n_grid = grid.total_size();
    let sigma = config.gaussian_sigma.unwrap_or(DEFAULT_GAUSSIAN_SIGMA);

    // Build ρ_init(G) on the FFT grid in G-space
    let mut rho_g = vec![Complex64::new(0.0, 0.0); n_grid];

    for atom in &crystal.atoms {
        // SAFETY: ScfContext::new validates all atoms have matching PPs before
        // this function is called. A missing PP here would be a programming error.
        let pp = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials)
            .expect("BUG: atom has no matching pseudopotential (should have been validated at startup)");

        let tau = atom.cart_position(&crystal.lattice);
        let z_val = pp.z_valence;

        if pp.has_rho_atom() {
            // Use atomic density from pseudopotential
            add_atomic_density_from_pp(pp, &tau, grid, omega, &mut rho_g);
        } else {
            // Gaussian model: ρ_atom(G) = (Z_val/Ω) × exp(-|G|²σ²/2) × S(G)
            add_gaussian_density(&tau, z_val, sigma, grid, omega, &mut rho_g);
        }
    }

    // Inverse FFT to get ρ(r) in real space
    // rho_g currently has Fourier convention: ρ(G) = (1/N) Σ_r ρ(r) e^{-iGr}
    // So ρ(r) = Σ_G ρ(G) e^{iGr} = N × IFFT(ρ(G))
    let mut rho_r_complex = rho_g;
    grid.fft.inverse(&mut rho_r_complex);

    // Extract real part (imaginary should be negligible)
    let mut rho_r: Vec<f64> = rho_r_complex.iter().map(|c| c.re).collect();

    // Ensure non-negative density (Gaussian model can produce tiny negative values
    // from FFT aliasing)
    for v in &mut rho_r {
        if *v < 0.0 {
            *v = 0.0;
        }
    }

    // Normalize to N_electrons
    let dvol = omega / n_grid as f64;
    let integral: f64 = rho_r.iter().sum::<f64>() * dvol;
    if integral.abs() > 1e-15 {
        let scale = n_electrons / integral;
        for v in &mut rho_r {
            *v *= scale;
        }
    }

    rho_r
}

/// Add Gaussian model charge for one atom to ρ(G).
///
/// ρ_atom(G) = (Z_val/Ω) × exp(-|G|²σ²/2) × exp(-iG·τ)
fn add_gaussian_density(
    tau: &nalgebra::Vector3<f64>,
    z_val: f64,
    sigma: f64,
    grid: &FftGrid,
    omega: f64,
    rho_g: &mut [Complex64],
) {
    let n_grid = grid.total_size();
    let half_sigma2 = 0.5 * sigma * sigma;
    let prefactor = z_val / omega;

    for (idx, rho_g_val) in rho_g.iter_mut().enumerate().take(n_grid) {
        let g = grid.g_vector_at(idx);
        let g2 = g.norm_squared();
        let gauss = (-g2 * half_sigma2).exp();
        let phase = -g.dot(tau);
        let sf = Complex64::cis(phase);
        *rho_g_val += sf * (prefactor * gauss);
    }
}

/// Add atomic density from PP_RHOATOM for one atom to ρ(G).
///
/// Uses the radial ρ_atom(r) from the pseudopotential, Bessel-transformed
/// to G-space: ρ_atom(G) = (1/Ω) × 4π ∫ ρ_atom(r) j_0(|G|r) r² dr × S(G)
///
/// PP_RHOATOM stores 4πr²ρ(r), so the integral is:
/// ρ_atom(G) = (1/Ω) × ∫ [4πr²ρ(r)] × j_0(|G|r) × dr × S(G)
fn add_atomic_density_from_pp(
    pp: &PseudopotentialData,
    tau: &nalgebra::Vector3<f64>,
    grid: &FftGrid,
    omega: f64,
    rho_g: &mut [Complex64],
) {
    let n_grid = grid.total_size();
    let r_grid = &pp.r_grid;
    let rab = &pp.rab;
    let rho_at = &pp.rho_atom;

    for (idx, rho_g_val) in rho_g.iter_mut().enumerate().take(n_grid) {
        let g = grid.g_vector_at(idx);
        let g_norm = g.norm();

        // Bessel transform of radial atomic density (Simpson's rule)
        let n = r_grid.len();
        let mut integrand = vec![0.0; n];
        for i in 0..n {
            let r = r_grid[i];
            let gr = g_norm * r;
            let j0 = if gr < 1e-10 {
                1.0 - gr * gr / 6.0
            } else {
                gr.sin() / gr
            };
            integrand[i] = rho_at[i] * j0;
        }
        let integral = crate::numerics::simpson_integrate(&integrand, rab);

        let phase = -g.dot(tau);
        let sf = Complex64::cis(phase);
        *rho_g_val += sf * (integral / omega);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::{Atom, Crystal, Lattice};

    use approx::relative_eq;
    use nalgebra::Vector3;

    fn si_crystal() -> Crystal {
        let a = 5.431;
        Crystal {
            lattice: Lattice::new(
                a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
            ),
            atoms: vec![
                Atom::new(14, [0.0, 0.0, 0.0]),
                Atom::new(14, [0.25, 0.25, 0.25]),
            ],
        }
    }

    fn make_grid(crystal: &Crystal) -> FftGrid {
        let basis = crate::basis::BasisSet::new(&crystal.lattice, 100.0);
        FftGrid::new(&basis, &crystal.lattice, 4, None)
    }

    #[test]
    fn test_gaussian_density_integrates_to_n_electrons() {
        let crystal = si_crystal();
        let mut grid = make_grid(&crystal);
        let pp = crate::pseudopotential::load(
            &std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap();

        let config = InitialDensityConfig::non_magnetic(2);
        let rho = generate_initial_density(&crystal, &mut grid, &[&pp], 8.0, &config);

        let omega = crystal.lattice.volume();
        let dvol = omega / rho.len() as f64;
        let integral: f64 = rho.iter().sum::<f64>() * dvol;
        assert!(
            relative_eq!(integral, 8.0, epsilon = 0.01),
            "density integrates to {integral}, expected 8.0"
        );
    }

    #[test]
    fn test_gaussian_density_non_negative() {
        let crystal = si_crystal();
        let mut grid = make_grid(&crystal);
        let pp = crate::pseudopotential::load(
            &std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap();

        let config = InitialDensityConfig::non_magnetic(2);
        let rho = generate_initial_density(&crystal, &mut grid, &[&pp], 8.0, &config);

        for &v in &rho {
            assert!(v >= 0.0, "negative density: {v}");
        }
    }

    #[test]
    fn test_gaussian_density_peaked_at_atoms() {
        let crystal = si_crystal();
        let mut grid = make_grid(&crystal);
        let pp = crate::pseudopotential::load(
            &std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap();

        let config = InitialDensityConfig::non_magnetic(2);
        let rho = generate_initial_density(&crystal, &mut grid, &[&pp], 8.0, &config);

        // Density should be positive and have structure (not uniform)
        let rho_max = rho.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let rho_min = rho.iter().copied().fold(f64::INFINITY, f64::min);
        let rho_mean = rho.iter().sum::<f64>() / rho.len() as f64;
        assert!(rho_mean > 0.0, "mean density should be positive");
        assert!(
            rho_max > rho_min * 1.1,
            "density should have spatial variation: max={rho_max:.4e}, min={rho_min:.4e}"
        );
    }

    #[test]
    fn test_gaussian_density_has_diamond_symmetry() {
        let crystal = si_crystal();
        let mut grid = make_grid(&crystal);
        let pp = crate::pseudopotential::load(
            &std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap();

        let config = InitialDensityConfig::non_magnetic(2);
        let rho = generate_initial_density(&crystal, &mut grid, &[&pp], 8.0, &config);

        let [_nx, ny, nz] = grid.dims;
        // For FCC Si with equal grid dims, permuting (x,y,z) should give same density
        // (cubic symmetry). Check a few points.
        let rho_at = |ix: usize, iy: usize, iz: usize| rho[ix * ny * nz + iy * nz + iz];

        // (1,0,0) vs (0,1,0) vs (0,0,1)
        let v100 = rho_at(1, 0, 0);
        let v010 = rho_at(0, 1, 0);
        let v001 = rho_at(0, 0, 1);
        assert!(
            relative_eq!(v100, v010, epsilon = 1e-10),
            "ρ(1,0,0)={v100} ≠ ρ(0,1,0)={v010}"
        );
        assert!(
            relative_eq!(v100, v001, epsilon = 1e-10),
            "ρ(1,0,0)={v100} ≠ ρ(0,0,1)={v001}"
        );
    }
}
