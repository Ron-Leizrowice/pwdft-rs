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

/// Gaussian width parameter for model atomic charge (Å).
/// Controls the spatial extent of the initial guess.
/// A value of ~1.0 Å gives a physically reasonable charge distribution
/// for most elements. Doesn't need to be precise — the SCF will refine.
const DEFAULT_GAUSSIAN_SIGMA: f64 = 1.0;

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
    let omega = crystal.lattice.volume().abs();
    let n_grid = grid.total_size();
    let sigma = config.gaussian_sigma.unwrap_or(DEFAULT_GAUSSIAN_SIGMA);

    // Build ρ_init(G) on the FFT grid in G-space
    let mut rho_g = vec![Complex64::new(0.0, 0.0); n_grid];

    for atom in &crystal.atoms {
        let pp = crate::pseudopotential::find_for_atom(atom.z, pseudopotentials);

        let tau = atom.cart_position(&crystal.lattice);
        let z_val = pp.z_valence;

        if pp.has_rho_atom() {
            // Use atomic density from pseudopotential
            add_atomic_density_from_pp(pp, &tau, z_val, grid, omega, &mut rho_g);
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
    _z_val: f64,
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

        // Bessel transform of radial atomic density
        let mut integral = 0.0;
        for (r, (&dr, &rho_r_at)) in r_grid.iter().zip(rab.iter().zip(rho_at.iter())) {
            let gr = g_norm * r;
            let j0 = if gr < 1e-10 {
                1.0 - gr * gr / 6.0
            } else {
                gr.sin() / gr
            };

            integral += rho_r_at * j0 * dr;
        }

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
            &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Si.UPF"),
        )
        .unwrap();

        let config = InitialDensityConfig::non_magnetic(2);
        let rho = generate_initial_density(&crystal, &mut grid, &[&pp], 8.0, &config);

        let omega = crystal.lattice.volume().abs();
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
            &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Si.UPF"),
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
            &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Si.UPF"),
        )
        .unwrap();

        let config = InitialDensityConfig::non_magnetic(2);
        let rho = generate_initial_density(&crystal, &mut grid, &[&pp], 8.0, &config);

        // Origin (0,0,0) is atom position — density should be higher there than average
        let rho_origin = rho[0];
        let rho_mean = rho.iter().sum::<f64>() / rho.len() as f64;
        assert!(
            rho_origin > rho_mean * 1.3,
            "density at origin ({rho_origin}) should be higher than mean ({rho_mean})"
        );
    }

    #[test]
    fn test_gaussian_density_has_diamond_symmetry() {
        let crystal = si_crystal();
        let mut grid = make_grid(&crystal);
        let pp = crate::pseudopotential::load(
            &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Si.UPF"),
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
