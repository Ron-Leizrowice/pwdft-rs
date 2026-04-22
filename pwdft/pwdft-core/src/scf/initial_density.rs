//! Initial electron density for SCF convergence.
//!
//! Provides Superposition of Atomic Densities (SAD) using either:
//! - Atomic charge density from the pseudopotential file (PP_RHOATOM)
//! - Gaussian model charges when PP_RHOATOM is unavailable
//!
//! The initial density is constructed in G-space on the FFT grid:
//! ρ(G) = Σ_i S_i(G) ρ_atom,i(|G|) / Ω
//! where S_i(G) = exp(-iG·τ_i) is the structure factor.

use std::collections::HashMap;

use elements_rs::Element;
use num_complex::Complex64;
use rayon::prelude::*;

use super::FftGrid;
use crate::{basis::BasisSet, crystal::Crystal, pseudopotential::UpfPseudoPotential};

/// Default Gaussian width parameter for the SAD model charge (Å).
pub(crate) const DEFAULT_GAUSSIAN_SIGMA: f64 = 1.0;

/// Configuration for initial density generation.
pub struct InitialDensityConfig {
    /// Initial magnetic moment fraction per atom, in [-1, 1].
    pub magnetic_moments: Vec<f64>,
    /// Gaussian width override (Å). None uses the default.
    pub gaussian_sigma: Option<f64>,
}

impl InitialDensityConfig {
    pub fn non_magnetic(n_atoms: usize) -> Self {
        Self {
            magnetic_moments: vec![0.0; n_atoms],
            gaussian_sigma: None,
        }
    }
}

/// Statistics for diagnosing the SAD generation pipeline.
#[derive(Debug, Clone, Copy, Default)]
pub struct SadDiagnosticStats {
    pub integrated_pre_clamp: f64,
    pub negative_mass_clamped: f64,
    pub renorm_scale: f64,
}

/// Computes the SAD density on the FFT grid.
///
/// Returns the real-space density ρ(r) and diagnostic statistics.
pub(super) fn generate_initial_density(
    crystal: &Crystal,
    grid: &mut FftGrid,
    pseudopotentials: &HashMap<Element, UpfPseudoPotential>,
    n_electrons: f64,
    config: &InitialDensityConfig,
) -> Vec<f64> {
    let (_, _, rho_final, _) = assemble_sad_pipeline(crystal, grid, pseudopotentials, n_electrons, config);
    rho_final
}

/// Master pipeline for SAD assembly, shared by production and diagnostic entry
/// points.
fn assemble_sad_pipeline(
    crystal: &Crystal,
    grid: &mut FftGrid,
    pseudopotentials: &HashMap<Element, UpfPseudoPotential>,
    n_electrons: f64,
    config: &InitialDensityConfig,
) -> ([usize; 3], Vec<f64>, Vec<f64>, SadDiagnosticStats) {
    let omega = crystal.lattice.volume();
    let n_grid = grid.total_size();
    let sigma = config.gaussian_sigma.unwrap_or(DEFAULT_GAUSSIAN_SIGMA);

    // 1. Pre-calculate atomic positions and associated PPs to avoid lookups in the hot loop
    let atom_data: Vec<_> = crystal
        .atoms
        .iter()
        .map(|atom| {
            let pp = pseudopotentials
                .get(&atom.symbol)
                .expect("BUG: Missing PP during initial density generation");
            (atom.cart_position(&crystal.lattice), pp)
        })
        .collect();

    // 2. Assemble ρ(G) in parallel
    // ρ(G) = (1/Ω) Σ_i exp(-iG·τ_i) * ρ_at,i(|G|)
    let mut rho_g: Vec<Complex64> = (0..n_grid)
        .into_par_iter()
        .map(|idx| {
            let g = grid.g_vector_at(idx);
            let g_norm = g.norm();

            atom_data.iter().fold(Complex64::default(), |acc, (tau, pp)| {
                let rho_at_g = if pp.has_rho_atom() {
                    compute_bessel_rho_at(pp, g_norm)
                } else {
                    // Gaussian: ρ(G) = Z * exp(-|G|²σ²/2)
                    pp.z_valence * (-0.5 * g_norm * g_norm * sigma * sigma).exp()
                };

                let phase = -g.dot(tau);
                acc + Complex64::cis(phase) * (rho_at_g / omega)
            })
        })
        .collect();

    // 3. Inverse FFT to real space: ρ(r) = Σ_G ρ(G) exp(iG·r)
    grid.fft.inverse(&mut rho_g);
    let rho_pre_clamp: Vec<f64> = rho_g.into_iter().map(|c| c.re).collect();

    // 4. Post-processing: Clamp and Renormalize
    let dvol = omega / n_grid as f64;
    let integrated_pre_clamp = rho_pre_clamp.iter().sum::<f64>() * dvol;
    let negative_mass_clamped = rho_pre_clamp.iter().filter(|&&v| v < 0.0).sum::<f64>() * dvol;

    let mut rho_final = rho_pre_clamp.clone();
    for v in &mut rho_final {
        if *v < 0.0 {
            *v = 0.0;
        }
    }

    let post_clamp_integral = rho_final.iter().sum::<f64>() * dvol;
    let renorm_scale = if post_clamp_integral > 1e-15 {
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

    (grid.dims, rho_pre_clamp, rho_final, stats)
}

/// Compute the radial Fourier transform of the atomic density:
/// ρ_at(G) = ∫ [4πr²ρ(r)] j₀(|G|r) dr
///
/// Note: UPF's PP_RHOATOM already stores the 4πr²ρ(r) factor.
fn compute_bessel_rho_at(pp: &UpfPseudoPotential, g_norm: f64) -> f64 {
    let integrand: Vec<f64> = pp
        .rho_atom
        .iter()
        .zip(&pp.r_grid)
        .map(|(&rho_at, &r)| {
            let gr = g_norm * r;
            let j0 = if gr < 1e-10 {
                1.0 - gr * gr / 6.0 // Sinc expansion
            } else {
                gr.sin() / gr
            };
            rho_at * j0
        })
        .collect();

    crate::numerics::simpson_integrate(&integrand, &pp.rab)
}

// --- Diagnostic Wrappers ---

pub fn build_sad_density_for_diagnostic(
    crystal: &Crystal,
    pseudopotentials: &HashMap<Element, UpfPseudoPotential>,
    n_electrons: f64,
    ecutwfc_ev: f64,
    ecutrho_ratio: u32,
    explicit_dims: Option<[usize; 3]>,
    config: &InitialDensityConfig,
) -> ([usize; 3], Vec<f64>) {
    let basis = BasisSet::new(&crystal.lattice, ecutwfc_ev);
    let mut grid = FftGrid::new(&basis, &crystal.lattice, ecutrho_ratio, explicit_dims);
    let (_, _, rho_final, _) = assemble_sad_pipeline(crystal, &mut grid, pseudopotentials, n_electrons, config);
    (grid.dims, rho_final)
}

pub fn build_sad_density_for_diagnostic_verbose(
    crystal: &Crystal,
    pseudopotentials: &HashMap<Element, UpfPseudoPotential>,
    n_electrons: f64,
    ecutwfc_ev: f64,
    ecutrho_ratio: u32,
    explicit_dims: Option<[usize; 3]>,
    config: &InitialDensityConfig,
) -> ([usize; 3], Vec<f64>, Vec<f64>, SadDiagnosticStats) {
    let basis = BasisSet::new(&crystal.lattice, ecutwfc_ev);
    let mut grid = FftGrid::new(&basis, &crystal.lattice, ecutrho_ratio, explicit_dims);
    assemble_sad_pipeline(crystal, &mut grid, pseudopotentials, n_electrons, config)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    use super::*;
    use crate::crystal::{Atom, Crystal, Lattice};

    fn si_test_env() -> (Crystal, HashMap<Element, UpfPseudoPotential>) {
        let a = 5.431;
        let crystal = Crystal {
            lattice: Lattice::new(
                a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
            ),
            atoms: vec![
                Atom::new(Element::Si, [0.0, 0.0, 0.0]),
                Atom::new(Element::Si, [0.25, 0.25, 0.25]),
            ],
        };
        let mut pps = HashMap::new();
        pps.insert(Element::Si, UpfPseudoPotential::load("Si").unwrap());
        (crystal, pps)
    }

    #[test]
    fn test_sad_normalization() {
        let (crystal, pps) = si_test_env();
        let basis = BasisSet::new(&crystal.lattice, 100.0);
        let mut grid = FftGrid::new(&basis, &crystal.lattice, 4, None);
        let n_electrons = 8.0;

        let rho = generate_initial_density(
            &crystal,
            &mut grid,
            &pps,
            n_electrons,
            &InitialDensityConfig::non_magnetic(2),
        );

        let dvol = crystal.lattice.volume() / rho.len() as f64;
        let integral: f64 = rho.iter().sum::<f64>() * dvol;

        assert_relative_eq!(integral, n_electrons, epsilon = 1e-7);
    }

    #[test]
    fn test_sad_positivity() {
        let (crystal, pps) = si_test_env();
        let basis = BasisSet::new(&crystal.lattice, 100.0);
        let mut grid = FftGrid::new(&basis, &crystal.lattice, 4, None);

        let rho = generate_initial_density(&crystal, &mut grid, &pps, 8.0, &InitialDensityConfig::non_magnetic(2));

        assert!(
            rho.iter().all(|&v| v >= 0.0),
            "Density contained negative values post-clamp"
        );
    }

    #[test]
    fn test_sad_symmetry() {
        let (crystal, pps) = si_test_env();
        let basis = BasisSet::new(&crystal.lattice, 100.0);
        let mut grid = FftGrid::new(&basis, &crystal.lattice, 4, Some([16, 16, 16]));

        let rho = generate_initial_density(&crystal, &mut grid, &pps, 8.0, &InitialDensityConfig::non_magnetic(2));

        let [nx, ny, nz] = grid.dims;
        let get_rho = |ix, iy, iz| rho[(ix % nx) * ny * nz + (iy % ny) * nz + (iz % nz)];

        // Check cubic symmetry points (1,0,0) vs (0,1,0)
        assert_relative_eq!(get_rho(1, 0, 0), get_rho(0, 1, 0), epsilon = 1e-10);
        assert_relative_eq!(get_rho(1, 0, 0), get_rho(0, 0, 1), epsilon = 1e-10);
    }
}
