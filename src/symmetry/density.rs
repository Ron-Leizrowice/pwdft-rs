//! Charge density symmetrization on the FFT grid.
//!
//! Applies all crystal symmetry operations to the density and averages:
//!   ρ_sym(r) = (1/N_ops) Σ_S ρ(S⁻¹ r)

use super::SymmetryInfo;

/// Symmetrize a real-space charge density on the FFT grid.
///
/// For each symmetry operation S = {R|τ}, the inverse S⁻¹ = {R⁻¹|-R⁻¹τ}
/// maps each grid point to another grid point (for compatible grids).
/// The symmetrized density is the average over all operations.
///
/// **Important**: The FFT grid dimensions must be compatible with the symmetry
/// operations. Use [`check_grid_compatibility`] before calling this.
pub fn symmetrize_density(rho: &mut [f64], dims: [usize; 3], symmetry: &SymmetryInfo) {
    let [nx, ny, nz] = dims;
    let n_grid = nx * ny * nz;
    assert_eq!(rho.len(), n_grid);

    if symmetry.n_ops <= 1 {
        return; // nothing to symmetrize
    }

    let n_ops = symmetry.n_ops as f64;
    let rho_orig = rho.to_vec();

    // Zero out and accumulate
    for v in rho.iter_mut() {
        *v = 0.0;
    }

    for op in &symmetry.operations {
        let s_inv = op.inverse();
        let r_inv = s_inv.rotation;
        let tau_inv = s_inv.translation;

        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    // Fractional coordinates of this grid point
                    let f = [
                        ix as f64 / nx as f64,
                        iy as f64 / ny as f64,
                        iz as f64 / nz as f64,
                    ];

                    // Apply S⁻¹: f' = R⁻¹·f + τ_inv
                    let fp = [
                        r_inv[0][0] as f64 * f[0]
                            + r_inv[0][1] as f64 * f[1]
                            + r_inv[0][2] as f64 * f[2]
                            + tau_inv[0],
                        r_inv[1][0] as f64 * f[0]
                            + r_inv[1][1] as f64 * f[1]
                            + r_inv[1][2] as f64 * f[2]
                            + tau_inv[1],
                        r_inv[2][0] as f64 * f[0]
                            + r_inv[2][1] as f64 * f[1]
                            + r_inv[2][2] as f64 * f[2]
                            + tau_inv[2],
                    ];

                    // Map to grid indices with periodic boundary conditions
                    let jx = ((fp[0] * nx as f64).round() as i64 % nx as i64 + nx as i64) as usize
                        % nx;
                    let jy = ((fp[1] * ny as f64).round() as i64 % ny as i64 + ny as i64) as usize
                        % ny;
                    let jz = ((fp[2] * nz as f64).round() as i64 % nz as i64 + nz as i64) as usize
                        % nz;

                    let src_idx = jx * ny * nz + jy * nz + jz;
                    let dst_idx = ix * ny * nz + iy * nz + iz;
                    rho[dst_idx] += rho_orig[src_idx];
                }
            }
        }
    }

    // Normalize by number of operations
    for v in rho.iter_mut() {
        *v /= n_ops;
    }
}

/// Check if the FFT grid dimensions are compatible with all symmetry operations.
///
/// For integer rotation R in fractional coords, the grid point (ix, iy, iz) maps
/// to another exact grid point if and only if R_{ij} × n_j ≡ 0 (mod n_i) for all i, j.
///
/// Returns true if all operations are compatible.
pub fn check_grid_compatibility(dims: [usize; 3], symmetry: &SymmetryInfo) -> bool {
    let ns = [dims[0] as i32, dims[1] as i32, dims[2] as i32];
    for op in &symmetry.operations {
        let r = &op.rotation;
        for i in 0..3 {
            for j in 0..3 {
                if (r[i][j] * ns[j]) % ns[i] != 0 {
                    return false;
                }
            }
        }
    }
    true
}

/// Find the smallest FFT-friendly grid dimensions compatible with the symmetry.
///
/// Starts from the given minimum dimensions and increases until compatibility
/// is achieved. Returns adjusted dimensions.
pub fn compatible_grid_dims(
    min_dims: [usize; 3],
    symmetry: &SymmetryInfo,
) -> [usize; 3] {
    // For cubic symmetry, making all dimensions equal is usually sufficient
    let max_dim = *min_dims.iter().max().unwrap();
    let mut dims = [max_dim; 3];

    // Try increasing until compatible
    for _ in 0..100 {
        let candidate = [
            crate::fft::fft_grid_size(dims[0] as i32 / 2),
            crate::fft::fft_grid_size(dims[1] as i32 / 2),
            crate::fft::fft_grid_size(dims[2] as i32 / 2),
        ];
        // Make all equal to the max for safety
        let m = *candidate.iter().max().unwrap();
        let candidate = [m, m, m];
        if check_grid_compatibility(candidate, symmetry) {
            return candidate;
        }
        dims = [dims[0] + 1, dims[1] + 1, dims[2] + 1];
    }

    // Fallback: just use the input dims (may not be compatible)
    min_dims
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::{Atom, Crystal, Lattice};
    use approx::relative_eq;
    use nalgebra::Vector3;

    fn si_fcc() -> Crystal {
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

    #[test]
    fn test_symmetrize_preserves_integral() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12]; // compatible with O_h
        let n = dims[0] * dims[1] * dims[2];

        // Random-ish density
        let mut rho: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin().abs() + 0.1).collect();
        let integral_before: f64 = rho.iter().sum();

        symmetrize_density(&mut rho, dims, &symmetry);

        let integral_after: f64 = rho.iter().sum();
        assert!(
            relative_eq!(integral_before, integral_after, epsilon = 1e-10),
            "integral changed: {integral_before} → {integral_after}"
        );
    }

    #[test]
    fn test_symmetrize_uniform_unchanged() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12];
        let n = dims[0] * dims[1] * dims[2];

        // Uniform density should be unchanged
        let mut rho = vec![1.0; n];
        symmetrize_density(&mut rho, dims, &symmetry);
        for &v in &rho {
            assert!(
                relative_eq!(v, 1.0, epsilon = 1e-14),
                "uniform density changed to {v}"
            );
        }
    }

    #[test]
    fn test_symmetrize_single_point_orbit() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12];
        let n = dims[0] * dims[1] * dims[2];

        // Set density at a single non-special point
        let mut rho = vec![0.0; n];
        let test_idx = 1 * dims[1] * dims[2] + 2 * dims[2] + 3; // (1,2,3)
        rho[test_idx] = 48.0; // put n_ops worth of charge so average = 1.0 per orbit point

        symmetrize_density(&mut rho, dims, &symmetry);

        // After symmetrization, all orbit points should have equal density
        let nonzero: Vec<f64> = rho.iter().filter(|&&v| v > 1e-10).cloned().collect();
        assert!(
            !nonzero.is_empty(),
            "symmetrized density has no nonzero points"
        );
        let ref_val = nonzero[0];
        for &v in &nonzero {
            assert!(
                relative_eq!(v, ref_val, epsilon = 1e-10),
                "orbit points not equal: {v} vs {ref_val}"
            );
        }
    }

    #[test]
    fn test_grid_compatibility_cubic() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);

        // Equal cubic grid should be compatible with O_h
        assert!(check_grid_compatibility([12, 12, 12], &symmetry));
        assert!(check_grid_compatibility([18, 18, 18], &symmetry));
        assert!(check_grid_compatibility([20, 20, 20], &symmetry));

        // Unequal grid is not compatible with full cubic symmetry
        assert!(!check_grid_compatibility([12, 12, 15], &symmetry));
    }

    #[test]
    fn test_idempotent() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12];
        let n = dims[0] * dims[1] * dims[2];

        let mut rho: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin().abs()).collect();
        symmetrize_density(&mut rho, dims, &symmetry);

        let rho_once = rho.clone();
        symmetrize_density(&mut rho, dims, &symmetry);

        // Applying symmetrization twice should give the same result
        for (a, b) in rho.iter().zip(rho_once.iter()) {
            assert!(
                (a - b).abs() < 1e-12,
                "symmetrization not idempotent: {a} vs {b}"
            );
        }
    }
}
