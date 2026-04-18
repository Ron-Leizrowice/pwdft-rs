//! Charge density symmetrization on the FFT grid.
//!
//! Two forms are implemented:
//!
//! * [`symmetrize_density`] — real-space form
//!   `ρ_sym(r) = (1/N_ops) Σ_S ρ(S⁻¹ r)`, implemented by rotating each
//!   grid point and rounding to the nearest neighbour. Exact only when
//!   `τ_{S,i} · n_i ∈ ℤ` for every operation S and axis i — i.e. the
//!   fractional translation lands on an integer grid point.
//!
//! * [`symmetrize_density_g`] — G-space form (preferred for SCF)
//!   `ρ_sym(G) = (1/N_ops) Σ_S exp(-i G · τ_S) · ρ(R_S⁻¹ G)`,
//!   exact for any fractional translation because the phase factor is
//!   analytic. Matches QE 7.5 `PW/src/symme.f90::sym_rho_serial`.
//!
//! The real-space form silently smears density across the wrong grid
//! points for non-symmorphic space groups whose τ doesn't land on the
//! grid (e.g. Fd-3m with τ=(¼,¼,¼) on an 18³ grid). Always use the
//! G-space form in SCF; keep the real-space form for direct-grid unit
//! tests where translation exactness is guaranteed by construction.
//!
//! See `proposals/completed/PCFX-symmetrize-rho-g-space.md`.

mod g_space;
mod real_space;

use super::SymmetryInfo;

#[allow(deprecated)]
pub use real_space::symmetrize_density;

pub use g_space::symmetrize_density_g;

/// Check if the FFT grid dimensions are compatible with all symmetry operations.
///
/// For integer rotation R in fractional coords, the grid point (ix, iy, iz) maps
/// to another exact grid point if and only if R_{ij} × n_j ≡ 0 (mod n_i) for all i, j.
///
/// Returns true if all operations are compatible.
#[must_use]
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "FFT grid dims are asserted <= MAX_FFT_DIM (1024) at `scf::grid::FftGrid::new`; fit in i32 trivially"
)]
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
#[must_use]
pub fn compatible_grid_dims(min_dims: [usize; 3], symmetry: &SymmetryInfo) -> [usize; 3] {
    // For cubic symmetry, making all dimensions equal is usually sufficient
    // SAFETY: min_dims is [usize; 3], always has 3 elements -- max() cannot be None.
    let max_dim = *min_dims.iter().max().expect("BUG: empty fixed-size array");
    let mut dims = [max_dim; 3];

    // Try increasing until compatible
    for _ in 0..100 {
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_possible_wrap,
            reason = "FFT grid dims are asserted <= MAX_FFT_DIM (1024) at `scf::grid::FftGrid::new`; fit in i32 trivially"
        )]
        let candidate = [
            crate::fft::fft_grid_size(dims[0] as i32 / 2),
            crate::fft::fft_grid_size(dims[1] as i32 / 2),
            crate::fft::fft_grid_size(dims[2] as i32 / 2),
        ];
        // Make all equal to the max for safety
        // SAFETY: candidate is [usize; 3], always has 3 elements.
        let m = *candidate.iter().max().expect("BUG: empty fixed-size array");
        let candidate = [m, m, m];
        if check_grid_compatibility(candidate, symmetry) {
            return candidate;
        }
        dims = [dims[0] + 1, dims[1] + 1, dims[2] + 1];
    }

    // Fallback: just use the input dims (may not be compatible)
    min_dims
}
