//! Charge density symmetrization on the FFT grid.
//!
//! The SCF uses the G-space form:
//!
//! ```text
//! ρ_sym(G) = (1/N_ops) Σ_S exp(-i G · τ_S) · ρ(R_S^T · G)
//! ```
//!
//! implemented by [`symmetrize_density_g`]. The fractional translation
//! enters analytically as an `exp(i·G·τ)` phase, so the result is exact
//! for any `τ` on any sufficiently band-limited density — including the
//! non-symmorphic case `τ=(¼,¼,¼)` on an 18³ grid that defeats a
//! real-space `nint`-based averager. Matches QE 7.5
//! `PW/src/symme.f90::sym_rho_serial`.

mod g_space;

use super::SymmetryInfo;

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
///
/// # Panics
///
/// Panics with a `BUG:` message if the internal fixed-size `[usize; 3]`
/// dims array is somehow empty when `iter().max()` is called. That is
/// structurally unreachable — a `[usize; 3]` always has three elements —
/// and the expect is present purely as a safety-belt against a future
/// refactor that changes the array shape.
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
            reason = "FFT grid dims are asserted <= MAX_FFT_DIM (1024) at `scf::grid::FftGrid::new`; fit in u32 trivially"
        )]
        let candidate = [
            crate::fft::fft_grid_size((dims[0] / 2) as u32),
            crate::fft::fft_grid_size((dims[1] / 2) as u32),
            crate::fft::fft_grid_size((dims[2] / 2) as u32),
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
