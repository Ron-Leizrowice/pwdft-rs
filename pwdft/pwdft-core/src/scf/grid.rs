//! FFT grid management for SCF calculations.
//!
//! Handles the mapping between plane-wave basis, Miller indices,
//! and the real-space FFT grid used for potentials and density.

use nalgebra::Vector3;

use crate::{
    basis::BasisSet,
    fft::{FFT3D, fft_grid_size},
};

/// Hard upper bound on per-axis FFT grid dimensions.
///
/// Load-bearing for the `cast_possible_truncation` / `cast_possible_wrap`
/// annotations in this module and in `src/symmetry/density/`: those sites
/// rely on grid dims fitting in `i32` (for signed Miller arithmetic).
/// 1024 is well above any realistic physical grid — pwdft-core SCF runs
/// typically use 18-200 per axis, a 100 Ry ecut on a tight cell reaches
/// ~256, and the jump to 1024 leaves ~4× head-room before silent
/// truncation would occur at `as i32`.
///
/// Enforced by [`FftGrid::new`] via a runtime `assert!` so violating
/// inputs abort cleanly rather than producing wrap-around corruption
/// in Miller-to-flat index arithmetic.
pub(crate) const MAX_FFT_DIM: usize = 1024;

/// FFT grid with index mapping utilities.
pub(crate) struct FftGrid {
    pub dims: [usize; 3],
    pub fft: FFT3D,
    pub recip: crate::crystal::Lattice,
}

impl FftGrid {
    /// Create an FFT grid sized for the charge density.
    ///
    /// # Panics
    ///
    /// Panics if any resulting FFT dimension exceeds [`MAX_FFT_DIM`].
    /// This is the runtime enforcement of the `grid ≤ 1024` invariant
    /// that the `cast_*` `#[allow]` sites in this module and in
    /// `src/symmetry/density/` rely on for correctness.
    pub fn new(
        basis: &BasisSet,
        lattice: &crate::crystal::Lattice,
        ecutrho_ratio: u32,
        explicit_dims: Option<[usize; 3]>,
    ) -> Self {
        let dims = if let Some(d) = explicit_dims {
            d
        } else {
            let miller = basis.miller_indices();
            // Maximum absolute Miller index per axis. `unsigned_abs` returns
            // `u32` directly, which matches `fft_grid_size`'s signature.
            let n_max: Vec<u32> = (0..3)
                .map(|dim| miller.iter().map(|m| m[dim].unsigned_abs()).max().unwrap_or(0))
                .collect();
            // `ecutrho_ratio` is a small integer (typically 4). The grid
            // scaling factor is `ceil(sqrt(ecutrho_ratio))`; compute via
            // integer `isqrt` so the value stays in the type system.
            let isqrt = ecutrho_ratio.isqrt();
            let scale = if isqrt * isqrt == ecutrho_ratio {
                isqrt
            } else {
                isqrt + 1
            };
            [
                fft_grid_size(scale * n_max[0]),
                fft_grid_size(scale * n_max[1]),
                fft_grid_size(scale * n_max[2]),
            ]
        };
        assert!(
            dims[0] <= MAX_FFT_DIM && dims[1] <= MAX_FFT_DIM && dims[2] <= MAX_FFT_DIM,
            "FFT grid too large: {}x{}x{} exceeds MAX_FFT_DIM={} \
             (per-axis bound is load-bearing for i32/u32 casts in grid / symmetry / GPU code)",
            dims[0],
            dims[1],
            dims[2],
            MAX_FFT_DIM,
        );
        let fft = FFT3D::new(dims[0], dims[1], dims[2]);
        let recip = lattice.reciprocal();
        Self { dims, fft, recip }
    }

    pub fn total_size(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    pub fn miller_to_idx(&self, n1: i32, n2: i32, n3: i32) -> usize {
        miller_to_idx(self.dims, n1, n2, n3)
    }

    pub fn g_vector_at(&self, idx: usize) -> Vector3<f64> {
        g_vector_at_dims(idx, self.dims, &self.recip)
    }

    pub fn basis_to_fft(&self, basis: &BasisSet) -> Vec<usize> {
        basis
            .miller_indices()
            .iter()
            .map(|&[n1, n2, n3]| self.miller_to_idx(n1, n2, n3))
            .collect()
    }
}

/// Compute G-vector from FFT grid index (standalone, safe for parallel
/// contexts).
///
/// All `usize -> i32` casts below are safe because FFT grid dimensions
/// are bounded by [`MAX_FFT_DIM`] = 1024, asserted at [`FftGrid::new`]
/// construction time; i32::MAX is 2^31 ≈ 2.1e9.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "FFT grid dims nx,ny,nz are asserted <= MAX_FFT_DIM (1024) at FftGrid::new; i1,i2,i3 < nx,ny,nz so all fit in i32. Miller indices n1,n2,n3 = i - n*bool are in [-n/2, n/2] and fit in i32 likewise."
)]
pub(crate) fn g_vector_at_dims(idx: usize, dims: [usize; 3], recip: &crate::crystal::Lattice) -> Vector3<f64> {
    let [nx, ny, nz] = dims;
    let i1 = idx / (ny * nz);
    let i2 = (idx / nz) % ny;
    let i3 = idx % nz;
    let n1 = if i1 > nx / 2 { i1 as i32 - nx as i32 } else { i1 as i32 };
    let n2 = if i2 > ny / 2 { i2 as i32 - ny as i32 } else { i2 as i32 };
    let n3 = if i3 > nz / 2 { i3 as i32 - nz as i32 } else { i3 as i32 };
    f64::from(n1) * recip.a + f64::from(n2) * recip.b + f64::from(n3) * recip.c
}

/// Map Miller indices to a flat FFT grid index.
///
/// The `((x % d) + d) as usize % d` idiom computes a non-negative
/// remainder: `(n % d)` is in `[-(d-1), d-1]`, so `(n % d) + d` is in
/// `[1, 2d-1]` — always non-negative, so `as usize` loses no sign.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "dims are FFT grid sizes asserted <= MAX_FFT_DIM (1024) at FftGrid::new; Miller indices n1,n2,n3 are bounded by ecut (fit in i32). The sum `(n % d) + d` is mathematically non-negative so `as usize` loses no sign."
)]
pub(crate) fn miller_to_idx(dims: [usize; 3], n1: i32, n2: i32, n3: i32) -> usize {
    let i1 = ((n1 % dims[0] as i32) + dims[0] as i32) as usize % dims[0];
    let i2 = ((n2 % dims[1] as i32) + dims[1] as i32) as usize % dims[1];
    let i3 = ((n3 % dims[2] as i32) + dims[2] as i32) as usize % dims[2];
    i1 * dims[1] * dims[2] + i2 * dims[2] + i3
}
