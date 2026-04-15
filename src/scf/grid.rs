//! FFT grid management for SCF calculations.
//!
//! Handles the mapping between plane-wave basis, Miller indices,
//! and the real-space FFT grid used for potentials and density.

use nalgebra::Vector3;

use crate::{
    basis::BasisSet,
    fft::{fft_grid_size, FFT3D},
};

/// FFT grid with index mapping utilities.
pub(crate) struct FftGrid {
    pub dims: [usize; 3],
    pub fft: FFT3D,
    pub recip: crate::crystal::Lattice,
}

impl FftGrid {
    /// Create an FFT grid sized for the charge density.
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
            let n_max: Vec<i32> = (0..3)
                .map(|dim| miller.iter().map(|m| m[dim].abs()).max().unwrap_or(0))
                .collect();
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

/// Compute G-vector from FFT grid index (standalone, safe for parallel contexts).
pub(crate) fn g_vector_at_dims(
    idx: usize,
    dims: [usize; 3],
    recip: &crate::crystal::Lattice,
) -> Vector3<f64> {
    let [nx, ny, nz] = dims;
    let i1 = idx / (ny * nz);
    let i2 = (idx / nz) % ny;
    let i3 = idx % nz;
    let n1 = if i1 > nx / 2 { i1 as i32 - nx as i32 } else { i1 as i32 };
    let n2 = if i2 > ny / 2 { i2 as i32 - ny as i32 } else { i2 as i32 };
    let n3 = if i3 > nz / 2 { i3 as i32 - nz as i32 } else { i3 as i32 };
    n1 as f64 * recip.a + n2 as f64 * recip.b + n3 as f64 * recip.c
}

/// Map Miller indices to a flat FFT grid index.
pub(crate) fn miller_to_idx(dims: [usize; 3], n1: i32, n2: i32, n3: i32) -> usize {
    let i1 = ((n1 % dims[0] as i32) + dims[0] as i32) as usize % dims[0];
    let i2 = ((n2 % dims[1] as i32) + dims[1] as i32) as usize % dims[1];
    let i3 = ((n3 % dims[2] as i32) + dims[2] as i32) as usize % dims[2];
    i1 * dims[1] * dims[2] + i2 * dims[2] + i3
}
