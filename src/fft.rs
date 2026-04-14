use num_complex::Complex64;

// ============================================================================
// FFTW3 backend (feature = "fftw")
// ============================================================================

#[cfg(feature = "fftw")]
mod backend {
    use super::Complex64;
    use fftw::array::AlignedVec;
    use fftw::plan::*;
    use fftw::types::*;

    /// 3D FFT backed by FFTW3.
    ///
    /// Uses native multi-dimensional C2C transforms — no manual gather/scatter.
    /// Plans are created with MEASURE for optimal performance (amortized over
    /// many SCF iterations).
    pub struct FFT3D {
        dims: [usize; 3],
        plan_fwd: C2CPlan64,
        plan_inv: C2CPlan64,
        /// Aligned scratch buffer for FFTW execution
        buf: AlignedVec<c64>,
    }

    impl FFT3D {
        pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
            let n = nx * ny * nz;
            let plan_fwd =
                C2CPlan64::aligned(&[nx, ny, nz], Sign::Forward, Flag::MEASURE).unwrap();
            let plan_inv =
                C2CPlan64::aligned(&[nx, ny, nz], Sign::Backward, Flag::MEASURE).unwrap();
            let buf = AlignedVec::new(n);
            Self {
                dims: [nx, ny, nz],
                plan_fwd,
                plan_inv,
                buf,
            }
        }

        pub fn dims(&self) -> [usize; 3] {
            self.dims
        }

        pub fn total_size(&self) -> usize {
            self.dims[0] * self.dims[1] * self.dims[2]
        }

        /// Forward FFT: real-space → reciprocal-space (unnormalized).
        pub fn forward(&mut self, data: &mut [Complex64]) {
            assert_eq!(data.len(), self.total_size());
            let n = self.total_size();
            // Copy into aligned buffer, execute, copy back.
            // fftw::c64 and Complex64 are the same type.
            let in_buf = unsafe { std::slice::from_raw_parts_mut(data.as_mut_ptr().cast(), n) };
            let out_buf: &mut [c64] = &mut self.buf;
            self.plan_fwd.c2c(in_buf, out_buf).unwrap();
            // Copy result back
            let out_slice =
                unsafe { std::slice::from_raw_parts(out_buf.as_ptr().cast::<Complex64>(), n) };
            data.copy_from_slice(out_slice);
        }

        /// Inverse FFT: reciprocal-space → real-space (unnormalized).
        pub fn inverse(&mut self, data: &mut [Complex64]) {
            assert_eq!(data.len(), self.total_size());
            let n = self.total_size();
            let in_buf = unsafe { std::slice::from_raw_parts_mut(data.as_mut_ptr().cast(), n) };
            let out_buf: &mut [c64] = &mut self.buf;
            self.plan_inv.c2c(in_buf, out_buf).unwrap();
            let out_slice =
                unsafe { std::slice::from_raw_parts(out_buf.as_ptr().cast::<Complex64>(), n) };
            data.copy_from_slice(out_slice);
        }

        /// Inverse FFT with normalization (divides by N).
        pub fn inverse_normalized(&mut self, data: &mut [Complex64]) {
            self.inverse(data);
            let norm = 1.0 / self.total_size() as f64;
            for v in data.iter_mut() {
                *v *= norm;
            }
        }
    }
}

// ============================================================================
// rustfft backend (default, no feature flag)
// ============================================================================

#[cfg(not(feature = "fftw"))]
mod backend {
    use super::Complex64;
    use rayon::prelude::*;
    use rustfft::{FftPlanner, num_complex::Complex};
    use std::sync::Arc;

    /// 3D FFT backed by rustfft with rayon parallelism.
    ///
    /// Performs forward and inverse transforms by batching 1D FFTs
    /// along each dimension (z → y → x). The z-dimension transforms
    /// operate on contiguous memory. The y and x transforms use
    /// gather/scatter with thread-local buffers.
    pub struct FFT3D {
        dims: [usize; 3],
        fwd: [Arc<dyn rustfft::Fft<f64>>; 3],
        inv: [Arc<dyn rustfft::Fft<f64>>; 3],
    }

    impl FFT3D {
        pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
            let mut planner = FftPlanner::<f64>::new();
            let fwd = [
                planner.plan_fft_forward(nx),
                planner.plan_fft_forward(ny),
                planner.plan_fft_forward(nz),
            ];
            let inv = [
                planner.plan_fft_inverse(nx),
                planner.plan_fft_inverse(ny),
                planner.plan_fft_inverse(nz),
            ];
            Self { dims: [nx, ny, nz], fwd, inv }
        }

        pub fn dims(&self) -> [usize; 3] {
            self.dims
        }

        pub fn total_size(&self) -> usize {
            self.dims[0] * self.dims[1] * self.dims[2]
        }

        /// Forward FFT: real-space → reciprocal-space (unnormalized).
        pub fn forward(&mut self, data: &mut [Complex64]) {
            assert_eq!(data.len(), self.total_size());
            self.transform_all_dims(data, &self.fwd.clone());
        }

        /// Inverse FFT: reciprocal-space → real-space (unnormalized).
        pub fn inverse(&mut self, data: &mut [Complex64]) {
            assert_eq!(data.len(), self.total_size());
            self.transform_all_dims(data, &self.inv.clone());
        }

        /// Inverse FFT with normalization (divides by N).
        pub fn inverse_normalized(&mut self, data: &mut [Complex64]) {
            self.inverse(data);
            let norm = 1.0 / self.total_size() as f64;
            for v in data.iter_mut() {
                *v *= norm;
            }
        }

        fn transform_all_dims(
            &self,
            data: &mut [Complex64],
            plans: &[Arc<dyn rustfft::Fft<f64>>; 3],
        ) {
            let [nx, ny, nz] = self.dims;

            // Transform along z (innermost, contiguous)
            {
                let plan = &plans[2];
                data.par_chunks_mut(nz).for_each(|slab| {
                    let mut scratch =
                        vec![Complex::new(0.0, 0.0); plan.get_inplace_scratch_len()];
                    let slab = unsafe {
                        std::slice::from_raw_parts_mut(slab.as_mut_ptr().cast(), nz)
                    };
                    plan.process_with_scratch(slab, &mut scratch);
                });
            }

            // Transform along y — parallel over ix slabs
            {
                let plan = &plans[1];
                let slab_ny_nz = ny * nz;
                let slabs: Vec<&mut [Complex64]> = data.chunks_mut(slab_ny_nz).collect();
                slabs.into_par_iter().for_each(|slab| {
                    let mut scratch =
                        vec![Complex::new(0.0, 0.0); plan.get_inplace_scratch_len()];
                    let mut buf = vec![Complex::new(0.0, 0.0); ny];
                    for iz in 0..nz {
                        for iy in 0..ny {
                            let v = slab[iy * nz + iz];
                            buf[iy] = Complex::new(v.re, v.im);
                        }
                        plan.process_with_scratch(&mut buf, &mut scratch);
                        for iy in 0..ny {
                            slab[iy * nz + iz] = Complex64::new(buf[iy].re, buf[iy].im);
                        }
                    }
                });
            }

            // Transform along x — parallel over (iy, iz) pairs
            {
                let plan = &plans[0];
                let stride = ny * nz;
                let base = data.as_mut_ptr() as usize;
                let data_len = data.len();

                (0..ny * nz).into_par_iter().for_each(|yz| {
                    let ptr = base as *mut Complex64;
                    let mut scratch =
                        vec![Complex::new(0.0, 0.0); plan.get_inplace_scratch_len()];
                    let mut buf = vec![Complex::new(0.0, 0.0); nx];
                    for ix in 0..nx {
                        let idx = ix * stride + yz;
                        debug_assert!(idx < data_len);
                        let v = unsafe { *ptr.add(idx) };
                        buf[ix] = Complex::new(v.re, v.im);
                    }
                    plan.process_with_scratch(&mut buf, &mut scratch);
                    for ix in 0..nx {
                        let idx = ix * stride + yz;
                        unsafe {
                            *ptr.add(idx) = Complex64::new(buf[ix].re, buf[ix].im)
                        };
                    }
                });
            }
        }
    }
}

// Re-export the active backend
pub use backend::FFT3D;

/// Choose FFT-friendly grid dimensions for a given basis.
///
/// Returns the smallest `n >= 2*n_max + 1` that is a product of small primes (2,3,5).
pub fn fft_grid_size(n_max: i32) -> usize {
    let min_n = (2 * n_max + 1) as usize;
    let mut n = min_n;
    loop {
        if is_fft_friendly(n) {
            return n;
        }
        n += 1;
    }
}

/// Check if n is a product of 2, 3, 5 only (FFT-friendly size).
fn is_fft_friendly(mut n: usize) -> bool {
    if n == 0 {
        return false;
    }
    for &p in &[2, 3, 5] {
        while n.is_multiple_of(p) {
            n /= p;
        }
    }
    n == 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fft_roundtrip() {
        let mut fft = FFT3D::new(4, 4, 4);
        let n = fft.total_size();
        let mut data: Vec<Complex64> = (0..n)
            .map(|i| Complex64::new(i as f64, 0.0))
            .collect();
        let original = data.clone();

        fft.forward(&mut data);
        fft.inverse_normalized(&mut data);

        for (i, (got, want)) in data.iter().zip(original.iter()).enumerate() {
            assert!(
                (got - want).norm() < 1e-10,
                "roundtrip failed at index {i}: got {got}, expected {want}"
            );
        }
    }

    #[test]
    fn test_fft_parseval() {
        // Parseval's theorem: sum|f(x)|^2 = (1/N) sum|F(k)|^2
        let mut fft = FFT3D::new(8, 8, 8);
        let n = fft.total_size();
        let mut data: Vec<Complex64> = (0..n)
            .map(|i| Complex64::new((i as f64 * 0.1).sin(), (i as f64 * 0.2).cos()))
            .collect();
        let real_sum: f64 = data.iter().map(|c| c.norm_sqr()).sum();

        fft.forward(&mut data);
        let recip_sum: f64 = data.iter().map(|c| c.norm_sqr()).sum();

        let ratio = recip_sum / real_sum;
        assert!(
            (ratio - n as f64).abs() < 1e-8,
            "Parseval: recip/real = {ratio}, expected {n}"
        );
    }

    #[test]
    fn test_fft_grid_size() {
        assert_eq!(fft_grid_size(3), 8);
        assert_eq!(fft_grid_size(4), 9);
        assert_eq!(fft_grid_size(5), 12);
    }

    #[test]
    fn test_is_fft_friendly() {
        assert!(is_fft_friendly(1));
        assert!(is_fft_friendly(2));
        assert!(is_fft_friendly(4));
        assert!(is_fft_friendly(8));
        assert!(is_fft_friendly(30));
        assert!(!is_fft_friendly(7));
        assert!(!is_fft_friendly(11));
    }
}
