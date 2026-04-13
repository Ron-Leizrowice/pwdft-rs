use num_complex::Complex64;
use rustfft::{FftPlanner, num_complex::Complex};
use std::sync::Arc;

/// 3D FFT on a regular grid, wrapping `rustfft`.
///
/// Performs forward and inverse transforms by batching 1D FFTs
/// along each dimension sequentially (x → y → z).
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
        Self {
            dims: [nx, ny, nz],
            fwd,
            inv,
        }
    }

    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    pub fn total_size(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Forward FFT: real-space → reciprocal-space (unnormalized).
    pub fn forward(&self, data: &mut [Complex64]) {
        assert_eq!(data.len(), self.total_size());
        self.transform_all_dims(data, &self.fwd);
    }

    /// Inverse FFT: reciprocal-space → real-space (unnormalized).
    /// Divide by `total_size()` afterwards for proper normalization.
    pub fn inverse(&self, data: &mut [Complex64]) {
        assert_eq!(data.len(), self.total_size());
        self.transform_all_dims(data, &self.inv);
    }

    /// Inverse FFT with normalization (divides by N).
    pub fn inverse_normalized(&self, data: &mut [Complex64]) {
        self.inverse(data);
        let norm = 1.0 / self.total_size() as f64;
        for v in data.iter_mut() {
            *v *= norm;
        }
    }

    fn transform_all_dims(&self, data: &mut [Complex64], plans: &[Arc<dyn rustfft::Fft<f64>>; 3]) {
        let [nx, ny, nz] = self.dims;

        // Transform along z (innermost, contiguous)
        {
            let plan = &plans[2];
            let mut scratch = vec![Complex::new(0.0, 0.0); plan.get_inplace_scratch_len()];
            for ix in 0..nx {
                for iy in 0..ny {
                    let offset = (ix * ny + iy) * nz;
                    let slice = &mut data[offset..offset + nz];
                    // rustfft Complex and num_complex Complex64 have identical layout
                    let slice =
                        unsafe { std::slice::from_raw_parts_mut(slice.as_mut_ptr().cast(), nz) };
                    plan.process_with_scratch(slice, &mut scratch);
                }
            }
        }

        // Transform along y
        {
            let plan = &plans[1];
            let mut scratch = vec![Complex::new(0.0, 0.0); plan.get_inplace_scratch_len()];
            let mut buf = vec![Complex::new(0.0, 0.0); ny];
            for ix in 0..nx {
                for iz in 0..nz {
                    // Gather y-stride into contiguous buffer
                    for iy in 0..ny {
                        let idx = (ix * ny + iy) * nz + iz;
                        buf[iy] = Complex::new(data[idx].re, data[idx].im);
                    }
                    plan.process_with_scratch(&mut buf, &mut scratch);
                    // Scatter back
                    for iy in 0..ny {
                        let idx = (ix * ny + iy) * nz + iz;
                        data[idx] = Complex64::new(buf[iy].re, buf[iy].im);
                    }
                }
            }
        }

        // Transform along x
        {
            let plan = &plans[0];
            let mut scratch = vec![Complex::new(0.0, 0.0); plan.get_inplace_scratch_len()];
            let mut buf = vec![Complex::new(0.0, 0.0); nx];
            for iy in 0..ny {
                for iz in 0..nz {
                    // Gather x-stride
                    for ix in 0..nx {
                        let idx = (ix * ny + iy) * nz + iz;
                        buf[ix] = Complex::new(data[idx].re, data[idx].im);
                    }
                    plan.process_with_scratch(&mut buf, &mut scratch);
                    // Scatter back
                    for ix in 0..nx {
                        let idx = (ix * ny + iy) * nz + iz;
                        data[idx] = Complex64::new(buf[ix].re, buf[ix].im);
                    }
                }
            }
        }
    }
}

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
        while n % p == 0 {
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
        let fft = FFT3D::new(4, 4, 4);
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
        let fft = FFT3D::new(8, 8, 8);
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
        assert_eq!(fft_grid_size(3), 8); // 2*3+1=7 → next FFT-friendly is 8
        assert_eq!(fft_grid_size(4), 9); // 2*4+1=9 = 3^2
        assert_eq!(fft_grid_size(5), 12); // 2*5+1=11 → 12
    }

    #[test]
    fn test_is_fft_friendly() {
        assert!(is_fft_friendly(1));
        assert!(is_fft_friendly(2));
        assert!(is_fft_friendly(4));
        assert!(is_fft_friendly(8));
        assert!(is_fft_friendly(30)); // 2*3*5
        assert!(!is_fft_friendly(7));
        assert!(!is_fft_friendly(11));
    }
}
