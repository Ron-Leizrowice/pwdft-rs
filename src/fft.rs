//! 3D FFT via ndrustfft (safe wrapper over rustfft + ndarray).
//!
//! Zero unsafe code. Uses ndarray's safe strided access for all
//! three dimension transforms.

use ndarray::Array3;
use ndrustfft::{FftHandler, Normalization, ndfft, ndifft};
use num_complex::Complex64;

/// 3D FFT on a regular grid.
///
/// Performs forward and inverse complex-to-complex transforms using
/// ndrustfft's safe ndarray-based 1D FFTs along each axis.
pub struct FFT3D {
    dims: [usize; 3],
    fwd_handlers: [FftHandler<f64>; 3],
    inv_handlers: [FftHandler<f64>; 3],
}

impl FFT3D {
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        // Forward: no normalization (standard convention: unnormalized forward)
        // Inverse: no normalization (we apply 1/N manually in inverse_normalized)
        Self {
            dims: [nx, ny, nz],
            fwd_handlers: [
                FftHandler::new(nx),
                FftHandler::new(ny),
                FftHandler::new(nz),
            ],
            inv_handlers: [
                FftHandler::<f64>::new(nx).normalization(Normalization::None),
                FftHandler::<f64>::new(ny).normalization(Normalization::None),
                FftHandler::<f64>::new(nz).normalization(Normalization::None),
            ],
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
        let [nx, ny, nz] = self.dims;
        assert_eq!(data.len(), nx * ny * nz);

        let mut a = Array3::from_shape_vec((nx, ny, nz), data.to_vec()).unwrap();
        let mut b = Array3::zeros((nx, ny, nz));

        ndfft(&a, &mut b, &self.fwd_handlers[0], 0);
        ndfft(&b, &mut a, &self.fwd_handlers[1], 1);
        ndfft(&a, &mut b, &self.fwd_handlers[2], 2);

        data.copy_from_slice(b.as_slice().unwrap());
    }

    /// Inverse FFT: reciprocal-space → real-space (unnormalized).
    /// Divide by `total_size()` afterwards for proper normalization.
    pub fn inverse(&mut self, data: &mut [Complex64]) {
        let [nx, ny, nz] = self.dims;
        assert_eq!(data.len(), nx * ny * nz);

        let mut a = Array3::from_shape_vec((nx, ny, nz), data.to_vec()).unwrap();
        let mut b = Array3::zeros((nx, ny, nz));

        ndifft(&a, &mut b, &self.inv_handlers[0], 0);
        ndifft(&b, &mut a, &self.inv_handlers[1], 1);
        ndifft(&a, &mut b, &self.inv_handlers[2], 2);

        data.copy_from_slice(b.as_slice().unwrap());
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
