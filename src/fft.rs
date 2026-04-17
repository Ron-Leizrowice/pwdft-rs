//! 3D FFT via ndrustfft (safe wrapper over rustfft + ndarray).
//!
//! Zero unsafe code. Uses ndarray's safe strided access for all
//! three dimension transforms.

use ndarray::Array3;
use ndrustfft::{FftHandler, Normalization, ndfft, ndifft};
use num_complex::Complex64;

/// 3D FFT via batched 1D transforms (z → y → x).
///
/// Convention:
///   Forward:  f̃(G) = Σ_r f(r) e^{-iG·r}     (unnormalized)
///   Inverse:  f(r) = Σ_G f̃(G) e^{+iG·r}     (unnormalized)
///
/// The forward FFT is unnormalized; callers must divide by N = nx·ny·nz
/// to get Fourier coefficients. Use `inverse_normalized()` for the
/// convention f(r) = (1/N) Σ_G f̃(G) e^{+iG·r}.
pub struct FFT3D {
    dims: [usize; 3],
    fwd_handlers: [FftHandler<f64>; 3],
    inv_handlers: [FftHandler<f64>; 3],
    /// Scratch buffer A — reused across `forward`/`inverse` calls to avoid
    /// per-call `Array3<Complex64>` allocations (~1 MB at 32^3).
    buf_a: Array3<Complex64>,
    /// Scratch buffer B — paired with `buf_a` for ping-pong 1D transforms.
    buf_b: Array3<Complex64>,
}

impl FFT3D {
    #[must_use]
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
            buf_a: Array3::zeros((nx, ny, nz)),
            buf_b: Array3::zeros((nx, ny, nz)),
        }
    }

    #[must_use]
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    #[must_use]
    pub fn total_size(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Forward FFT: real-space → reciprocal-space (unnormalized).
    pub fn forward(&mut self, data: &mut [Complex64]) {
        let [nx, ny, nz] = self.dims;
        assert_eq!(data.len(), nx * ny * nz);

        // Reuse pre-allocated scratch buffers (allocated once in `new`).
        // `buf_a` is a row-major contiguous Array3, so `as_slice_mut` always
        // succeeds. Copy input into `buf_a` without reallocating.
        self.buf_a
            .as_slice_mut()
            .expect("BUG: Array3 should be contiguous")
            .copy_from_slice(data);

        ndfft(&self.buf_a, &mut self.buf_b, &self.fwd_handlers[0], 0);
        ndfft(&self.buf_b, &mut self.buf_a, &self.fwd_handlers[1], 1);
        ndfft(&self.buf_a, &mut self.buf_b, &self.fwd_handlers[2], 2);

        data.copy_from_slice(
            self.buf_b
                .as_slice()
                .expect("BUG: Array3 should be contiguous"),
        );
    }

    /// Inverse FFT: reciprocal-space → real-space (unnormalized).
    /// Divide by `total_size()` afterwards for proper normalization.
    pub fn inverse(&mut self, data: &mut [Complex64]) {
        let [nx, ny, nz] = self.dims;
        assert_eq!(data.len(), nx * ny * nz);

        // Reuse pre-allocated scratch buffers.
        self.buf_a
            .as_slice_mut()
            .expect("BUG: Array3 should be contiguous")
            .copy_from_slice(data);

        ndifft(&self.buf_a, &mut self.buf_b, &self.inv_handlers[0], 0);
        ndifft(&self.buf_b, &mut self.buf_a, &self.inv_handlers[1], 1);
        ndifft(&self.buf_a, &mut self.buf_b, &self.inv_handlers[2], 2);

        data.copy_from_slice(
            self.buf_b
                .as_slice()
                .expect("BUG: Array3 should be contiguous"),
        );
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
/// Find the smallest FFT-friendly grid size n ≥ 2·n_max + 1.
///
/// The factor 2·n_max + 1 is the Nyquist criterion: G-vectors range from
/// -n_max to +n_max, requiring at least 2·n_max + 1 grid points to avoid
/// aliasing when computing products like V(G-G') in the Hamiltonian.
///
/// Grid sizes that are products of small primes (2, 3, 5) give optimal FFT
/// performance; arbitrary sizes may be much slower.
#[must_use]
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
