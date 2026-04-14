use num_complex::Complex64;

// ============================================================================
// FFTW3 backend (feature = "fftw")
// ============================================================================

// Link the FFTW threads library for multi-threaded transforms.
// The search path is set via build.rs for the system FFTW installation.
#[cfg(feature = "fftw")]
#[link(name = "fftw3_threads")]
unsafe extern "C" {}

#[cfg(feature = "fftw")]
mod backend {
    use super::Complex64;
    use std::sync::{Mutex, Once};

    /// Global mutex for FFTW plan creation/destruction (not thread-safe in FFTW).
    static FFTW_PLANNER_LOCK: Mutex<()> = Mutex::new(());

    /// Initialize FFTW threading once.
    static INIT_THREADS: Once = Once::new();

    fn init_fftw_threads() {
        INIT_THREADS.call_once(|| {
            let ok = unsafe { fftw_sys::fftw_init_threads() };
            assert!(ok != 0, "fftw_init_threads failed");
        });
    }

    /// Choose thread count based on grid size.
    /// Threading overhead dominates for small grids; use 1 thread below 32³.
    fn threads_for_size(n: usize) -> i32 {
        if n < 32 * 32 * 32 {
            1
        } else {
            std::thread::available_parallelism()
                .map(|n| n.get() as i32)
                .unwrap_or(4)
        }
    }

    /// 3D FFT backed by FFTW3 with in-place transforms (zero copy).
    ///
    /// Uses fftw_sys directly to create in-place plans where in == out.
    /// Since fftw_complex = Complex64 (same type), we operate directly on
    /// the caller's data with no copies.
    pub struct FFT3D {
        dims: [usize; 3],
        plan_fwd: fftw_sys::fftw_plan,
        plan_inv: fftw_sys::fftw_plan,
        /// Planning buffer kept alive for FFTW's internal reference.
        /// FFTW_MEASURE profiles algorithms using this buffer; it may cache
        /// layout information, so the buffer must outlive the plans.
        _plan_buf: Vec<Complex64>,
    }

    // Safety: FFTW plan execution (fftw_execute_dft) is thread-safe for
    // distinct data pointers. Only plan creation/destruction is non-thread-safe,
    // and we serialize those with FFTW_PLANNER_LOCK.
    unsafe impl Send for FFT3D {}

    impl Drop for FFT3D {
        fn drop(&mut self) {
            let _lock = FFTW_PLANNER_LOCK.lock().unwrap();
            unsafe {
                fftw_sys::fftw_destroy_plan(self.plan_fwd);
                fftw_sys::fftw_destroy_plan(self.plan_inv);
            }
        }
    }

    impl FFT3D {
        pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
            init_fftw_threads();
            let _lock = FFTW_PLANNER_LOCK.lock().unwrap();
            let n = nx * ny * nz;
            unsafe { fftw_sys::fftw_plan_with_nthreads(threads_for_size(n)) };
            let mut plan_buf = vec![Complex64::new(0.0, 0.0); n];
            let ptr = plan_buf.as_mut_ptr();
            let dims = [nx as i32, ny as i32, nz as i32];
            let flags = fftw_sys::FFTW_MEASURE | fftw_sys::FFTW_UNALIGNED;

            let plan_fwd = unsafe {
                fftw_sys::fftw_plan_dft(
                    3,
                    dims.as_ptr(),
                    ptr, ptr,
                    fftw_sys::FFTW_FORWARD,
                    flags,
                )
            };
            assert!(!plan_fwd.is_null(), "FFTW forward plan creation failed");

            let plan_inv = unsafe {
                fftw_sys::fftw_plan_dft(
                    3,
                    dims.as_ptr(),
                    ptr, ptr,
                    fftw_sys::FFTW_BACKWARD as i32,
                    flags,
                )
            };
            assert!(!plan_inv.is_null(), "FFTW inverse plan creation failed");

            Self {
                dims: [nx, ny, nz],
                plan_fwd,
                plan_inv,
                _plan_buf: plan_buf,
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
            unsafe {
                fftw_sys::fftw_execute_dft(self.plan_fwd, data.as_mut_ptr(), data.as_mut_ptr());
            }
        }

        /// Inverse FFT: reciprocal-space → real-space (unnormalized).
        pub fn inverse(&mut self, data: &mut [Complex64]) {
            assert_eq!(data.len(), self.total_size());
            unsafe {
                fftw_sys::fftw_execute_dft(self.plan_inv, data.as_mut_ptr(), data.as_mut_ptr());
            }
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
// ndrustfft backend (default, no feature flag) — zero unsafe
// ============================================================================

#[cfg(not(feature = "fftw"))]
mod backend {
    use super::Complex64;
    use ndarray::Array3;
    use ndrustfft::{FftHandler, Normalization, ndfft, ndifft};

    /// 3D FFT backed by ndrustfft (safe wrapper over rustfft + ndarray).
    ///
    /// Uses ndarray's safe strided access for all three dimension transforms.
    /// No unsafe code, no raw pointer arithmetic, no provenance violations.
    pub struct FFT3D {
        dims: [usize; 3],
        fwd_handlers: [FftHandler<f64>; 3],
        inv_handlers: [FftHandler<f64>; 3],
    }

    impl FFT3D {
        pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
            // Forward: no normalization (matches FFTW/rustfft convention)
            // Inverse: no normalization (we do 1/N manually in inverse_normalized)
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

            // Wrap flat slice as 3D array (copies into owned array)
            let mut a = Array3::from_shape_vec((nx, ny, nz), data.to_vec()).unwrap();
            let mut b = Array3::zeros((nx, ny, nz));

            // Transform each axis: a → b → a → b
            ndfft(&a, &mut b, &self.fwd_handlers[0], 0);
            ndfft(&b, &mut a, &self.fwd_handlers[1], 1);
            ndfft(&a, &mut b, &self.fwd_handlers[2], 2);

            // Copy result back to caller's slice
            data.copy_from_slice(b.as_slice().unwrap());
        }

        /// Inverse FFT: reciprocal-space → real-space (unnormalized).
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
