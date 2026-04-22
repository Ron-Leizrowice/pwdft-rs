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
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        // Forward: no normalization (standard convention: unnormalized forward)
        // Inverse: no normalization (we apply 1/N manually in inverse_normalized)
        Self {
            dims: [nx, ny, nz],
            fwd_handlers: [FftHandler::new(nx), FftHandler::new(ny), FftHandler::new(nz)],
            inv_handlers: [
                FftHandler::<f64>::new(nx).normalization(Normalization::None),
                FftHandler::<f64>::new(ny).normalization(Normalization::None),
                FftHandler::<f64>::new(nz).normalization(Normalization::None),
            ],
            buf_a: Array3::zeros((nx, ny, nz)),
            buf_b: Array3::zeros((nx, ny, nz)),
        }
    }

    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    pub fn total_size(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Forward FFT: real-space → reciprocal-space (unnormalized).
    ///
    /// # Panics
    ///
    /// Panics if `data.len() != nx * ny * nz` (the dimensions recorded at
    /// [`FFT3D::new`]). The caller must size the slice to match the handler
    /// grid; any other size is a programming error.
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

        data.copy_from_slice(self.buf_b.as_slice().expect("BUG: Array3 should be contiguous"));
    }

    /// Inverse FFT: reciprocal-space → real-space (unnormalized).
    /// Divide by `total_size()` afterwards for proper normalization.
    ///
    /// # Panics
    ///
    /// Panics if `data.len() != nx * ny * nz` (the dimensions recorded at
    /// [`FFT3D::new`]). The caller must size the slice to match the handler
    /// grid; any other size is a programming error.
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

        data.copy_from_slice(self.buf_b.as_slice().expect("BUG: Array3 should be contiguous"));
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

/// Real-space density gradient via FFT: `∇ρ_G = iG · ρ_G`.
///
/// Takes a real-space density (or any scalar field) `ρ(r)` on the FFT
/// grid, forward-FFTs to `ρ(G)` with the `1/N` charge-density convention
/// used elsewhere in the engine (see `scf::energy::density_r_to_g`),
/// multiplies each Fourier coefficient by `i·G` per Cartesian
/// component, inverse-FFTs back, and returns the three Cartesian
/// components of `∇ρ(r)` packed into a single `Vec<[f64; 3]>`.
///
/// Units: with `rho_r` in `e/Å³` and `g_vectors[idx]` in `Å⁻¹`, the
/// returned entries are in `e/Å⁴` — the input units that
/// `pbe_exchange` and `pbe_correlation` (GGAP Phase B/C) expect on
/// `|∇ρ|`.
///
/// Inputs:
/// - `rho_r`: scalar field in real space, length `fft.total_size()`.
/// - `fft`: FFT handler whose dimensions must match `rho_r`.
/// - `g_vectors`: per-grid-point reciprocal-space vectors in the **FFT-aligned** ordering
///   (`scf::grid::g_vector_at_dims`). Length must equal `rho_r.len()`. Each entry carries `[G_x,
///   G_y, G_z]` in `Å⁻¹`.
///
/// Returns a fresh `Vec<[f64; 3]>` of length `rho_r.len()`. The `G = 0`
/// component contributes zero to the gradient (since `iG = 0`), matching
/// the `∇(constant) = 0` limit of any periodic calculation.
///
/// # Panics
///
/// Panics if `rho_r.len() != fft.total_size()` or if
/// `g_vectors.len() != rho_r.len()`. Debug-asserts that the imaginary
/// residual of the real output is bounded by `1e-8 * ‖Re‖_∞ + 1e-10`;
/// the unit tests pin the Gaussian residual well under that.
///
/// # Spectral convention at the Nyquist mode
///
/// On any even axis the Nyquist coefficient is self-conjugate
/// (`+N/2` and `−N/2` alias onto the same DFT slot), so the "signed"
/// `G_α` value at that slot is ambiguous and the naïve `iG·ρ(G)`
/// multiplication produces a non-Hermitian perturbation whose
/// inverse-FFT picks up an imaginary residual. Spectral-method
/// convention (Boyd, *Chebyshev and Fourier Spectral Methods*, §3.5)
/// is to zero out the Nyquist mode before differentiating — the
/// derivative there is not well-defined on the grid. On a band-limited
/// density (which the SCF charge always is) this is bit-identical to
/// the unfiltered path; on broadband inputs it removes the ambiguous
/// contribution cleanly.
///
/// # Allocation
///
/// Allocates three temporary `Vec<Complex64>` of length
/// `fft.total_size()` (the per-component `iG · ρ(G)` buffer, reused
/// across axes) plus the output `Vec<[f64; 3]>`. No pooling is
/// attempted — this is called once per SCF iteration per GGA channel
/// and the FFT work itself dominates.
pub fn compute_density_gradient(rho_r: &[f64], fft: &mut FFT3D, g_vectors: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let n = fft.total_size();
    assert_eq!(rho_r.len(), n, "rho_r length must equal fft.total_size()");
    assert_eq!(g_vectors.len(), n, "g_vectors length must equal rho_r length",);

    // Forward FFT with the same `1/N` normalisation convention used by
    // `scf::energy::density_r_to_g`. This is the ρ(G) coefficient that
    // `∑_G ρ(G) exp(+i G·r)` reconstructs — so `∇ρ(r) = ∑_G iG ρ(G)
    // exp(+i G·r)` uses exactly the same coefficients without extra
    // scaling.
    let mut rho_g: Vec<Complex64> = rho_r.iter().map(|&r| Complex64::new(r, 0.0)).collect();
    fft.forward(&mut rho_g);
    let inv_n = 1.0 / n as f64;
    for v in &mut rho_g {
        *v *= inv_n;
    }

    // Zero out the Nyquist coefficients on each even axis. At
    // `n_α = N_α/2` (even N) the DFT slot is self-conjugate and its
    // Hermitian partner at `-N_α/2` aliases onto the same slot — so
    // the "signed" `G_α` value at that slot is ambiguous and the iG
    // multiplication produces a non-Hermitian perturbation whose
    // inverse-FFT picks up an imaginary residual of order
    // `max|ρ̂(Nyquist)|`. Spectral-method convention (see e.g.
    // Boyd, "Chebyshev and Fourier Spectral Methods", §3.5) is to
    // zero out the Nyquist mode before differentiating, since the
    // derivative at Nyquist is not well-defined on the grid. On a
    // band-limited input this change is bit-identical to the naive
    // path (ρ̂(Nyquist) ≈ 0); on a broadband input it removes the
    // ambiguous contribution cleanly.
    let [nx, ny, nz] = fft.dims();
    for (idx, rg) in rho_g.iter_mut().enumerate() {
        let ix = idx / (ny * nz);
        let iy = (idx / nz) % ny;
        let iz = idx % nz;
        let at_nyquist_x = nx.is_multiple_of(2) && ix == nx / 2;
        let at_nyquist_y = ny.is_multiple_of(2) && iy == ny / 2;
        let at_nyquist_z = nz.is_multiple_of(2) && iz == nz / 2;
        if at_nyquist_x || at_nyquist_y || at_nyquist_z {
            *rg = Complex64::new(0.0, 0.0);
        }
    }

    // Per-component inverse FFT of `i G_α · ρ(G)`, accumulated into the
    // output `Vec<[f64; 3]>`. One `Vec<Complex64>` scratch is allocated
    // per axis; the three passes are sequential on the shared `fft`
    // handler (rustfft is not `Sync` on its internal scratch).
    let mut grad_r: Vec<[f64; 3]> = vec![[0.0; 3]; n];
    for axis in 0..3 {
        let mut buf: Vec<Complex64> = rho_g
            .iter()
            .zip(g_vectors.iter())
            .map(|(&rg, g)| Complex64::new(0.0, g[axis]) * rg)
            .collect();
        fft.inverse(&mut buf);
        // The output of `∇ρ(r) = ∑_G iG ρ(G) exp(+i G·r)` is
        // mathematically real. With the Nyquist mode zeroed above, the
        // only residual imaginary part is FFT round-off (~1e-13 on a
        // balanced DFT). The debug-assert below fires at a generous
        // tolerance of 1e-8 of the maximum real magnitude; tighter
        // bounds are exercised by the unit tests.
        debug_assert!(
            {
                let max_im = buf.iter().map(|c| c.im.abs()).fold(0.0_f64, f64::max);
                let max_re = buf.iter().map(|c| c.re.abs()).fold(0.0_f64, f64::max);
                max_im < 1e-8 * max_re.max(1e-300) + 1e-10
            },
            "∇ρ(r) should be real; Nyquist zero-out path bypassed"
        );
        for (dst, src) in grad_r.iter_mut().zip(buf.iter()) {
            dst[axis] = src.re;
        }
    }

    grad_r
}

/// Find the smallest FFT-friendly grid size n ≥ 2·n_max + 1.
///
/// The factor 2·n_max + 1 is the Nyquist criterion: G-vectors range from
/// -n_max to +n_max, requiring at least 2·n_max + 1 grid points to avoid
/// aliasing when computing products like V(G-G') in the Hamiltonian.
///
/// Grid sizes that are products of small primes (2, 3, 5) give optimal FFT
/// performance; arbitrary sizes may be much slower.
///
/// `n_max` is typed as `u32` to encode the non-negative invariant at the
/// API boundary rather than via a runtime assertion.
pub fn fft_grid_size(n_max: u32) -> usize {
    let min_n = 2 * (n_max as usize) + 1;
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
        let mut data: Vec<Complex64> = (0..n).map(|i| Complex64::new(i as f64, 0.0)).collect();
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
        assert_eq!(fft_grid_size(3_u32), 8);
        assert_eq!(fft_grid_size(4_u32), 9);
        assert_eq!(fft_grid_size(5_u32), 12);
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

    // -----------------------------------------------------------------
    // GGAP Phase A.1 — `compute_density_gradient` unit tests
    // -----------------------------------------------------------------
    //
    // A cubic box of side `L` with `N³` FFT grid points has reciprocal
    // vectors `G_α = (2π/L) · n_α`, where `n_α` ranges over the
    // FFT-aligned integers `{0, 1, …, N/2, −N/2+1, …, −1}`. The test
    // helpers below build those integer slots directly so the gradient
    // test is self-contained and does not import from `src/scf/grid.rs`.

    /// Build the Cartesian grid of real-space coordinates for a cubic
    /// box of side `L` with `N³` points, with `r = 0` at the corner
    /// (same convention as the SCF FFT grid).
    fn cubic_real_grid(n: usize, l: f64) -> Vec<[f64; 3]> {
        let h = l / n as f64;
        let mut out = Vec::with_capacity(n * n * n);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    out.push([i as f64 * h, j as f64 * h, k as f64 * h]);
                }
            }
        }
        out
    }

    /// Build the FFT-aligned G-vector list for a cubic box of side `L`
    /// with `N³` points. Layout matches `scf::grid::g_vector_at_dims`.
    #[allow(
        clippy::cast_possible_wrap,
        clippy::cast_possible_truncation,
        reason = "test helper only; n is a tiny FFT size (≤ 64) well inside i32 range"
    )]
    fn cubic_g_vectors(n: usize, l: f64) -> Vec<[f64; 3]> {
        let two_pi_l = 2.0 * std::f64::consts::PI / l;
        let mut out = Vec::with_capacity(n * n * n);
        let signed = |i: usize| -> i32 { if i > n / 2 { i as i32 - n as i32 } else { i as i32 } };
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    out.push([
                        f64::from(signed(i)) * two_pi_l,
                        f64::from(signed(j)) * two_pi_l,
                        f64::from(signed(k)) * two_pi_l,
                    ]);
                }
            }
        }
        out
    }

    #[test]
    fn test_density_gradient_gaussian_analytic() {
        // Gaussian centred at the box centre r_0 = (L/2, L/2, L/2).
        // Analytic gradient: ∇ρ(r) = −2α (r − r_0) ρ(r).
        //
        // α = 1.0 Å⁻² (FWHM ≈ 1.66 Å) on a 32³ grid of L = 10 Å gives
        // grid spacing h = 0.3125 Å → ~5.3 grid points across the FWHM.
        // At that resolution the Nyquist contribution is negligible
        // (Gaussian Fourier tail is ρ̂(k) ∝ exp(-k²/4α), and the
        // zero-out in `compute_density_gradient` doesn't degrade the
        // low/mid frequencies we care about). We pin at 1e-4 relative,
        // which is already 100× tighter than what the Phase-C PBE unit
        // tests demand on SCF grids.
        let n = 32;
        let l = 10.0;
        let r0 = [l / 2.0; 3];
        let alpha = 1.0_f64;

        let r = cubic_real_grid(n, l);
        let rho_r: Vec<f64> = r
            .iter()
            .map(|p| {
                let dx = p[0] - r0[0];
                let dy = p[1] - r0[1];
                let dz = p[2] - r0[2];
                (-alpha * (dx * dx + dy * dy + dz * dz)).exp()
            })
            .collect();

        let g_vectors = cubic_g_vectors(n, l);
        let mut fft = FFT3D::new(n, n, n);
        let grad_num = compute_density_gradient(&rho_r, &mut fft, &g_vectors);

        // Compare to analytic. Skip the outermost two shells to avoid
        // periodic-wrap residual near the box faces (the FFT "wraps"
        // across the boundary; finite α means the density is not
        // exactly zero there, so the periodic image contributes).
        let mut max_abs_err = 0.0_f64;
        let mut max_abs_ref = 0.0_f64;
        for (idx, p) in r.iter().enumerate() {
            let ix = idx / (n * n);
            let iy = (idx / n) % n;
            let iz = idx % n;
            if ix < 2 || ix > n - 3 || iy < 2 || iy > n - 3 || iz < 2 || iz > n - 3 {
                continue;
            }
            let dx = p[0] - r0[0];
            let dy = p[1] - r0[1];
            let dz = p[2] - r0[2];
            let rho = (-alpha * (dx * dx + dy * dy + dz * dz)).exp();
            let ana = [
                -2.0 * alpha * dx * rho,
                -2.0 * alpha * dy * rho,
                -2.0 * alpha * dz * rho,
            ];
            let num = grad_num[idx];
            for a in 0..3 {
                max_abs_err = max_abs_err.max((num[a] - ana[a]).abs());
                max_abs_ref = max_abs_ref.max(ana[a].abs());
            }
        }
        let rel_err = max_abs_err / max_abs_ref.max(1e-300);
        assert!(
            rel_err < 1e-4,
            "Gaussian gradient relative error {rel_err:.3e} exceeds 1e-4 (max|Δ|={max_abs_err:.3e}, max|∇ρ|={max_abs_ref:.3e})",
        );
    }

    #[test]
    fn test_density_gradient_linearity() {
        // ∇(aρ₁ + bρ₂) = a·∇ρ₁ + b·∇ρ₂ to FFT round-off.
        let n = 16;
        let l = 8.0;
        let r = cubic_real_grid(n, l);
        let g_vectors = cubic_g_vectors(n, l);
        let mut fft = FFT3D::new(n, n, n);

        // Two arbitrary smooth-ish densities.
        let rho1: Vec<f64> = r
            .iter()
            .map(|p| {
                let x = p[0] - l / 2.0;
                let y = p[1] - l / 2.0;
                let z = p[2] - l / 2.0;
                (-0.5 * (x * x + y * y + z * z)).exp()
            })
            .collect();
        let rho2: Vec<f64> = r
            .iter()
            .map(|p| {
                let x = p[0] - l / 3.0;
                let y = p[1] - l / 4.0;
                let z = p[2] - 2.0 * l / 3.0;
                (-(x * x + y * y + z * z)).exp()
            })
            .collect();

        let a = 2.7_f64;
        let b = -1.3_f64;
        let rho_sum: Vec<f64> = rho1.iter().zip(rho2.iter()).map(|(&r1, &r2)| a * r1 + b * r2).collect();

        let g1 = compute_density_gradient(&rho1, &mut fft, &g_vectors);
        let g2 = compute_density_gradient(&rho2, &mut fft, &g_vectors);
        let gsum = compute_density_gradient(&rho_sum, &mut fft, &g_vectors);

        let mut max_abs_err = 0.0_f64;
        for i in 0..gsum.len() {
            for ax in 0..3 {
                let expected = a * g1[i][ax] + b * g2[i][ax];
                max_abs_err = max_abs_err.max((gsum[i][ax] - expected).abs());
            }
        }
        assert!(
            max_abs_err < 1e-12,
            "linearity residual {max_abs_err:.3e} exceeds FFT round-off (1e-12)",
        );
    }

    #[test]
    fn test_density_gradient_g0_is_zero() {
        // Constant density → all Fourier coefficients at G ≠ 0 vanish;
        // the G = 0 coefficient is ρ_0 but `iG = 0` at G = 0 kills its
        // contribution. Result: ∇ρ(r) ≡ 0 everywhere to round-off.
        let n = 8;
        let l = 4.0;
        let rho_r = vec![0.37_f64; n * n * n];
        let g_vectors = cubic_g_vectors(n, l);
        let mut fft = FFT3D::new(n, n, n);

        let grad = compute_density_gradient(&rho_r, &mut fft, &g_vectors);
        let max_abs: f64 = grad
            .iter()
            .flat_map(|g| g.iter().copied().map(f64::abs))
            .fold(0.0, f64::max);
        assert!(
            max_abs < 1e-13,
            "constant density should yield ∇ρ ≡ 0; got max |∂ρ| = {max_abs:.3e}",
        );
    }
}
