//! Numerical integration utilities for radial integrals.
//!
//! These routines operate on logarithmic radial grids with `rab[i]` weights,
//! matching the conventions used by Quantum ESPRESSO's pseudopotential library.

/// Simpson's 1/3 rule integration on a radial grid with `rab` weights.
///
/// Computes: `sum_i c_i * func[i] * rab[i]`
///
/// where c_i alternate 2/3, 4/3, with endpoints at 1/3. Matches QE's
/// `simpsn.f90` exactly, including the even-mesh boundary correction
/// from DFTK.
///
/// For odd mesh (standard composite Simpson):
///   weights = [1, 4, 2, 4, 2, ..., 4, 1] / 3
///
/// For even mesh (Simpson + boundary correction):
///   Same as odd up to index n-4, then corrected at the boundary:
///   ..., 2, 15/12, 1, 5/12  (i.e., subtract 1/4 from c_{n-3},
///   keep c_{n-2}, add 5/4 to c_{n-1})
///
/// # Panics
///
/// Panics if `func` and `rab` have different lengths.
pub fn simpson_integrate(func: &[f64], rab: &[f64]) -> f64 {
    let n = func.len();
    assert_eq!(n, rab.len(), "func and rab must have the same length");

    if n < 3 {
        // Fallback to trapezoidal for tiny grids (shouldn't happen in practice)
        return func.iter().zip(rab).map(|(&f, &dr)| f * dr).sum();
    }

    // Interior sum: i = 1..n-2 (0-based), matching Fortran i = 2..mesh-1 (1-based)
    // Weight: 4 for 0-based odd index (Fortran even), 2 for 0-based even index (Fortran odd)
    let mut sum = 0.0;
    for i in 1..n - 1 {
        let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * func[i] * rab[i];
    }

    if n % 2 == 1 {
        // Odd mesh: standard Simpson's rule
        (sum + func[0] * rab[0] + func[n - 1] * rab[n - 1]) / 3.0
    } else {
        // Even mesh: boundary correction (matches QE/DFTK formula)
        (sum + func[0] * rab[0]
            - func[n - 3] * rab[n - 3] * 0.25
            + func[n - 2] * rab[n - 2]
            + func[n - 1] * rab[n - 1] * 1.25)
            / 3.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

    #[test]
    fn test_simpson_exact_for_cubic() {
        // Simpson's rule is exact for polynomials up to degree 3.
        // Integrate f(x) = x^3 from 0 to 1 on a uniform grid.
        // Exact answer: 1/4 = 0.25
        let n = 101; // odd
        let h = 1.0 / (n - 1) as f64;
        let func: Vec<f64> = (0..n).map(|i| {
            let x = i as f64 * h;
            x * x * x
        }).collect();
        let rab: Vec<f64> = vec![h; n];

        let result = simpson_integrate(&func, &rab);
        assert!(
            relative_eq!(result, 0.25, epsilon = 1e-12),
            "Simpson integral of x^3 = {result}, expected 0.25"
        );
    }

    #[test]
    fn test_simpson_exact_for_quadratic() {
        // Integrate f(x) = x^2 from 0 to 2. Exact answer: 8/3
        let n = 51; // odd
        let h = 2.0 / (n - 1) as f64;
        let func: Vec<f64> = (0..n).map(|i| {
            let x = i as f64 * h;
            x * x
        }).collect();
        let rab: Vec<f64> = vec![h; n];

        let result = simpson_integrate(&func, &rab);
        assert!(
            relative_eq!(result, 8.0 / 3.0, epsilon = 1e-12),
            "Simpson integral of x^2 = {result}, expected {}", 8.0 / 3.0
        );
    }

    #[test]
    fn test_simpson_even_mesh() {
        // Even mesh: integrate f(x) = x^2 from 0 to 2.
        let n = 50; // even
        let h = 2.0 / (n - 1) as f64;
        let func: Vec<f64> = (0..n).map(|i| {
            let x = i as f64 * h;
            x * x
        }).collect();
        let rab: Vec<f64> = vec![h; n];

        let result = simpson_integrate(&func, &rab);
        // Even mesh correction is less accurate but should still be close
        assert!(
            (result - 8.0 / 3.0).abs() < 1e-4,
            "Simpson (even mesh) integral of x^2 = {result}, expected ~{}", 8.0 / 3.0
        );
    }

    #[test]
    fn test_simpson_vs_trapezoidal_gaussian() {
        // On a radial grid, Simpson should be much more accurate than trapezoidal
        // for a smooth Gaussian integrand.
        // Integrate exp(-r^2) from 0 to 5. Exact: sqrt(pi)/2 * erf(5) ~ 0.886227
        let n = 101;
        let h = 5.0 / (n - 1) as f64;
        let func: Vec<f64> = (0..n).map(|i| {
            let r = i as f64 * h;
            (-r * r).exp()
        }).collect();
        let rab: Vec<f64> = vec![h; n];

        let exact = std::f64::consts::PI.sqrt() / 2.0; // 0.886226925...
        let simp = simpson_integrate(&func, &rab);
        let trap: f64 = func.iter().zip(rab.iter()).map(|(&f, &dr)| f * dr).sum();

        let simp_err = (simp - exact).abs();
        let trap_err = (trap - exact).abs();

        // Simpson should be orders of magnitude more accurate
        assert!(
            simp_err < trap_err * 0.01,
            "Simpson error ({simp_err:.4e}) should be << trapezoidal error ({trap_err:.4e})"
        );
    }

    #[test]
    fn test_simpson_tiny_grid() {
        // n < 3: falls back to trapezoidal
        let func = vec![1.0, 2.0];
        let rab = vec![0.5, 0.5];
        let result = simpson_integrate(&func, &rab);
        assert!(relative_eq!(result, 1.5, epsilon = 1e-15));
    }
}
