//! Density mixing schemes for SCF convergence.
//!
//! Anderson (Pulay) mixing: stores a history of input/output density pairs
//! and finds the optimal linear combination.

use ndarray::{Array1, ArrayView1};

/// Anderson/Pulay density mixer.
pub struct AndersonMixer {
    beta: f64,
    max_history: usize,
    history_in: Vec<Array1<f64>>,
    history_res: Vec<Array1<f64>>,
}

impl AndersonMixer {
    pub fn new(beta: f64, max_history: usize, _n_grid: usize) -> Self {
        Self {
            beta,
            max_history,
            history_in: Vec::new(),
            history_res: Vec::new(),
        }
    }

    /// Mix input density with output density.
    ///
    /// `rho_in`: current input density.
    /// `rho_out`: density computed from KS eigenstates.
    ///
    /// Returns the new input density for the next iteration.
    pub fn mix(&mut self, rho_in: &[f64], rho_out: &[f64]) -> Vec<f64> {
        let rho_in_arr = ArrayView1::from(rho_in);
        let rho_out_arr = ArrayView1::from(rho_out);
        let residual = &rho_out_arr - &rho_in_arr;

        self.history_in.push(rho_in_arr.to_owned());
        self.history_res.push(residual.clone());

        // Trim history
        if self.history_in.len() > self.max_history {
            self.history_in.remove(0);
            self.history_res.remove(0);
        }

        let m = self.history_in.len();
        if m < 2 {
            // Simple linear mixing for first iteration
            let rho_new = &rho_in_arr + &(self.beta * &residual);
            return rho_new.to_vec();
        }

        // Anderson mixing: find coefficients that minimize |Σ α_i R_i|²
        // subject to Σ α_i = 1
        // Solve: A α = b where A_{ij} = <R_i - R_m | R_j - R_m>, b_i = -<R_i - R_m | R_m>
        let last = m - 1;
        let mm = m - 1; // number of equations

        let r_last = &self.history_res[last];

        // Build the system using ndarray dot products
        let mut a_mat = vec![0.0; mm * mm];
        let mut b_vec = vec![0.0; mm];

        // Precompute delta residuals
        let dr: Vec<Array1<f64>> = (0..mm)
            .map(|i| &self.history_res[i] - r_last)
            .collect();

        for i in 0..mm {
            b_vec[i] = -dr[i].dot(r_last);
            for j in 0..mm {
                a_mat[i * mm + j] = dr[i].dot(&dr[j]);
            }
        }

        // Solve via simple Gauss elimination (mm is small, typically 2-8)
        let alpha_prev = solve_linear_system(&a_mat, &b_vec, mm);
        let alpha_last = 1.0 - alpha_prev.iter().sum::<f64>();

        // Construct mixed density: Σ α_i (ρ_in_i + β R_i)
        let mut rho_new = alpha_last * (&self.history_in[last] + &(self.beta * r_last));
        for j in 0..mm {
            rho_new += &(alpha_prev[j] * (&self.history_in[j] + &(self.beta * &self.history_res[j])));
        }

        rho_new.to_vec()
    }
}

/// Solve A x = b for small systems via Gauss elimination with partial pivoting.
fn solve_linear_system(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }

    // Forward elimination
    for col in 0..n {
        // Partial pivoting
        let mut max_row = col;
        let mut max_val = aug[col * (n + 1) + col].abs();
        for row in col + 1..n {
            let val = aug[row * (n + 1) + col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }
        if max_row != col {
            for j in 0..=n {
                aug.swap(col * (n + 1) + j, max_row * (n + 1) + j);
            }
        }

        let pivot = aug[col * (n + 1) + col];
        if pivot.abs() < 1e-15 {
            // Singular — return equal weights
            return vec![1.0 / (n + 1) as f64; n];
        }

        for row in col + 1..n {
            let factor = aug[row * (n + 1) + col] / pivot;
            for j in col..=n {
                aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i * (n + 1) + n];
        for j in i + 1..n {
            sum -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] = sum / aug[i * (n + 1) + i];
    }

    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_mixing() {
        let mut mixer = AndersonMixer::new(0.3, 4, 10);
        let rho_in = vec![1.0; 10];
        let rho_out = vec![2.0; 10];
        let result = mixer.mix(&rho_in, &rho_out);
        // First iteration: linear mixing ρ_new = ρ_in + β(ρ_out - ρ_in) = 1 + 0.3 = 1.3
        for &v in &result {
            assert!(
                (v - 1.3).abs() < 1e-10,
                "expected 1.3, got {v}"
            );
        }
    }

    #[test]
    fn test_solve_linear_system() {
        // 2x + y = 5, x + 3y = 7 → x = 1.6, y = 1.8
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let b = vec![5.0, 7.0];
        let x = solve_linear_system(&a, &b, 2);
        assert!((x[0] - 1.6).abs() < 1e-10);
        assert!((x[1] - 1.8).abs() < 1e-10);
    }
}
