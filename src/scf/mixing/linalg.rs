//! Small dense linear solver used by the DIIS/Broyden mixers.
//!
//! Both Anderson (DIIS) and modified Broyden mixing need to solve a small
//! (typically 2–8 × 2–8) linear system every iteration to extract the
//! mixing coefficients. Gauss elimination with partial pivoting is
//! sufficient at these sizes; a larger solver (faer LU) would dominate in
//! setup cost. See MODR's "Flagged for follow-up" for a future swap.

/// Solve A x = b for small systems via Gauss elimination with partial pivoting.
pub(super) fn solve_linear_system(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
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

    for col in 0..n {
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
            log::warn!(
                "Anderson mixer: singular overlap matrix (pivot={pivot:.2e}), \
                 falling back to uniform coefficients"
            );
            return vec![1.0 / (n + 1) as f64; n];
        }

        for row in col + 1..n {
            let factor = aug[row * (n + 1) + col] / pivot;
            for j in col..=n {
                aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }

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
    fn test_solve_linear_system() {
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let b = vec![5.0, 7.0];
        let x = solve_linear_system(&a, &b, 2);
        assert!((x[0] - 1.6).abs() < 1e-10);
        assert!((x[1] - 1.8).abs() < 1e-10);
    }
}
