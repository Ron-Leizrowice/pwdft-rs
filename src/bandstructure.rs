//! Non-self-consistent band structure along a k-path.
//!
//! Diagonalizes a fixed Hamiltonian at each k-point on a supplied path
//! (typically a Γ-X-W-L-Γ high-symmetry circuit produced by
//! [`crate::kpoints::high_symmetry_path`]) and collects the lowest
//! `n_bands` eigenvalues per point, together with the cumulative
//! path-distance coordinate used as the x-axis in band plots.
//!
//! The current entry point [`compute_band_structure`] runs the
//! free-electron (kinetic-only) spectrum — a diagnostic mode useful for
//! validating the basis and k-path. Self-consistent bands come from the
//! [`crate::scf`] driver; a post-SCF hook that reuses the converged
//! potential is tracked separately.

use crate::{
    basis::BasisSet,
    eigensolver::dense,
    error::Result,
    hamiltonian,
    kpoints::KPoint,
};

/// Result of a band structure calculation.
pub struct BandStructure {
    /// Cumulative k-path distance for each k-point.
    pub distances: Vec<f64>,
    /// Labels at high-symmetry points: (distance, label).
    pub labels: Vec<(f64, String)>,
    /// Eigenvalues `[k_index][band_index]` in eV.
    pub eigenvalues: Vec<Vec<f64>>,
}

/// Compute the free-electron (kinetic-only) band structure along a k-path.
///
/// For each k-point, builds the kinetic-energy Hamiltonian
/// `H_{G,G'}(k) = δ_{GG'} · (ℏ²/2m) |k + G|²` and diagonalizes it,
/// keeping the lowest `n_bands` eigenvalues. This is a diagnostic / validation
/// mode — the nearly-free-electron spectrum that a converged pseudopotential
/// calculation should reproduce at high `|k + G|`. Self-consistent bands come
/// from the SCF driver, not this routine.
///
/// # Errors
/// Returns `PwdftError::Eigensolver` if any eigendecomposition fails.
pub fn compute_band_structure(
    basis: &BasisSet,
    kpoints: &[KPoint],
    distances: &[f64],
    n_bands: usize,
) -> Result<BandStructure> {
    let mut eigenvalues = Vec::with_capacity(kpoints.len());
    let mut labels = Vec::new();

    for (i, kp) in kpoints.iter().enumerate() {
        let h = hamiltonian::build_kinetic(basis, &kp.k);
        let result = dense::diagonalize_lowest(&h, n_bands)?;
        eigenvalues.push(result.eigenvalues);

        if let Some(ref label) = kp.label {
            labels.push((distances[i], label.clone()));
        }
    }

    Ok(BandStructure {
        distances: distances.to_vec(),
        labels,
        eigenvalues,
    })
}

impl BandStructure {
    /// Write band structure as TSV suitable for plotting.
    ///
    /// Format: distance  band_0  band_1  band_2  ...
    ///
    /// # Errors
    ///
    /// Returns any `std::io::Error` produced by the underlying writer (e.g.
    /// disk full, broken pipe). No error is synthesized internally — the
    /// function only forwards I/O failures from `write!`/`writeln!`.
    pub fn write_tsv<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        if self.eigenvalues.is_empty() {
            return Ok(());
        }
        let n_bands = self.eigenvalues[0].len();

        // Header
        write!(writer, "# k_distance")?;
        for i in 0..n_bands {
            write!(writer, "\tband_{i}")?;
        }
        writeln!(writer)?;

        // High-symmetry point labels as comments
        for (dist, label) in &self.labels {
            writeln!(writer, "# {label} at k_distance = {dist:.6}")?;
        }

        // Data
        for (i, dist) in self.distances.iter().enumerate() {
            write!(writer, "{dist:.6}")?;
            for ev in &self.eigenvalues[i] {
                write!(writer, "\t{ev:.6}")?;
            }
            writeln!(writer)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use crate::{
        crystal::Lattice,
        consts::HBAR2_OVER_2M,
    };
    use approx::relative_eq;

    fn si_lattice() -> Lattice {
        let a = 5.431;
        Lattice::new(
            a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        )
    }

    #[test]
    fn test_free_electron_gamma() {
        let lattice = si_lattice();
        let basis = BasisSet::new(&lattice, 200.0);

        let gamma = KPoint {
            k: Vector3::zeros(),
            weight: 1.0,
            label: Some("Γ".into()),
        };

        let bs = compute_band_structure(&basis, &[gamma], &[0.0], 10).unwrap();

        // At Γ, lowest eigenvalue should be 0 (G=0, |k+G|=0)
        assert!(
            relative_eq!(bs.eigenvalues[0][0], 0.0, epsilon = 1e-10),
            "lowest band at Γ should be 0 eV, got {}",
            bs.eigenvalues[0][0]
        );

        // Eigenvalues should be non-decreasing
        for i in 1..bs.eigenvalues[0].len() {
            assert!(
                bs.eigenvalues[0][i] >= bs.eigenvalues[0][i - 1] - 1e-10,
                "eigenvalues not sorted at Γ"
            );
        }
    }

    #[test]
    fn test_free_electron_parabolic() {
        let lattice = si_lattice();
        let basis = BasisSet::new(&lattice, 200.0);

        // At a generic k, the free-electron energy for band with G-vector G_i is:
        // E = (ℏ²/2m)|k + G_i|²
        // The lowest band at any k should be (ℏ²/2m)|k|² (from G=0)
        let k = Vector3::new(0.1, 0.05, 0.02);
        let kp = KPoint {
            k,
            weight: 1.0,
            label: None,
        };

        let bs = compute_band_structure(&basis, &[kp], &[0.0], 5).unwrap();
        let expected_lowest = HBAR2_OVER_2M * k.norm_squared();

        assert!(
            relative_eq!(bs.eigenvalues[0][0], expected_lowest, epsilon = 1e-10),
            "lowest band = {}, expected {}",
            bs.eigenvalues[0][0],
            expected_lowest
        );
    }
}
