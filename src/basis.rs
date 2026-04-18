use nalgebra::Vector3;

use crate::{consts::HBAR2_OVER_2M, crystal::Lattice};

pub struct BasisSet {
    /// G-vectors in Cartesian coordinates (1/Å).
    pw: Vec<Vector3<f64>>,
    /// Integer Miller indices (n1, n2, n3) for each G-vector.
    ///
    /// Stored as `[i16; 3]` per TYPE-A. `i16`'s range `[-32 768, 32 767]`
    /// covers every physically reasonable `ecut · Ω` combination with
    /// three orders of magnitude to spare — reaching i16 overflow would
    /// require `ecut` beyond `10^7 Ry`. All arithmetic on Miller indices
    /// (rotation, FFT-index mapping, k+G) widens to `i32` at the read
    /// site to preserve overflow-free behavior.
    miller: Vec<[i16; 3]>,
    /// Plane-wave energy cutoff (eV).
    ecut: f64,
}

impl BasisSet {
    /// Construct the plane-wave basis set for a given energy cutoff.
    ///
    /// Includes all reciprocal lattice vectors G = n₁b₁ + n₂b₂ + n₃b₃
    /// satisfying the kinetic energy cutoff:
    ///   (ħ²/2m) |G|² ≤ E_cut
    ///
    /// where b_i are reciprocal lattice vectors (2π/V × a_j × a_k).
    /// The number of basis functions scales as N_pw ∝ E_cut^{3/2} × Ω.
    #[must_use]
    pub fn new(lattice: &Lattice, ecut: f64) -> Self {
        let recip = lattice.reciprocal();
        let g_max_sq = ecut / HBAR2_OVER_2M;

        let g_max = g_max_sq.sqrt();
        let (b1, b2, b3) = (recip.a, recip.b, recip.c);

        #[allow(
            clippy::cast_possible_truncation,
            reason = "n_i_max bounded by ecut and reciprocal lattice norms; exceeding i32::MAX would require ecut > 10^18 eV"
        )]
        let n1_max = (g_max / b1.norm()).ceil() as i32;
        #[allow(
            clippy::cast_possible_truncation,
            reason = "n_i_max bounded by ecut and reciprocal lattice norms; exceeding i32::MAX would require ecut > 10^18 eV"
        )]
        let n2_max = (g_max / b2.norm()).ceil() as i32;
        #[allow(
            clippy::cast_possible_truncation,
            reason = "n_i_max bounded by ecut and reciprocal lattice norms; exceeding i32::MAX would require ecut > 10^18 eV"
        )]
        let n3_max = (g_max / b3.norm()).ceil() as i32;

        let mut pw = Vec::new();
        let mut miller = Vec::new();

        // Miller indices are narrowed from i32 to i16 for storage; the
        // enumeration driver (n1/n2/n3) remains i32 to keep the loop
        // bounds and intermediate arithmetic overflow-safe for any
        // physically meaningful ecut. `i16::try_from` is infallible under
        // the proposal's argument (ecut > 10^7 Ry is non-physical); the
        // panic path is a guard against a caller who passes an absurdly
        // large lattice.
        let to_i16 = |n: i32| {
            i16::try_from(n).expect(
                "BasisSet::new: Miller index exceeds i16 range; ecut must be below 10^7 Ry",
            )
        };

        for n1 in -n1_max..=n1_max {
            for n2 in -n2_max..=n2_max {
                for n3 in -n3_max..=n3_max {
                    let g = n1 as f64 * b1 + n2 as f64 * b2 + n3 as f64 * b3;
                    if g.norm_squared() <= g_max_sq {
                        pw.push(g);
                        miller.push([to_i16(n1), to_i16(n2), to_i16(n3)]);
                    }
                }
            }
        }

        Self { pw, miller, ecut }
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.pw.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pw.is_empty()
    }

    #[must_use]
    pub fn ecut(&self) -> f64 {
        self.ecut
    }

    /// Access G-vectors (Cartesian, 1/Å).
    #[must_use]
    pub fn g_vectors(&self) -> &[Vector3<f64>] {
        &self.pw
    }

    /// Access integer Miller indices for each G-vector.
    ///
    /// Stored as `[i16; 3]` for cache density (TYPE-A). Callers that feed
    /// the values into `i32`-sized arithmetic (rotation, FFT index wrap,
    /// `miller_to_idx`) should widen with `i32::from(m[k])` — the
    /// widening is infallible and typically folds into the load.
    #[must_use]
    pub fn miller_indices(&self) -> &[[i16; 3]] {
        &self.miller
    }

    /// Linear lookup: given Miller indices, return the index into the basis.
    ///
    /// `O(n_pw)` scan. Retained as public API because it is used by a
    /// handful of unit and integration tests for G-vector pin assertions
    /// (`index_of(0, 0, 0)`, small-shell sanity probes); production SCF
    /// paths work in the opposite direction (`miller_to_idx` on an FFT
    /// grid), so the former `HashMap<(i32, i32, i32), usize>` cache was
    /// ~30 kB of dead memory at `n_pw = 725`. The linear scan's cost
    /// (~n_pw comparisons per call, a handful of calls per test run) is
    /// invisible next to SCF wall time.
    ///
    /// Accepts `i32` so out-of-range queries return `None` naturally
    /// (no narrowing panic at the boundary).
    #[must_use]
    pub fn index_of(&self, n1: i32, n2: i32, n3: i32) -> Option<usize> {
        // `try_into` folds out-of-i16-range queries into a `None` result,
        // matching the pre-TYPE-A `HashMap::get` semantics.
        let n1: i16 = n1.try_into().ok()?;
        let n2: i16 = n2.try_into().ok()?;
        let n3: i16 = n3.try_into().ok()?;
        self.miller.iter().position(|m| m == &[n1, n2, n3])
    }

    /// Kinetic energy (eV) for each G-vector: (ℏ²/2m)|G|².
    #[must_use]
    pub fn kinetic_energy(&self) -> Vec<f64> {
        self.pw
            .iter()
            .map(|g| HBAR2_OVER_2M * g.norm_squared())
            .collect()
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "ERR2 § Phase 0: in-src test modules are allowed to panic")]
mod tests {
    use super::*;

    fn si_basis(ecut: f64) -> BasisSet {
        let a = 5.431;
        let si = Lattice::new(
            a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        );
        BasisSet::new(&si, ecut)
    }

    #[test]
    fn test_si_basis_count() {
        let basis = si_basis(200.0);
        assert_eq!(
            basis.len(),
            259,
            "expected 259 G-vectors at 200 eV, got {}",
            basis.len()
        );
    }

    #[test]
    fn test_kinetic_energy_contains_zero() {
        let basis = si_basis(200.0);
        let ke = basis.kinetic_energy();
        let min_ke = ke.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            min_ke.abs() < 1e-12,
            "G=0 should have zero kinetic energy, got {min_ke}"
        );
    }

    #[test]
    fn test_kinetic_energy_within_cutoff() {
        let ecut = 200.0;
        let basis = si_basis(ecut);
        let ke = basis.kinetic_energy();
        let max_ke = ke.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_ke <= ecut + 1e-10, "max KE {max_ke} exceeds cutoff {ecut}");
    }

    #[test]
    fn test_index_of_gamma() {
        let basis = si_basis(200.0);
        // G=0 must exist
        let idx = basis.index_of(0, 0, 0);
        assert!(idx.is_some(), "G=0 should be in the basis");
        let g = basis.g_vectors()[idx.unwrap()];
        assert!(g.norm() < 1e-12, "G=0 should be the zero vector");
    }

    #[test]
    fn test_index_of_out_of_range() {
        let basis = si_basis(200.0);
        // Very large indices should not be in the basis
        assert!(basis.index_of(100, 100, 100).is_none());
        // Values outside i16 range also return None, not panic
        assert!(basis.index_of(100_000, 0, 0).is_none());
    }

    #[test]
    fn test_miller_indices_consistency() {
        let basis = si_basis(200.0);
        assert_eq!(basis.g_vectors().len(), basis.miller_indices().len());
        // Every miller index should roundtrip through the linear scan
        for (i, &[n1, n2, n3]) in basis.miller_indices().iter().enumerate() {
            assert_eq!(
                basis.index_of(i32::from(n1), i32::from(n2), i32::from(n3)),
                Some(i)
            );
        }
    }
}
