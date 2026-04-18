//! Crystal symmetry detection, k-point reduction, and density symmetrization.

pub mod density;
pub mod detect;
pub mod kpoints;
pub mod operations;

pub use operations::{SpaceGroupOp, SymmOp};

use crate::crystal::Crystal;

/// Complete symmetry information for a crystal.
#[derive(Debug, Clone)]
pub struct SymmetryInfo {
    /// All space group operations {R|τ}.
    pub operations: Vec<SpaceGroupOp>,
    /// Number of operations.
    pub n_ops: usize,
    /// Whether the crystal has spatial inversion symmetry.
    pub has_inversion: bool,
    /// Whether time-reversal symmetry applies (true for non-magnetic systems).
    pub has_time_reversal: bool,
    /// Tolerance used for symmetry detection (fractional coordinates).
    pub tolerance: f64,
}

impl SymmetryInfo {
    /// Detect all symmetry operations of a crystal.
    ///
    /// `tolerance`: maximum distance (in fractional coordinates) for atoms
    /// to be considered equivalent under a symmetry operation. Typical: 1e-5.
    #[must_use]
    pub fn from_crystal(crystal: &Crystal, tolerance: f64) -> Self {
        let ops = detect::find_symmetry_operations(crystal, tolerance);
        let n_ops = ops.len();

        let has_inversion = ops
            .iter()
            .any(|op| op.rotation == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]);

        Self {
            operations: ops,
            n_ops,
            has_inversion,
            has_time_reversal: true, // default for non-magnetic
            tolerance,
        }
    }

    /// Construct a trivial symmetry group containing only the identity
    /// operation, with time-reversal disabled.
    ///
    /// Used when the user sets `symmetry.enabled = false` in the input.
    /// Every crystal has at least the identity — this makes the "no
    /// symmetry" setting expressible as a real (if minimal) group rather
    /// than as a missing value.
    ///
    /// With this group:
    /// - K-point reduction leaves every input k-point unchanged (each
    ///   orbit has size 1, so weights remain `1/n_total`).
    /// - Density symmetrization is a no-op (see
    ///   [`density::symmetrize_density`], which short-circuits when
    ///   `n_ops <= 1`).
    ///
    /// The resulting behavior is bit-identical to the legacy path that
    /// skipped symmetrization via `Option::None`.
    #[must_use]
    pub fn identity_only() -> Self {
        Self {
            operations: vec![SpaceGroupOp::identity()],
            n_ops: 1,
            has_inversion: false,
            has_time_reversal: false,
            tolerance: 1e-5,
        }
    }

    /// Whether this group is trivial in the sense that applying it to any
    /// quantity (density, k-point mesh, …) is guaranteed to be a no-op.
    ///
    /// Returns `true` iff the group contains **only** the identity operation
    /// **and** time-reversal symmetry is disabled. Both conditions are
    /// required: a P1 crystal (triclinic, single identity op) with
    /// time-reversal on still folds k ↔ −k under [`kpoints::reduce_kpoints`],
    /// so treating it as "trivial" would skip real symmetry work.
    ///
    /// Cheap predicate — an identity check on the single op plus two scalar
    /// comparisons. Intended for callers that want to short-circuit work
    /// which would genuinely be a no-op under this group; callers needing a
    /// weaker check (e.g. "is this the spatial-symmetry identity?") should
    /// inspect `n_ops` and `operations` directly rather than generalize
    /// this predicate.
    #[must_use]
    pub fn is_trivial(&self) -> bool {
        self.n_ops == 1
            && self.operations.len() == 1
            && !self.has_time_reversal
            && self.operations[0].is_identity(self.tolerance.max(1e-12))
    }

    /// Verify group closure: every composition of two operations is in the set,
    /// and every operation has an inverse in the set.
    ///
    /// Panics with a descriptive message if the group is not closed.
    pub fn verify_group_closure(&self) {
        for (i, a) in self.operations.iter().enumerate() {
            // Check inverse exists
            let a_inv = a.inverse();
            let has_inv = self
                .operations
                .iter()
                .any(|b| b.approx_eq(&a_inv, self.tolerance));
            assert!(
                has_inv,
                "operation {i} has no inverse in the group: R={:?} τ={:?}",
                a.rotation, a.translation
            );

            // Check closure under composition
            for (j, b) in self.operations.iter().enumerate() {
                let ab = a.compose(b);
                let in_set = self
                    .operations
                    .iter()
                    .any(|c| c.approx_eq(&ab, self.tolerance));
                assert!(
                    in_set,
                    "composition of ops {i} and {j} not in group: R={:?} τ={:?}",
                    ab.rotation, ab.translation
                );
            }
        }
    }
}

#[cfg(test)]
// `symmetrize_with_identity_only_is_noop` exercises the deprecated legacy
// real-space `density::symmetrize_density` on purpose — pinning the
// n_ops ≤ 1 short-circuit so the SCF identity-only fallback stays
// bit-identical to the no-symmetrization path.
#[allow(deprecated)]
mod tests {
    use super::*;

    #[test]
    fn identity_only_has_single_identity_op() {
        let s = SymmetryInfo::identity_only();
        assert_eq!(s.n_ops, 1);
        assert_eq!(s.operations.len(), 1);
        assert!(s.operations[0].is_identity(1e-12));
        assert!(!s.has_inversion);
        assert!(!s.has_time_reversal);
    }

    #[test]
    fn identity_only_is_trivial() {
        let s = SymmetryInfo::identity_only();
        assert!(s.is_trivial());
    }

    #[test]
    fn identity_plus_time_reversal_is_not_trivial() {
        // A P1 crystal has only the identity space-group op but still admits
        // time-reversal symmetry (k ↔ −k), which folds the MP grid. Such a
        // group must NOT be reported as trivial — otherwise callers that
        // short-circuit on `is_trivial()` will skip real folding work. This
        // is the SOPT regression the predicate tightening guards against.
        let mut s = SymmetryInfo::identity_only();
        s.has_time_reversal = true;
        assert!(
            !s.is_trivial(),
            "identity-only with time-reversal must not be considered trivial: \
             k ↔ −k folding is still meaningful"
        );
    }

    #[test]
    fn symmetrize_with_identity_only_is_noop() {
        // Verify the "identity-only SymmetryInfo" is numerically indistinguishable
        // from the legacy "no symmetrization" path: the density must be left
        // bit-identical after a call to `symmetrize_density` because the
        // short-circuit `n_ops <= 1` applies.
        let s = SymmetryInfo::identity_only();
        let dims = [8, 8, 8];
        let n = dims[0] * dims[1] * dims[2];
        let rho_original: Vec<f64> = (0..n).map(|i| (i as f64 * 0.13).sin() + 1.0).collect();
        let mut rho = rho_original.clone();
        density::symmetrize_density(&mut rho, dims, &s);
        for (orig, sym) in rho_original.iter().zip(rho.iter()) {
            // Bit-identical, not just numerically close. The short-circuit
            // returns before touching `rho`, so no rounding occurs.
            assert_eq!(
                orig.to_bits(),
                sym.to_bits(),
                "identity-only symmetrization must leave density bit-identical"
            );
        }
    }

    #[test]
    fn reduce_kpoints_with_identity_only_preserves_full_grid() {
        // With only identity + no time-reversal, every orbit has size 1, so
        // `reduce_kpoints` must return the full input grid unchanged in both
        // k-points (Cartesian) and weights. This reproduces exactly the
        // "no symmetry" path that main.rs takes when `symmetry.enabled=false`.
        use crate::crystal::{Atom, Crystal, Lattice};
        use nalgebra::Vector3;

        let a = 5.431;
        let crystal = Crystal {
            lattice: Lattice::new(
                a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
            ),
            atoms: vec![
                Atom::new(14, [0.0, 0.0, 0.0]),
                Atom::new(14, [0.25, 0.25, 0.25]),
            ],
        };

        let grid = [4_u32, 4, 4];
        let full = crate::kpoints::monkhorst_pack(grid[0], grid[1], grid[2], &crystal.lattice);
        let s = SymmetryInfo::identity_only();
        let reduced = kpoints::reduce_kpoints(&full, grid, &s, &crystal.lattice);

        assert_eq!(
            reduced.len(),
            full.len(),
            "identity-only reduction must return the full grid (no folding)"
        );
        for (f, r) in full.iter().zip(reduced.iter()) {
            assert!(
                (f.k - r.k).norm() < 1e-12,
                "k-point position changed under identity-only reduction: {:?} vs {:?}",
                f.k,
                r.k
            );
            assert!(
                (f.weight - r.weight).abs() < 1e-12,
                "k-point weight changed: {} vs {}",
                f.weight,
                r.weight
            );
        }
        let w_sum: f64 = reduced.iter().map(|kp| kp.weight).sum();
        assert!(
            (w_sum - 1.0).abs() < 1e-12,
            "weights must sum to 1.0, got {w_sum}"
        );
    }
}
