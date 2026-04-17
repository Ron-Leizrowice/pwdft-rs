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

        let has_inversion = ops.iter().any(|op| {
            op.rotation == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
        });

        Self {
            operations: ops,
            n_ops,
            has_inversion,
            has_time_reversal: true, // default for non-magnetic
            tolerance,
        }
    }

    /// Verify group closure: every composition of two operations is in the set,
    /// and every operation has an inverse in the set.
    ///
    /// Panics with a descriptive message if the group is not closed.
    pub fn verify_group_closure(&self) {
        for (i, a) in self.operations.iter().enumerate() {
            // Check inverse exists
            let a_inv = a.inverse();
            let has_inv = self.operations.iter().any(|b| b.approx_eq(&a_inv, self.tolerance));
            assert!(
                has_inv,
                "operation {i} has no inverse in the group: R={:?} τ={:?}",
                a.rotation, a.translation
            );

            // Check closure under composition
            for (j, b) in self.operations.iter().enumerate() {
                let ab = a.compose(b);
                let in_set = self.operations.iter().any(|c| c.approx_eq(&ab, self.tolerance));
                assert!(
                    in_set,
                    "composition of ops {i} and {j} not in group: R={:?} τ={:?}",
                    ab.rotation, ab.translation
                );
            }
        }
    }
}
