//! Crystal symmetry detection algorithm.
//!
//! Finds all space group operations {R|τ} compatible with a given crystal structure.

use crate::crystal::Crystal;

use super::operations::{frac_distance, wrap_to_unit_cell, SpaceGroupOp, SymmOp};

/// Find all symmetry operations of a crystal.
///
/// Algorithm:
/// 1. Compute metric tensor M = Lᵀ L
/// 2. Enumerate all 3×3 integer matrices R with det(R) = ±1 and Rᵀ M R = M
/// 3. For each R, find fractional translation τ such that {R|τ} maps all atoms
///    to equivalent atoms (same species at equivalent positions mod 1)
#[must_use]
pub fn find_symmetry_operations(crystal: &Crystal, tolerance: f64) -> Vec<SpaceGroupOp> {
    let lattice_matrix = crystal.lattice.matrix();
    let metric = lattice_matrix.transpose() * lattice_matrix;

    // Step 1: Find all point group rotations preserving the metric tensor
    let rotations = find_metric_preserving_rotations(&metric, tolerance);

    // Step 2: For each rotation, find compatible translations via atom mapping.
    //
    // The internal enumeration works in `i32` (unchanged) but the public
    // `SpaceGroupOp` now stores rotations as `i8`. `narrow_rotation` performs
    // the `i32 → i8` narrow via `try_from` + `expect`: crystallographic
    // rotations never exceed `|R_ij| ≤ 3`, so the narrowing is infallible
    // in practice; a failure would indicate a bug in
    // `find_metric_preserving_rotations`.
    let mut ops = Vec::new();
    for rot in &rotations {
        if let Some(tau) = find_translation(crystal, rot, tolerance) {
            ops.push(SpaceGroupOp::new(narrow_rotation(rot), tau));
        }
    }

    // Sanity: identity must always be present
    assert!(
        ops.iter().any(|op| op.is_identity(tolerance)),
        "identity operation not found — bug in symmetry detection"
    );

    ops
}

/// Find all 3×3 integer matrices R with det(R) = ±1 and Rᵀ M R = M.
///
/// For each column j of R, enumerate integer vectors c satisfying cᵀ M c = M_{jj}.
/// Then filter valid 3-column combinations by the full metric condition and determinant.
fn find_metric_preserving_rotations(
    metric: &nalgebra::Matrix3<f64>,
    tolerance: f64,
) -> Vec<[[i32; 3]; 3]> {
    // For each column, find candidate integer vectors
    let candidates: Vec<Vec<[i32; 3]>> = (0..3)
        .map(|j| find_candidate_vectors(metric, metric[(j, j)], tolerance))
        .collect();

    let mut rotations = Vec::new();

    for c0 in &candidates[0] {
        // Check c0ᵀ M c0 = M_00 (already guaranteed by candidate generation)
        for c1 in &candidates[1] {
            // Early exit: check c0ᵀ M c1 = M_01
            if (dot_metric(c0, c1, metric) - metric[(0, 1)]).abs() > tolerance {
                continue;
            }
            for c2 in &candidates[2] {
                // Check determinant first (cheapest full-matrix check)
                let det = c0[0] * (c1[1] * c2[2] - c1[2] * c2[1])
                    - c0[1] * (c1[0] * c2[2] - c1[2] * c2[0])
                    + c0[2] * (c1[0] * c2[1] - c1[1] * c2[0]);
                if det != 1 && det != -1 {
                    continue;
                }

                // Check remaining metric conditions
                if (dot_metric(c0, c2, metric) - metric[(0, 2)]).abs() > tolerance {
                    continue;
                }
                if (dot_metric(c1, c2, metric) - metric[(1, 2)]).abs() > tolerance {
                    continue;
                }

                // Valid rotation found
                rotations.push([
                    [c0[0], c1[0], c2[0]],
                    [c0[1], c1[1], c2[1]],
                    [c0[2], c1[2], c2[2]],
                ]);
            }
        }
    }

    rotations
}

/// Find all integer vectors c such that cᵀ M c ≈ target_norm_sq.
fn find_candidate_vectors(
    metric: &nalgebra::Matrix3<f64>,
    target_norm_sq: f64,
    tolerance: f64,
) -> Vec<[i32; 3]> {
    // Search radius: |c_i| ≤ sqrt(target / M_ii) + 1
    #[allow(
        clippy::cast_possible_truncation,
        reason = "target_norm_sq comes from the metric matrix of the unit cell (Å²); sqrt is a handful of integer units and trivially fits in i32"
    )]
    let max_range: Vec<i32> = (0..3)
        .map(|i| ((target_norm_sq / metric[(i, i)]).sqrt() + 1.5) as i32)
        .collect();

    let mut candidates = Vec::new();
    for n0 in -max_range[0]..=max_range[0] {
        for n1 in -max_range[1]..=max_range[1] {
            for n2 in -max_range[2]..=max_range[2] {
                let c = [n0, n1, n2];
                let norm_sq = dot_metric(&c, &c, metric);
                if (norm_sq - target_norm_sq).abs() < tolerance {
                    candidates.push(c);
                }
            }
        }
    }
    candidates
}

/// Compute cᵀ M d for integer vectors c, d and real metric M.
fn dot_metric(c: &[i32; 3], d: &[i32; 3], metric: &nalgebra::Matrix3<f64>) -> f64 {
    let mut sum = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            sum += c[i] as f64 * metric[(i, j)] * d[j] as f64;
        }
    }
    sum
}

/// For a given rotation R, find a fractional translation τ such that {R|τ}
/// maps every atom to an equivalent atom (same species, same position mod 1).
///
/// Returns None if no valid translation exists.
fn find_translation(
    crystal: &Crystal,
    rotation: &[[i32; 3]; 3],
    tolerance: f64,
) -> Option<[f64; 3]> {
    if crystal.atoms.is_empty() {
        return Some([0.0, 0.0, 0.0]);
    }

    let ref_atom = &crystal.atoms[0];
    let ref_z = ref_atom.z;
    let ref_pos = ref_atom.position;

    // Rotated position of reference atom
    let r_pos = SymmOp { rotation: narrow_rotation(rotation) }.apply(&ref_pos);

    // For each atom of the same species, compute candidate translation
    for target in crystal.atoms.iter().filter(|a| a.z == ref_z) {
        let tau = [
            target.position[0] - r_pos[0],
            target.position[1] - r_pos[1],
            target.position[2] - r_pos[2],
        ];
        let tau = wrap_to_unit_cell(tau);

        // Verify this translation works for ALL atoms
        if all_atoms_map(crystal, rotation, &tau, tolerance) {
            return Some(tau);
        }
    }

    None
}

/// Narrow a `[[i32; 3]; 3]` rotation matrix to the `[[i8; 3]; 3]` storage
/// used by `SymmOp` / `SpaceGroupOp`. Panics if any entry overflows `i8`,
/// which would indicate a non-crystallographic rotation escaping from
/// `find_metric_preserving_rotations`.
fn narrow_rotation(r: &[[i32; 3]; 3]) -> [[i8; 3]; 3] {
    #[expect(
        clippy::expect_used,
        reason = "BUG: TYPE-A narrowing; rotations returned by find_metric_preserving_rotations have crystallographic entries in {-2..=2}, well within i8 range. A panic here indicates a non-crystallographic rotation escaped that routine."
    )]
    let to_i8 = |v: i32| {
        i8::try_from(v).expect("symmetry::detect: rotation entry exceeds i8 range")
    };
    [
        [to_i8(r[0][0]), to_i8(r[0][1]), to_i8(r[0][2])],
        [to_i8(r[1][0]), to_i8(r[1][1]), to_i8(r[1][2])],
        [to_i8(r[2][0]), to_i8(r[2][1]), to_i8(r[2][2])],
    ]
}

/// Check if {R|τ} maps every atom to an equivalent atom.
fn all_atoms_map(
    crystal: &Crystal,
    rotation: &[[i32; 3]; 3],
    tau: &[f64; 3],
    tolerance: f64,
) -> bool {
    let r_i8 = narrow_rotation(rotation);
    for atom in &crystal.atoms {
        let r_pos = SymmOp { rotation: r_i8 }.apply(&atom.position);
        let mapped = wrap_to_unit_cell([r_pos[0] + tau[0], r_pos[1] + tau[1], r_pos[2] + tau[2]]);

        let found = crystal.atoms.iter().any(|other| {
            other.z == atom.z && frac_distance(&mapped, &other.position) < tolerance
        });

        if !found {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::{Atom, Crystal, Lattice};
    use nalgebra::Vector3;

    fn si_fcc() -> Crystal {
        let a = 5.431;
        Crystal {
            lattice: Lattice::new(
                a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
            ),
            atoms: vec![
                Atom::new(14, [0.0, 0.0, 0.0]),
                Atom::new(14, [0.25, 0.25, 0.25]),
            ],
        }
    }

    fn bcc_fe() -> Crystal {
        let a = 2.87;
        Crystal {
            lattice: Lattice::new(
                a / 2.0 * Vector3::new(-1.0, 1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, -1.0, 1.0),
                a / 2.0 * Vector3::new(1.0, 1.0, -1.0),
            ),
            atoms: vec![Atom::new(26, [0.0, 0.0, 0.0])],
        }
    }

    fn triclinic_p1() -> Crystal {
        // Two atoms of DIFFERENT species at general positions → no symmetry possible
        Crystal {
            lattice: Lattice::new(
                Vector3::new(3.1, 0.2, 0.1),
                Vector3::new(0.3, 4.2, -0.1),
                Vector3::new(-0.1, 0.15, 5.3),
            ),
            atoms: vec![
                Atom::new(6, [0.13, 0.27, 0.41]),
                Atom::new(7, [0.6, 0.3, 0.1]),
            ],
        }
    }

    fn triclinic_p1bar() -> Crystal {
        // Centrosymmetric: atom at general position + its inversion image
        Crystal {
            lattice: Lattice::new(
                Vector3::new(3.1, 0.2, 0.1),
                Vector3::new(0.3, 4.2, -0.1),
                Vector3::new(-0.1, 0.15, 5.3),
            ),
            atoms: vec![
                Atom::new(6, [0.13, 0.27, 0.41]),
                Atom::new(6, [0.87, 0.73, 0.59]), // inversion of (0.13, 0.27, 0.41)
            ],
        }
    }

    #[test]
    fn test_si_fcc_48_operations() {
        let crystal = si_fcc();
        let ops = find_symmetry_operations(&crystal, 1e-5);
        assert_eq!(
            ops.len(),
            48,
            "FCC Si (Fd-3m) should have 48 operations, found {}",
            ops.len()
        );
    }

    #[test]
    fn test_bcc_fe_48_operations() {
        let crystal = bcc_fe();
        let ops = find_symmetry_operations(&crystal, 1e-5);
        assert_eq!(
            ops.len(),
            48,
            "BCC Fe (Im-3m) should have 48 operations, found {}",
            ops.len()
        );
    }

    #[test]
    fn test_triclinic_p1() {
        let crystal = triclinic_p1();
        let ops = find_symmetry_operations(&crystal, 1e-5);
        assert_eq!(ops.len(), 1, "P1 should have only identity, found {}", ops.len());
    }

    #[test]
    fn test_triclinic_p1bar() {
        let crystal = triclinic_p1bar();
        let ops = find_symmetry_operations(&crystal, 1e-5);
        assert_eq!(
            ops.len(),
            2,
            "P-1 should have 2 operations (identity + inversion), found {}",
            ops.len()
        );
    }

    #[test]
    fn test_si_has_identity() {
        let crystal = si_fcc();
        let ops = find_symmetry_operations(&crystal, 1e-5);
        assert!(ops.iter().any(|op| op.is_identity(1e-5)));
    }

    #[test]
    fn test_si_has_inversion() {
        let crystal = si_fcc();
        let ops = find_symmetry_operations(&crystal, 1e-5);
        let has_inv = ops
            .iter()
            .any(|op| op.rotation == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]);
        assert!(has_inv, "FCC Si should have inversion symmetry");
    }

    #[test]
    fn test_si_group_closure() {
        let crystal = si_fcc();
        let info = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        info.verify_group_closure();
    }

    #[test]
    fn test_bcc_group_closure() {
        let crystal = bcc_fe();
        let info = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        info.verify_group_closure();
    }

    #[test]
    fn test_all_dets_valid() {
        let crystal = si_fcc();
        let ops = find_symmetry_operations(&crystal, 1e-5);
        for op in &ops {
            let d = op.det();
            assert!(d == 1 || d == -1, "invalid det {d} for {:?}", op.rotation);
        }
    }

    #[test]
    fn test_si_proper_improper_count() {
        let crystal = si_fcc();
        let ops = find_symmetry_operations(&crystal, 1e-5);
        let n_proper = ops.iter().filter(|op| op.det() == 1).count();
        let n_improper = ops.iter().filter(|op| op.det() == -1).count();
        // O_h has 24 proper + 24 improper = 48
        assert_eq!(n_proper, 24, "expected 24 proper rotations");
        assert_eq!(n_improper, 24, "expected 24 improper rotations");
    }
}
