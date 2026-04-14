//! Crystallographic symmetry operations in fractional coordinates.
//!
//! A symmetry operation {R|τ} transforms a point in fractional coordinates as:
//!   r' = R · r + τ  (mod 1)
//!
//! R is always a 3×3 integer matrix with det(R) = ±1.
//! τ is a fractional translation with components in [0, 1).

/// A crystallographic symmetry operation {R|τ} in fractional coordinates.
///
/// The rotation matrix R operates on column vectors of fractional coordinates.
/// To convert to Cartesian: x' = L R L⁻¹ x, where L is the lattice matrix
/// (columns = lattice vectors).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SymmOp {
    /// 3×3 integer rotation/reflection matrix in fractional coordinates.
    /// Row-major: rotation[i][j] is element (i,j).
    pub rotation: [[i32; 3]; 3],
}

/// A symmetry operation with an associated fractional translation.
/// Separated from SymmOp because the translation requires floating-point
/// comparison and cannot participate in Eq/Hash.
#[derive(Debug, Clone)]
pub struct SpaceGroupOp {
    /// The point group rotation part.
    pub rotation: [[i32; 3]; 3],
    /// Fractional translation, each component in [0, 1).
    pub translation: [f64; 3],
}

impl SymmOp {
    /// The identity operation.
    pub fn identity() -> Self {
        Self {
            rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        }
    }

    /// Inversion operation.
    pub fn inversion() -> Self {
        Self {
            rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
        }
    }

    /// Construct from a flat row-major array.
    pub fn from_flat(m: [i32; 9]) -> Self {
        Self {
            rotation: [
                [m[0], m[1], m[2]],
                [m[3], m[4], m[5]],
                [m[6], m[7], m[8]],
            ],
        }
    }

    /// Determinant of the rotation matrix.
    /// +1 for proper rotations, -1 for improper (inversion, mirrors, rotoinversion).
    pub fn det(&self) -> i32 {
        let r = &self.rotation;
        r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
    }

    /// Trace of the rotation matrix.
    pub fn trace(&self) -> i32 {
        self.rotation[0][0] + self.rotation[1][1] + self.rotation[2][2]
    }

    /// Whether this is the identity operation.
    pub fn is_identity(&self) -> bool {
        self.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    }

    /// Inverse rotation matrix.
    /// Since det = ±1 and all entries are integers, the inverse is the adjugate
    /// divided by the determinant (which is also integer).
    pub fn inverse(&self) -> Self {
        let r = &self.rotation;
        let d = self.det();
        assert!(d == 1 || d == -1, "invalid rotation: det = {d}");

        // Adjugate (cofactor matrix transposed)
        let adj = [
            [
                (r[1][1] * r[2][2] - r[1][2] * r[2][1]) * d,
                (r[0][2] * r[2][1] - r[0][1] * r[2][2]) * d,
                (r[0][1] * r[1][2] - r[0][2] * r[1][1]) * d,
            ],
            [
                (r[1][2] * r[2][0] - r[1][0] * r[2][2]) * d,
                (r[0][0] * r[2][2] - r[0][2] * r[2][0]) * d,
                (r[0][2] * r[1][0] - r[0][0] * r[1][2]) * d,
            ],
            [
                (r[1][0] * r[2][1] - r[1][1] * r[2][0]) * d,
                (r[0][1] * r[2][0] - r[0][0] * r[2][1]) * d,
                (r[0][0] * r[1][1] - r[0][1] * r[1][0]) * d,
            ],
        ];

        Self { rotation: adj }
    }

    /// Compose two operations: self ∘ other = R_self · R_other.
    pub fn compose(&self, other: &Self) -> Self {
        let a = &self.rotation;
        let b = &other.rotation;
        let mut r = [[0i32; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
            }
        }
        Self { rotation: r }
    }

    /// Apply rotation to a fractional coordinate vector (no translation).
    pub fn apply(&self, f: &[f64; 3]) -> [f64; 3] {
        let r = &self.rotation;
        [
            r[0][0] as f64 * f[0] + r[0][1] as f64 * f[1] + r[0][2] as f64 * f[2],
            r[1][0] as f64 * f[0] + r[1][1] as f64 * f[1] + r[1][2] as f64 * f[2],
            r[2][0] as f64 * f[0] + r[2][1] as f64 * f[1] + r[2][2] as f64 * f[2],
        ]
    }

    /// Transpose of the inverse: (R⁻¹)ᵀ.
    /// This is the transformation matrix for reciprocal-space vectors.
    pub fn inverse_transpose(&self) -> Self {
        let inv = self.inverse();
        let r = &inv.rotation;
        Self {
            rotation: [
                [r[0][0], r[1][0], r[2][0]],
                [r[0][1], r[1][1], r[2][1]],
                [r[0][2], r[1][2], r[2][2]],
            ],
        }
    }
}

impl SpaceGroupOp {
    /// Construct from rotation and translation.
    pub fn new(rotation: [[i32; 3]; 3], translation: [f64; 3]) -> Self {
        Self {
            rotation,
            translation: wrap_to_unit_cell(translation),
        }
    }

    /// The identity operation (no rotation, no translation).
    pub fn identity() -> Self {
        Self::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0])
    }

    /// Compute the inverse: {R|τ}⁻¹ = {R⁻¹ | -R⁻¹τ}.
    pub fn inverse(&self) -> Self {
        let r_op = SymmOp {
            rotation: self.rotation,
        };
        let r_inv = r_op.inverse();
        let neg_rinv_tau = r_inv.apply(&self.translation);
        Self::new(
            r_inv.rotation,
            [-neg_rinv_tau[0], -neg_rinv_tau[1], -neg_rinv_tau[2]],
        )
    }

    /// Compose: {R₁|τ₁} ∘ {R₂|τ₂} = {R₁R₂ | R₁τ₂ + τ₁}.
    pub fn compose(&self, other: &Self) -> Self {
        let r1 = SymmOp {
            rotation: self.rotation,
        };
        let r2 = SymmOp {
            rotation: other.rotation,
        };
        let r12 = r1.compose(&r2);
        let r1_tau2 = r1.apply(&other.translation);
        Self::new(
            r12.rotation,
            [
                r1_tau2[0] + self.translation[0],
                r1_tau2[1] + self.translation[1],
                r1_tau2[2] + self.translation[2],
            ],
        )
    }

    /// Apply to a fractional coordinate: r' = R·r + τ (mod 1).
    pub fn apply_to_fractional(&self, f: &[f64; 3]) -> [f64; 3] {
        let r = &self.rotation;
        wrap_to_unit_cell([
            r[0][0] as f64 * f[0] + r[0][1] as f64 * f[1] + r[0][2] as f64 * f[2]
                + self.translation[0],
            r[1][0] as f64 * f[0] + r[1][1] as f64 * f[1] + r[1][2] as f64 * f[2]
                + self.translation[1],
            r[2][0] as f64 * f[0] + r[2][1] as f64 * f[1] + r[2][2] as f64 * f[2]
                + self.translation[2],
        ])
    }

    /// Determinant of the rotation part.
    pub fn det(&self) -> i32 {
        SymmOp {
            rotation: self.rotation,
        }
        .det()
    }

    /// Check if this is approximately the identity operation.
    pub fn is_identity(&self, tol: f64) -> bool {
        let r = &self.rotation;
        r == &[[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            && self.translation.iter().all(|&t| t.abs() < tol || (1.0 - t).abs() < tol)
    }

    /// Check approximate equality with another operation.
    pub fn approx_eq(&self, other: &Self, tol: f64) -> bool {
        self.rotation == other.rotation && frac_distance(&self.translation, &other.translation) < tol
    }
}

/// Wrap fractional coordinates to [0, 1).
pub fn wrap_to_unit_cell(f: [f64; 3]) -> [f64; 3] {
    [
        f[0] - f[0].floor(),
        f[1] - f[1].floor(),
        f[2] - f[2].floor(),
    ]
}

/// Distance between two fractional coordinate vectors with minimum image convention.
pub fn frac_distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let mut d2 = 0.0;
    for i in 0..3 {
        let mut dx = a[i] - b[i];
        dx -= dx.round();
        d2 += dx * dx;
    }
    d2.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity() {
        let id = SymmOp::identity();
        assert_eq!(id.det(), 1);
        assert_eq!(id.trace(), 3);
        assert!(id.is_identity());
    }

    #[test]
    fn test_inversion() {
        let inv = SymmOp::inversion();
        assert_eq!(inv.det(), -1);
        assert_eq!(inv.trace(), -3);
        assert!(!inv.is_identity());
    }

    #[test]
    fn test_inverse_of_identity() {
        let id = SymmOp::identity();
        let inv = id.inverse();
        assert!(inv.is_identity());
    }

    #[test]
    fn test_inverse_of_inversion() {
        let inv = SymmOp::inversion();
        let inv2 = inv.inverse();
        assert_eq!(inv, inv2); // inversion is its own inverse
    }

    #[test]
    fn test_compose_identity() {
        let id = SymmOp::identity();
        let rot = SymmOp::from_flat([0, -1, 0, 1, 0, 0, 0, 0, 1]); // 90° around z
        assert_eq!(id.compose(&rot), rot);
        assert_eq!(rot.compose(&id), rot);
    }

    #[test]
    fn test_90_degree_rotation_order_4() {
        // 90° rotation about z-axis in cubic fractional coords
        let rot = SymmOp::from_flat([0, -1, 0, 1, 0, 0, 0, 0, 1]);
        assert_eq!(rot.det(), 1);
        assert_eq!(rot.trace(), 1); // trace=1 → 90° rotation

        // rot^4 = identity
        let rot2 = rot.compose(&rot);
        let rot3 = rot2.compose(&rot);
        let rot4 = rot3.compose(&rot);
        assert!(rot4.is_identity());
        assert!(!rot2.is_identity());
    }

    #[test]
    fn test_inverse_compose_is_identity() {
        let rot = SymmOp::from_flat([0, 1, 0, 0, 0, 1, 1, 0, 0]); // 3-fold rotation
        let inv = rot.inverse();
        let prod = rot.compose(&inv);
        assert!(prod.is_identity());
        let prod2 = inv.compose(&rot);
        assert!(prod2.is_identity());
    }

    #[test]
    fn test_apply_identity() {
        let id = SymmOp::identity();
        let f = [0.25, 0.5, 0.75];
        let result = id.apply(&f);
        assert!((result[0] - 0.25).abs() < 1e-15);
        assert!((result[1] - 0.5).abs() < 1e-15);
        assert!((result[2] - 0.75).abs() < 1e-15);
    }

    #[test]
    fn test_apply_inversion() {
        let inv = SymmOp::inversion();
        let f = [0.25, 0.3, 0.1];
        let result = inv.apply(&f);
        assert!((result[0] + 0.25).abs() < 1e-15);
        assert!((result[1] + 0.3).abs() < 1e-15);
        assert!((result[2] + 0.1).abs() < 1e-15);
    }

    #[test]
    fn test_inverse_transpose() {
        // For a proper rotation R, (R⁻¹)ᵀ = R for orthogonal matrices in Cartesian.
        // But in fractional coords this is NOT generally true.
        let rot = SymmOp::from_flat([0, -1, 0, 1, 0, 0, 0, 0, 1]);
        let rit = rot.inverse_transpose();
        // (R⁻¹)ᵀ should also have det ±1
        assert!(rit.det().abs() == 1);
        // Applying rit then rot should give identity for the transpose product
        // R (R⁻¹)ᵀ·v is NOT identity in general, but (R⁻¹)ᵀ · R^T = I
    }

    #[test]
    fn test_space_group_op_compose() {
        let op1 = SpaceGroupOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.5, 0.0, 0.0]);
        let op2 = SpaceGroupOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.5, 0.0, 0.0]);
        let prod = op1.compose(&op2);
        // Two half-translations compose to full translation = identity (mod 1)
        assert!(prod.is_identity(1e-10));
    }

    #[test]
    fn test_space_group_op_inverse() {
        let op = SpaceGroupOp::new(
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            [0.25, 0.25, 0.25],
        );
        let inv = op.inverse();
        let prod = op.compose(&inv);
        assert!(
            prod.is_identity(1e-10),
            "op · op⁻¹ should be identity: rotation={:?}, translation={:?}",
            prod.rotation,
            prod.translation
        );
    }

    #[test]
    fn test_space_group_op_apply() {
        let op = SpaceGroupOp::new(
            [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], // inversion
            [0.5, 0.5, 0.5],                        // + half translation
        );
        // Apply to (0.25, 0.25, 0.25): -0.25 + 0.5 = 0.25 → (0.25, 0.25, 0.25)
        let result = op.apply_to_fractional(&[0.25, 0.25, 0.25]);
        assert!((result[0] - 0.25).abs() < 1e-10);
        assert!((result[1] - 0.25).abs() < 1e-10);
        assert!((result[2] - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_wrap_to_unit_cell() {
        assert!((wrap_to_unit_cell([1.3, -0.2, 0.5])[0] - 0.3).abs() < 1e-10);
        assert!((wrap_to_unit_cell([1.3, -0.2, 0.5])[1] - 0.8).abs() < 1e-10);
        assert!((wrap_to_unit_cell([1.3, -0.2, 0.5])[2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_frac_distance() {
        // Same point
        assert!(frac_distance(&[0.0, 0.0, 0.0], &[0.0, 0.0, 0.0]) < 1e-15);
        // Across boundary: (0.01, 0, 0) and (0.99, 0, 0) are distance 0.02
        assert!((frac_distance(&[0.01, 0.0, 0.0], &[0.99, 0.0, 0.0]) - 0.02).abs() < 1e-10);
        // Diagonal
        let d = frac_distance(&[0.0, 0.0, 0.0], &[0.5, 0.5, 0.5]);
        assert!((d - 0.5 * 3.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_all_cubic_rotations_have_valid_det() {
        // Generate some known cubic rotations and verify properties
        let rots = vec![
            SymmOp::identity(),
            SymmOp::inversion(),
            SymmOp::from_flat([0, 1, 0, 0, 0, 1, 1, 0, 0]),   // 3-fold [111]
            SymmOp::from_flat([0, -1, 0, 1, 0, 0, 0, 0, 1]),  // 4-fold [001]
            SymmOp::from_flat([-1, 0, 0, 0, 1, 0, 0, 0, 1]),  // mirror perp to [100]
        ];
        for r in &rots {
            let d = r.det();
            assert!(d == 1 || d == -1, "det = {d} for {:?}", r.rotation);
            // R · R⁻¹ = I
            let prod = r.compose(&r.inverse());
            assert!(prod.is_identity(), "R·R⁻¹ ≠ I for {:?}", r.rotation);
        }
    }
}
