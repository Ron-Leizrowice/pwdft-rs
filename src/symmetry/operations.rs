//! Crystallographic symmetry operations in fractional coordinates.
//!
//! A symmetry operation {R|τ} transforms a point in fractional coordinates as:
//!   r' = R · r + τ  (mod 1)
//!
//! R is always a 3×3 integer matrix with det(R) = ±1.
//! τ is a fractional translation with components in [0, 1).
//!
//! Storage note: the rotation matrix is held as `[[i8; 3]; 3]`. Entries in
//! the fractional basis are `{-1, 0, 1}` for cubic groups and at most
//! `{-2, -1, 0, 1, 2}` for rare hexagonal settings — `i8`'s range
//! `[-128, 127]` gives ~40× headroom over any value `detect.rs`
//! produces. Arithmetic (det, compose, adjugate) widens entries to `i32`
//! at method entry to avoid any risk of overflow in the worst-case
//! triple product (max ~27 in the cubic case). This narrower storage
//! shaves 27 bytes per op vs the prior `[[i32; 3]; 3]` layout and
//! packs all 48 ops of `Fd-3m` into ~2 cache lines.

/// A crystallographic symmetry operation {R|τ} in fractional coordinates.
///
/// The rotation matrix R operates on column vectors of fractional coordinates.
/// To convert to Cartesian: x' = L R L⁻¹ x, where L is the lattice matrix
/// (columns = lattice vectors).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SymmOp {
    /// 3×3 integer rotation/reflection matrix in fractional coordinates.
    /// Row-major: `rotation[i][j]` is element (i,j). Entries are in the
    /// bounded set produced by `detect::find_symmetry_operations` (at most
    /// `|R_ij| ≤ 3` for crystallographic point groups); all arithmetic on
    /// the matrix widens to `i32` to preserve overflow-free behavior.
    pub rotation: [[i8; 3]; 3],
}

/// A symmetry operation with an associated fractional translation.
/// Separated from SymmOp because the translation requires floating-point
/// comparison and cannot participate in Eq/Hash.
#[derive(Debug, Clone)]
pub struct SpaceGroupOp {
    /// The point group rotation part. Stored as `i8` (see `SymmOp::rotation`
    /// for the range argument).
    pub rotation: [[i8; 3]; 3],
    /// Fractional translation, each component in [0, 1).
    pub translation: [f64; 3],
}

impl SymmOp {
    /// The identity operation.
    #[must_use]
    pub fn identity() -> Self {
        Self {
            rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        }
    }

    /// Inversion operation.
    #[must_use]
    pub fn inversion() -> Self {
        Self {
            rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
        }
    }

    /// Construct from a flat row-major array.
    ///
    /// Accepts `i32` for backward-compatible call sites; entries are
    /// narrowed to the internal `i8` storage. Values outside `i8::MIN..=i8::MAX`
    /// are a programming error and will panic (checked via `i8::try_from`).
    #[must_use]
    pub fn from_flat(m: [i32; 9]) -> Self {
        #[expect(
            clippy::expect_used,
            reason = "BUG: TYPE-A narrowing; crystallographic rotation entries are in {-2..=2}, well within i8 range. The try_from here is structurally infallible — a panic indicates a caller passing non-crystallographic input."
        )]
        let to_i8 = |v: i32| {
            i8::try_from(v).expect("SymmOp::from_flat: rotation entry out of i8 range")
        };
        Self {
            rotation: [
                [to_i8(m[0]), to_i8(m[1]), to_i8(m[2])],
                [to_i8(m[3]), to_i8(m[4]), to_i8(m[5])],
                [to_i8(m[6]), to_i8(m[7]), to_i8(m[8])],
            ],
        }
    }

    /// Determinant of the rotation matrix.
    /// +1 for proper rotations, -1 for improper (inversion, mirrors, rotoinversion).
    ///
    /// Widens each entry to `i32` before multiplying; the worst-case triple
    /// product for crystallographic entries (`|R_ij| ≤ 3`) is 81, far below
    /// `i32::MAX`.
    #[must_use]
    pub fn det(&self) -> i32 {
        let r = self.rotation_i32();
        r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
    }

    /// Trace of the rotation matrix.
    #[must_use]
    pub fn trace(&self) -> i32 {
        i32::from(self.rotation[0][0])
            + i32::from(self.rotation[1][1])
            + i32::from(self.rotation[2][2])
    }

    /// Whether this is the identity operation.
    #[must_use]
    pub fn is_identity(&self) -> bool {
        self.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    }

    /// Inverse rotation matrix.
    /// Since det = ±1 and all entries are integers, the inverse is the adjugate
    /// divided by the determinant (which is also integer).
    #[must_use]
    pub fn inverse(&self) -> Self {
        let r = self.rotation_i32();
        let d = self.det();
        assert!(d == 1 || d == -1, "invalid rotation: det = {d}");

        #[expect(
            clippy::expect_used,
            reason = "BUG: TYPE-A narrowing; adjugate entries of a crystallographic rotation (|R_ij| <= 2) are bounded by 2*2 - (-2)*(-2) = 0..=8 times det=±1, well within i8 range. Infallible by construction."
        )]
        let to_i8 = |v: i32| {
            i8::try_from(v).expect("SymmOp::inverse: adjugate entry out of i8 range")
        };

        // Adjugate (cofactor matrix transposed)
        let adj = [
            [
                to_i8((r[1][1] * r[2][2] - r[1][2] * r[2][1]) * d),
                to_i8((r[0][2] * r[2][1] - r[0][1] * r[2][2]) * d),
                to_i8((r[0][1] * r[1][2] - r[0][2] * r[1][1]) * d),
            ],
            [
                to_i8((r[1][2] * r[2][0] - r[1][0] * r[2][2]) * d),
                to_i8((r[0][0] * r[2][2] - r[0][2] * r[2][0]) * d),
                to_i8((r[0][2] * r[1][0] - r[0][0] * r[1][2]) * d),
            ],
            [
                to_i8((r[1][0] * r[2][1] - r[1][1] * r[2][0]) * d),
                to_i8((r[0][1] * r[2][0] - r[0][0] * r[2][1]) * d),
                to_i8((r[0][0] * r[1][1] - r[0][1] * r[1][0]) * d),
            ],
        ];

        Self { rotation: adj }
    }

    /// Compose two operations: self ∘ other = R_self · R_other.
    #[must_use]
    pub fn compose(&self, other: &Self) -> Self {
        let a = self.rotation_i32();
        let b = other.rotation_i32();
        #[expect(
            clippy::expect_used,
            reason = "BUG: TYPE-A narrowing; product of two crystallographic rotations has entries bounded by 3*|R_ij|_max^2 <= 3*2*2 = 12 (cubic case; hexagonal worst case <= 27), well within i8 range. Infallible by construction."
        )]
        let to_i8 = |v: i32| {
            i8::try_from(v).expect("SymmOp::compose: product entry out of i8 range")
        };
        let mut r = [[0i8; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = to_i8(a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j]);
            }
        }
        Self { rotation: r }
    }

    /// Apply rotation to a fractional coordinate vector (no translation).
    #[must_use]
    pub fn apply(&self, f: &[f64; 3]) -> [f64; 3] {
        let r = &self.rotation;
        [
            f64::from(r[0][0]) * f[0] + f64::from(r[0][1]) * f[1] + f64::from(r[0][2]) * f[2],
            f64::from(r[1][0]) * f[0] + f64::from(r[1][1]) * f[1] + f64::from(r[1][2]) * f[2],
            f64::from(r[2][0]) * f[0] + f64::from(r[2][1]) * f[1] + f64::from(r[2][2]) * f[2],
        ]
    }

    /// Transpose of the inverse: (R⁻¹)ᵀ.
    /// This is the transformation matrix for reciprocal-space vectors.
    #[must_use]
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

    /// Widen the rotation matrix to `i32` for arithmetic that can overflow
    /// `i8` (adjugate, compose, determinant). Inlined; the compiler emits a
    /// handful of sign-extension moves for a 9-element `i8` read.
    #[inline]
    fn rotation_i32(&self) -> [[i32; 3]; 3] {
        let r = &self.rotation;
        [
            [i32::from(r[0][0]), i32::from(r[0][1]), i32::from(r[0][2])],
            [i32::from(r[1][0]), i32::from(r[1][1]), i32::from(r[1][2])],
            [i32::from(r[2][0]), i32::from(r[2][1]), i32::from(r[2][2])],
        ]
    }
}

impl SpaceGroupOp {
    /// Construct from rotation and translation.
    #[must_use]
    pub fn new(rotation: [[i8; 3]; 3], translation: [f64; 3]) -> Self {
        Self {
            rotation,
            translation: wrap_to_unit_cell(translation),
        }
    }

    /// The identity operation (no rotation, no translation).
    #[must_use]
    pub fn identity() -> Self {
        Self::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0])
    }

    /// Compute the inverse: {R|τ}⁻¹ = {R⁻¹ | -R⁻¹τ}.
    #[must_use]
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
    #[must_use]
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
    #[must_use]
    pub fn apply_to_fractional(&self, f: &[f64; 3]) -> [f64; 3] {
        let r = &self.rotation;
        wrap_to_unit_cell([
            f64::from(r[0][0]) * f[0]
                + f64::from(r[0][1]) * f[1]
                + f64::from(r[0][2]) * f[2]
                + self.translation[0],
            f64::from(r[1][0]) * f[0]
                + f64::from(r[1][1]) * f[1]
                + f64::from(r[1][2]) * f[2]
                + self.translation[1],
            f64::from(r[2][0]) * f[0]
                + f64::from(r[2][1]) * f[1]
                + f64::from(r[2][2]) * f[2]
                + self.translation[2],
        ])
    }

    /// Determinant of the rotation part.
    #[must_use]
    pub fn det(&self) -> i32 {
        SymmOp {
            rotation: self.rotation,
        }
        .det()
    }

    /// Check if this is approximately the identity operation.
    #[must_use]
    pub fn is_identity(&self, tol: f64) -> bool {
        let r = &self.rotation;
        r == &[[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            && self.translation.iter().all(|&t| t.abs() < tol || (1.0 - t).abs() < tol)
    }

    /// Check approximate equality with another operation.
    #[must_use]
    pub fn approx_eq(&self, other: &Self, tol: f64) -> bool {
        self.rotation == other.rotation && frac_distance(&self.translation, &other.translation) < tol
    }
}

/// Wrap fractional coordinates to [0, 1).
#[must_use]
pub fn wrap_to_unit_cell(f: [f64; 3]) -> [f64; 3] {
    [
        f[0] - f[0].floor(),
        f[1] - f[1].floor(),
        f[2] - f[2].floor(),
    ]
}

/// Distance between two fractional coordinate vectors with minimum image convention.
#[must_use]
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
