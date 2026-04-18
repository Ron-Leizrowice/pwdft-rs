use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::crystal::Lattice;

/// A single k-point with its Cartesian coordinates and BZ integration weight.
#[derive(Debug, Clone)]
pub struct KPoint {
    /// k-vector in Cartesian reciprocal space (1/Å).
    pub k: Vector3<f64>,
    /// Integration weight (sums to 1 over the full grid).
    pub weight: f64,
    /// Optional label for high-symmetry points.
    pub label: Option<String>,
}

/// Shift convention for a Monkhorst-Pack k-grid.
///
/// Selects between the two common formulas for a uniform k-mesh:
///
/// - [`GammaCentered`](KGridShift::GammaCentered) — `f_j = (i_j − 1)/N_j`
///   (QE's default `K_POINTS automatic / nk1 nk2 nk3 0 0 0`). For N=4
///   this yields `{0, 1/4, 1/2, 3/4}`, i.e. **includes** Γ and the
///   high-symmetry BZ-boundary points (X, L on FCC).
/// - [`MP1976`](KGridShift::MP1976) — `f_j = (2·i_j − N_j + 1)/(2·N_j)`
///   (the original Monkhorst & Pack, *Phys. Rev. B* **13**, 5188 (1976)
///   Eq. 4; equivalent to QE's `1 1 1` half-shift). For N=4 this yields
///   `{−3/8, −1/8, 1/8, 3/8}` — **no** point at Γ, symmetric about Γ.
/// - [`Custom`](KGridShift::Custom) — per-axis half-shift flags
///   `k_α ∈ {0, 1}` matching QE's free-form third line: `f_j = (i_j − 1)/N_j
///   + k_j/(2·N_j)`. `Custom(\[0,0,0\])` is identical to `GammaCentered`.
///
/// For odd N both `GammaCentered` and `MP1976` produce the same mesh
/// modulo a cyclic permutation; for even N they are genuinely distinct
/// physical samples of the BZ.
///
/// Default: `GammaCentered` — matches Quantum ESPRESSO's default k-mesh
/// for byte-identical validation comparisons.
///
/// ## YAML syntax
///
/// ```yaml
/// shift: gamma_centered     # or the symbols below
/// shift: mp1976             # aliases: shifted, mp_1976
/// shift:
///   custom: [1, 0, 1]       # per-axis half-shift flags
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum KGridShift {
    /// Γ-centered grid (QE default, `0 0 0`).
    #[default]
    GammaCentered,
    /// Original Monkhorst & Pack 1976 shifted grid (QE `1 1 1`).
    #[serde(rename = "mp1976", alias = "mp_1976", alias = "shifted")]
    MP1976,
    /// Per-axis half-shift flags (QE's third line on `K_POINTS automatic`).
    Custom([u32; 3]),
}

impl KGridShift {
    /// Return the per-axis half-shift integers `k_α ∈ {0, 1}` that
    /// parameterise this shift under the QE convention.
    #[must_use]
    pub fn axis_flags(self) -> [u32; 3] {
        match self {
            Self::GammaCentered => [0, 0, 0],
            Self::MP1976 => [1, 1, 1],
            Self::Custom(flags) => flags,
        }
    }
}

/// Compute the fractional reciprocal coordinate of grid point `(i1,i2,i3)`
/// on an `N1×N2×N3` Monkhorst-Pack mesh with the given shift.
///
/// Formula (QE `kpoint_grid.f90`): `f_α = (i_α)/N_α + k_α/(2·N_α)` for
/// `i_α ∈ {0, …, N_α−1}` and `k_α ∈ {0, 1}`. The raw result lies in
/// `[0, 1)` for `k_α = 0` and in `[1/(2N), 1)` for `k_α = 1`; this
/// function wraps into the first Brillouin zone `[-1/2, 1/2)` so the
/// SCF sees k-points in canonical form. (Mathematically k and k+G give
/// identical physics, but at finite `ecut` the shared plane-wave basis
/// adapts unequally to them — wrapping keeps kinetic energies minimized
/// and keeps the result bit-compatible with the pre-MPSH mesh.)
#[must_use]
pub fn mp_fractional_coord(
    i1: u32,
    i2: u32,
    i3: u32,
    grid: [u32; 3],
    shift: KGridShift,
) -> [f64; 3] {
    let [k1, k2, k3] = shift.axis_flags();
    let raw = [
        f64::from(i1) / f64::from(grid[0]) + f64::from(k1) / (2.0 * f64::from(grid[0])),
        f64::from(i2) / f64::from(grid[1]) + f64::from(k2) / (2.0 * f64::from(grid[1])),
        f64::from(i3) / f64::from(grid[2]) + f64::from(k3) / (2.0 * f64::from(grid[2])),
    ];
    // Wrap each axis into the first BZ (−1/2, 1/2]: subtract floor(f + 1/2).
    [
        raw[0] - (raw[0] + 0.5).floor(),
        raw[1] - (raw[1] + 0.5).floor(),
        raw[2] - (raw[2] + 0.5).floor(),
    ]
}

/// Generate a uniform Monkhorst-Pack k-point grid.
///
/// Produces the full `N1·N2·N3` grid (no symmetry reduction applied). The
/// `shift` argument selects between the Γ-centered and MP-1976 (shifted)
/// conventions — see [`KGridShift`] for the formulas. Fractional
/// coordinates are converted to Cartesian reciprocal space via `lattice`.
///
/// ## Convention vs QE
///
/// - `KGridShift::GammaCentered` ↔ QE `K_POINTS automatic / nk1 nk2 nk3 0 0 0`.
/// - `KGridShift::MP1976` ↔ QE `K_POINTS automatic / nk1 nk2 nk3 1 1 1`.
/// - `KGridShift::Custom([k1,k2,k3])` ↔ QE `K_POINTS automatic / nk1 nk2 nk3 k1 k2 k3`.
#[must_use]
pub fn monkhorst_pack(
    n1: u32,
    n2: u32,
    n3: u32,
    shift: KGridShift,
    lattice: &Lattice,
) -> Vec<KPoint> {
    let recip = lattice.reciprocal();
    let ntotal = f64::from(n1 * n2 * n3);
    let weight = 1.0 / ntotal;

    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "ntotal = n1*n2*n3 with each n_i a small MP mesh count (typically <= 64); f64->usize via truncation is exact for integer-valued f64s below 2^52"
    )]
    let mut kpoints = Vec::with_capacity(ntotal as usize);
    for i1 in 0..n1 {
        for i2 in 0..n2 {
            for i3 in 0..n3 {
                let [f1, f2, f3] = mp_fractional_coord(i1, i2, i3, [n1, n2, n3], shift);

                let k = f1 * recip.a + f2 * recip.b + f3 * recip.c;
                kpoints.push(KPoint {
                    k,
                    weight,
                    label: None,
                });
            }
        }
    }
    kpoints
}

/// A named high-symmetry point in fractional reciprocal coordinates.
#[derive(Debug, Clone)]
pub struct HighSymPoint {
    pub label: String,
    pub frac: [f64; 3],
}

/// Standard high-symmetry points for FCC Brillouin zone.
#[must_use]
pub fn fcc_high_sym_points() -> Vec<HighSymPoint> {
    vec![
        HighSymPoint {
            label: "Γ".into(),
            frac: [0.0, 0.0, 0.0],
        },
        HighSymPoint {
            label: "X".into(),
            frac: [0.5, 0.0, 0.5],
        },
        HighSymPoint {
            label: "W".into(),
            frac: [0.5, 0.25, 0.75],
        },
        HighSymPoint {
            label: "K".into(),
            frac: [0.375, 0.375, 0.75],
        },
        HighSymPoint {
            label: "L".into(),
            frac: [0.5, 0.5, 0.5],
        },
        HighSymPoint {
            label: "U".into(),
            frac: [0.625, 0.25, 0.625],
        },
    ]
}

/// Generate a k-point path through high-symmetry points for band structure.
///
/// `segments` is a list of (label, fractional_coords) pairs defining the path vertices.
/// `npoints_per_segment` controls the density of sampling between each pair.
/// Returns k-points with cumulative distance for plotting.
#[must_use]
pub fn high_symmetry_path(
    segments: &[HighSymPoint],
    npoints_per_segment: usize,
    lattice: &Lattice,
) -> (Vec<KPoint>, Vec<f64>) {
    let recip = lattice.reciprocal();
    let mut kpoints = Vec::new();
    let mut distances = Vec::new();
    let mut cumulative_dist = 0.0;

    for seg_idx in 0..segments.len() - 1 {
        let start = &segments[seg_idx];
        let end = &segments[seg_idx + 1];

        let k_start =
            start.frac[0] * recip.a + start.frac[1] * recip.b + start.frac[2] * recip.c;
        let k_end = end.frac[0] * recip.a + end.frac[1] * recip.b + end.frac[2] * recip.c;
        let dk = k_end - k_start;
        let seg_len = dk.norm();

        let n = if seg_idx == segments.len() - 2 {
            npoints_per_segment + 1 // include endpoint on last segment
        } else {
            npoints_per_segment // exclude endpoint to avoid duplication
        };

        for i in 0..n {
            let t = i as f64 / npoints_per_segment as f64;
            let k = k_start + t * dk;

            let label = if i == 0 {
                Some(start.label.clone())
            } else if seg_idx == segments.len() - 2 && i == n - 1 {
                Some(end.label.clone())
            } else {
                None
            };

            kpoints.push(KPoint {
                k,
                weight: 0.0, // band structure points have zero weight for integration
                label,
            });
            distances.push(cumulative_dist + t * seg_len);
        }
        cumulative_dist += seg_len;
    }

    (kpoints, distances)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "ERR2 § Phase 0: in-src test modules are allowed to panic")]
mod tests {
    use super::*;
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
    fn test_monkhorst_pack_count() {
        let kpts = monkhorst_pack(4, 4, 4, KGridShift::GammaCentered, &si_lattice());
        assert_eq!(kpts.len(), 64);
    }

    #[test]
    fn test_monkhorst_pack_weights_sum_to_one() {
        let kpts = monkhorst_pack(4, 4, 4, KGridShift::GammaCentered, &si_lattice());
        let total: f64 = kpts.iter().map(|k| k.weight).sum();
        assert!(
            relative_eq!(total, 1.0, epsilon = 1e-12),
            "weights sum to {total}, expected 1.0"
        );
    }

    #[test]
    fn test_monkhorst_pack_gamma_centered_even_includes_gamma() {
        // Γ-centered 4×4×4 includes Γ and BZ-boundary points exactly.
        let kpts = monkhorst_pack(4, 4, 4, KGridShift::GammaCentered, &si_lattice());
        let has_gamma = kpts.iter().any(|kp| kp.k.norm() < 1e-10);
        assert!(has_gamma, "Γ-centered 4×4×4 MP grid must include Γ");
    }

    #[test]
    fn test_monkhorst_pack_mp1976_even_excludes_gamma() {
        // Original MP-1976 shifted grid: {−3/8, −1/8, 1/8, 3/8} — no Γ.
        let kpts = monkhorst_pack(4, 4, 4, KGridShift::MP1976, &si_lattice());
        let has_gamma = kpts.iter().any(|kp| kp.k.norm() < 1e-10);
        assert!(!has_gamma, "MP-1976 4×4×4 grid must NOT include Γ");
    }

    #[test]
    fn test_monkhorst_pack_gamma_centered_odd() {
        // Odd grids include Γ point under either convention.
        let kpts = monkhorst_pack(3, 3, 3, KGridShift::GammaCentered, &si_lattice());
        let has_gamma = kpts.iter().any(|kp| kp.k.norm() < 1e-10);
        assert!(has_gamma, "3×3×3 Γ-centered MP grid should include Γ point");
    }

    #[test]
    fn test_kgrid_shift_default_is_gamma() {
        // Project-wide default follows QE: Γ-centered.
        assert_eq!(KGridShift::default(), KGridShift::GammaCentered);
        assert_eq!(KGridShift::default().axis_flags(), [0, 0, 0]);
        assert_eq!(KGridShift::MP1976.axis_flags(), [1, 1, 1]);
        assert_eq!(KGridShift::Custom([1, 0, 1]).axis_flags(), [1, 0, 1]);
    }

    #[test]
    fn test_kgrid_shift_custom_equals_gamma_for_all_zero() {
        // Custom([0,0,0]) must produce the same grid as GammaCentered.
        let a = monkhorst_pack(4, 4, 4, KGridShift::GammaCentered, &si_lattice());
        let b = monkhorst_pack(4, 4, 4, KGridShift::Custom([0, 0, 0]), &si_lattice());
        assert_eq!(a.len(), b.len());
        for (ka, kb) in a.iter().zip(b.iter()) {
            assert!(
                (ka.k - kb.k).norm() < 1e-14,
                "Custom([0,0,0]) should match GammaCentered"
            );
        }
    }

    #[test]
    fn test_kgrid_shift_custom_equals_mp1976_for_all_one() {
        // Custom([1,1,1]) must produce the same grid as MP1976.
        let a = monkhorst_pack(4, 4, 4, KGridShift::MP1976, &si_lattice());
        let b = monkhorst_pack(4, 4, 4, KGridShift::Custom([1, 1, 1]), &si_lattice());
        assert_eq!(a.len(), b.len());
        for (ka, kb) in a.iter().zip(b.iter()) {
            assert!(
                (ka.k - kb.k).norm() < 1e-14,
                "Custom([1,1,1]) should match MP1976"
            );
        }
    }

    #[test]
    fn test_high_symmetry_path() {
        let lattice = si_lattice();
        let points = vec![
            HighSymPoint {
                label: "Γ".into(),
                frac: [0.0, 0.0, 0.0],
            },
            HighSymPoint {
                label: "X".into(),
                frac: [0.5, 0.0, 0.5],
            },
            HighSymPoint {
                label: "Γ".into(),
                frac: [0.0, 0.0, 0.0],
            },
        ];
        let (kpts, dists) = high_symmetry_path(&points, 10, &lattice);
        // 10 points on first segment + 11 on last (including endpoint) = 21
        assert_eq!(kpts.len(), 21);
        assert_eq!(dists.len(), 21);
        // First point should be Γ
        assert!(kpts[0].k.norm() < 1e-10);
        // Last point should also be Γ
        assert!(kpts.last().unwrap().k.norm() < 1e-10);
        // Distances should be monotonically non-decreasing
        for i in 1..dists.len() {
            assert!(dists[i] >= dists[i - 1] - 1e-12);
        }
    }
}
