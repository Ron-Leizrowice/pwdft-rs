use nalgebra::Vector3;

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

/// Generate a uniform Monkhorst-Pack k-point grid (shifted convention).
///
/// Produces fractional coordinates
/// `f_j = (2·i_j − N_j + 1) / (2·N_j)` for `i_j ∈ {0, …, N_j−1}`, then
/// converted to Cartesian reciprocal space. For N=4 this yields
/// `{−3/8, −1/8, 1/8, 3/8}` — the original **shifted** Monkhorst-Pack
/// grid (Monkhorst & Pack, *Phys. Rev. B* **13**, 5188 (1976), Eq. 4).
/// No symmetry reduction is applied (full grid).
///
/// ## Convention caveat vs QE
///
/// This is equivalent to Quantum ESPRESSO's
/// `K_POINTS automatic / nk1 nk2 nk3 1 1 1` (half-shift along every
/// axis), **not** the default `0 0 0` Γ-centred grid. The Γ-centred
/// version would produce `{0, 1/4, 1/2, 3/4}` and contain Γ, X, L
/// exactly; the shifted version avoids all BZ-boundary points.
///
/// For even N these are physically distinct k-meshes of the same
/// density, and their IBZ reductions generally differ in the number
/// of irreducible points. This is a convention choice, not a bug —
/// see `proposals/completed/SYKP-symmetry-ibz-audit.md`.
#[must_use]
pub fn monkhorst_pack(n1: u32, n2: u32, n3: u32, lattice: &Lattice) -> Vec<KPoint> {
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
                // MP mesh counts (n_i, i_i) are bounded by O(100) in practice;
                // i32 casts here cannot wrap for any physical input.
                #[allow(
                    clippy::cast_possible_wrap,
                    reason = "MP mesh counts bounded by O(100); u32 < i32::MAX trivially"
                )]
                let f1 = f64::from(2 * i1 as i32 - n1 as i32 + 1) / (2.0 * f64::from(n1));
                #[allow(
                    clippy::cast_possible_wrap,
                    reason = "MP mesh counts bounded by O(100); u32 < i32::MAX trivially"
                )]
                let f2 = f64::from(2 * i2 as i32 - n2 as i32 + 1) / (2.0 * f64::from(n2));
                #[allow(
                    clippy::cast_possible_wrap,
                    reason = "MP mesh counts bounded by O(100); u32 < i32::MAX trivially"
                )]
                let f3 = f64::from(2 * i3 as i32 - n3 as i32 + 1) / (2.0 * f64::from(n3));

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
        let kpts = monkhorst_pack(4, 4, 4, &si_lattice());
        assert_eq!(kpts.len(), 64);
    }

    #[test]
    fn test_monkhorst_pack_weights_sum_to_one() {
        let kpts = monkhorst_pack(4, 4, 4, &si_lattice());
        let total: f64 = kpts.iter().map(|k| k.weight).sum();
        assert!(
            relative_eq!(total, 1.0, epsilon = 1e-12),
            "weights sum to {total}, expected 1.0"
        );
    }

    #[test]
    fn test_monkhorst_pack_gamma_centered_odd() {
        // Odd grids include Γ point
        let kpts = monkhorst_pack(3, 3, 3, &si_lattice());
        let has_gamma = kpts.iter().any(|kp| kp.k.norm() < 1e-10);
        assert!(has_gamma, "3×3×3 MP grid should include Γ point");
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
