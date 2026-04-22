//! Brillouin-zone sampling: Monkhorst-Pack grids and high-symmetry paths.
//!
//! A [`KPoint`] is a Cartesian reciprocal-space vector (1/Å) paired with
//! an integration weight. The two standard constructors are
//! [`monkhorst_pack`] (uniform `N₁×N₂×N₃` grids with a selectable
//! [`KGridShift`] — gamma-centered or the 1976 MP offset) and
//! [`high_symmetry_path`] (piecewise-linear k-paths through named BZ
//! corners, used for band-structure plots).
//!
//! Symmetry reduction of the full grid to the irreducible BZ lives in
//! [`crate::symmetry::kpoints`]; this module generates the unreduced
//! sampling.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::crystal::Lattice;

/// A collection of k-points with enough provenance to symmetrize densities
/// and interpret band-path distances. Construct via the inherent
/// constructors; do not build one by hand outside this module.
#[derive(Debug, Clone)]
pub struct KPointSet {
    points: Vec<KPoint>,
    kind: SamplingKind,
}

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

impl KPointSet {
    /// Γ-only sampling. Single k-point at the origin with weight 1.
    ///
    /// Used by molecular tests, free-electron Γ spot-checks, and any
    /// mixing-algorithm regression that wants to bypass a Brillouin-zone
    /// integration. Tagged as `MonkhorstPack { grid: [1,1,1], shift: GammaCentered }`
    /// so downstream consumers treating this as a degenerate 1×1×1 grid
    /// get consistent metadata.
    pub fn gamma_only() -> Self {
        Self {
            points: vec![KPoint {
                k: Vector3::zeros(),
                weight: 1.0,
                label: Some("Γ".into()),
            }],
            kind: SamplingKind::MonkhorstPack {
                grid: [1, 1, 1],
                shift: KGridShift::GammaCentered,
            },
        }
    }

    /// Full (unreduced) Monkhorst-Pack grid. See [`monkhorst_pack`] for
    /// the underlying formula.
    pub fn monkhorst_pack(grid: [u32; 3], shift: KGridShift, lattice: &Lattice) -> Self {
        Self {
            points: monkhorst_pack(grid[0], grid[1], grid[2], shift, lattice),
            kind: SamplingKind::MonkhorstPack { grid, shift },
        }
    }

    /// Piecewise-linear path through high-symmetry points for band
    /// structure. See [`high_symmetry_path`] for the interpolation rule.
    pub fn band_path(segments: &[HighSymPoint], npoints_per_segment: usize, lattice: &Lattice) -> Self {
        let (points, distances) = high_symmetry_path(segments, npoints_per_segment, lattice);
        Self {
            points,
            kind: SamplingKind::BandPath { distances },
        }
    }

    /// Wrap the output of an IBZ reduction into a `KPointSet`. Called by
    /// [`crate::symmetry::kpoints::reduce_kpoints`]; callers outside the
    /// symmetry module should not need this.
    ///
    /// `grid` and `shift` must describe the *parent* full MP grid — they
    /// are retained so density symmetrization can unfold the IBZ back to
    /// the full grid.
    pub fn from_reduced(points: Vec<KPoint>, grid: [u32; 3], shift: KGridShift) -> Self {
        Self {
            points,
            kind: SamplingKind::Irreducible { grid, shift },
        }
    }
}

impl KPointSet {
    /// Underlying slice of k-points.
    pub fn as_slice(&self) -> &[KPoint] {
        &self.points
    }

    /// Iterator over k-points.
    pub fn iter(&self) -> std::slice::Iter<'_, KPoint> {
        self.points.iter()
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Provenance tag. Match on this when the caller needs to know how
    /// the set was generated (e.g. density symmetrization reaches for the
    /// parent MP grid via [`Self::grid_and_shift`]).
    pub fn sampling(&self) -> &SamplingKind {
        &self.kind
    }

    /// Parent MP grid size and shift, for `MonkhorstPack` and
    /// `Irreducible` sets. Returns `None` for band paths.
    pub fn grid_and_shift(&self) -> Option<([u32; 3], KGridShift)> {
        match self.kind {
            SamplingKind::MonkhorstPack { grid, shift } | SamplingKind::Irreducible { grid, shift } => {
                Some((grid, shift))
            },
            SamplingKind::BandPath { .. } => None,
        }
    }

    /// Cumulative along-path distances for band-structure sets. Returns
    /// `None` for MP grids.
    pub fn distances(&self) -> Option<&[f64]> {
        match &self.kind {
            SamplingKind::BandPath { distances } => Some(distances),
            _ => None,
        }
    }

    /// Sum of integration weights. Invariants:
    /// - `MonkhorstPack` and `Irreducible` sets: 1.0 (modulo rounding).
    /// - `BandPath`: 0.0 (band-path points do not participate in BZ integration).
    pub fn total_weight(&self) -> f64 {
        self.points.iter().map(|kp| kp.weight).sum()
    }
}

impl<'a> IntoIterator for &'a KPointSet {
    type IntoIter = std::slice::Iter<'a, KPoint>;
    type Item = &'a KPoint;

    fn into_iter(self) -> Self::IntoIter {
        self.points.iter()
    }
}

impl std::ops::Index<usize> for KPointSet {
    type Output = KPoint;

    fn index(&self, i: usize) -> &KPoint {
        &self.points[i]
    }
}

/// Provenance tag for a [`KPointSet`].
#[derive(Debug, Clone)]
pub enum SamplingKind {
    /// Full unreduced MP grid. Rarely consumed directly — usually the
    /// input to an IBZ reduction.
    MonkhorstPack { grid: [u32; 3], shift: KGridShift },
    /// MP grid reduced to the irreducible BZ under a space group.
    /// `grid` and `shift` describe the *parent* full grid; density
    /// symmetrization needs them to unfold the IBZ contribution.
    Irreducible { grid: [u32; 3], shift: KGridShift },
    /// Piecewise-linear band path through high-symmetry points.
    /// `distances` has the same length as the point list and gives
    /// cumulative reciprocal-space distance (1/Å) along the path.
    BandPath { distances: Vec<f64> },
}

/// Shift convention for a Monkhorst-Pack k-grid.
///
/// Selects between the two common formulas for a uniform k-mesh:
///
/// - [`GammaCentered`](KGridShift::GammaCentered) — `f_j = (i_j − 1)/N_j` (half-shift flags `0 0
///   0`). For N=4 this yields `{0, 1/4, 1/2, 3/4}`, i.e. **includes** Γ and the high-symmetry
///   BZ-boundary points (X, L on FCC).
/// - [`MP1976`](KGridShift::MP1976) — `f_j = (2·i_j − N_j + 1)/(2·N_j)` (the original Monkhorst &
///   Pack, *Phys. Rev. B* **13**, 5188 (1976) Eq. 4; half-shift flags `1 1 1`). For N=4 this yields
///   `{−3/8, −1/8, 1/8, 3/8}` — **no** point at Γ, symmetric about Γ.
/// - [`Custom`](KGridShift::Custom) — per-axis half-shift flags `k_α ∈ {0, 1}` giving `f_j = (i_j −
///   1)/N_j + k_j/(2·N_j)`. `Custom([0,0,0])` is identical to `GammaCentered`.
///
/// For odd N both `GammaCentered` and `MP1976` produce the same mesh
/// modulo a cyclic permutation; for even N they are genuinely distinct
/// physical samples of the BZ.
///
/// Default: `GammaCentered` — the mesh includes Γ and the
/// high-symmetry BZ-boundary points.
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
    /// Γ-centered grid (half-shift flags `0 0 0`).
    #[default]
    GammaCentered,
    /// Original Monkhorst & Pack 1976 shifted grid (half-shift flags `1 1 1`).
    #[serde(rename = "mp1976", alias = "mp_1976", alias = "shifted")]
    MP1976,
    /// Per-axis half-shift flags.
    Custom([u32; 3]),
}

impl KGridShift {
    /// Return the per-axis half-shift integers `k_α ∈ {0, 1}` that
    /// parameterise this shift.
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
/// Formula: `f_α = (i_α)/N_α + k_α/(2·N_α)` for `i_α ∈ {0, …, N_α−1}` and
/// `k_α ∈ {0, 1}`. The raw result lies in `[0, 1)` for `k_α = 0` and in
/// `[1/(2N), 1)` for `k_α = 1`; this function wraps into the first
/// Brillouin zone `[-1/2, 1/2)` so the SCF sees k-points in canonical
/// form. (Mathematically k and k+G give identical physics, but at finite
/// `ecut` the shared plane-wave basis adapts unequally to them — wrapping
/// keeps kinetic energies minimized.)
pub fn mp_fractional_coord(i1: u32, i2: u32, i3: u32, grid: [u32; 3], shift: KGridShift) -> [f64; 3] {
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
pub fn monkhorst_pack(n1: u32, n2: u32, n3: u32, shift: KGridShift, lattice: &Lattice) -> Vec<KPoint> {
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
                kpoints.push(KPoint { k, weight, label: None });
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
/// `segments` is a list of (label, fractional_coords) pairs defining the path
/// vertices. `npoints_per_segment` controls the density of sampling between
/// each pair. Returns k-points with cumulative distance for plotting.
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

        let k_start = start.frac[0] * recip.a + start.frac[1] * recip.b + start.frac[2] * recip.c;
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
mod tests {
    use approx::relative_eq;

    use super::*;

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
            assert!((ka.k - kb.k).norm() < 1e-14, "Custom([1,1,1]) should match MP1976");
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

    // ─── KPointSet coverage ──────────────────────────────────────────

    #[test]
    fn kpointset_gamma_only_is_single_origin_weight_one() {
        let set = KPointSet::gamma_only();
        assert_eq!(set.len(), 1);
        assert!(!set.is_empty());
        assert!(set[0].k.norm() < 1e-14);
        assert!(relative_eq!(set[0].weight, 1.0, epsilon = 1e-14));
        assert!(relative_eq!(set.total_weight(), 1.0, epsilon = 1e-14));
        assert_eq!(set[0].label.as_deref(), Some("Γ"));
        assert!(set.distances().is_none());
        assert_eq!(set.grid_and_shift(), Some(([1, 1, 1], KGridShift::GammaCentered)));
    }

    #[test]
    fn kpointset_monkhorst_pack_matches_free_function() {
        let lattice = si_lattice();
        let set = KPointSet::monkhorst_pack([4, 4, 4], KGridShift::GammaCentered, &lattice);
        let raw = monkhorst_pack(4, 4, 4, KGridShift::GammaCentered, &lattice);

        assert_eq!(set.len(), raw.len());
        for (a, b) in set.iter().zip(raw.iter()) {
            assert!((a.k - b.k).norm() < 1e-14);
            assert!(relative_eq!(a.weight, b.weight, epsilon = 1e-14));
        }
        assert!(matches!(
            set.sampling(),
            SamplingKind::MonkhorstPack {
                grid: [4, 4, 4],
                shift: KGridShift::GammaCentered
            }
        ));
        assert!(relative_eq!(set.total_weight(), 1.0, epsilon = 1e-12));
    }

    #[test]
    fn kpointset_band_path_stores_distances() {
        let lattice = si_lattice();
        let segments = vec![
            HighSymPoint {
                label: "Γ".into(),
                frac: [0.0, 0.0, 0.0],
            },
            HighSymPoint {
                label: "X".into(),
                frac: [0.5, 0.0, 0.5],
            },
        ];
        let set = KPointSet::band_path(&segments, 10, &lattice);

        assert_eq!(set.len(), 11); // 10 + endpoint on the last segment
        let dists = set.distances().expect("band path must expose distances");
        assert_eq!(dists.len(), set.len());
        assert!(dists[0].abs() < 1e-14);
        assert!(dists.windows(2).all(|w| w[1] >= w[0] - 1e-12));

        assert!(set.grid_and_shift().is_none());
        assert!(matches!(set.sampling(), SamplingKind::BandPath { .. }));
        assert!(set.total_weight().abs() < 1e-14);
    }

    #[test]
    fn kpointset_from_reduced_tags_irreducible_and_preserves_parent_grid() {
        let lattice = si_lattice();
        // Fake "reduced" input: three arbitrary points with hand-picked weights
        // summing to 1.0. We're not testing the reduction math here — only
        // that `from_reduced` wraps the data with the right provenance.
        let pts = vec![
            KPoint {
                k: Vector3::new(0.1, 0.0, 0.0),
                weight: 0.25,
                label: None,
            },
            KPoint {
                k: Vector3::new(0.2, 0.1, 0.0),
                weight: 0.5,
                label: None,
            },
            KPoint {
                k: Vector3::new(0.3, 0.2, 0.1),
                weight: 0.25,
                label: None,
            },
        ];
        let set = KPointSet::from_reduced(pts, [4, 4, 4], KGridShift::GammaCentered);

        assert_eq!(set.len(), 3);
        assert!(matches!(
            set.sampling(),
            SamplingKind::Irreducible {
                grid: [4, 4, 4],
                shift: KGridShift::GammaCentered
            }
        ));
        assert_eq!(set.grid_and_shift(), Some(([4, 4, 4], KGridShift::GammaCentered)));
        assert!(relative_eq!(set.total_weight(), 1.0, epsilon = 1e-12));
        let _ = lattice; // silence unused if future edits drop the binding
    }

    #[test]
    fn kpointset_iterates_and_indexes() {
        let set = KPointSet::monkhorst_pack([2, 2, 2], KGridShift::GammaCentered, &si_lattice());
        let via_index: Vec<_> = (0..set.len()).map(|i| set[i].k).collect();
        let via_iter: Vec<_> = (&set).into_iter().map(|kp| kp.k).collect();
        assert_eq!(via_index.len(), via_iter.len());
        for (a, b) in via_index.iter().zip(via_iter.iter()) {
            assert!((a - b).norm() < 1e-14);
        }
    }
}
