//! Charge density symmetrization on the FFT grid.
//!
//! Two forms are implemented:
//!
//! * [`symmetrize_density`] — real-space form
//!   `ρ_sym(r) = (1/N_ops) Σ_S ρ(S⁻¹ r)`, implemented by rotating each
//!   grid point and rounding to the nearest neighbour. Exact only when
//!   `τ_{S,i} · n_i ∈ ℤ` for every operation S and axis i — i.e. the
//!   fractional translation lands on an integer grid point.
//!
//! * [`symmetrize_density_g`] — G-space form (preferred for SCF)
//!   `ρ_sym(G) = (1/N_ops) Σ_S exp(-i G · τ_S) · ρ(R_S⁻¹ G)`,
//!   exact for any fractional translation because the phase factor is
//!   analytic. Matches QE 7.5 `PW/src/symme.f90::sym_rho_serial`.
//!
//! The real-space form silently smears density across the wrong grid
//! points for non-symmorphic space groups whose τ doesn't land on the
//! grid (e.g. Fd-3m with τ=(¼,¼,¼) on an 18³ grid). Always use the
//! G-space form in SCF; keep the real-space form for direct-grid unit
//! tests where translation exactness is guaranteed by construction.
//!
//! See `proposals/completed/PCFX-symmetrize-rho-g-space.md`.

use num_complex::Complex64;

use super::SymmetryInfo;
use crate::fft::FFT3D;

/// Map a fractional coordinate to the nearest grid index in [0, n).
///
/// Uses nearest-integer (round) mapping, matching QE's `nint()` convention.
/// Handles negative coordinates and periodic wrapping via double-modulo.
fn frac_to_grid_idx(frac: f64, n: usize) -> usize {
    let ni = n as i64;
    let idx = (frac * n as f64).round() as i64;
    ((idx % ni) + ni) as usize % n
}

/// Symmetrize a real-space charge density on the FFT grid.
///
/// For each symmetry operation S = {R|τ}, the inverse S⁻¹ = {R⁻¹|-R⁻¹τ}
/// maps each grid point to another grid point (for compatible grids).
/// The symmetrized density is the average over all operations.
///
/// **Important**: The FFT grid dimensions must be compatible with the symmetry
/// operations. Use [`check_grid_compatibility`] before calling this.
pub fn symmetrize_density(rho: &mut [f64], dims: [usize; 3], symmetry: &SymmetryInfo) {
    let [nx, ny, nz] = dims;
    let n_grid = nx * ny * nz;
    assert_eq!(rho.len(), n_grid);

    if symmetry.n_ops <= 1 {
        return; // nothing to symmetrize
    }

    let n_ops = symmetry.n_ops as f64;
    let rho_orig = rho.to_vec();

    // Zero out and accumulate
    for v in rho.iter_mut() {
        *v = 0.0;
    }

    for op in &symmetry.operations {
        let s_inv = op.inverse();
        let r_inv = s_inv.rotation;
        let tau_inv = s_inv.translation;

        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    // Fractional coordinates of this grid point
                    let f = [
                        ix as f64 / nx as f64,
                        iy as f64 / ny as f64,
                        iz as f64 / nz as f64,
                    ];

                    // Apply S⁻¹: f' = R⁻¹·f + τ_inv
                    let fp = [
                        r_inv[0][0] as f64 * f[0]
                            + r_inv[0][1] as f64 * f[1]
                            + r_inv[0][2] as f64 * f[2]
                            + tau_inv[0],
                        r_inv[1][0] as f64 * f[0]
                            + r_inv[1][1] as f64 * f[1]
                            + r_inv[1][2] as f64 * f[2]
                            + tau_inv[1],
                        r_inv[2][0] as f64 * f[0]
                            + r_inv[2][1] as f64 * f[1]
                            + r_inv[2][2] as f64 * f[2]
                            + tau_inv[2],
                    ];

                    // Map to grid indices with periodic boundary conditions
                    let jx = frac_to_grid_idx(fp[0], nx);
                    let jy = frac_to_grid_idx(fp[1], ny);
                    let jz = frac_to_grid_idx(fp[2], nz);

                    let src_idx = jx * ny * nz + jy * nz + jz;
                    let dst_idx = ix * ny * nz + iy * nz + iz;
                    rho[dst_idx] += rho_orig[src_idx];
                }
            }
        }
    }

    // Normalize by number of operations
    for v in rho.iter_mut() {
        *v /= n_ops;
    }
}

/// Check if the FFT grid dimensions are compatible with all symmetry operations.
///
/// For integer rotation R in fractional coords, the grid point (ix, iy, iz) maps
/// to another exact grid point if and only if R_{ij} × n_j ≡ 0 (mod n_i) for all i, j.
///
/// Returns true if all operations are compatible.
#[must_use]
pub fn check_grid_compatibility(dims: [usize; 3], symmetry: &SymmetryInfo) -> bool {
    let ns = [dims[0] as i32, dims[1] as i32, dims[2] as i32];
    for op in &symmetry.operations {
        let r = &op.rotation;
        for i in 0..3 {
            for j in 0..3 {
                if (r[i][j] * ns[j]) % ns[i] != 0 {
                    return false;
                }
            }
        }
    }
    true
}

/// Find the smallest FFT-friendly grid dimensions compatible with the symmetry.
///
/// Starts from the given minimum dimensions and increases until compatibility
/// is achieved. Returns adjusted dimensions.
#[must_use]
pub fn compatible_grid_dims(min_dims: [usize; 3], symmetry: &SymmetryInfo) -> [usize; 3] {
    // For cubic symmetry, making all dimensions equal is usually sufficient
    // SAFETY: min_dims is [usize; 3], always has 3 elements -- max() cannot be None.
    let max_dim = *min_dims.iter().max().expect("BUG: empty fixed-size array");
    let mut dims = [max_dim; 3];

    // Try increasing until compatible
    for _ in 0..100 {
        let candidate = [
            crate::fft::fft_grid_size(dims[0] as i32 / 2),
            crate::fft::fft_grid_size(dims[1] as i32 / 2),
            crate::fft::fft_grid_size(dims[2] as i32 / 2),
        ];
        // Make all equal to the max for safety
        // SAFETY: candidate is [usize; 3], always has 3 elements.
        let m = *candidate.iter().max().expect("BUG: empty fixed-size array");
        let candidate = [m, m, m];
        if check_grid_compatibility(candidate, symmetry) {
            return candidate;
        }
        dims = [dims[0] + 1, dims[1] + 1, dims[2] + 1];
    }

    // Fallback: just use the input dims (may not be compatible)
    min_dims
}

/// Apply an integer rotation matrix (Miller-index basis) to a Miller index triple.
///
/// Fractional-basis rotations R are integer matrices with det = ±1, so
/// Miller indices transform as integer vectors: `n' = R · n`.
#[inline]
fn rotate_miller(r: &[[i32; 3]; 3], n: [i32; 3]) -> [i32; 3] {
    [
        r[0][0] * n[0] + r[0][1] * n[1] + r[0][2] * n[2],
        r[1][0] * n[0] + r[1][1] * n[1] + r[1][2] * n[2],
        r[2][0] * n[0] + r[2][1] * n[1] + r[2][2] * n[2],
    ]
}

/// Transpose of an integer 3×3 rotation matrix.
#[inline]
fn transpose_rotation(r: &[[i32; 3]; 3]) -> [[i32; 3]; 3] {
    [
        [r[0][0], r[1][0], r[2][0]],
        [r[0][1], r[1][1], r[2][1]],
        [r[0][2], r[1][2], r[2][2]],
    ]
}

/// Map signed Miller indices to the FFT flat index with periodic wrap.
///
/// Mirrors `scf::grid::miller_to_idx` so the G-space symmetrizer agrees
/// with the rest of the codebase on how negative `n_i` aliases to
/// `n_i + N_i`. Kept local because `scf::grid` is `pub(crate)` and
/// `symmetry::density` lives in a different module subtree.
#[inline]
fn miller_to_flat(dims: [usize; 3], n: [i32; 3]) -> usize {
    let wrap = |ni: i32, dim: usize| -> usize {
        let d = dim as i32;
        (((ni % d) + d) as usize) % dim
    };
    let i1 = wrap(n[0], dims[0]);
    let i2 = wrap(n[1], dims[1]);
    let i3 = wrap(n[2], dims[2]);
    i1 * dims[1] * dims[2] + i2 * dims[2] + i3
}

/// Recover the signed Miller indices (n1, n2, n3) of the G-vector stored at
/// flat FFT index `idx`. For a dimension `N`, positive frequencies occupy
/// indices `0..=N/2` and negative frequencies occupy `N/2+1..N`, matching
/// the standard FFT layout and `scf::grid::g_vector_at_dims`.
#[inline]
fn flat_to_miller(dims: [usize; 3], idx: usize) -> [i32; 3] {
    let [nx, ny, nz] = dims;
    let i1 = idx / (ny * nz);
    let i2 = (idx / nz) % ny;
    let i3 = idx % nz;
    let n1 = if i1 > nx / 2 {
        i1 as i32 - nx as i32
    } else {
        i1 as i32
    };
    let n2 = if i2 > ny / 2 {
        i2 as i32 - ny as i32
    } else {
        i2 as i32
    };
    let n3 = if i3 > nz / 2 {
        i3 as i32 - nz as i32
    } else {
        i3 as i32
    };
    [n1, n2, n3]
}

/// Symmetrize the charge density in G-space via phase factors.
///
/// Computes, on the full FFT grid,
///
/// ```text
/// ρ_sym(G) = (1/N_ops) Σ_S  exp(-i G · τ_S) · ρ(R_S^T · G)
/// ```
///
/// where each symmetry operation `S = {R_S | τ_S}` contributes an
/// analytic phase. Our `SpaceGroupOp::rotation` is the integer rotation
/// that acts on direct-space fractional coords as `r' = R·r + τ`, and
/// G-vectors in the Miller basis rotate as `n → R^T · n` under the
/// induced Fourier action. Unlike the real-space form
/// ([`symmetrize_density`]) this is **exact** for any fractional
/// translation on any sufficiently-band-limited density: a glide of
/// `τ=(¼,¼,¼)` on an 18³ grid incurs no rounding error (the real-space
/// form does, because `18·¼ = 4.5 ∉ ℤ`).
///
/// The formula matches QE 7.5 `PW/src/symme.f90::sym_rho_serial`
/// up to a group-level `S → S⁻¹` relabelling: QE stores its `s(:,:,ns)`
/// as the *transpose* of the direct-space fractional rotation (atoms
/// rotate as `rau = s^T · xau`, see `symm_base.f90:533`), so its
/// `s(:,:,invs(ns))·g0 = R^{-T}·g0` matches our `R^T·n` under the
/// relabelling, and the phase conventions are equivalent.
///
/// Work per call (leading order): one forward FFT, one inverse FFT,
/// plus `N_grid · N_ops` complex multiplies/accumulates. At Si 18³ with
/// 48 operations that is ~280k complex MACs plus two FFTs, fully
/// negligible next to eigen-decomposition cost.
///
/// # Band-limitation requirement
///
/// Because DFT coefficients are periodic in the Miller index with
/// period `N`, a rotated Miller `R^T · m` that falls outside the
/// representable range `[−N/2, N/2)` aliases to a different canonical
/// representative. In general this aliasing introduces a residual
/// phase `exp(−i·2π·N·δ·τ)` (with `δ` an integer vector) that is
/// only a group-unit when `N·τ ∈ ℤ` — the same grid-compatibility
/// condition as the real-space form. To avoid this, the input density
/// must be band-limited to `|G|² < |G|²_Nyquist / |R|_op,max²` so that
/// rotations never wrap. In the pwdft-rs SCF this is automatic: the
/// density `ρ(r) = Σ_nk |ψ_nk(r)|²` has Fourier support on
/// `|G|² ≤ 4 · ecutwfc = ecutrho`, and the FFT grid is chosen via
/// `scf::grid::FftGrid::new` (with `ecutrho_ratio ≥ 4`) to be strictly
/// larger, so the density cutoff sits inside the representable Miller
/// range with enough margin for the largest rotation coefficient
/// (cubic groups have |R| ≤ 3).
///
/// # Parameters
///
/// * `rho_r` — density on the FFT grid in row-major `(ix, iy, iz)` order,
///   matching `FftGrid::total_size()`. Modified in place.
/// * `dims` — FFT grid dimensions `[nx, ny, nz]`.
/// * `fft` — a reusable [`FFT3D`] instance sized for `dims`.
/// * `symmetry` — space-group operations. If `n_ops <= 1` the call is a
///   no-op (identity-only group produces the original density).
pub fn symmetrize_density_g(
    rho_r: &mut [f64],
    dims: [usize; 3],
    fft: &mut FFT3D,
    symmetry: &SymmetryInfo,
) {
    let n_grid = dims[0] * dims[1] * dims[2];
    assert_eq!(rho_r.len(), n_grid);
    assert_eq!(fft.dims(), dims, "FFT3D dims must match density dims");

    if symmetry.n_ops <= 1 {
        return;
    }

    // 1. Forward FFT: ρ(r) → unnormalized ρ̃(G).
    //    We keep the unnormalized coefficients throughout; normalization
    //    by 1/N cancels between forward and inverse FFT.
    let mut rho_g: Vec<Complex64> = rho_r.iter().map(|&r| Complex64::new(r, 0.0)).collect();
    fft.forward(&mut rho_g);

    // 2. For every destination Miller index n_dst, accumulate contributions
    //    from all operations.
    //
    //    Convention. Our `SpaceGroupOp::rotation` (= `op.rotation`) is the
    //    integer matrix `R` that acts on **direct-space fractional
    //    coordinates** as `r' = R·r + τ`. Under the pull-back action
    //    `(S·ρ)(r) = ρ(S⁻¹ r)`, which is the standard left group action
    //    on functions, the Fourier coefficient of `S·ρ` is
    //
    //        (S·ρ̃)(n) = e^{-i·2π·n·τ_S} · ρ̃(R_S^T · n)
    //
    //    (derivation: substitute `r = S r'` in `ρ̃(G) = ∫ e^{-iG·r} ρ(r) dr`
    //    and pull the translation phase out). The group-average projector
    //    is therefore
    //
    //        ρ̃_sym(n) = (1/N_ops) Σ_S e^{-i·2π·n·τ_S} · ρ̃(R_S^T · n)
    //
    //    For every destination Miller `n` the source is `R_S^T · n` and
    //    the phase uses the *destination* Miller `n` dotted with τ_S.
    //    The map `S → (S·)` is a left group homomorphism — one can
    //    check that `(S₁·S₂)·ρ̃ = S₁·(S₂·ρ̃)` because the Seitz
    //    composition `τ_{S₁·S₂} = τ_{S₁} + R_{S₁}·τ_{S₂}` cancels the
    //    `(R_{S₁}^T n)·τ_{S₂} = n·(R_{S₁} τ_{S₂})` cross term exactly.
    //    This guarantees `P² = P`.
    //
    //    Matches QE 7.5 `PW/src/symme.f90::sym_rho_serial` after the
    //    group-level relabelling `S → S⁻¹`: QE stores its `s(:,:,ns)` as
    //    the transpose of our direct-space `R` (atoms rotate as
    //    `rau = s^T xau`, see `symm_base.f90:533`), so QE's forward
    //    accumulation `sg = s(:,:,invs(ns)) · g0 = R^{-T}·g0` and its
    //    phase `exp(-i·2π·sg·τ_ns)` correspond to the `S → S⁻¹`
    //    relabelling of the projector above.
    let n_ops = symmetry.n_ops;
    let two_pi = std::f64::consts::TAU;
    let inv_n_ops = 1.0 / n_ops as f64;
    let mut rho_g_sym: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n_grid];

    let r_ts: Vec<[[i32; 3]; 3]> = symmetry
        .operations
        .iter()
        .map(|op| transpose_rotation(&op.rotation))
        .collect();

    for (idx_dst, dst_slot) in rho_g_sym.iter_mut().enumerate() {
        let n_dst = flat_to_miller(dims, idx_dst);
        let mut acc = Complex64::new(0.0, 0.0);
        for (op, r_t) in symmetry.operations.iter().zip(r_ts.iter()) {
            // Source Miller: n_src = R_direct^T · n_dst.
            let n_src = rotate_miller(r_t, n_dst);
            let idx_src = miller_to_flat(dims, n_src);
            let rho_src = rho_g[idx_src];

            // Phase: exp(-i · 2π · (n_dst · τ)). The phase argument
            // uses n_dst (destination), not n_src — required for
            // `P² = P` by the Seitz composition identity
            // `m·τ_{S₁·S₂} = m·τ_{S₁} + (R_{S₁}^T m)·τ_{S₂}`.
            let arg = two_pi
                * (n_dst[0] as f64 * op.translation[0]
                    + n_dst[1] as f64 * op.translation[1]
                    + n_dst[2] as f64 * op.translation[2]);
            let (s, c) = arg.sin_cos();
            let phase = Complex64::new(c, -s);

            acc += phase * rho_src;
        }
        *dst_slot = acc * inv_n_ops;
    }

    // 4. Inverse FFT (normalized) back to real space. The real part of the
    //    result is the symmetrized density; the imaginary part is noise at
    //    the ~1e-16 · ||ρ||_∞ level for a true real-valued input density.
    fft.inverse_normalized(&mut rho_g_sym);
    for (out, c) in rho_r.iter_mut().zip(rho_g_sym.iter()) {
        *out = c.re;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::{Atom, Crystal, Lattice};
    use approx::relative_eq;
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

    #[test]
    fn test_frac_to_grid_idx_basics() {
        assert_eq!(frac_to_grid_idx(0.0, 10), 0);
        assert_eq!(frac_to_grid_idx(0.5, 10), 5);
        assert_eq!(frac_to_grid_idx(1.0, 10), 0); // wraps
        assert_eq!(frac_to_grid_idx(0.3, 10), 3);
        assert_eq!(frac_to_grid_idx(0.95, 10), 10 % 10); // rounds to 10, wraps to 0
    }

    #[test]
    fn test_frac_to_grid_idx_negative() {
        assert_eq!(frac_to_grid_idx(-0.1, 10), 9); // -1 + 10 = 9
        assert_eq!(frac_to_grid_idx(-0.5, 10), 5); // -5 + 10 = 5
        assert_eq!(frac_to_grid_idx(-1.0, 10), 0); // -10 + 10 = 0
        // -0.05 * 10 = -0.5: round(-0.5) is implementation-defined
        let idx = frac_to_grid_idx(-0.05, 10);
        assert!(idx == 0 || idx == 9, "round(-0.5) should give 0 or 9, got {idx}");
    }

    #[test]
    fn test_frac_to_grid_idx_boundary() {
        // At half-integer: round(0.5) = 0 or 1 depending on banker's rounding
        // Either is acceptable as long as it's consistent
        let idx = frac_to_grid_idx(0.05, 10); // 0.05 * 10 = 0.5
        assert!(idx == 0 || idx == 1, "half-integer should map to 0 or 1, got {idx}");
    }

    #[test]
    fn test_symmetrize_preserves_integral() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12];
        let n = dims[0] * dims[1] * dims[2];

        let mut rho: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin().abs() + 0.1).collect();
        let integral_before: f64 = rho.iter().sum();

        symmetrize_density(&mut rho, dims, &symmetry);

        let integral_after: f64 = rho.iter().sum();
        assert!(
            relative_eq!(integral_before, integral_after, epsilon = 1e-10),
            "integral changed: {integral_before} → {integral_after}"
        );
    }

    #[test]
    fn test_symmetrize_uniform_unchanged() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12];
        let n = dims[0] * dims[1] * dims[2];

        let mut rho = vec![1.0; n];
        symmetrize_density(&mut rho, dims, &symmetry);
        for &v in &rho {
            assert!(
                relative_eq!(v, 1.0, epsilon = 1e-14),
                "uniform density changed to {v}"
            );
        }
    }

    #[test]
    fn test_symmetrize_single_point_orbit() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12];
        let n = dims[0] * dims[1] * dims[2];

        let mut rho = vec![0.0; n];
        let test_idx = dims[1] * dims[2] + 2 * dims[2] + 3; // (1,2,3)
        rho[test_idx] = 48.0;

        symmetrize_density(&mut rho, dims, &symmetry);

        let nonzero: Vec<f64> = rho.iter().filter(|&&v| v > 1e-10).copied().collect();
        assert!(!nonzero.is_empty());
        let ref_val = nonzero[0];
        for &v in &nonzero {
            assert!(
                relative_eq!(v, ref_val, epsilon = 1e-10),
                "orbit points not equal: {v} vs {ref_val}"
            );
        }
    }

    #[test]
    fn test_grid_compatibility_cubic() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);

        assert!(check_grid_compatibility([12, 12, 12], &symmetry));
        assert!(check_grid_compatibility([18, 18, 18], &symmetry));
        assert!(check_grid_compatibility([20, 20, 20], &symmetry));
        assert!(!check_grid_compatibility([12, 12, 15], &symmetry));
    }

    #[test]
    fn test_idempotent() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12];
        let n = dims[0] * dims[1] * dims[2];

        let mut rho: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin().abs()).collect();
        symmetrize_density(&mut rho, dims, &symmetry);

        let rho_once = rho.clone();
        symmetrize_density(&mut rho, dims, &symmetry);

        for (a, b) in rho.iter().zip(rho_once.iter()) {
            assert!(
                (a - b).abs() < 1e-12,
                "symmetrization not idempotent: {a} vs {b}"
            );
        }
    }

    // ------------------------------------------------------------------
    // G-space symmetrization tests (PCFX)
    // ------------------------------------------------------------------

    #[test]
    fn test_symmetrize_g_preserves_integral_and_electron_count() {
        // The DC Fourier coefficient is the volume-averaged density, so
        // ρ(G=0) must be invariant under symmetrization and the real-space
        // integral (= ρ(G=0)·N) must match before and after.
        //
        // Si Fd-3m on an 18³ grid is the canonical failing case for
        // real-space symmetrization — the fractional translation
        // τ=(¼,¼,¼) lands between grid points since 18·¼=4.5. The G-space
        // form must recover both invariants exactly regardless.
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [18, 18, 18];
        let n = dims[0] * dims[1] * dims[2];
        let mut fft = FFT3D::new(dims[0], dims[1], dims[2]);

        let mut rho: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin().abs() + 0.1).collect();
        let sum_before: f64 = rho.iter().sum();

        symmetrize_density_g(&mut rho, dims, &mut fft, &symmetry);

        let sum_after: f64 = rho.iter().sum();
        assert!(
            relative_eq!(sum_before, sum_after, epsilon = 1e-9),
            "G-space symmetrization changed total integral: {sum_before} → {sum_after}"
        );
    }

    #[test]
    fn test_symmetrize_g_uniform_unchanged() {
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [18, 18, 18];
        let n = dims[0] * dims[1] * dims[2];
        let mut fft = FFT3D::new(dims[0], dims[1], dims[2]);

        let mut rho = vec![1.0; n];
        symmetrize_density_g(&mut rho, dims, &mut fft, &symmetry);
        for &v in &rho {
            assert!(
                relative_eq!(v, 1.0, epsilon = 1e-12),
                "uniform density changed to {v}"
            );
        }
    }

    /// Build a band-limited test density with Fourier support in the
    /// open ball `|m_i| ≤ m_max`, sized so cubic rotations (max row sum
    /// 3 for Fd-3m) keep the orbit inside `(−N/2, N/2)` even for a
    /// non-symmorphic grid — the regime in which the G-space projector
    /// is guaranteed to satisfy P² = P.
    fn band_limited_density(dims: [usize; 3], m_max: i32) -> Vec<f64> {
        let [nx, ny, nz] = dims;
        let n = nx * ny * nz;
        let mut rho = vec![0.0_f64; n];
        // Sum a handful of cosine modes. Real input preserves
        // Hermitian symmetry of the FFT without any extra effort.
        let modes = [
            ([0_i32, 0, 0], 1.0_f64),
            ([1, 0, 0], 0.17),
            ([0, 1, 0], 0.11),
            ([0, 0, 1], 0.13),
            ([1, 1, 0], 0.07),
            ([1, 0, 1], 0.05),
            ([1, 1, 1], 0.04),
            ([2, 0, 0], 0.03),
        ];
        assert!(
            modes.iter().all(|(k, _)| k[0].abs() <= m_max
                && k[1].abs() <= m_max
                && k[2].abs() <= m_max),
            "band-limited helper: a mode exceeds |m|≤{m_max}"
        );
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let fx = ix as f64 / nx as f64;
                    let fy = iy as f64 / ny as f64;
                    let fz = iz as f64 / nz as f64;
                    let mut s = 0.0;
                    for (k, amp) in &modes {
                        s += amp
                            * (std::f64::consts::TAU
                                * (k[0] as f64 * fx + k[1] as f64 * fy + k[2] as f64 * fz))
                                .cos();
                    }
                    rho[ix * ny * nz + iy * nz + iz] = s.abs() + 0.5;
                }
            }
        }
        rho
    }

    #[test]
    fn test_symmetrize_g_idempotent_on_18_cubed_bandlimited() {
        // P² = P for a band-limited density on Fd-3m's canonical
        // "incompatible" grid (18 is not divisible by 4, so the
        // real-space form smears τ=(¼,¼,¼) between neighbours). This
        // is the PCFX regression guard: the failing case for the
        // real-space form is a passing case for the G-space form,
        // because the phase factor handles τ analytically.
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [18, 18, 18];
        let mut fft = FFT3D::new(dims[0], dims[1], dims[2]);

        // Band-limit to |m|≤2: max orbit is 3·2 = 6, well inside [-9, 9]
        // for Nyquist 9 on an 18³ grid.
        let mut rho = band_limited_density(dims, 2);

        symmetrize_density_g(&mut rho, dims, &mut fft, &symmetry);
        let rho_once = rho.clone();
        symmetrize_density_g(&mut rho, dims, &mut fft, &symmetry);

        let max_diff = rho
            .iter()
            .zip(rho_once.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        let max_rho = rho_once.iter().copied().fold(0.0f64, f64::max);
        assert!(
            max_diff < 1e-12 * max_rho.max(1.0),
            "G-space symmetrization not idempotent: max |Δ| = {max_diff:.2e}, \
             max |ρ| = {max_rho:.2e}"
        );
    }

    #[test]
    fn test_symmetrize_g_matches_real_space_on_compatible_grid() {
        // When the grid captures every fractional translation on an
        // integer grid point, the real-space and G-space forms produce
        // the same density. For Si Fd-3m with τ=(¼,¼,¼), a 12³ grid
        // (12 divisible by 4) satisfies this. Any disagreement would
        // indicate a convention mismatch (wrong sign of τ, wrong
        // direction of R, wrong FFT layout, etc.).
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12];
        let n = dims[0] * dims[1] * dims[2];
        let mut fft = FFT3D::new(dims[0], dims[1], dims[2]);

        let rho0: Vec<f64> = (0..n)
            .map(|i| (i as f64 * 0.23).sin().abs() + 0.05)
            .collect();

        let mut rho_real = rho0.clone();
        symmetrize_density(&mut rho_real, dims, &symmetry);

        let mut rho_g = rho0;
        symmetrize_density_g(&mut rho_g, dims, &mut fft, &symmetry);

        let max_diff = rho_real
            .iter()
            .zip(rho_g.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        let max_rho = rho_real.iter().copied().fold(0.0f64, f64::max);
        assert!(
            max_diff < 1e-10 * max_rho.max(1.0),
            "real-space and G-space disagree on compatible grid: max |Δ| = {max_diff:.2e}"
        );
    }

    #[test]
    fn test_symmetrize_g_identity_only_is_noop() {
        // Identity-only group (SOPT semantics for "symmetry off"): the
        // call must leave the density bit-identical because of the
        // `n_ops <= 1` short-circuit.
        let s = crate::symmetry::SymmetryInfo::identity_only();
        let dims = [8, 8, 8];
        let n = dims[0] * dims[1] * dims[2];
        let mut fft = FFT3D::new(dims[0], dims[1], dims[2]);

        let original: Vec<f64> = (0..n).map(|i| (i as f64 * 0.13).sin() + 1.0).collect();
        let mut rho = original.clone();
        symmetrize_density_g(&mut rho, dims, &mut fft, &s);
        for (a, b) in original.iter().zip(rho.iter()) {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "identity-only G-space symmetrization must leave density bit-identical"
            );
        }
    }

    #[test]
    fn test_symmetrize_g_projects_pre_symmetric_density() {
        // A density that is already space-group symmetric — here produced
        // by the *real-space* symmetrizer on a compatible 12³ grid where
        // τ=(¼,¼,¼) lands exactly on grid points — must be a fixed point
        // of `symmetrize_density_g`. This is the projector property
        // `P · (P · ρ) = P · ρ` in a form that ties the two
        // implementations together: whatever G-space convention we use
        // for R and τ, it must fix every function that the real-space
        // form on a compatible grid produces. If this test fails, the
        // G-space formula's rotation / phase sign is wrong.
        let crystal = si_fcc();
        let symmetry = crate::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
        let dims = [12, 12, 12];
        let n = dims[0] * dims[1] * dims[2];
        let mut fft = FFT3D::new(dims[0], dims[1], dims[2]);

        let mut rho: Vec<f64> = (0..n)
            .map(|i| (i as f64 * 0.17).sin().abs() + 0.01)
            .collect();
        symmetrize_density(&mut rho, dims, &symmetry);
        let rho_pre = rho.clone();

        symmetrize_density_g(&mut rho, dims, &mut fft, &symmetry);

        let max_diff = rho
            .iter()
            .zip(rho_pre.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(
            max_diff < 1e-11,
            "G-space symmetrization altered a pre-symmetric density: max |Δ| = {max_diff:.2e}"
        );
    }
}
