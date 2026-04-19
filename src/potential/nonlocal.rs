//! Non-local pseudopotential in the Kleinman-Bylander separable form.
//!
//! V_NL = Σ_{atom} Σ_{i,j} |β_i⟩ D_{ij} ⟨β_j| × S_atom
//!
//! In reciprocal space at k-point k:
//! V_NL(k)_{G,G'} = Σ_{atom} S_atom(G-G') × Σ_{i,j} β_i(k+G) D_{ij} β_j*(k+G')
//!
//! where β_i(k+G) is the Fourier transform of the real-space projector:
//! β_{l}(q) = 4π i^l ∫ r·β_l(r) j_l(qr) r dr × Y_lm(q̂)
//!
//! (UPF stores r·β(r), so the integral is ∫ [r·β(r)] j_l(qr) r dr)
//!
//! # Implementation (VNLM, 2026-04-18)
//!
//! The Hamiltonian contribution is assembled via a single BLAS-3 GEMM:
//!
//! ```text
//! H_NL = B · D · B^H
//! ```
//!
//! where the "expanded" projector matrix
//! `B ∈ ℂ^(n_pw × n_channels)` (one channel per (atom, radial projector, m))
//! has entries
//!
//! ```text
//! B[G, α] = (1/√Ω) · e^{−iG·τ_α} · F_{i_α}(|k+G|) · Y_{l_α m_α}(q̂_{k+G})
//! ```
//!
//! and `D` is block-diagonal: only atoms of the same type and channels with
//! the same `(l, m)` contribute. The identity
//! `Σ_m Y_lm(q̂) Y*_lm(q̂') = (2l+1)/(4π) P_l(cos θ)` (spherical-harmonic
//! addition theorem) makes this algebraically equivalent to the compact
//! Legendre form used historically.
//!
//! The previous implementation evaluated the Legendre-folded form with a
//! nested `(ig, jg, i, j)` scalar loop of cost O(n_pw² · n_proj²). The GEMM
//! path is O(n_pw² · n_channels) with a hot kernel dispatched to faer/gemm
//! SIMD — ~10× faster at n_pw = 725 on Apple M2.
//!
//! Phase factor `i^{l_α}` cancels exactly: it appears as `i^{l_α} · (i^{l_β})*`
//! in H_NL, and the block-diagonal structure of D forces `l_α = l_β`, so the
//! combined phase is 1. We omit it from B entirely.

use faer::Mat;
use faer::linalg::matmul::matmul;
use nalgebra::Vector3;
use num_complex::Complex64;
use std::f64::consts::PI;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    error::{PwdftError, Result},
    pseudopotential::PseudopotentialData,
};

/// Precomputed non-local KB projector data at a single k-point.
///
/// On construction we build the expanded complex projector matrix
/// `B ∈ ℂ^(n_pw × n_channels)` and its right-acting counterpart
/// `D_over_omega · B^H ∈ ℂ^(n_channels × n_pw)` so that
/// `add_to_hamiltonian` reduces to a single GEMM.
pub struct NonlocalPotential {
    /// Number of plane waves at this k-point (rows of B).
    n_pw: usize,
    /// Expanded KB projector matrix:
    ///   B[G, α] = (1/√Ω) · exp(−iG·τ_α) · F_{i_α}(|k+G|) · Y_{l_α m_α}(q̂_{k+G})
    /// indexed as B[ig, channel]. Size n_pw × n_channels.
    b: Mat<Complex64>,
    /// Pre-applied operator D · B^H of shape n_channels × n_pw.
    ///
    /// Because `D` is block-diagonal (only same-atom, same-(l,m) entries
    /// survive), we can compute this in O(n_pw · Σ_atom n_proj_l²) once
    /// per k-point and then issue a single GEMM in `add_to_hamiltonian`.
    ///
    /// (The 1/Ω normalization is carried entirely by the 1/√Ω prefactor
    /// in each leg of B, so D enters unscaled.)
    d_bh: Mat<Complex64>,
}

/// Per-element cache used while building the KB projector matrix: the
/// radial form factors F_i(q) at the current k-point's G-vectors, plus
/// the (l, D_ij, n_proj) metadata shared by every atom of that type.
struct TypeCache {
    /// form_factors[proj][g_index] = F_i(|k+G|)
    form_factors: Vec<Vec<f64>>,
    /// angular momentum of each radial projector
    ls: Vec<i32>,
    /// D_ij matrix (row-major, n_proj × n_proj)
    dij: Vec<f64>,
    n_proj: usize,
}

impl NonlocalPotential {
    /// Precompute projector form factors for a given k-point and assemble
    /// the cached `(B, D·B^H)` pair used by [`Self::add_to_hamiltonian`].
    ///
    /// For each projector i with angular momentum l:
    /// F_i(|k+G|) = 4π ∫ [r·β_i(r)] j_l(|k+G|·r) r dr
    ///
    /// The full projector in G-space is:
    /// β_i(k+G) = F_i(|k+G|) × i^l × Y_lm(k̂+G)
    ///
    /// # Errors
    /// Returns `PwdftError::MissingPseudopotential` if any atom type lacks a loaded PP.
    pub fn new(
        crystal: &Crystal,
        basis: &BasisSet,
        k: &Vector3<f64>,
        pseudopotentials: &[&PseudopotentialData],
    ) -> Result<Self> {
        let n_pw = basis.len();
        let omega = crystal.lattice.volume();
        let inv_sqrt_omega = 1.0 / omega.sqrt();

        // Precompute q-vectors, |q| values, and the real-spherical-harmonic
        // table Y_lm(q̂) for every (G, l, m) up to the maximum l of any
        // projector in the calculation.
        let g_vecs = basis.g_vectors();
        let q_vecs: Vec<Vector3<f64>> = g_vecs.iter().map(|g| k + g).collect();
        let q_norms: Vec<f64> = q_vecs.iter().map(|q| q.norm()).collect();

        // Angular momentum `l` is a non-negative quantum number by physics
        // (l ∈ {0, 1, 2, 3} for s/p/d/f). A negative `proj.l` would be a
        // malformed pseudopotential (the UPF parser takes whatever integer
        // is in `angular_momentum=` verbatim). Reaching this code with any
        // negative `l` would silently blow up an allocation via sign-loss;
        // we assert here to turn it into a clean panic.
        let mut lmax: i32 = 0;
        for pp in pseudopotentials {
            for proj in &pp.beta_projectors {
                assert!(
                    proj.l >= 0,
                    "NonlocalPotential::new: beta projector has negative angular momentum l={}",
                    proj.l
                );
                lmax = lmax.max(proj.l);
            }
        }
        #[allow(
            clippy::cast_sign_loss,
            reason = "lmax is the maximum of all beta-projector `l` values, each asserted non-negative above"
        )]
        let ylm_stride = ((lmax + 1) * (lmax + 1)) as usize; // (lmax+1)^2 entries per G
        let mut ylm = vec![0.0_f64; n_pw * ylm_stride];
        for (ig, q) in q_vecs.iter().enumerate() {
            real_sph_harmonics(q, lmax, &mut ylm[ig * ylm_stride..(ig + 1) * ylm_stride]);
        }

        // Precompute per-atom G-phase: exp(-iG·τ).
        // (k·τ cancels between bra and ket in the matrix element, so we omit it.)
        let atom_phases: Vec<Vec<Complex64>> = crystal
            .atoms
            .iter()
            .map(|atom| {
                let tau = atom.cart_position(&crystal.lattice);
                g_vecs
                    .iter()
                    .map(|g| Complex64::cis(-g.dot(&tau)))
                    .collect()
            })
            .collect();

        // Per-atom radial form-factor cache: form_factor[iatom][projector][ig]
        // (keyed by atom so each site reuses its type's F_i(q) table).
        let mut form_factor_by_atom: Vec<Vec<Vec<f64>>> = Vec::with_capacity(crystal.atoms.len());
        // Keep parallel list of (l,) per-atom projector list.
        let mut proj_l_by_atom: Vec<Vec<i32>> = Vec::with_capacity(crystal.atoms.len());
        // Per-atom D_ij table (row-major, n_proj×n_proj).
        let mut dij_by_atom: Vec<Vec<f64>> = Vec::with_capacity(crystal.atoms.len());
        let mut n_proj_by_atom: Vec<usize> = Vec::with_capacity(crystal.atoms.len());

        // Cache F_i(q) per *type* so atoms of the same element reuse the table.
        let mut type_ff: std::collections::HashMap<u32, TypeCache> =
            std::collections::HashMap::new();

        for atom in &crystal.atoms {
            let z = atom.z;
            if let std::collections::hash_map::Entry::Vacant(e) = type_ff.entry(z) {
                let pp = crate::pseudopotential::find_for_atom(z, pseudopotentials)
                    .ok_or_else(|| {
                        PwdftError::MissingPseudopotential(format!(
                            "Z={z} not found in loaded pseudopotentials"
                        ))
                    })?;
                let mut ff: Vec<Vec<f64>> = Vec::with_capacity(pp.beta_projectors.len());
                let mut ls: Vec<i32> = Vec::with_capacity(pp.beta_projectors.len());
                for proj in &pp.beta_projectors {
                    let l = proj.l;
                    ls.push(l);
                    let mut fi = Vec::with_capacity(n_pw);
                    for &q in &q_norms {
                        fi.push(bessel_transform_projector(
                            &pp.r_grid,
                            &pp.rab,
                            &proj.values,
                            l,
                            q,
                        ));
                    }
                    ff.push(fi);
                }
                e.insert(TypeCache {
                    form_factors: ff,
                    ls,
                    dij: pp.dij.clone(),
                    n_proj: pp.n_projectors(),
                });
            }
            let entry = &type_ff[&z];
            form_factor_by_atom.push(entry.form_factors.clone());
            proj_l_by_atom.push(entry.ls.clone());
            dij_by_atom.push(entry.dij.clone());
            n_proj_by_atom.push(entry.n_proj);
        }

        // Compute the channel layout. A "channel" = (atom, radial projector, m)
        // where m runs over 2l+1 real spherical-harmonic components of that
        // projector's angular momentum.
        let mut n_channels: usize = 0;
        // For each atom, record the starting channel index of each radial
        // projector (first m=-l of that projector). We use this to apply D
        // in a per-(atom, l) block.
        let mut atom_channel_starts: Vec<Vec<usize>> = Vec::with_capacity(crystal.atoms.len());
        for ls in &proj_l_by_atom {
            let mut starts = Vec::with_capacity(ls.len());
            for &l in ls {
                starts.push(n_channels);
                #[allow(
                    clippy::cast_sign_loss,
                    reason = "each projector's l asserted non-negative at the start of new()"
                )]
                let l_channels = (2 * l + 1) as usize;
                n_channels += l_channels;
            }
            atom_channel_starts.push(starts);
        }

        // Build the expanded B matrix.
        //   B[G, α] = (1/√Ω) · e^{-iG·τ_a(α)} · F_{i(α)}(|k+G|) · Y_{l(α) m(α)}(q̂_{k+G})
        let mut b: Mat<Complex64> = Mat::zeros(n_pw, n_channels);
        for (iatom, atom_phase) in atom_phases.iter().enumerate() {
            let ls = &proj_l_by_atom[iatom];
            let ff = &form_factor_by_atom[iatom];
            let starts = &atom_channel_starts[iatom];
            for (iproj, &l) in ls.iter().enumerate() {
                #[allow(
                    clippy::cast_sign_loss,
                    reason = "each projector's l asserted non-negative at the start of new()"
                )]
                let l_usize = l as usize;
                let base_lm = l_usize * l_usize; // starting index of (l, m=-l) in ylm row
                let ff_row = &ff[iproj];
                let ch0 = starts[iproj];
                for ig in 0..n_pw {
                    let base = ig * ylm_stride;
                    let phase = atom_phase[ig];
                    let f = ff_row[ig];
                    let scale = inv_sqrt_omega * f;
                    // m index within the (2l+1) block, using the same ordering
                    // as real_sph_harmonics: m = 0, +1, -1, +2, -2, ...
                    for m_off in 0..(2 * l_usize + 1) {
                        let y = ylm[base + base_lm + m_off];
                        b[(ig, ch0 + m_off)] = phase * Complex64::new(scale * y, 0.0);
                    }
                }
            }
        }

        // Build (D/Ω) · B^H directly into a (n_channels × n_pw) matrix.
        //
        // D is block-diagonal: D[(a,i,m), (a',j,m')] nonzero only when
        // a = a', l_i = l_j, m = m'. So for each atom, for each (l_i, l_j)
        // pair with l_i == l_j, for each m in -l..=+l:
        //   DB_H[channel(a,i,m), :] += (D_ij / Ω) · conj(B[:, channel(a,j,m)])
        let mut d_bh: Mat<Complex64> = Mat::zeros(n_channels, n_pw);
        for iatom in 0..crystal.atoms.len() {
            let n_proj = n_proj_by_atom[iatom];
            let dij = &dij_by_atom[iatom];
            let ls = &proj_l_by_atom[iatom];
            let starts = &atom_channel_starts[iatom];
            for i in 0..n_proj {
                for j in 0..n_proj {
                    if ls[i] != ls[j] {
                        continue;
                    }
                    // The 1/Ω normalization is carried by the two √Ω factors
                    // in B (one per projector leg); D enters unscaled.
                    let d_scaled = dij[i * n_proj + j];
                    if d_scaled.abs() < 1e-20 {
                        continue;
                    }
                    #[allow(
                        clippy::cast_sign_loss,
                        reason = "each projector's l asserted non-negative at the start of new()"
                    )]
                    let l_usize = ls[i] as usize;
                    for m_off in 0..(2 * l_usize + 1) {
                        let ci = starts[i] + m_off;
                        let cj = starts[j] + m_off;
                        // DB_H[ci, :] += d_scaled · conj(B[:, cj])
                        let d_c = Complex64::new(d_scaled, 0.0);
                        for ig in 0..n_pw {
                            d_bh[(ci, ig)] += d_c * b[(ig, cj)].conj();
                        }
                    }
                }
            }
        }

        Ok(Self { n_pw, b, d_bh })
    }

    /// Add V_NL to the Hamiltonian matrix at the k-point this potential was
    /// constructed for.
    ///
    /// `H += B · D · B^H`, via a single BLAS-3 GEMM. Memory access is
    /// cache-friendly; on Apple M2 this is ~10× faster than the previous
    /// scalar-loop path at n_pw = 725. (The 1/Ω normalization is already
    /// folded into B at construction time.)
    ///
    /// The `crystal`, `basis`, and `k` arguments are kept for API
    /// compatibility; the heavy lifting was done at construction.
    pub fn add_to_hamiltonian(
        &self,
        h: &mut faer::Mat<Complex64>,
        _crystal: &Crystal,
        _basis: &BasisSet,
        _k: &Vector3<f64>,
    ) {
        debug_assert_eq!(h.nrows(), self.n_pw);
        debug_assert_eq!(h.ncols(), self.n_pw);
        // H += 1 · B · d_bh
        matmul(
            h.as_mut(),
            faer::Accum::Add,
            self.b.as_ref(),
            self.d_bh.as_ref(),
            Complex64::new(1.0, 0.0),
            faer::Par::Seq,
        );
    }
}

/// Real spherical harmonics Y_lm(q̂) for 0 ≤ l ≤ lmax, written into `out`
/// in the order (l, m) = (0,0), (1,0), (1,+1), (1,-1), (2,0), (2,+1), (2,-1),
/// (2,+2), (2,-2), …, matching QE's `ylmr2` layout.
///
/// The algorithm is adapted from QE's `ylmr2_gpu.f90`:
///   1. Compute cos θ, sin θ, φ from the direction of q.
///   2. Build Q(l, m) := sqrt((l−m)! / (l+m)!) · P_l^m(cos θ) for 0 ≤ m ≤ l
///      using the standard associated-Legendre recurrence.
///   3. Multiply by the normalization and cos(mφ) / sin(mφ) for real harmonics.
///
/// At |q| = 0, q̂ is undefined; we set Y_lm = 0 for l > 0 (physical: F_i(0) = 0
/// for l > 0 so the product is zero anyway) and Y_00 = 1/√(4π).
#[allow(
    clippy::cast_sign_loss,
    reason = "lmax is a non-negative angular-momentum bound; the release-mode assert! on the next line enforces lmax >= 0 even when caller-side guards are absent"
)]
fn real_sph_harmonics(q: &Vector3<f64>, lmax: i32, out: &mut [f64]) {
    // Release-mode guard (not debug_assert!): this function is private today
    // but future callers could bypass the NonlocalPotential::new check, and
    // the `(lmax + 1) * (lmax + 1) as usize` at the next line wraps silently
    // on a negative lmax. One branch per call vs. hundreds of FLOPs is cheap.
    assert!(lmax >= 0, "real_sph_harmonics: lmax must be non-negative, got {lmax}");
    debug_assert_eq!(out.len(), ((lmax + 1) * (lmax + 1)) as usize);
    let fpi = 4.0 * PI;
    let inv_sqrt_fpi = (1.0 / fpi).sqrt();

    if lmax == 0 {
        out[0] = inv_sqrt_fpi;
        return;
    }

    let gmod = q.norm();
    let eps = 1e-9;

    // At |q|=0 Y_lm is angular-indeterminate; zero the l>0 block.
    // Y_00 is still 1/√(4π).
    for slot in out.iter_mut() {
        *slot = 0.0;
    }
    out[0] = inv_sqrt_fpi;
    if gmod < eps {
        return;
    }

    let cost = q.z / gmod;
    let sint = (1.0 - cost * cost).max(0.0).sqrt();
    // φ = atan2(qy, qx) — QE uses a custom branch; atan2 is equivalent and
    // numerically cleaner at the axes.
    let phi = q.y.atan2(q.x);

    // Q[l][m] for m in 0..=l. We only need the last two l rows during
    // recurrence, but lmax is tiny (≤ ~6 for any real PP) so a flat
    // (lmax+1)×(lmax+1) buffer is fine.
    let lm1 = (lmax + 1) as usize;
    let mut q_lm = vec![0.0_f64; lm1 * lm1];
    let qidx = |l: usize, m: usize| -> usize { l * lm1 + m };
    q_lm[qidx(0, 0)] = 1.0;
    q_lm[qidx(1, 0)] = cost;
    q_lm[qidx(1, 1)] = -sint / 2.0_f64.sqrt();

    // l=0: Y_00 already written at line 347 (the |q|=0 guard path initializes
    // it to the same value and returns before reaching here, so the invariant
    // holds whether we took that branch or not).
    // l=1: Y_10, Y_1,+1, Y_1,-1
    let c1 = (3.0 / fpi).sqrt();
    out[1] = c1 * q_lm[qidx(1, 0)];
    out[2] = c1 * 2.0_f64.sqrt() * q_lm[qidx(1, 1)] * phi.cos();
    out[3] = c1 * 2.0_f64.sqrt() * q_lm[qidx(1, 1)] * phi.sin();

    let mut lm_idx: usize = 4;
    for l in 2..=lmax as usize {
        let c = ((2 * l + 1) as f64 / fpi).sqrt();
        let l_f = l as f64;
        // Recurrence on l for Q(l, m), m = 0 ..= l-2.
        // `saturating_sub` keeps the range empty for l < 2 (unreachable here
        // because the outer `for l in 2..=lmax`, but belt+braces).
        for m in 0..=l.saturating_sub(2) {
            let m_f = m as f64;
            let llmm = (l_f * l_f - m_f * m_f).sqrt();
            let llm1 = ((l_f - 1.0) * (l_f - 1.0) - m_f * m_f).sqrt();
            q_lm[qidx(l, m)] = (cost * (2.0 * l_f - 1.0) * q_lm[qidx(l - 1, m)]
                - llm1 * q_lm[qidx(l - 2, m)])
                / llmm;
        }
        // m = l-1
        q_lm[qidx(l, l - 1)] = cost * (2.0 * l_f - 1.0).sqrt() * q_lm[qidx(l - 1, l - 1)];
        // m = l
        q_lm[qidx(l, l)] = -((2.0 * l_f - 1.0) / (2.0 * l_f)).sqrt() * sint * q_lm[qidx(l - 1, l - 1)];

        // Y_l,0 at lm_idx
        out[lm_idx] = c * q_lm[qidx(l, 0)];
        // Y_l,m (cos), Y_l,-m (sin) for m = 1..=l
        for m in 1..=l {
            let m_f = m as f64;
            let cosmphi = (m_f * phi).cos();
            let sinmphi = (m_f * phi).sin();
            out[lm_idx + 2 * m - 1] = c * 2.0_f64.sqrt() * q_lm[qidx(l, m)] * cosmphi;
            out[lm_idx + 2 * m] = c * 2.0_f64.sqrt() * q_lm[qidx(l, m)] * sinmphi;
        }
        lm_idx += 2 * l + 1;
    }
}

/// Spherical Bessel transform of a projector:
/// F(q) = 4π ∫₀^∞ [r·β(r)] j_l(qr) r dr
///
/// where `r_beta` stores r·β(r) (the UPF convention).
/// Spherical Bessel transform of a projector:
///   F(q) = 4π ∫₀^∞ [r·β(r)] j_l(qr) r dr
///
/// `r_grid`: radial grid points (Å).
/// `rab`: integration weights dr (Å). For log grids, rab[i] = r[i] × log_step.
/// `r_beta`: r·β(r) in Å^{-1/2} (UPF convention: projectors stored as r×β).
/// `l`: angular momentum quantum number.
/// `q`: wavevector magnitude |k+G| (Å⁻¹).
fn bessel_transform_projector(
    r_grid: &[f64],
    rab: &[f64],
    r_beta: &[f64],
    l: i32,
    q: f64,
) -> f64 {
    use crate::numerics::simpson_integrate;

    let n = r_grid.len();
    let mut integrand = vec![0.0; n];

    for i in 0..n {
        let r = r_grid[i];
        let rb = r_beta[i]; // r·β(r)
        let qr = q * r;
        let jl = spherical_bessel_j(l, qr);

        // Integrand: [r·β(r)] × j_l(qr) × r
        integrand[i] = rb * jl * r;
    }

    4.0 * PI * simpson_integrate(&integrand, rab)
}

/// Spherical Bessel function j_l(x) for arbitrary l >= 0.
///
/// Uses explicit formulas for l = 0, 1 and upward recurrence for l >= 2:
///   j_{l+1}(x) = (2l+1)/x × j_l(x) − j_{l-1}(x)
///
/// Note: upward recurrence is stable for l < x. For DFT pseudopotentials
/// l <= 6 is typical, and qr values are always moderate, so this is safe.
fn spherical_bessel_j(l: i32, x: f64) -> f64 {
    assert!(l >= 0, "spherical_bessel_j: l must be non-negative, got {l}");
    if x.abs() < 1e-10 {
        return if l == 0 { 1.0 } else { 0.0 };
    }
    if l == 0 {
        return x.sin() / x;
    }
    if l == 1 {
        return x.sin() / (x * x) - x.cos() / x;
    }
    // Upward recurrence from j_0, j_1
    let mut jlm1 = x.sin() / x;
    let mut jl = x.sin() / (x * x) - x.cos() / x;
    for n in 1..l {
        let jlp1 = ((2 * n + 1) as f64 / x).mul_add(jl, -jlm1);
        jlm1 = jl;
        jl = jlp1;
    }
    jl
}

/// Legendre polynomial P_l(x) for arbitrary l >= 0.
///
/// Uses Bonnet's recurrence relation:
///   (n+1) P_{n+1}(x) = (2n+1) x P_n(x) − n P_{n-1}(x)
///
/// Kept for unit-test validation of the spherical-harmonic addition theorem
/// (see `test_ylm_addition_theorem`); the hot-path V_NL assembly no longer
/// calls this (the angular sum is done via a GEMM over expanded real Y_lm
/// channels instead — see the module-level docs).
#[cfg(test)]
fn legendre_p(l: i32, x: f64) -> f64 {
    assert!(l >= 0, "legendre_p: l must be non-negative, got {l}");
    if l == 0 {
        return 1.0;
    }
    if l == 1 {
        return x;
    }
    let mut plm1 = 1.0;
    let mut pl = x;
    for n in 1..l {
        let plp1 = ((2 * n + 1) as f64 * x).mul_add(pl, -(n as f64 * plm1)) / (n + 1) as f64;
        plm1 = pl;
        pl = plp1;
    }
    pl
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

    #[test]
    fn test_spherical_bessel_j0() {
        assert!(relative_eq!(spherical_bessel_j(0, 0.0), 1.0, epsilon = 1e-10));
        // j_0(π) = sin(π)/π = 0
        assert!(spherical_bessel_j(0, PI).abs() < 1e-10);
        // j_0(1) = sin(1)/1 ≈ 0.8415
        assert!(relative_eq!(
            spherical_bessel_j(0, 1.0),
            1.0_f64.sin(),
            epsilon = 1e-10
        ));
    }

    #[test]
    fn test_spherical_bessel_j1() {
        assert!(relative_eq!(spherical_bessel_j(1, 0.0), 0.0, epsilon = 1e-10));
        // j_1(x) = sin(x)/x² - cos(x)/x
        let x: f64 = 2.0;
        let expected = x.sin() / (x * x) - x.cos() / x;
        assert!(relative_eq!(
            spherical_bessel_j(1, x),
            expected,
            epsilon = 1e-10
        ));
    }

    #[test]
    fn test_spherical_bessel_higher_l() {
        // j_2(x) = (3/x² - 1) sin(x)/x - 3 cos(x)/x²
        let x: f64 = 2.5;
        let j2_exact = (3.0 / (x * x) - 1.0) * x.sin() / x - 3.0 * x.cos() / (x * x);
        assert!(relative_eq!(spherical_bessel_j(2, x), j2_exact, epsilon = 1e-10));

        // j_3(x) = (15/x³ - 6/x) sin(x)/x - (15/x² - 1) cos(x)/x
        let j3_exact = (15.0 / (x * x * x) - 6.0 / x) * x.sin() / x
            - (15.0 / (x * x) - 1.0) * x.cos() / x;
        assert!(relative_eq!(spherical_bessel_j(3, x), j3_exact, epsilon = 1e-10));

        // j_4(3) ≈ 0.05615 (computed via recurrence from j_0, j_1)
        let j4_3 = spherical_bessel_j(4, 3.0);
        assert!(
            relative_eq!(j4_3, 0.056_149_714_328_844, epsilon = 1e-10),
            "j_4(3) = {j4_3}"
        );

        // j_5(5) ≈ 0.10681 (computed via recurrence)
        let j5_5 = spherical_bessel_j(5, 5.0);
        assert!(
            relative_eq!(j5_5, 0.106_811_161_456_505, epsilon = 1e-10),
            "j_5(5) = {j5_5}"
        );

        // j_4(5) ≈ 0.18702
        let j4_5 = spherical_bessel_j(4, 5.0);
        assert!(
            relative_eq!(j4_5, 0.187_017_655_344_889, epsilon = 1e-10),
            "j_4(5) = {j4_5}"
        );

        // j_l(0) = 0 for all l > 0
        for l in 2..=6 {
            assert!(
                spherical_bessel_j(l, 0.0).abs() < 1e-10,
                "j_{l}(0) should be 0"
            );
        }
    }

    #[test]
    fn test_legendre_higher_l() {
        // P_4(x) = (35x⁴ - 30x² + 3) / 8
        let x: f64 = 0.6;
        let p4_exact = (35.0 * x.powi(4) - 30.0 * x * x + 3.0) / 8.0;
        assert!(relative_eq!(legendre_p(4, x), p4_exact, epsilon = 1e-12));

        // P_5(x) = (63x⁵ - 70x³ + 15x) / 8
        let p5_exact = (63.0 * x.powi(5) - 70.0 * x.powi(3) + 15.0 * x) / 8.0;
        assert!(relative_eq!(legendre_p(5, x), p5_exact, epsilon = 1e-12));

        // P_6(x) = (231x⁶ - 315x⁴ + 105x² - 5) / 16
        let p6_exact =
            (231.0 * x.powi(6) - 315.0 * x.powi(4) + 105.0 * x * x - 5.0) / 16.0;
        assert!(relative_eq!(legendre_p(6, x), p6_exact, epsilon = 1e-12));

        // P_l(1) = 1 for all l
        for l in 0..=10 {
            assert!(
                relative_eq!(legendre_p(l, 1.0), 1.0, epsilon = 1e-12),
                "P_{l}(1) should be 1"
            );
        }

        // P_l(-1) = (-1)^l
        for l in 0..=10 {
            let expected = if l % 2 == 0 { 1.0 } else { -1.0 };
            assert!(
                relative_eq!(legendre_p(l, -1.0), expected, epsilon = 1e-12),
                "P_{l}(-1) should be {expected}"
            );
        }
    }

    #[test]
    fn test_legendre_orthogonality() {
        // ∫₋₁¹ P_l(x) P_m(x) dx = 2/(2l+1) δ_{lm}
        // Extended to l,m up to 6 (validates recurrence for higher l)
        let n = 1000;
        for l in 0..=6 {
            for m in 0..=6 {
                let mut integral = 0.0;
                for i in 0..=n {
                    let x = -1.0 + 2.0 * i as f64 / n as f64;
                    let w = if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
                    integral += w * legendre_p(l, x) * legendre_p(m, x);
                }
                integral *= 2.0 / (3.0 * n as f64);
                let expected = if l == m {
                    2.0 / (2 * l + 1) as f64
                } else {
                    0.0
                };
                assert!(
                    (integral - expected).abs() < 0.01,
                    "P_{l} · P_{m} integral: {integral:.6}, expected {expected:.6}"
                );
            }
        }
    }

    #[test]
    fn test_bessel_transform_delta() {
        // For a delta-like projector at r=0, the Bessel transform is constant
        // Not easily testable, but ensure it returns finite values
        let r = vec![0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0];
        let rab = vec![0.01, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.0];
        let rbeta: Vec<f64> = r.iter().map(|&ri: &f64| (-ri * ri).exp()).collect();

        let f0 = bessel_transform_projector(&r, &rab, &rbeta, 0, 0.0);
        let f1 = bessel_transform_projector(&r, &rab, &rbeta, 0, 1.0);
        let f5 = bessel_transform_projector(&r, &rab, &rbeta, 0, 5.0);

        assert!(f0.is_finite());
        assert!(f1.is_finite());
        assert!(f5.is_finite());
        // Should decay with q
        assert!(f5.abs() < f0.abs());
    }

    /// Defense-in-depth m-channel pin on V_NL (VNMT, 2026-04-18).
    ///
    /// The complementary check to `test_ylm_addition_theorem`: addition
    /// theorem only pins the *sum* `Σ_m Y_{l,m}(q̂₁)·Y_{l,m}(q̂₂)` against
    /// a Legendre identity, so it would still pass if a future edit
    /// shuffled normalization between m-channels in a way that leaves
    /// the sum (or the sum's Legendre equivalent) invariant under the
    /// `q̂₁ = q̂₂` special case. This test instead pins
    /// `H_NL[G₁, G₂]` at a chosen (non-degenerate) G-vector pair where
    /// three m-channels (m=0, +1, +2) each contribute a *distinct,
    /// individually hand-computed* amount — a normalization error on
    /// any single m-slot (including the √2-sensitive `(l=2, m=+2)`
    /// case called out by the VNLM review) shifts the total detectably.
    ///
    /// ## Setup
    ///
    /// * Simple cubic lattice with `a = 2π` Å, so the reciprocal
    ///   lattice is `b_i = ê_i` and integer Miller indices coincide
    ///   with Cartesian G-vectors in 1/Å.
    /// * Single Si atom at the origin → structure factor
    ///   `exp(-iG·τ) = 1` for every G.
    /// * `D_ij` zeroed everywhere except the diagonal entry on the
    ///   first `l=2` projector (row-major index `4 · 6 + 4 = 28`
    ///   of Si ONCV's 6-projector layout `[s,s,p,p,d,d]`). This
    ///   leaves exactly one `l=2` radial channel alive; its
    ///   contribution sums over all five `m ∈ {-2,-1,0,+1,+2}`.
    /// * k = Γ, so q = G.
    /// * G₁ = (1, 0, 1) · 1/Å, G₂ = (2, 0, 1) · 1/Å. Both have
    ///   φ = 0 and positive z-component, so Y_{2,-1} and Y_{2,-2}
    ///   vanish on both (leaving three live m-channels).
    ///
    /// ## Hand-computed reference
    ///
    /// With `ch = (atom=0, proj=4, m)`:
    ///
    /// ```text
    /// B[G, ch(m)]  = (1/√Ω) · F₄(|G|) · Y_{2,m}(Ĝ)
    /// (D·B^H)[ch(m), G] = D_{44} · conj(B[G, ch(m)]) = D_{44} · B[G, ch(m)]
    ///   (B is real here — τ=0 → phase=1 and Y_{l,m} is real)
    /// H_NL[G₁, G₂] = Σ_m B[G₁, ch(m)] · (D·B^H)[ch(m), G₂]
    ///              = (D_{44}/Ω) · F₄(|G₁|) · F₄(|G₂|) · Σ_m Y_{2,m}(Ĝ₁)·Y_{2,m}(Ĝ₂)
    /// ```
    ///
    /// Per-m breakdown (QE `ylmr2` convention — see `real_sph_harmonics` above):
    ///   * Y_{2,0}    = ½ · √(5 / 4π) · (3cos²θ − 1)
    ///   * Y_{2,+1}   = −√(15/4π)    · cosθ sinθ cos φ
    ///   * Y_{2,-1}   = −√(15/4π)    · cosθ sinθ sin φ
    ///   * Y_{2,+2}   =  √(15/16π)   · sin²θ cos(2φ)
    ///   * Y_{2,-2}   =  √(15/16π)   · sin²θ sin(2φ)
    ///
    /// G₁ = (1,0,1):  |G₁| = √2, cosθ = 1/√2, sinθ = 1/√2, φ = 0
    ///   * Y_{2,0}(Ĝ₁)  = ½ · √(5/4π) · (3/2 − 1)   = ¼ · √(5/4π)
    ///   * Y_{2,+1}(Ĝ₁) = −√(15/4π)   · ½ · 1       = −½ · √(15/4π)
    ///   * Y_{2,-1}(Ĝ₁) = 0                           (sin φ = 0)
    ///   * Y_{2,+2}(Ĝ₁) =  √(15/16π)  · ½ · 1       =  ½ · √(15/16π)
    ///   * Y_{2,-2}(Ĝ₁) = 0                           (sin 2φ = 0)
    ///
    /// G₂ = (2,0,1):  |G₂| = √5, cosθ = 1/√5, sinθ = 2/√5, φ = 0
    ///   * Y_{2,0}(Ĝ₂)  = ½ · √(5/4π) · (3/5 − 1)   = −1/5 · √(5/4π)
    ///   * Y_{2,+1}(Ĝ₂) = −√(15/4π)   · 2/5 · 1     = −2/5 · √(15/4π)
    ///   * Y_{2,-1}(Ĝ₂) = 0
    ///   * Y_{2,+2}(Ĝ₂) =  √(15/16π)  · 4/5 · 1     =  4/5 · √(15/16π)
    ///   * Y_{2,-2}(Ĝ₂) = 0
    ///
    /// Products and sum:
    ///   * m=0:   ¼ · √(5/4π) · (−1/5) · √(5/4π)           = −(1/20) · 5/(4π)  = −1/(16π)
    ///   * m=+1:  (−½) · √(15/4π) · (−2/5) · √(15/4π)      =  (1/5)  · 15/(4π) =  12/(16π) = 3/(4π)
    ///   * m=-1:  0
    ///   * m=+2:  ½ · √(15/16π) · (4/5) · √(15/16π)        =  (2/5)  · 15/(16π) = 6/(16π) = 3/(8π)
    ///   * m=-2:  0
    ///   * Σ_m = (−1 + 12 + 6) / (16π) = 17/(16π)
    ///
    /// Sanity via addition theorem:
    ///   cosθ₁₂ = Ĝ₁·Ĝ₂ = (1·2 + 0·0 + 1·1) / (√2·√5) = 3/√10
    ///   (2l+1)/(4π) · P_2(cosθ₁₂) = 5/(4π) · (3·(9/10) − 1)/2
    ///                             = 5/(4π) · 17/20 = 17/(16π)  ✓
    ///
    /// Expected:
    ///   H_NL[G₁, G₂] = (D_{44}/Ω) · F₄(√2) · F₄(√5) · 17/(16π)
    ///
    /// `F_4` is a Bessel transform of the Si UPF radial projector;
    /// we compute it via the same `bessel_transform_projector` helper
    /// the production path uses (NOT part of the "defense" — any bug
    /// in that function would also affect the production V_NL), so
    /// the angular Y_{2,m} pin is what this test actually guards.
    ///
    /// Tolerance 1e-10 (ULP headroom; the whole pipeline is scalar
    /// f64 on a handful of values).
    #[test]
    fn test_single_channel_l2_m_isolation() {
        use crate::pseudopotential::load;
        use std::path::PathBuf;

        // 1. Load Si ONCV PP and zero D_ij except diagonal entry for
        //    the first l=2 projector (projector index 4, row-major
        //    offset 4·6 + 4 = 28 in the 6×6 D matrix).
        let si_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("pseudopotentials/nc/lda/Si.upf");
        let mut pp = load(&si_path).expect("Si UPF must load");
        assert_eq!(pp.n_projectors(), 6, "Si ONCV PP expected 6 projectors");
        assert_eq!(pp.beta_projectors[4].l, 2, "projector 4 must be l=2");

        // Preserve the original D_{44} value before zeroing the rest.
        let d44 = pp.dij[4 * 6 + 4];
        assert!(d44.abs() > 1e-6, "Si ONCV D_{{44}} must be nontrivial, got {d44}");
        pp.dij.iter_mut().for_each(|d| *d = 0.0);
        pp.dij[4 * 6 + 4] = d44;

        // 2. Cubic lattice a = 2π → b_i = ê_i. G-vectors are integer
        //    Cartesian triples in 1/Å.
        let a = 2.0 * PI;
        let lattice = crate::crystal::Lattice::new(
            Vector3::new(a, 0.0, 0.0),
            Vector3::new(0.0, a, 0.0),
            Vector3::new(0.0, 0.0, a),
        );
        let omega = lattice.volume(); // (2π)³

        // ecut = 50 eV > HBAR2_OVER_2M · |(2,0,1)|² = 3.81 · 5 ≈ 19 eV
        // comfortably includes G₁=(1,0,1), G₂=(2,0,1).
        let basis = BasisSet::new(&lattice, 50.0);
        let g1_idx = basis
            .index_of(1, 0, 1)
            .expect("G=(1,0,1) must be in basis");
        let g2_idx = basis
            .index_of(2, 0, 1)
            .expect("G=(2,0,1) must be in basis");
        // Sanity: confirm the Cartesian coordinates (the assumption that
        // the cubic reciprocal lattice is the identity on Miller indices).
        assert!(
            relative_eq!(basis.g_vectors()[g1_idx], Vector3::new(1.0, 0.0, 1.0), epsilon = 1e-12)
        );
        assert!(
            relative_eq!(basis.g_vectors()[g2_idx], Vector3::new(2.0, 0.0, 1.0), epsilon = 1e-12)
        );

        // Single Si atom at origin: τ = 0 → phase = 1 for all G.
        let crystal = crate::crystal::Crystal {
            atoms: vec![crate::crystal::Atom::new(14, [0.0, 0.0, 0.0])],
            lattice,
        };

        // 3. Assemble V_NL at Γ via the production path.
        let k = Vector3::new(0.0, 0.0, 0.0);
        let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp])
            .expect("NonlocalPotential::new");

        let n_pw = basis.len();
        let mut h: Mat<Complex64> = Mat::zeros(n_pw, n_pw);
        vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k);
        let h_g1_g2 = h[(g1_idx, g2_idx)];

        // 4. Compute the expected value by hand — Σ_m Y_{2,m}(Ĝ₁)·Y_{2,m}(Ĝ₂) = 17/(16π).
        //    The two radial form factors F₄(|G|) are evaluated via the
        //    same Bessel-transform helper the production code calls.
        let fpi = 4.0 * PI;
        // Per-m products, re-computed explicitly from the QE convention (see doc comment).
        //   m=0:   ¼ · √(5/4π) · (−1/5) · √(5/4π) = −(1/20) · (5/4π) = −1/(16π)
        let y_m0 = -1.0 / (16.0 * PI);
        //   m=+1:  (−½) · √(15/4π) · (−2/5) · √(15/4π) = (1/5) · (15/4π) = 3/(4π) = 12/(16π)
        let y_mp1 = 3.0 / fpi;
        //   m=+2:  ½ · √(15/16π) · (4/5) · √(15/16π) = (2/5) · (15/16π) = 3/(8π) = 6/(16π)
        let y_mp2 = 3.0 / (8.0 * PI);
        let ang_sum = y_m0 + y_mp1 + y_mp2;
        // Consistency check with the closed-form addition theorem sum.
        assert!(
            relative_eq!(ang_sum, 17.0 / (16.0 * PI), epsilon = 1e-14),
            "per-m sum mismatch: got {ang_sum}, closed form {}",
            17.0 / (16.0 * PI)
        );

        let g1_norm = (2.0_f64).sqrt();
        let g2_norm = (5.0_f64).sqrt();
        let f4_g1 = bessel_transform_projector(
            &pp.r_grid,
            &pp.rab,
            &pp.beta_projectors[4].values,
            2,
            g1_norm,
        );
        let f4_g2 = bessel_transform_projector(
            &pp.r_grid,
            &pp.rab,
            &pp.beta_projectors[4].values,
            2,
            g2_norm,
        );

        let expected = Complex64::new(d44 / omega * f4_g1 * f4_g2 * ang_sum, 0.0);

        // 5. The computed and expected matrix elements must agree to
        //    ULP-ish precision; any per-m Y_{l,m} normalization error
        //    throws this off by ≥ ~1e-2 · |expected|.
        let diff = (h_g1_g2 - expected).norm();
        let scale = expected.norm().max(1.0);
        assert!(
            diff < 1e-10 * scale,
            "H_NL[G1,G2] mismatch: got {h_g1_g2}, expected {expected}, diff={diff}, scale={scale}"
        );
    }

    /// Verify the spherical-harmonic addition theorem on a realistic
    /// q-vector pair:
    ///   Σ_m Y_lm(q̂₁) Y_lm(q̂₂) = (2l+1)/(4π) · P_l(q̂₁·q̂₂)
    /// This is the identity that makes the GEMM-lifted KB assembly exact.
    #[test]
    #[allow(
        clippy::cast_sign_loss,
        reason = "test-only: lmax ranges over 0..=5 by construction"
    )]
    fn test_ylm_addition_theorem() {
        let qs = [
            Vector3::new(0.3, 0.7, -0.5),
            Vector3::new(-1.1, 0.2, 0.4),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        ];
        for lmax in 0..=5_i32 {
            let stride = ((lmax + 1) * (lmax + 1)) as usize;
            for q1 in &qs {
                for q2 in &qs {
                    let mut y1 = vec![0.0; stride];
                    let mut y2 = vec![0.0; stride];
                    real_sph_harmonics(q1, lmax, &mut y1);
                    real_sph_harmonics(q2, lmax, &mut y2);
                    for l in 0..=lmax {
                        let l_us = l as usize;
                        let base = l_us * l_us;
                        let mut lhs = 0.0;
                        for m in 0..(2 * l_us + 1) {
                            lhs += y1[base + m] * y2[base + m];
                        }
                        let cos_theta = q1.dot(q2) / (q1.norm() * q2.norm());
                        let rhs = (2 * l + 1) as f64 / (4.0 * PI) * legendre_p(l, cos_theta);
                        assert!(
                            (lhs - rhs).abs() < 1e-12,
                            "addition theorem l={l}: lhs={lhs} rhs={rhs}"
                        );
                    }
                }
            }
        }
    }
}
