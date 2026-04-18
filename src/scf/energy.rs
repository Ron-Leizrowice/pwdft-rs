//! Total energy computation and density utilities for SCF.
//!
//! ## Nonlinear core correction (NLCC)
//!
//! Several helpers here — `add_core_density`, `xc_energy_corrected` — are
//! shared between the standard Kohn-Sham path and the NLCC path. NLCC
//! (Louie, Froyen, Cohen, *Phys. Rev. B* **26**, 1738 (1982)) restores
//! the nonlinear coupling
//! ```text
//!     E_xc[ρ_val + ρ_core] − E_xc[ρ_val]
//! ```
//! that is dropped when the core is frozen and orthogonalized out of the
//! valence problem. The implementation keeps the core *only* inside the
//! XC functional:
//!
//! - `ρ_val + ρ_core` enters `ε_xc[·]` and `v_xc[·]` (see
//!   [`xc_energy_corrected`] and QE `PW/src/v_of_rho.f90:511`).
//! - `ρ_val` alone enters the Hartree source, the electron count, and the
//!   double-counting integral `∫ ρ_val · v_xc dr`.
//! - In LSDA, `ρ_core` is spin-unpolarized and split evenly as
//!   `ρ_core/2` between the two spin channels before being added to each
//!   `ρ_σ` (see `scf::driver_spin::run_scf_spin`).
//!
//! `ρ_core` itself is built on the FFT grid by
//! [`scf::potentials::compute_core_density`](super::potentials::compute_core_density)
//! from the PP's `PP_NLCC` block (see `src/pseudopotential/upf/convert.rs`
//! for the storage-unit convention — bare ρ_core(r) in e/Å³, *not* the
//! 4πr²·ρ convention used by `PP_RHOATOM`).

use nalgebra::Vector3;
use num_complex::Complex64;

use crate::{
    basis::BasisSet,
    consts::HBAR2_OVER_2M,
    fft::FFT3D,
    potential::xc,
};

// ---------------------------------------------------------------------------
// Energy computation
// ---------------------------------------------------------------------------

/// Sum of occupied Kohn-Sham eigenvalues weighted by k-point and occupation.
///
/// ```text
///     E_band = Σ_{n,k} f_{n,k} · w_k · ε_{n,k}
/// ```
/// - `ε_{n,k}` in eV (band index `n`, k-point index `k`);
/// - `f_{n,k}` dimensionless occupation in `[0, spin_factor]`
///   (`spin_factor = 2` for nspin=1, `1` for nspin=2);
/// - `w_k` k-point weight with `Σ_k w_k = 1` (IBZ-reduced; see
///   [`crate::symmetry::kpoints::reduce_kpoints`]).
///
/// E_band is **not** the total KS energy — it double-counts Hartree and
/// XC. See [`total_energy`] for the corrected expression. Returned in eV.
pub(crate) fn band_energy(
    eigenvalues: &[Vec<f64>],
    occupations: &[Vec<f64>],
    kpoint_weights: &[f64],
) -> f64 {
    eigenvalues
        .iter()
        .zip(occupations.iter())
        .zip(kpoint_weights.iter())
        .map(|((evs, occs), &w)| {
            evs.iter()
                .zip(occs.iter())
                .map(|(&e, &f)| f * w * e)
                .sum::<f64>()
        })
        .sum()
}

/// Classical electrostatic (Hartree) energy of the electron density,
/// summed in reciprocal space.
///
/// ```text
///     E_H = (Ω/2) · Σ_{G ≠ 0} |ρ(G)|² · 4πe² / |G|²
/// ```
/// Derivation: Parseval on `E_H = (1/2) ∫∫ ρ(r) ρ(r')/|r−r'| d³r d³r'` with
/// the convention `ρ(r) = (1/Ω) Σ_G ρ(G) e^{iG·r}` and `v_C(G) = 4πe²/|G|²`.
/// - `rho_g[ig]` complex Fourier coefficient `ρ(G)` in e/Å³ on the FFT
///   grid; `rho_g[0]` is the G=0 component (average density = N_el / Ω);
/// - `g_squared[ig]` = `|G|²` in Å⁻² for the same index;
/// - `omega` cell volume Ω in Å³;
/// - `4πe² = 4π · E2_COULOMB` with `E2_COULOMB ≈ 14.3996 eV·Å` (see
///   [`crate::consts::E2_COULOMB`]).
///
/// The G=0 divergence is excised: a neutral compensating background from
/// the ion lattice makes the full electrostatic sum finite. The
/// compensation is paid back in [`with_g0_shift`] and [`crate::ewald::ewald_energy`].
/// Returns E_H in eV.
pub(crate) fn hartree_energy(rho_g: &[Complex64], g_squared: &[f64], omega: f64) -> f64 {
    let fourpi_e2 = 4.0 * std::f64::consts::PI * crate::consts::E2_COULOMB;
    rho_g
        .iter()
        .zip(g_squared.iter())
        .map(|(rho, &g2)| {
            if g2 > crate::consts::G2_ZERO_THRESHOLD {
                rho.norm_sqr() * fourpi_e2 / g2
            } else {
                0.0
            }
        })
        .sum::<f64>()
        * 0.5
        * omega
}

/// Exchange-correlation energy with the Kohn-Sham double-counting subtraction.
///
/// ```text
///     E_xc − E_dc = E_xc[ρ_xc] − ∫ ρ_val(r) · v_xc[ρ_xc](r) d³r
/// ```
/// where `ρ_xc = ρ_val + ρ_core` under NLCC and `ρ_xc = ρ_val` otherwise.
///
/// The double-counting subtraction removes the `∫ρ·v_xc` piece that is
/// implicitly present in the band sum `Σ f w ε` (because `v_xc` enters
/// the Hamiltonian whose eigenvalues are summed). With NLCC, `v_xc` is
/// evaluated on the total density but the subtrahend integrates against
/// `ρ_val` only — the core is frozen out of the valence problem and must
/// not contribute to the band sum
/// (Louie, Froyen, Cohen, *Phys. Rev. B* **26**, 1738 (1982); see
/// QE `qe-7.5/PW/src/v_of_rho.f90:511`).
///
/// - `rho_xc`, `rho_val`: densities on the FFT grid in e/Å³;
/// - `exc_r[i]`: energy density per electron `ε_xc(ρ(r_i))` in eV;
/// - `vxc_r[i]`: XC potential `v_xc(r_i) = δE_xc/δρ(r_i)` in eV;
/// - `omega`: cell volume Ω in Å³.
///
/// The `E_xc[ρ_xc]` piece is `Σ_r ρ_xc(r) · ε_xc(r) · dV` via
/// [`xc::lda_xc_energy`]. Returns (E_xc − E_dc) in eV.
pub(crate) fn xc_energy_corrected(
    rho_xc: &[f64],
    rho_val: &[f64],
    exc_r: &[f64],
    vxc_r: &[f64],
    omega: f64,
) -> f64 {
    let n_grid = rho_xc.len();
    let dvol = omega / n_grid as f64;

    let e_xc = xc::lda_xc_energy(rho_xc, exc_r, omega);
    let e_vxc: f64 = rho_val
        .iter()
        .zip(vxc_r.iter())
        .map(|(&rho, &vxc)| rho * vxc * dvol)
        .sum();

    e_xc - e_vxc
}

/// Kohn-Sham total energy assembled from the band sum plus double-counting
/// corrections, evaluated on the **output** density.
///
/// ```text
///     E_KS = E_band − E_H[ρ_out] + (E_xc[ρ_out] − E_vxc[ρ_out]) + E_ion-ion
/// ```
/// Derivation: `E_band = Σ f w ε` contains `E_H + E_vxc + E_local + E_nl +
/// E_kin` (each band eigenvalue equals `⟨ψ|T + V_ext + V_H + V_xc|ψ⟩`),
/// so `E_H` and `E_vxc` are double-counted and must be removed, leaving
/// `E_xc` as the genuine functional value plus the ion-ion Ewald sum.
///
/// All arguments in eV. Caller is responsible for also applying the
/// `V_local(G=0) · N_el` compensating shift via [`with_g0_shift`] (the
/// G=0 of the local PP is zeroed to keep the Hamiltonian diagonal
/// finite; see [`crate::scf::context::ScfContext::new`] and NCFX).
///
/// Using the **output** density for the double-counting terms gives the
/// variationally exact KS energy once SCF is converged. Away from self-
/// consistency this estimator is only linear in (ρ_out − ρ_in);
/// [`harris_foulkes_energy`] gives a better early-iteration estimate.
pub(crate) fn total_energy(
    e_band: f64,
    e_hartree: f64,
    e_xc_corrected: f64,
    e_ewald: f64,
) -> f64 {
    e_band - e_hartree + e_xc_corrected + e_ewald
}

/// Harris-Foulkes non-variational energy estimator, evaluated on the
/// **input** density.
///
/// ```text
///     E_HF = E_band − E_H[ρ_in] + (E_xc[ρ_in] − E_vxc[ρ_in]) + E_ion-ion
/// ```
/// Unlike [`total_energy`], which mixes output eigenvalues with output
/// density, E_HF pairs the output eigenvalues (from diagonalizing
/// `H[ρ_in]`) with double-counting terms built from `ρ_in`. The functional
/// `E_KS[ρ]` is stationary at the self-consistent density, so the first
/// variation vanishes and
/// ```text
///     E_HF − E_KS = O(‖ρ_out − ρ_in‖²).
/// ```
/// That quadratic convergence makes `|E_HF − E_total|` a sensitive
/// self-consistency diagnostic — the driver warns when it stays large
/// after the density threshold is met.
///
/// Harris, *Phys. Rev. B* **31**, 1770 (1985);
/// Foulkes & Haydock, *Phys. Rev. B* **39**, 12520 (1989).
/// All arguments and return value in eV. Caller applies [`with_g0_shift`].
pub(crate) fn harris_foulkes_energy(
    e_band: f64,
    e_hartree_in: f64,
    e_xc_corrected_in: f64,
    e_ewald: f64,
) -> f64 {
    e_band - e_hartree_in + e_xc_corrected_in + e_ewald
}

/// Re-add the G=0 uniform-background piece of the local pseudopotential
/// that was subtracted to keep the Hamiltonian diagonal finite.
///
/// ```text
///     E_corrected = E + V_local(G=0) · N_el
/// ```
/// where `V_local(G=0) = (1/Ω) ∫ V_local(r) d³r` is the spatial average
/// of the local PP in eV, and `N_el` is the total valence electron
/// count (dimensionless). The G=0 component of every local PP diverges
/// as `−4πZ_α / |G|²` as G→0, so the bare sum is ill-defined. Zeroing
/// `V_local(G=0)` gauge-shifts the one-body Hamiltonian by a constant;
/// the constant re-enters here, with the Ewald sum
/// ([`crate::ewald::ewald_energy`]) providing the matching divergent
/// ion-ion piece so the total electrostatic energy is cutoff-
/// independent.
///
/// See NCFX / VGCMP proposals for the derivation;
/// [`crate::scf::context::ScfContext::new`] does the G=0 zeroing.
/// `energy` and the return value in eV.
pub(crate) fn with_g0_shift(energy: f64, ctx: &super::context::ScfContext<'_>) -> f64 {
    energy + ctx.v_local_g0 * ctx.n_electrons
}

// ---------------------------------------------------------------------------
// Density utilities
// ---------------------------------------------------------------------------

/// Volume-averaged RMS density difference, used as the primary SCF
/// convergence scalar.
///
/// ```text
///     Δρ_rms = sqrt( (1/Ω) · Σ_r |ρ_new(r) − ρ_old(r)|² · dV )
///            = sqrt( (1/N) · Σ_r |ρ_new(r) − ρ_old(r)|² )
/// ```
/// where `dV = Ω / N_grid` is the volume element and `N = N_grid` is
/// the total number of real-space grid points. Units: e/Å³. The Ω
/// factor cancels algebraically, so the result is a pure per-grid-point
/// RMS; it is kept in the signature for dimensional clarity.
///
/// Compared element-wise against `ScfParams::conv_threshold`. For
/// nspin=2 the driver takes `max(Δρ_up, Δρ_down)` rather than the
/// total-density difference — see SPNC comment in
/// `scf::driver_spin::run_scf_spin`.
pub(crate) fn density_diff(rho_old: &[f64], rho_new: &[f64], omega: f64, n_grid: usize) -> f64 {
    let dvol = omega / n_grid as f64;
    let sum_sq: f64 = rho_old
        .iter()
        .zip(rho_new.iter())
        .map(|(&a, &b)| (a - b).powi(2) * dvol)
        .sum();
    (sum_sq / omega).sqrt()
}

/// Forward FFT with the `1/N` normalization convention used throughout
/// the SCF pipeline.
///
/// ```text
///     ρ(G) = (1/N) · Σ_r ρ(r) · exp(−i G · r)
/// ```
/// where `N = n_x · n_y · n_z = fft.total_size()`. This normalization
/// is the reciprocal of [`crate::fft::FFT3D::inverse`], so a forward
/// followed by inverse recovers the input exactly. With this
/// convention `ρ(G=0) = (1/N) Σ_r ρ(r) = ⟨ρ⟩ = N_el / Ω`.
///
/// `rho_r` in e/Å³; `rho_g[ig]` returned as a complex coefficient in
/// e/Å³. Length must equal `fft.total_size()`.
pub(crate) fn density_r_to_g(fft: &mut FFT3D, rho_r: &[f64], rho_g: &mut [Complex64]) {
    for (i, &r) in rho_r.iter().enumerate() {
        rho_g[i] = Complex64::new(r, 0.0);
    }
    fft.forward(rho_g);
    let norm = 1.0 / fft.total_size() as f64;
    for v in rho_g.iter_mut() {
        *v *= norm;
    }
}

/// Owning-allocation convenience wrapper around [`density_r_to_g`].
///
/// Identical `1/N`-normalized forward FFT convention; returns a freshly
/// allocated `Vec<Complex64>` of length `data_r.len()`. Any real-space
/// field (density, potential, core charge) in the FFT grid's natural
/// units passes through unchanged — see [`density_r_to_g`] for the
/// explicit transform and normalization.
pub(crate) fn real_to_g_space(data_r: &[f64], fft: &mut FFT3D) -> Vec<Complex64> {
    let mut data_g = vec![Complex64::new(0.0, 0.0); data_r.len()];
    density_r_to_g(fft, data_r, &mut data_g);
    data_g
}

/// Element-wise sum of the three local-potential contributions to the
/// Kohn-Sham effective potential in reciprocal space.
///
/// ```text
///     V_eff(G) = V_local(G) + V_H(G) + V_xc(G)
/// ```
/// - `V_local(G)`: ionic local PP (G=0 zeroed; compensated by
///   [`with_g0_shift`]);
/// - `V_H(G) = 4πe² · ρ(G) / |G|²`: classical electron repulsion;
/// - `V_xc(G)`: Fourier transform of the LDA XC potential `v_xc(r)
///   = δE_xc/δρ`.
///
/// All three arrays live on the FFT grid with consistent index ordering;
/// all entries in eV. Returned `V_eff(G)` is the reciprocal-space
/// representation folded back onto the wavefunction basis by
/// [`crate::scf::potentials::build_hamiltonian_with_v_eff`]
/// to form the Hamiltonian matrix elements
/// `⟨G | V_eff | G'⟩ = V_eff(G − G')`. The non-local KB term is added
/// separately; see [`crate::potential::nonlocal::NonlocalPotential::add_to_hamiltonian`].
pub(crate) fn assemble_v_eff(
    v_local: &[Complex64],
    v_h: &[Complex64],
    v_xc: &[Complex64],
) -> Vec<Complex64> {
    use rayon::prelude::*;
    v_local
        .par_iter()
        .zip(v_h.par_iter())
        .zip(v_xc.par_iter())
        .map(|((&vl, &vh), &vxc)| vl + vh + vxc)
        .collect()
}

/// Add NLCC core density to valence density for XC evaluation.
///
/// Implements `ρ_xc(r) = max(ρ_val(r) + ρ_core(r), 0)`. If `rho_core` is
/// empty (no NLCC), returns `rho_val` unchanged.
///
/// The clamp protects the LDA XC functional from spurious negative
/// densities that can arise from FFT-wrap round-off in ρ_core or from
/// density mixing; without it, `ρ^(1/3)` in the exchange term would
/// produce NaN. Mirrors QE `PW/src/v_of_rho.f90:511` (adds `rho_core` to
/// `rho%of_r(ir,1)` before calling `xc_lda`).
///
/// Reference: Louie, Froyen, Cohen, *Phys. Rev. B* **26**, 1738 (1982).
pub(crate) fn add_core_density(rho_val: &[f64], rho_core: &[f64]) -> Vec<f64> {
    if rho_core.is_empty() {
        rho_val.to_vec()
    } else {
        rho_val
            .iter()
            .zip(rho_core.iter())
            .map(|(&v, &c)| (v + c).max(0.0))
            .collect()
    }
}

/// Solve Poisson's equation in reciprocal space on the full FFT grid.
///
/// ```text
///     V_H(G) = 4πe² · ρ(G) / |G|²     (G ≠ 0)
///     V_H(G=0) = 0
/// ```
/// The real-space form `∇² V_H = −4πe² ρ` diagonalizes on plane waves
/// (`∇² → −|G|²`). The `G=0` component is set to zero: in a neutral
/// periodic system the electron–electron and electron–ion `G=0` pieces
/// cancel, and `E_ion-ion` (Ewald) supplies the finite remainder.
///
/// - `rho_g`: total electron density in G-space (e/Å³);
/// - `g_squared[ig]` = |G|² in Å⁻², thresholded by
///   [`crate::consts::G2_ZERO_THRESHOLD`] to identify G=0;
/// - `4πe² = 4π · E2_COULOMB` with `E2_COULOMB ≈ 14.3996 eV·Å`.
///
/// Returns `V_H(G)` in eV on the same FFT grid.
pub(crate) fn hartree_on_fft_grid(rho_g: &[Complex64], g_squared: &[f64]) -> Vec<Complex64> {
    use rayon::prelude::*;
    let fourpi_e2 = 4.0 * std::f64::consts::PI * crate::consts::E2_COULOMB;

    rho_g
        .par_iter()
        .zip(g_squared.par_iter())
        .map(|(&rho, &g2)| {
            if g2 > crate::consts::G2_ZERO_THRESHOLD {
                rho * fourpi_e2 / g2
            } else {
                Complex64::new(0.0, 0.0)
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Per-component diagnostics (VGC5)
// ---------------------------------------------------------------------------

/// Kinetic-energy expectation value summed over occupied bands and
/// k-points.
///
/// ```text
///     E_kin = Σ_{n,k} f_{n,k} · w_k · ⟨ψ_{n,k}| T |ψ_{n,k}⟩
///           = Σ_{n,k} f_{n,k} · w_k · Σ_G |c_{n,k}(G)|² · (ℏ²/2m) · |k + G|²
/// ```
/// The kinetic operator `T = −(ℏ²/2m) ∇²` is diagonal in the plane-wave
/// basis, giving `⟨k+G|T|k+G'⟩ = (ℏ²/2m) |k+G|² · δ_{G,G'}`.
/// `ℏ²/(2m) = HBAR2_OVER_2M ≈ 3.810 eV·Å²` ([`crate::consts::HBAR2_OVER_2M`]).
///
/// - `c_{n,k}(G)` = `wavefunctions[ik][(ig, nb)]` plane-wave coefficients,
///   assumed orthonormal in the Bloch sense `Σ_G |c_{n,k}(G)|² = 1`;
/// - `k + G` in Å⁻¹; `|k+G|²` in Å⁻²;
/// - band occupation `f_{n,k}` dimensionless (Fermi-Dirac or other;
///   [`crate::scf::smearing`]);
/// - k-point weight `w_k` with `Σ_k w_k = 1`.
///
/// Returns E_kin in eV. Parallelized over k-points via rayon.
pub(crate) fn kinetic_expectation(
    basis: &BasisSet,
    k_points: &[Vector3<f64>],
    kpoint_weights: &[f64],
    wavefunctions: &[faer::Mat<Complex64>],
    occupations: &[Vec<f64>],
) -> f64 {
    use rayon::prelude::*;
    let g_vecs = basis.g_vectors();
    (0..k_points.len())
        .into_par_iter()
        .map(|ik| {
            let k = &k_points[ik];
            let wfn = &wavefunctions[ik];
            let occ = &occupations[ik];
            let w = kpoint_weights[ik];
            let n_pw = wfn.nrows();
            let n_bands = wfn.ncols();
            let mut e = 0.0f64;
            for nb in 0..n_bands {
                let f = occ[nb];
                if f == 0.0 {
                    continue;
                }
                let mut t_band = 0.0f64;
                for ig in 0..n_pw {
                    let c = wfn[(ig, nb)];
                    let ke = HBAR2_OVER_2M * (k + g_vecs[ig]).norm_squared();
                    t_band += c.norm_sqr() * ke;
                }
                e += f * w * t_band;
            }
            e
        })
        .sum()
}

/// Electron–ion local-PP energy, integrated on the real-space FFT grid
/// with the G=0 piece excluded.
///
/// ```text
///     E_local(G ≠ 0) = ∫ ρ(r) · V_local(r) d³r  ≈  Σ_r ρ(r) · V_local(r) · dV
/// ```
/// with `dV = Ω / N_grid`. `v_local_fft_r` is the inverse FFT of the
/// reciprocal-space local PP with its G=0 component pre-zeroed by
/// [`crate::scf::context::ScfContext::new`] (so `∫V_local dV = 0` by
/// construction). The compensating uniform background is accounted for
/// separately in `EnergyComponents::e_local_g0_shift` via
/// [`with_g0_shift`].
///
/// - `rho_r`: valence density in e/Å³ (total density in spin-
///   polarized runs — the local PP is spin-independent);
/// - `v_local_fft_r`: `V_local(r)` in eV, G=0 removed;
/// - `omega`: Ω in Å³.
///
/// Returns the G ≠ 0 piece of E_local in eV.
pub(crate) fn local_pp_energy_grid(
    rho_r: &[f64],
    v_local_fft_r: &[f64],
    omega: f64,
    n_grid: usize,
) -> f64 {
    let dvol = omega / n_grid as f64;
    rho_r
        .iter()
        .zip(v_local_fft_r.iter())
        .map(|(&r, &v)| r * v * dvol)
        .sum()
}

/// Non-local (Kleinman-Bylander) pseudopotential energy summed over
/// occupied states.
///
/// ```text
///     E_nl = Σ_{n,k} f_{n,k} · w_k · ⟨ψ_{n,k}| V_NL |ψ_{n,k}⟩
///          = Σ_{n,k} f_{n,k} · w_k · Σ_{G,G'} c*_{n,k}(G) · H_NL(G,G') · c_{n,k}(G')
/// ```
/// with the separable KB form
/// `V_NL = Σ_α Σ_{lm} D_{lm}^{(α)} |β_{lm}^{(α)}⟩⟨β_{lm}^{(α)}|`
/// (Kleinman & Bylander, *Phys. Rev. Lett.* **48**, 1425 (1982)).
/// `H_NL(G, G')` is built per k-point by the cached
/// [`crate::potential::nonlocal::NonlocalPotential::add_to_hamiltonian`];
/// see the VNLM single-GEMM assembly for the reciprocal-space form.
///
/// - `wavefunctions[ik]`: n_pw × n_bands column-major coefficient matrix;
/// - occupations, k-point weights: same conventions as
///   [`kinetic_expectation`].
///
/// Returns E_nl in eV. Parallelized over k-points via rayon; each
/// k-point pays an O(n_pw²) matrix-vector cost per band.
pub(crate) fn nonlocal_expectation(
    basis: &BasisSet,
    crystal: &crate::crystal::Crystal,
    k_points: &[Vector3<f64>],
    kpoint_weights: &[f64],
    wavefunctions: &[faer::Mat<Complex64>],
    occupations: &[Vec<f64>],
    vnl_cache: &[crate::potential::nonlocal::NonlocalPotential],
) -> f64 {
    use rayon::prelude::*;
    (0..k_points.len())
        .into_par_iter()
        .map(|ik| {
            let k = &k_points[ik];
            let wfn = &wavefunctions[ik];
            let occ = &occupations[ik];
            let w = kpoint_weights[ik];
            let n_pw = wfn.nrows();
            let n_bands = wfn.ncols();

            let mut h_nl = faer::Mat::<Complex64>::zeros(n_pw, n_pw);
            vnl_cache[ik].add_to_hamiltonian(&mut h_nl, crystal, basis, k);

            let mut e = 0.0f64;
            for nb in 0..n_bands {
                let f = occ[nb];
                if f == 0.0 {
                    continue;
                }
                // ⟨ψ|H_NL|ψ⟩ = Σ_{G,G'} c*(G) H_NL[G,G'] c(G')
                let mut acc = Complex64::new(0.0, 0.0);
                for ig in 0..n_pw {
                    let mut row_sum = Complex64::new(0.0, 0.0);
                    for jg in 0..n_pw {
                        row_sum += h_nl[(ig, jg)] * wfn[(jg, nb)];
                    }
                    acc += wfn[(ig, nb)].conj() * row_sum;
                }
                e += f * w * acc.re;
            }
            e
        })
        .sum()
}

/// Bare (uncorrected) exchange-correlation energy.
///
/// ```text
///     E_xc = ∫ ρ_xc(r) · ε_xc(ρ_xc(r)) d³r  ≈  Σ_r ρ_xc(r) · ε_xc(r) · dV
/// ```
/// with `ρ_xc = ρ_val + ρ_core` under NLCC or `ρ_val` otherwise;
/// `ε_xc[ρ]` is the LDA energy density per electron from
/// [`crate::potential::xc`] (Perdew-Zunger parametrization of the
/// Ceperley-Alder Monte Carlo data,
/// *Phys. Rev. B* **23**, 5048 (1981)).
///
/// No Kohn-Sham double-counting subtraction — this is the raw
/// functional value used in `EnergyComponents::e_xc` for validation
/// against QE's "xc contribution" line. See [`xc_energy_corrected`] for
/// the version that enters the total energy.
/// All arguments in eV / e·Å⁻³ / Å³; returns E_xc in eV.
pub(crate) fn xc_energy_bare(rho_xc: &[f64], exc_r: &[f64], omega: f64) -> f64 {
    xc::lda_xc_energy(rho_xc, exc_r, omega)
}

// ---------------------------------------------------------------------------
// Per-component decomposition (VGC5 diagnostic)
// ---------------------------------------------------------------------------

/// Per-term decomposition of the Kohn-Sham total energy (VGC5
/// diagnostic).
///
/// Intended as a validation handle: each field is an independent direct
/// evaluation of one term in the KS functional, so their sum is
/// identical to the result from the double-counting assembly in
/// `total_energy` only at self-consistency. All fields in eV.
///
/// ## Direct-sum identity (converged density)
///
/// ```text
///     E_total = e_kinetic
///             + e_local
///             + e_local_g0_shift       (= V_local(G=0) · N_el)
///             + e_nonlocal
///             + e_hartree
///             + e_xc
///             + e_ewald
/// ```
///
/// ## Kohn-Sham double-counting identity
///
/// Each band eigenvalue satisfies
/// `ε_{n,k} = ⟨ψ_{n,k}| T + V_ext + V_H + V_xc |ψ_{n,k}⟩` (with
/// `V_ext = V_local + V_NL`), so
/// ```text
///     E_band = e_kinetic + e_local + e_nonlocal + 2·e_hartree + e_vxc,
/// ```
/// where `e_vxc = ∫ρ(r)·V_xc(r)dr` appears because V_H is self-linear
/// in ρ (factor of 2) while V_xc is not. The `total_energy` assembly
/// subtracts the double-counted pieces (`− E_H`, `− E_vxc`) and re-adds
/// the true functional value `E_xc`; agreement of the two routes
/// confirms the one-body / two-body accounting.
///
/// ## Correspondence with QE `pw.x` output
///
/// ```text
///     one-electron contribution = e_kinetic + e_local + e_nonlocal + e_local_g0_shift
///     hartree    contribution = e_hartree
///     xc         contribution = e_xc
///     ewald      contribution = e_ewald
/// ```
///
/// Computed once on the final (converged) iteration in both driver
/// paths with one extra pass over wavefunctions, the local PP on the
/// FFT grid, and the non-local operator. Not used inside the SCF hot
/// loop.
#[derive(Debug, Clone)]
pub struct EnergyComponents {
    /// Weighted sum of occupied Kohn-Sham eigenvalues (eV):
    /// `Σ_{n,k} f_{n,k} · w_k · ε_{n,k}`.
    pub e_band: f64,
    /// Kinetic-energy expectation `Σ f·w·⟨ψ|T|ψ⟩` in eV, with
    /// `⟨ψ|T|ψ⟩ = Σ_G |c_{n,k}(G)|² · (ℏ²/2m) · |k+G|²`.
    pub e_kinetic: f64,
    /// Electron–ion local-PP energy `∫ρ(r)·V_local(r)d³r` (G ≠ 0 piece,
    /// eV). The G=0 component of `V_local` is zeroed at setup to keep
    /// the Hamiltonian diagonal finite; the compensating uniform
    /// background is stored separately in `e_local_g0_shift`.
    pub e_local: f64,
    /// Uniform-background restoration `V_local(G=0) · N_el` (eV), with
    /// `V_local(G=0) = (1/Ω) ∫ V_local(r) d³r`. Present in every
    /// neutral periodic pseudopotential calculation; paired with the
    /// matching ion-ion Ewald divergence so the total electrostatic
    /// energy is cutoff-independent.
    pub e_local_g0_shift: f64,
    /// Kleinman-Bylander separable non-local PP energy
    /// `Σ f·w·⟨ψ|V_NL|ψ⟩` in eV.
    pub e_nonlocal: f64,
    /// Classical electrostatic self-energy of the electrons
    /// `(Ω/2) Σ_{G ≠ 0} |ρ(G)|²·4πe²/|G|²` in eV, evaluated on the
    /// **output** density (the one produced by the last diagonalization).
    pub e_hartree: f64,
    /// Bare exchange-correlation energy `∫ρ_xc(r)·ε_xc(r)d³r` in eV from
    /// the output density, with `ρ_xc = ρ_val + ρ_core` under NLCC and
    /// `ρ_xc = ρ_val` otherwise. Sign and normalization match QE's
    /// "xc contribution" line.
    pub e_xc: f64,
    /// XC double-counting integral `∫ρ_val(r)·V_xc(r)d³r` in eV from the
    /// output density. This is the `E_vxc` quantity that appears in the
    /// Kohn-Sham double-counting subtraction `E_xc − E_vxc` inside the
    /// total-energy assembly and in the band-sum identity
    /// `E_band = e_kinetic + e_local + e_nonlocal + 2·e_hartree + e_vxc`.
    ///
    /// LSDA: computed as `∫ρ↑·V_xc↑ d³r + ∫ρ↓·V_xc↓ d³r`, i.e. the two
    /// spin channels are summed into a single scalar. The core density is
    /// NOT included here (the core is frozen and does not contribute to
    /// the band sum), but `V_xc` itself is evaluated on
    /// `ρ_val + ρ_core` when NLCC is active.
    pub e_vxc: f64,
    /// Ewald ion-ion electrostatic energy in eV; independent of the
    /// electron density and of spin, so computed once in
    /// `ScfContext::new` and copied through every iteration. See
    /// [`crate::ewald::ewald_energy`].
    pub e_ewald: f64,
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "ERR2 § Phase 0: in-src test modules are allowed to panic")]
mod tests {
    use super::*;

    use approx::relative_eq;
    use crate::fft::FFT3D;

    #[test]
    fn test_real_to_g_space_dc_component() {
        // A constant real-space function f(r) = C should give
        // F(G=0) = C and F(G≠0) = 0.
        let mut fft = FFT3D::new(8, 8, 8);
        let n = fft.total_size();
        let c = 3.5;
        let data_r = vec![c; n];
        let data_g = real_to_g_space(&data_r, &mut fft);

        // G=0 component (index 0) should be C
        assert!(
            relative_eq!(data_g[0].re, c, epsilon = 1e-10),
            "DC component: expected {c}, got {}", data_g[0].re
        );
        assert!(data_g[0].im.abs() < 1e-10);

        // All other G-components should be ~0
        for (i, &v) in data_g.iter().enumerate().skip(1) {
            assert!(
                v.norm() < 1e-10,
                "G≠0 component at {i}: expected ~0, got {v}"
            );
        }
    }

    #[test]
    fn test_real_to_g_space_roundtrip() {
        let mut fft = FFT3D::new(8, 8, 8);
        let n = fft.total_size();
        let data_r: Vec<f64> = (0..n).map(|i| (i as f64 * 0.1).sin()).collect();
        let data_g = real_to_g_space(&data_r, &mut fft);

        // Inverse FFT should recover original (unnormalized → need N factor)
        let mut data_back = data_g;
        fft.inverse(&mut data_back);
        // real_to_g_space divides by N, inverse multiplies by N → should recover original
        for (i, (&orig, &back)) in data_r.iter().zip(data_back.iter()).enumerate() {
            assert!(
                relative_eq!(orig, back.re, epsilon = 1e-10),
                "Roundtrip failed at {i}: original={orig}, recovered={}", back.re
            );
            assert!(back.im.abs() < 1e-10, "Imaginary part at {i}: {}", back.im);
        }
    }

    #[test]
    fn test_assemble_v_eff_adds_correctly() {
        let n = 100;
        let v1: Vec<Complex64> = (0..n).map(|i| Complex64::new(i as f64, 0.0)).collect();
        let v2: Vec<Complex64> = (0..n).map(|i| Complex64::new(0.0, i as f64 * 0.1)).collect();
        let v3: Vec<Complex64> = (0..n).map(|i| Complex64::new(-(i as f64) * 0.5, 0.0)).collect();

        let result = assemble_v_eff(&v1, &v2, &v3);

        for i in 0..n {
            let expected = v1[i] + v2[i] + v3[i];
            assert!(
                (result[i] - expected).norm() < 1e-14,
                "V_eff mismatch at {i}: expected {expected}, got {}", result[i]
            );
        }
    }

    #[test]
    fn test_hartree_on_fft_grid_g0_zero() {
        // V_H(G=0) should be zero (no divergence)
        let rho_g = vec![Complex64::new(1.0, 0.0); 10];
        let g_squared = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let v_h = hartree_on_fft_grid(&rho_g, &g_squared);
        assert!(v_h[0].norm() < 1e-15, "V_H(G=0) should be zero, got {}", v_h[0]);
        // V_H(G≠0) should be finite and positive real for positive ρ
        for &v in &v_h[1..] {
            assert!(v.re > 0.0, "V_H should be positive for positive ρ: {v}");
        }
    }

    #[test]
    fn test_density_diff_identical() {
        let rho = vec![1.0; 100];
        let diff = density_diff(&rho, &rho, 40.0, 100);
        assert!(diff < 1e-15, "Identical densities should give zero diff: {diff}");
    }

    #[test]
    fn test_density_diff_known() {
        let omega = 40.0;
        let n = 100;
        let rho_a = vec![1.0; n];
        let rho_b = vec![2.0; n];
        // diff = sqrt(Σ(1.0)² × dvol / omega) = sqrt(n × dvol / omega) = sqrt(dvol × n / omega)
        // dvol = omega / n = 0.4
        // diff = sqrt(0.4 * 100 / 40) = sqrt(1.0) = 1.0
        let diff = density_diff(&rho_a, &rho_b, omega, n);
        assert!(
            relative_eq!(diff, 1.0, epsilon = 1e-10),
            "Expected diff=1.0, got {diff}"
        );
    }
}
