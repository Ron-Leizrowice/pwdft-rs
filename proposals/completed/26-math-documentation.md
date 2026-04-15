# Proposal 26: Mathematical Documentation — Reference Docs, Audit, and Docstrings

**Status:** Partially done. `docs/theory.md` and `docs/math-audit.md` are written (uncommitted). All docstring changes and the `(-1)^l` comment fix are still outstanding.

## Problem

The codebase implements ~50 physics-critical functions across 12 modules but lacks:

1. **Reference documentation** — no single document collecting the canonical formulas and references (Payne 1992, Martin textbook, PZ 1981, etc.) for debugging and onboarding.
2. **Formula verification** — no record of which formulas were audited against references. The Fe NLCC bug (Proposal 25) might have been caught earlier with a systematic audit.
3. **Docstrings on core functions** — many functions implementing non-trivial math lack formula documentation, unit specifications, or both.
4. **One incorrect comment** — `nonlocal.rs:50-52` claims `i^l × (i*)^l = (-1)^l` when the correct identity is `|i|^{2l} = 1`. The code is correct (does NOT apply the factor); only the comment is wrong.

## Part A: Reference Documentation

### `docs/theory.md` (~18 KB, already drafted)

Comprehensive reference covering all formulas used in the code:

1. **Total energy** — KS decomposition, double-counting correction `E = E_band - E_H + E_xc - E_vxc + E_ewald`, NLCC modification (`E_xc` uses `ρ_val + ρ_core`, `E_vxc` uses `ρ_val` only)
2. **Reciprocal-space formulation** — Bloch waves, PW expansion, KS matrix equation `H_{G,G'} c = ε c`, FFT convention (1/N on forward)
3. **Potentials** — Local PP with G=0 convention, Hartree `V_H(G) = 4πe²ρ(G)/|G|²`, LDA XC (Slater + PZ with all 7 parameters), NLCC
4. **Non-local pseudopotential** — KB separable form, PW matrix elements with angular sum `(2l+1)/(4π) P_l(cos θ)`, Bessel transform form factor, D_ij for UPF/PSP8
5. **Ewald summation** — All four terms with formulas, screening parameter `η = (Nπ/Ω)^{1/3}`
6. **Electron density** — From wavefunctions `ρ = Σ f w |ψ|²`, SAD initial density
7. **SCF convergence** — Linear mixing, Anderson/Pulay DIIS with overlap matrix, Kerker preconditioning `P(G) = |G|²/(|G|² + q_TF²)`
8. **Smearing and entropy** — FD, Gaussian, MP, Cold occupation and entropy formulas, Mermin free energy `F = E - TS`, sigma→0 extrapolation `E₀ = (E+F)/2`
9. **Units convention** — eV/Å table with all conversion factors
10. **Common pitfalls** — NLCC omission, V_local(G=0), FFT normalization, PSP8 D_ij, spin factor

Each section cites canonical references and links to code locations.

### `docs/math-audit.md` (~11 KB, already drafted)

Line-by-line verification of 17 formulas against the code. All verified **CORRECT**:

| # | Formula | Code location |
|---|---------|---------------|
| 1 | Total energy | `scf/energy.rs` |
| 2 | Band energy | `scf/energy.rs` |
| 3 | Hartree energy | `scf/energy.rs` |
| 4 | Hartree potential | `scf/potentials.rs` |
| 5 | XC energy | `potential/xc.rs` |
| 6 | Slater exchange | `potential/xc.rs` |
| 7 | PZ correlation (7 params) | `potential/xc.rs` |
| 8 | Kinetic energy | `scf/potentials.rs` |
| 9 | Local pseudopotential | `scf/potentials.rs` |
| 10 | Non-local KB | `potential/nonlocal.rs` |
| 11 | Ewald (4 terms) | `ewald.rs` |
| 12 | Density | `scf/density.rs` |
| 13 | Fermi energy search | `scf/smearing.rs` |
| 14 | Smearing (4 schemes) | `scf/smearing.rs` |
| 15 | Sigma→0 extrapolation | `scf/mod.rs` |
| 16 | NLCC | `scf/mod.rs` |
| 17 | Spin-polarized LSDA | `potential/xc.rs` |

## Part B: Docstring Improvements

### B1. Functions with no docstring on non-trivial math

**`scf/smearing.rs` — 4 occupation functions (lines 107-149):**

```rust
/// Fermi-Dirac occupation (before spin factor).
///
/// f(ε) = 1 / (1 + exp(x))  where x = (ε - E_F) / σ
///
/// At T=0 (σ→0): step function θ(E_F - ε), with f(E_F) = 1/2.
/// Overflow-protected for |x| > 40.
fn fermi_dirac_01(...) -> f64 {

/// Gaussian smearing occupation (before spin factor).
///
/// f(ε) = erfc(x) / 2  where x = (ε - E_F) / σ
///
/// Equivalent to f = [1 + erf((E_F - ε)/σ)] / 2.
fn gaussian_01(...) -> f64 {

/// Methfessel-Paxton order-1 occupation (before spin factor).
///
/// f(ε) = erfc(x)/2 - (x/2) exp(-x²) / √π
///
/// Reference: Methfessel & Paxton, Phys. Rev. B 40, 3616 (1989).
fn methfessel_paxton_01(...) -> f64 {

/// Marzari-Vanderbilt "cold" smearing occupation (before spin factor).
///
/// f(ε) = (1/2) erfc(x + 1/√2) + exp(-(x + 1/√2)²) / √(2π)
///
/// Designed to give positive-definite entropy. Argument shifted by 1/√2
/// so that f(E_F) = 1/2 exactly.
///
/// Reference: Marzari, Vanderbilt, De Vita, Payne, Phys. Rev. Lett. 82, 3296 (1999).
fn cold_01(...) -> f64 {
```

**`scf/smearing.rs` — entropy functions:**

```rust
/// Per-state entropy weight s(x) for reduced variable x = (ε - E_F)/σ.
///
/// The total entropy is TS = σ × spin_factor × Σ_{n,k} w_k s(x_{n,k}).
///
/// Formulas by scheme:
/// - Fermi-Dirac: s = -[f ln f + (1-f) ln(1-f)]
/// - Gaussian:    s = exp(-x²) / √π
/// - Methfessel-Paxton: s = (1/2 - x²) exp(-x²) / √π
/// - Cold:        s = (x + 1/√2) exp(-(x + 1/√2)²) / √π
fn entropy_weight(...) -> f64 {
```

### B2. Struct docstrings

**`scf/mod.rs` — `ScfParams`:**

```rust
/// Parameters controlling the self-consistent field iteration.
///
/// The SCF loop solves the Kohn-Sham equations iteratively:
/// 1. Construct V_eff = V_local + V_Hartree[ρ] + V_xc[ρ]
/// 2. Diagonalize H = T + V_eff + V_NL at each k-point
/// 3. Compute occupations from eigenvalues (Fermi-Dirac or other smearing)
/// 4. Reconstruct density ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²
/// 5. Mix input and output densities (Anderson/Pulay) and repeat
///
/// Convergence requires both density (Δρ < conv_threshold) and
/// energy (ΔE < energy_threshold) criteria to be met.
pub struct ScfParams {
```

**`scf/mod.rs` — `ScfResult`:**

```rust
/// Output of a converged SCF calculation.
///
/// All energies are in eV. The three energy quantities are:
/// - `total_energy`: E = E_band - E_H + E_xc - E_vxc + E_ewald + V_local(G=0)·N_el
/// - `free_energy`: F = E - TS (Mermin functional, variational at finite σ)
/// - `energy_sigma0`: E₀ = (E + F)/2 (best estimate of T=0 energy)
pub struct ScfResult {
```

**`scf/mixing.rs` — `AndersonMixer`:**

```rust
/// Anderson/Pulay (DIIS) density mixer with optional Kerker preconditioning.
///
/// Stores a history of input densities and residuals R^(n) = ρ_out^(n) - ρ_in^(n).
/// At each step, finds coefficients c_i (summing to 1) that minimize |Σ c_i R^(i)|²
/// by solving the DIIS linear system, then constructs the new density as:
///
///   ρ_in^{n+1} = Σ_i c_i [ρ_in^(i) + β R^(i)]
///
/// where β is the mixing parameter.
///
/// With Kerker preconditioning, the residual is modified in G-space before mixing:
///   R̃(G) = [|G|² / (|G|² + q_TF²)] R(G)
///
/// This damps long-wavelength charge sloshing, which is the dominant source of
/// SCF instability in metals and large-gap systems.
pub struct AndersonMixer {
```

### B3. Functions with docstrings missing formulas

**`ewald.rs` — `ewald_energy`:**

```rust
/// Compute the Ewald ion-ion energy for a crystal. Returns energy in eV.
///
/// Decomposes the Coulomb sum of periodic point charges into four terms:
///
/// E_real  = (e²/2) Σ'_{i,j,T} Z_i Z_j erfc(η|r_ij+T|) / |r_ij+T|
/// E_recip = (2πe²/Ω) Σ_{G≠0} |S(G)|² exp(-|G|²/(4η²)) / |G|²
/// E_self  = -(η/√π) e² Σ_i Z_i²
/// E_bg    = -πe² (Σ Z_i)² / (2Ωη²)
///
/// where S(G) = Σ_i Z_i exp(iG·r_i) is the charge-weighted structure factor,
/// η = (N_atoms π/Ω)^{1/3} balances real/reciprocal cost, and primed sum
/// excludes i=j when T=0 (self-interaction).
///
/// Cutoffs: g_max = 10η (reciprocal), r_max = 10/η (real).
pub fn ewald_energy(...) -> f64 {
```

**`fft.rs` — `FFT3D` struct:**

```rust
/// 3D FFT via batched 1D transforms (z → y → x).
///
/// Convention:
///   Forward:  f̃(G) = Σ_r f(r) e^{-iG·r}     (unnormalized)
///   Inverse:  f(r) = Σ_G f̃(G) e^{+iG·r}     (unnormalized)
///
/// The forward FFT is unnormalized; callers must divide by N = nx·ny·nz
/// to get Fourier coefficients. Use `inverse_normalized()` for the
/// convention f(r) = (1/N) Σ_G f̃(G) e^{+iG·r}.
pub struct FFT3D {
```

**`fft.rs` — `fft_grid_size`:**

```rust
/// Find the smallest FFT-friendly grid size n ≥ 2·n_max + 1.
///
/// The factor 2·n_max + 1 is the Nyquist criterion: G-vectors range from
/// -n_max to +n_max, requiring at least 2·n_max + 1 grid points to avoid
/// aliasing when computing products like V(G-G') in the Hamiltonian.
///
/// Grid sizes that are products of small primes (2, 3, 5) give optimal FFT
/// performance; arbitrary sizes may be much slower.
pub fn fft_grid_size(n_max: usize) -> usize {
```

**`eigensolver/dense.rs` — `diagonalize_hermitian`:**

```rust
/// Full Hermitian eigendecomposition of H via faer.
///
/// Solves Hψ = εψ for all eigenvalues and eigenvectors.
/// Returns eigenvalues in ascending order (ε₁ ≤ ε₂ ≤ ... ≤ εₙ).
///
/// Uses faer's `self_adjoint_eigen` (dense, O(n³) LAPACK-equivalent).
/// Only the lower triangle of H is read.
pub fn diagonalize_hermitian(h: &faer::Mat<Complex64>) -> EigenResult {
```

**`potential/xc.rs` — `lda_xc_spin`:**

```rust
/// Evaluate spin-polarized LDA XC at a single point.
///
/// Exchange: ε_x^σ = -(3/4)(6ρ_σ/π)^{1/3} (fully polarized gas per channel).
/// Correlation: interpolated between unpolarized (ζ=0) and fully polarized (ζ=1)
/// using the von Barth-Hedin interpolation function:
///   f(ζ) = [(1+ζ)^{4/3} + (1-ζ)^{4/3} - 2] / [2^{4/3} - 2]
///   ε_c(r_s, ζ) = ε_c^unpol + f(ζ) [ε_c^pol - ε_c^unpol]
///
/// `rho_up`, `rho_down` in e/ų. Returns (ε_xc, V_xc↑, V_xc↓) in eV.
pub fn lda_xc_spin(rho_up: f64, rho_down: f64) -> SpinXcPoint {
```

**`basis.rs` — `BasisSet::new`:**

```rust
/// Construct the plane-wave basis set for a given energy cutoff.
///
/// Includes all reciprocal lattice vectors G = n₁b₁ + n₂b₂ + n₃b₃
/// satisfying the kinetic energy cutoff:
///   (ħ²/2m) |G|² ≤ E_cut
///
/// where b_i are reciprocal lattice vectors (2π/V × a_j × a_k).
/// The number of basis functions scales as N_pw ∝ E_cut^{3/2} × Ω.
pub fn new(lattice: &Lattice, ecut: f64) -> Self {
```

### B4. Unit conventions

**`scf/density.rs` — `compute_density`:** Add `Output density is in e/ų.`

**`crystal.rs` — `Lattice::volume`:**

```rust
/// Cell volume Ω = a · (b × c).
///
/// Returns the signed scalar triple product. Positive for right-handed
/// lattice vectors, negative for left-handed. Use `.abs()` when a
/// positive volume is needed (e.g., normalization).
pub fn volume(&self) -> f64 {
```

**`potential/nonlocal.rs` — `bessel_transform_projector`:**

```rust
/// Spherical Bessel transform of a projector:
///   F(q) = 4π ∫₀^∞ [r·β(r)] j_l(qr) r dr
///
/// `r_grid`: radial grid points in Å.
/// `rab`: integration weights dr (spacing between grid points) in Å.
///        For logarithmic grids, rab[i] = r[i] × log_step.
/// `r_beta`: r·β(r) in Å^{-1/2} (UPF convention: projectors stored as r×β).
/// `l`: angular momentum quantum number.
/// `q`: wavevector magnitude |k+G| in Å⁻¹.
fn bessel_transform_projector(...) -> f64 {
```

### B5. Comment and docstring fixes

**`potential/nonlocal.rs` lines 49-52** — fix incorrect `(-1)^l` comment:

```rust
// Before:
///     F_i(|k+G|) D_{ij} F_j(|k+G'|) × (2l+1)/(4π) P_l(cos θ) × (-1)^l
///
/// (The i^l × (i*)^l = (-1)^l factor)

// After:
///     F_i(|k+G|) D_{ij} F_j(|k+G'|) × (2l+1)/(4π) P_l(cos θ)
///
/// (Phase factors i^l from bra and (i*)^l from ket give |i|^{2l} = 1.)
```

The code correctly does NOT apply `(-1)^l`. Proof: `i^l × (i*)^l = |i|^{2l} = 1` for all `l`. Si (l=0,1 projectors) matches QE, confirming the code is correct.

**`scf/mod.rs` — `run_scf_spin` misaligned first line:**

```rust
// Before:
/// Compute local pseudopotential V_local(G) on the FULL FFT grid.
/// Spin-polarized SCF loop (nspin=2).

// After:
/// Spin-polarized SCF loop (nspin=2).
///
/// Two spin channels with independent densities, XC potentials, and
/// Hamiltonians. Hartree and V_local are spin-independent (computed from
/// total density ρ↑ + ρ↓). V_xc is spin-dependent via LSDA.
/// NLCC core charge is split equally between channels: ρ_core/2 per spin.
```

## Files Modified

| File | Changes |
|------|---------|
| `docs/theory.md` | New — comprehensive theory reference (~18 KB) |
| `docs/math-audit.md` | New — formula-by-formula verification (~11 KB) |
| `src/scf/smearing.rs` | Docstrings on 4 occupation fns + 2 entropy fns |
| `src/scf/mod.rs` | Docstrings on `ScfParams`, `ScfResult`; `run_scf_spin` fix |
| `src/scf/mixing.rs` | Docstring on `AndersonMixer` struct |
| `src/scf/density.rs` | Unit specification on `compute_density` |
| `src/scf/initial_density.rs` | Docstring on `InitialDensityConfig` |
| `src/potential/xc.rs` | Formula on `lda_xc_spin` |
| `src/potential/nonlocal.rs` | `(-1)^l` comment fix, `rab` parameter docs |
| `src/ewald.rs` | Expanded formulas on `ewald_energy` |
| `src/fft.rs` | FFT convention on `FFT3D`, Nyquist on `fft_grid_size` |
| `src/eigensolver/dense.rs` | Algorithm and complexity on `diagonalize_hermitian` |
| `src/basis.rs` | Cutoff formula on `BasisSet::new` |
| `src/crystal.rs` | Volume sign convention on `Lattice::volume` |
| `src/hamiltonian.rs` | V_eff indexing clarification on `build_hamiltonian` |

## Acceptance Criteria

1. **`docs/theory.md`** exists with all 10 sections citing canonical references.
2. **`docs/math-audit.md`** exists with all 17 formula verifications.
3. **Every function implementing a physics formula** has the formula in its docstring.
4. **Every function with dimensional inputs/outputs** states units (eV, Å, e/ų).
5. **Key algorithms** (Anderson mixing, Ewald, FFT) described at struct/function level.
6. **References cited** where non-obvious: PZ (1981), MP (1989), MV (1999).
7. **`(-1)^l` comment fixed** in `nonlocal.rs`.
8. **`run_scf_spin` docstring fixed** — no longer says "Compute local pseudopotential".
9. **No code logic changes** — documentation only, `cargo test` passes unchanged.
