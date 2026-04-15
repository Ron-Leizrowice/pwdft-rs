# Proposal 27: Mathematical Documentation in Code Docstrings

## Problem

The codebase implements ~50 physics-critical functions across 12 modules. While some have excellent docstrings (e.g., `add_atomic_density_from_pp` shows the full Bessel transform formula, `ScfResult` fields document the Mermin functional), many functions that implement non-trivial mathematics lack formula documentation, unit specifications, or both. This makes the code harder to audit, debug, and extend — the Fe NLCC bug (Proposal 25) might have been caught earlier if the XC evaluation docstring had explicitly stated what density it operates on.

A full audit was performed across every physics function in the codebase. The findings are organized by severity.

## Audit Results

### Functions With No Docstring on Non-Trivial Math

| Function | File | Line | What it computes |
|----------|------|------|-----------------|
| `fermi_dirac_01` | `scf/smearing.rs` | 107 | f(ε) = 1/(1 + exp((ε-E_F)/σ)) |
| `gaussian_01` | `scf/smearing.rs` | 117 | f(ε) = erfc((ε-E_F)/σ)/2 |
| `methfessel_paxton_01` | `scf/smearing.rs` | 127 | f(ε) = erfc(x)/2 - x exp(-x²)/(2√π) |
| `cold_01` | `scf/smearing.rs` | 139 | f(ε) = erfc(x+1/√2)/2 + exp(-(x+1/√2)²)/√(2π) |
| `ScfParams` (struct) | `scf/mod.rs` | 31 | SCF iteration parameters |
| `ScfResult` (struct) | `scf/mod.rs` | 81 | Converged SCF output |
| `AndersonMixer` (struct) | `scf/mixing.rs` | 24 | Pulay/DIIS density mixer |
| `InitialDensityConfig` (struct) | `scf/initial_density.rs` | 31 | SAD configuration |

### Functions With Docstrings But Missing Math Formulas

| Function | File | Line | What's missing |
|----------|------|------|---------------|
| `run_scf` | `scf/mod.rs` | 180 | SCF algorithm overview (steps, convergence criteria) |
| `AndersonMixer::mix` | `scf/mixing.rs` | 82 | Anderson/Pulay mixing formula |
| `entropy_ts` | `scf/smearing.rs` | 158 | TS = σ × spin_factor × Σ w_k s(x_nk) |
| `entropy_weight` | `scf/smearing.rs` | 182 | Per-scheme entropy formulas |
| `lda_xc_spin` | `potential/xc.rs` | 166 | Spin-polarized XC formula |
| `ewald_energy` | `ewald.rs` | 21 | Individual component formulas (real, recip, self, bg) |
| `diagonalize_hermitian` | `eigensolver/dense.rs` | 18 | Algorithm and complexity |
| `FFT3D::forward` | `fft.rs` | 48 | FFT definition f̃(G) = Σ f(r) e^{-iG·r} |
| `FFT3D::inverse` | `fft.rs` | 64 | IFFT definition and normalization |
| `fft_grid_size` | `fft.rs` | 91 | Nyquist criterion: n ≥ 2n_max + 1 |
| `BasisSet::new` | `basis.rs` | 19 | Cutoff: (ħ²/2m)\|k+G\|² ≤ E_cut |
| `build_hamiltonian` | `hamiltonian.rs` | 33 | V_eff indexing convention |

### Functions With Docstrings But Missing Unit Conventions

| Function | File | What's missing |
|----------|------|---------------|
| `compute_density` | `scf/density.rs` | Output units (e/ų) |
| `generate_initial_density` | `scf/initial_density.rs` | Output units (e/ų) |
| `precondition_residual` | `scf/mixing.rs` | FFT normalization convention |
| `slater_exchange` | `potential/xc.rs` | Input/output unit conversion chain |
| `pz_correlation_spin` | `potential/xc.rs` | ζ clamp behavior |
| `Lattice::volume` | `crystal.rs` | Sign convention (can be negative) |

### Misaligned Docstring

| Function | File | Line | Issue |
|----------|------|------|-------|
| `run_scf_spin` | `scf/mod.rs` | 456 | First line says "Compute local pseudopotential" — should say "Spin-polarized SCF loop" |

### Incorrect Comment (from Proposal 26)

| Location | File | Line | Issue |
|----------|------|------|-------|
| KB phase factor | `potential/nonlocal.rs` | 50-52 | Claims `(-1)^l`, should be `1` |

## Implementation

The changes are documentation-only — no logic changes. Each section below shows the exact docstring to add or modify.

### 1. Smearing functions (`scf/smearing.rs`)

Add docstrings with formulas and references to each private occupation function:

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

Add formulas to entropy functions:

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

### 2. Structs (`scf/mod.rs`)

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

/// Output of a converged SCF calculation.
///
/// All energies are in eV. The three energy quantities are:
/// - `total_energy`: E = E_band - E_H + E_xc - E_vxc + E_ewald + V_local(G=0)·N_el
/// - `free_energy`: F = E - TS (Mermin functional, variational at finite σ)
/// - `energy_sigma0`: E₀ = (E + F)/2 (best estimate of T=0 energy)
pub struct ScfResult {
```

### 3. Anderson mixer (`scf/mixing.rs`)

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

### 4. Ewald summation (`ewald.rs`)

Expand the existing docstring on `ewald_energy`:

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

### 5. FFT (`fft.rs`)

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

Add to `fft_grid_size`:

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

### 6. Eigensolver (`eigensolver/dense.rs`)

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

### 7. XC spin-polarized (`potential/xc.rs`)

Add formula to `lda_xc_spin`:

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

### 8. Basis set (`basis.rs`)

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

### 9. Nonlocal comment fix (from Proposal 26)

In `potential/nonlocal.rs`, lines 49-52:

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

### 10. Nonlocal `rab` explanation

In `bessel_transform_projector` docstring, add:

```rust
/// Spherical Bessel transform of a projector:
///   F(q) = 4π ∫₀^∞ [r·β(r)] j_l(qr) r dr
///
/// `r_grid`: radial grid points in Å.
/// `rab`: integration weights dr/di (spacing between consecutive grid points) in Å.
///        For logarithmic grids, rab[i] = r[i] × log_step.
/// `r_beta`: r·β(r) in Å^{-1/2} (UPF convention: projectors stored as r×β).
/// `l`: angular momentum quantum number.
/// `q`: wavevector magnitude |k+G| in Å⁻¹.
fn bessel_transform_projector(...) -> f64 {
```

### 11. Crystal volume sign (`crystal.rs`)

```rust
/// Cell volume Ω = a · (b × c).
///
/// Returns the signed scalar triple product. Positive for right-handed
/// lattice vectors, negative for left-handed. Use `.abs()` when a
/// positive volume is needed (e.g., normalization).
pub fn volume(&self) -> f64 {
```

### 12. Density output units (`scf/density.rs`)

```rust
/// Compute the charge density on the real-space FFT grid from wavefunctions.
///
/// ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²
///
/// Steps per (k-point, band):
/// 1. Place PW coefficients c_{n,k}(G) onto FFT grid
/// 2. Inverse FFT → ψ_{n,k}(r) (unnormalized)
/// 3. Accumulate f × w × |ψ|²
///
/// The result is normalized so that ∫ρ(r)dr = N_electrons.
/// Output density is in e/ų (electrons per cubic Ångström).
///
/// K-point contributions are computed in parallel via rayon fold/reduce.
pub fn compute_density(...) -> Vec<f64> {
```

### 13. `run_scf_spin` misaligned docstring fix

```rust
// Before (line 456):
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
| `src/scf/smearing.rs` | Docstrings on 4 occupation functions + 2 entropy functions |
| `src/scf/mod.rs` | Docstrings on `ScfParams`, `ScfResult`, `run_scf_spin` fix |
| `src/scf/mixing.rs` | Docstring on `AndersonMixer` struct |
| `src/scf/density.rs` | Unit specification on `compute_density` |
| `src/scf/initial_density.rs` | Docstring on `InitialDensityConfig` |
| `src/potential/xc.rs` | Formula on `lda_xc_spin` |
| `src/potential/nonlocal.rs` | Comment fix ((-1)^l → 1), `rab` docs |
| `src/ewald.rs` | Expanded formulas on `ewald_energy` |
| `src/fft.rs` | FFT convention on `FFT3D`, Nyquist on `fft_grid_size` |
| `src/eigensolver/dense.rs` | Algorithm and complexity on `diagonalize_hermitian` |
| `src/basis.rs` | Cutoff formula on `BasisSet::new` |
| `src/crystal.rs` | Volume sign convention on `Lattice::volume` |
| `src/hamiltonian.rs` | V_eff indexing clarification on `build_hamiltonian` |

## Acceptance Criteria

1. **Every function implementing a physics formula** has the formula in its docstring.
2. **Every function with dimensional inputs/outputs** states units (eV, Å, e/ų, etc.).
3. **Key algorithms** (Anderson mixing, Ewald, FFT convention) are described at the struct or function level, not just in scattered inline comments.
4. **References cited** where non-obvious: Perdew-Zunger (1981), Methfessel-Paxton (1989), Marzari-Vanderbilt (1999).
5. **No code changes** — documentation only, so `cargo test` passes unchanged.
6. **Incorrect comment fixed** — `(-1)^l` removed from nonlocal.rs.
7. **Misaligned docstring fixed** — `run_scf_spin` title corrected.
