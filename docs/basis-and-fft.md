# Basis Set and FFT

## Plane-Wave Basis

Wavefunctions are expanded in plane waves:

```
ψ_{n,k}(r) = (1/√Ω) Σ_G c_{n,k}(G) e^{i(k+G)·r}
```

where G are reciprocal lattice vectors `G = n₁b₁ + n₂b₂ + n₃b₃`.

### Energy cutoff

Include all G-vectors satisfying:

```
(ℏ²/2m) |k+G|² ≤ E_cut
```

At k=0: `|G|² ≤ E_cut / HBAR2_OVER_2M` where `HBAR2_OVER_2M ≈ 3.81 eV·Å²`.

**Code:** `src/basis.rs:29` (g_max_sq), lines 34-50 (enumeration)

### G-vector storage

Stored in a Vec with HashMap index for O(1) lookup by Miller indices `(n1, n2, n3)`.

**Code:** `src/basis.rs` (BasisSet struct)

## Reciprocal Lattice

```
b₁ = 2π(a₂ × a₃) / Ω
b₂ = 2π(a₃ × a₁) / Ω
b₃ = 2π(a₁ × a₂) / Ω
```

Uses the **signed** triple product `Ω = a₁·(a₂ × a₃)` to get correct
reciprocal vector orientation. Volume function returns `|Ω|` (unsigned).

**Code:** `src/crystal.rs:49-58` (reciprocal), line 46 (volume)

## Kohn-Sham Hamiltonian

```
H_{G,G'}(k) = T_{G,G'}(k) + V_eff(G-G')
```

- Kinetic: `T_{G,G'} = (ℏ²/2m)|k+G|² δ_{G,G'}` (diagonal)
- Potential: `V_eff(G-G')` looked up from FFT grid via Miller index difference

**Code:** `src/scf/potentials.rs:110-136` (Hamiltonian), `src/hamiltonian.rs` (kinetic-only)

### Miller index lookup

`miller_to_idx(dims, n1, n2, n3)` maps (possibly negative) Miller indices
to FFT grid indices using modular wrapping.

**Code:** `src/scf/grid.rs:85-90`

## FFT Convention

```
Forward:   f̃(G) = Σ_r f(r) e^{-iG·r}        (unnormalized)
Inverse:   f(r) = Σ_G f̃(G) e^{+iG·r}        (unnormalized)
Normalized: f̃(G) = (1/N) × Forward[f(r)]     (divide by N_grid)
```

The 1/N normalization on the forward FFT means:
- `ρ̃(G=0) = (1/N) Σ_r ρ(r) = ⟨ρ⟩` (spatial average)
- `V(G-G') = (1/N) Σ_r V(r) e^{-i(G-G')·r}` (Fourier coefficient in eV)

With this convention, V_eff(G-G') in the Hamiltonian has the correct units
(eV) without extra volume factors.

**Code:** `src/fft.rs`

### Grid size selection

```
n ≥ 2·n_max + 1
```

where n_max is the maximum Miller index. Grid sizes are restricted to
products of {2, 3, 5} for FFT efficiency. The `ecutrho_ratio` parameter
(default 4) scales the grid relative to the wavefunction cutoff.

**Code:** `src/fft.rs` (fft_grid_size), `src/scf/grid.rs:27-41` (FftGrid::new)

### Parseval's theorem

Verified in tests: `Σ|F(G)|² = N × Σ|f(r)|²`

## Eigensolver

Dense Hermitian eigendecomposition via `faer::self_adjoint_eigen(Side::Lower)`.
Returns eigenvalues in ascending order. `diagonalize_lowest(n_bands)` keeps
only the first n_bands eigenpairs.

**Code:** `src/eigensolver/dense.rs`

## Audit Status

| Item | Status |
|------|--------|
| Energy cutoff criterion | CORRECT |
| G-vector enumeration | CORRECT |
| Reciprocal lattice formula | CORRECT |
| Volume (signed triple product) | CORRECT |
| Hamiltonian assembly (T + V_eff) | CORRECT |
| Miller index lookup | CORRECT |
| FFT forward/inverse convention | CORRECT |
| FFT normalization (1/N) | CORRECT, consistent throughout |
| Grid size selection (Nyquist) | CORRECT |
| Parseval's theorem | CORRECT |
| Eigensolver (faer) | CORRECT |
