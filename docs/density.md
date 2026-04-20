# Electron Density

## Density from Wavefunctions

```text
ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²
```

where `f_{n,k}` is the occupation (includes spin factor: 0-2 for nspin=1,
0-1 for nspin=2) and `w_k` is the k-point weight (sum to 1).

### Procedure

1. Place PW coefficients `c_{n,k}(G)` onto FFT grid at positions `G`
2. Inverse FFT (unnormalized) to get `ψ_{n,k}(r)`
3. Accumulate `ρ(r) += f_{n,k} × w_k × |ψ(r)|²`
4. Normalize so `∫ ρ(r) dr = (Ω/N_grid) Σ_r ρ(r) = N_electrons`

**Code:** `src/scf/density.rs:30-100`

### Parallel reduction

Uses rayon fold-reduce: each thread accumulates into a local buffer, then
buffers are summed. Avoids allocating n_kpoints intermediate arrays.

### Normalization note

The inverse FFT is called unnormalized (`fft.inverse()`, not `inverse_normalized()`),
so intermediate |ψ|² values are scaled by N_grid². The final normalization
(scaling to N_electrons) corrects this. Mathematically equivalent to using
`inverse_normalized()` followed by the same rescaling, but intermediate values
are larger.

## Initial Density (Superposition of Atomic Densities)

For each atom with pseudopotential data `4πr²ρ_atom(r)`:

```text
ρ_atom(G) = (1/Ω) ∫ [4πr²ρ(r)] j₀(|G|r) dr × S(G)
```

UPF `PP_RHOATOM` stores `4πr²ρ(r)`, so no extra factor is needed.

If PP_RHOATOM is unavailable, falls back to a Gaussian model:

```text
ρ_atom(G) = (Z_val/Ω) exp(-|G|²σ²/2) × S(G)
```

with σ = 1.0 Å.

**Code:** `src/scf/initial_density.rs`

### Spin-polarized initialization

```text
ρ_up = (1+m)/2 × ρ
ρ_down = (1-m)/2 × ρ
```

where m is the starting magnetization fraction per atom type.

## Density Clamping

After FFT, small negative values can appear due to aliasing.
These are clamped to 0.0 before XC evaluation.

**Code:** `src/scf/initial_density.rs:94-97`

## Audit Status

| Item | Status |
|------|--------|
| ρ = Σ f w | ψ |
| FFT + normalization | CORRECT (rescaled at end) |
| Parallel reduction | CORRECT |
| SAD Bessel transform | CORRECT (formula) |
| SAD Gaussian fallback | CORRECT |
| Spin initialization | CORRECT |
| **SAD quadrature** | **ISSUE** — see [Radial Integration](radial-integration.md) |
