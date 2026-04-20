# Ewald Summation (Ion-Ion Energy)

## Four-Term Decomposition

```text
E_ewald = E_recip + E_real + E_self + E_bg
```

### Reciprocal space

```text
E_recip = (2πe²/Ω) Σ_{G≠0} |S(G)|² exp(-|G|²/(4η²)) / |G|²
```

where `S(G) = Σ_i Z_i exp(iG·r_i)` is the charge-weighted structure factor.

**Code:** `src/ewald.rs:56-85`
Cutoff: `g_max = 10η`

### Real space

```text
E_real = (e²/2) Σ_T Σ'_{i,j} Z_i Z_j erfc(η|r_ij + T|) / |r_ij + T|
```

Prime excludes i=j when T=0.

**Code:** `src/ewald.rs:87-114`
Cutoff: `r_max = 10/η`. Self-interaction excluded via `r_norm < 1e-10`.

### Self-energy correction

```text
E_self = -(η/√π) e² Σ_i Z_i²
```

**Code:** `src/ewald.rs:118`

### Background charge

```text
E_bg = -πe² (Σ_i Z_i)² / (2Ωη²)
```

Neutralizes the divergent G=0 term. Nonzero even for neutral cells.

**Code:** `src/ewald.rs:122`

## Screening Parameter

```text
η = (N_atoms × π / Ω)^{1/3}
```

Balances computational cost between real and reciprocal sums.

**Code:** `src/ewald.rs:53`

## erfc Implementation

Uses `puruspe::erfc()`. Validated against NIST DLMF to 1e-7 relative precision.

## Audit Status

| Item | Status | Verified against |
|------|--------|-----------------|
| Reciprocal term | CORRECT | Standard Ewald |
| Real-space term | CORRECT | Standard Ewald |
| Self-energy | CORRECT | Standard Ewald |
| Background | CORRECT | Standard Ewald |
| η parameter | CORRECT | Standard choice |
| erfc values | CORRECT | NIST DLMF |
| NaCl Madelung test | CORRECT | M=1.747565 within 1% |
| Fe energy vs QE | CORRECT | Match to 0.006 eV |

### Notes

- Cutoffs (10η, 10/η) are standard and converge well for typical bulk systems.
- For highly anisotropic cells (slabs), the scalar cutoff may not be optimal.
  Consider aspect-ratio-weighted cutoffs for future 2D/slab support.
