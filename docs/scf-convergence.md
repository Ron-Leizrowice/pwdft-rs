# SCF Convergence

## The Fixed-Point Problem

The SCF cycle maps input density to output density: `ρ_out = F[ρ_in]`.
Convergence requires finding `ρ* = F[ρ*]`.

## Convergence Criteria

Dual criterion — **both** must be satisfied:

1. **Density:** `Δρ_RMS = √[(1/Ω) ∫(ρ_new - ρ_old)² dr] < conv_threshold`
2. **Energy:** `|E_new - E_old| < energy_threshold`

**Code:** `src/scf/mod.rs:288-289`

## Linear Mixing

```
ρ_in^{n+1} = ρ_in^n + α R^n
```

where `R^n = ρ_out^n - ρ_in^n` and `α = mixing_beta`.

## Anderson/Pulay (DIIS) Mixing

Construct optimal linear combination of history:

```
ρ̄_in = Σ_i c_i ρ_in^{n-m+i},   Σ c_i = 1
```

Minimize `|R̄|² = Σ_{ij} c_i c_j ⟨R^i|R^j⟩` subject to the constraint.

### Constraint embedding

The constraint `Σ c_i = 1` is embedded by eliminating the last coefficient:
`α_last = 1 - Σ α_prev`. The system solves for m-1 coefficients via:

```
A[i,j] = ΔR_i · ΔR_j     (where ΔR_i = R_i - R_last)
b[i] = -ΔR_i · R_last
```

This is equivalent to minimizing `|Σ c_i R_i|²` subject to `Σ c_i = 1`.

**Code:** `src/scf/mixing.rs:140-160`

### Singular fallback

If the DIIS matrix is singular (pivot < 1e-15), uniform coefficients
`c_i = 1/(n+1)` are returned. This masks convergence problems without
failing.

**Code:** `src/scf/mixing.rs:217-218`

## Kerker Preconditioning

Damps long-wavelength charge sloshing:

```
R̃(G) = [|G|² / (|G|² + q_TF²)] × R(G)
```

- G=0 completely suppressed (set to 0.0)
- `q_TF² = 4(3π²ρ)^{1/3}/π` in a.u., converted to Å⁻²

This is the Thomas-Fermi screening length from Ashcroft & Mermin (1976), Ch. 17.

**Code:** `src/scf/mixing.rs:60-80` (preconditioner), lines 187-198 (q_TF auto-estimate)

## Audit Status

| Item | Status |
|------|--------|
| DIIS coefficient solving | CORRECT |
| Constraint embedding | CORRECT |
| Kerker P(G) = G²/(G²+q_TF²) | CORRECT |
| q_TF auto-estimate | CORRECT (free-electron formula) |
| G=0 suppression | CORRECT |
| Singular matrix fallback | CORRECT (logs warning) |
