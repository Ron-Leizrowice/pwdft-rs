# SCF Convergence

## The Fixed-Point Problem

The SCF cycle maps input density to output density: `ρ_out = F[ρ_in]`.
Convergence requires finding `ρ* = F[ρ*]`.

## Convergence Criteria

Dual criterion — **both** must be satisfied:

1. **Density:** `Δρ_RMS = √[(1/Ω) ∫(ρ_new - ρ_old)² dr] < conv_threshold`
2. **Energy:** `|E_new - E_old| < energy_threshold`

**Code:** `src/scf/driver.rs` (non-spin) / `src/scf/driver_spin.rs` (spin) — per-iteration `converged` check.

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

**Code:** `src/scf/mixing/anderson.rs` (DIIS coefficient assembly + constraint embedding).

### Singular fallback

If the DIIS matrix is singular (pivot < 1e-15), uniform coefficients
`c_i = 1/(n+1)` are returned. This masks convergence problems without
failing.

**Code:** `src/scf/mixing/linalg.rs` (`solve_linear_system` — uniform-coefficient fallback on singular pivot).

## Modified Broyden Mixing

Alternative to DIIS. Builds an approximate inverse Jacobian `J^{-1}` from the
history of density residuals and applies Newton-style updates:

```
ρ_in^{n+1} = ρ_in^n - β J_n^{-1} R^n
```

Implements Johnson's modified Broyden scheme (*Phys. Rev. B* **38**, 12807,
1988) — the same algorithm as QE's `mix_rho.f90` and VASP's `IMIX=4`. Often
more robust than Anderson for difficult systems (metals, large cells,
charge sloshing). Selectable via `mixing_mode: broyden` or
`broyden_kerker`.

**Code:** `src/scf/mixing/broyden.rs`.

## Periodic Pulay Mixing

Banerjee, Suryanarayana, Pask, *J. Chem. Theory Comput.* **12**, 3053
(2016). Plain linear mixing on every iteration *except* every k-th, where a
DIIS extrapolation is performed using the accumulated history. Avoids
divergence from early-iteration DIIS (when history is too short to be
reliable) while keeping the acceleration once enough residuals have
accumulated. Paper reports 30–50% iteration-count reduction on
transition-metal oxides versus continuous Anderson. Selectable via
`mixing_mode: periodic_pulay` or `periodic_pulay_kerker`; period k is set
by `pulay_period` (default 3).

**Code:** `src/scf/mixing/anderson.rs` (`PeriodicPulayMixer` wraps
`AndersonMixer`).

## Kerker Preconditioning

Damps long-wavelength charge sloshing:

```
R̃(G) = [|G|² / (|G|² + q_TF²)] × R(G)
```

- G=0 completely suppressed (set to 0.0)
- `q_TF² = 4(3π²ρ)^{1/3}/π` in a.u., converted to Å⁻²

This is the Thomas-Fermi screening length from Ashcroft & Mermin (1976), Ch. 17.

**Code:** `src/scf/mixing/kerker.rs` (`precondition_residual` and `auto_q_tf_squared`).

## Audit Status

| Item | Status |
|------|--------|
| DIIS coefficient solving | CORRECT |
| Constraint embedding | CORRECT |
| Kerker P(G) = G²/(G²+q_TF²) | CORRECT |
| q_TF auto-estimate | CORRECT (free-electron formula) |
| G=0 suppression | CORRECT |
| Singular matrix fallback | CORRECT (logs warning) |
