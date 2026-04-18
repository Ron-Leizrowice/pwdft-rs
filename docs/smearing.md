# Smearing and Entropy

## Occupation Functions

All return occupation in [0, 1] for a single state. Multiply by `spin_factor`
(2/nspin) for the full occupation. Argument: `x = (ε - E_F) / σ`.

| Scheme | f(x) |
|--------|------|
| **Fermi-Dirac** | `1 / (1 + e^x)` |
| **Gaussian** | `erfc(x) / 2` |
| **Methfessel-Paxton (order 1)** | `erfc(x)/2 - (1/2) x exp(-x²)/√π` |
| **Marzari-Vanderbilt (cold)** | `(1/2) erfc(x + 1/√2) + exp(-(x+1/√2)²)/√(2π)` |

Overflow protection: `|x| > 40` is clamped.

**Code:** `src/scf/smearing.rs:100-150`

## Entropy Weights

Contribution per state: `S_n = σ × s(x)`.

| Scheme | s(x) |
|--------|------|
| **Fermi-Dirac** | `-(f ln f + (1-f) ln(1-f))` |
| **Gaussian** | `exp(-x²)/√π` |
| **Methfessel-Paxton** | `(1/2 - x²) exp(-x²)/√π` |
| **Cold** | `(x + 1/√2) exp(-(x+1/√2)²)/√π` |

**Code:** `src/scf/smearing.rs:180-205`

## Free Energy and σ→0 Extrapolation

```
T×S = σ × spin_factor × Σ_{n,k} w_k × s(x_{n,k})
F = E - T×S                    (Mermin functional, variational at finite σ)
E₀ = (E + F) / 2 = E - T×S/2  (best estimate of T=0 energy)
```

The σ→0 formula (Gillan-De Vita-Caro) is exact for Methfessel-Paxton and a
good approximation for Fermi-Dirac and Gaussian.

**Code:** `src/scf/driver.rs` / `src/scf/driver_spin.rs` (per-driver final summary assembles `free_energy` and `energy_sigma0` into the returned `ScfResult`).

## Fermi Energy Search

Bisection: find `E_F` such that `N_el = Σ_{n,k} w_k f(ε, E_F, σ) × spin_factor`.

- Bounds: `E_min = ε_min - 10σ`, `E_max = ε_max + 10σ`
- 200 iterations → precision far beyond machine epsilon
- Tolerance: `(E_max - E_min) < 1e-14`

**Code:** `src/scf/smearing.rs:50-91`

## Audit Status

| Item | Status |
|------|--------|
| Fermi-Dirac occupation | CORRECT |
| Gaussian occupation | CORRECT |
| Methfessel-Paxton occupation | CORRECT |
| Cold (Marzari-Vanderbilt) occupation | CORRECT |
| All entropy formulas | CORRECT |
| Fermi energy bisection | CORRECT |
| σ→0 extrapolation | CORRECT |
| Overflow protection | CORRECT |
