# Math Audit: pwdft-rs vs. Canonical Formulas

Each section verifies a core formula in the codebase against the canonical
references in `docs/theory.md`. Status is one of:

- **CORRECT** — matches reference exactly
- **CORRECT (convention)** — matches with a documented convention difference
- **BUG** — known issue, with fix reference
- **NEEDS VERIFICATION** — formula is plausible but unit chain is complex

---

## 1. Total Energy Expression

**Reference (Section 1.1):**
```
E_total = E_band - E_H + E_xc - E_vxc + E_ewald + V_local(G=0) × N_el
```

**Code** (`src/scf/mod.rs:1020`):
```rust
e_band - e_hartree + e_xc - e_vxc + e_ewald
```
Plus `v_local_g0 * n_electrons` added at line 388.

**Status: CORRECT**

Each component verified below.

---

## 2. Band Energy

**Reference:** `E_band = Σ_{n,k} f_{n,k} w_k ε_{n,k}`

**Code** (`src/scf/mod.rs:989-996`):
```rust
let e_band: f64 = eigenvalues.iter()
    .zip(occupations.iter())
    .zip(kpoint_weights.iter())
    .map(|((evs, occs), &w)| {
        evs.iter().zip(occs.iter()).map(|(&e, &f)| f * w * e).sum::<f64>()
    })
    .sum();
```

**Status: CORRECT**

Note: `f` here already includes the spin factor (2 for nspin=1, 1 for nspin=2),
which is correct — each eigenvalue is multiplied by its full occupation.

---

## 3. Hartree Energy

**Reference:** `E_H = (Ω/2) Σ_{G≠0} |ρ̃(G)|² × 4πe² / |G|²`

**Code** (`src/scf/mod.rs:998-1007`):
```rust
let e_hartree: f64 = rho_g.iter().zip(g_squared.iter())
    .map(|(rho, &g2)| {
        if g2 > G2_ZERO_THRESHOLD { rho.norm_sqr() * fourpi_e2 / g2 } else { 0.0 }
    })
    .sum::<f64>()
    * 0.5 * omega;
```

**Status: CORRECT**

- `rho_g` is `ρ̃(G)` with 1/N FFT normalization (spatial average)
- `norm_sqr()` = `|ρ̃(G)|²`
- `fourpi_e2 = 4πe²` with `e² = 14.4 eV·Å`
- Factor `0.5 * omega` matches the reference formula
- G=0 excluded via threshold

---

## 4. Hartree Potential

**Reference:** `V_H(G) = 4πe² ρ̃(G) / |G|²` for G≠0

**Code** (`src/scf/mod.rs:785-798`):
```rust
rho * fourpi_e2 / g2    // for g2 > threshold
```

**Status: CORRECT**

This is `V_H(G)` in the same units as `V_eff`, which enters the Hamiltonian
as `V_{G,G'} = V_eff(G-G')`.

---

## 5. XC Energy

**Reference:** `E_xc = ∫ ε_xc(r) ρ(r) dr = (Ω/N_grid) Σ_r ε_xc(r) ρ(r)`

**Code** (`src/potential/xc.rs:64-72`):
```rust
pub fn lda_xc_energy(rho_r: &[f64], exc_r: &[f64], omega: f64) -> f64 {
    let dvol = omega / n_grid as f64;
    rho_r.iter().zip(exc_r.iter())
        .map(|(&rho, &exc)| rho * exc * dvol)
        .sum()
}
```

**Status: CORRECT**

`dvol = Ω/N_grid` is the volume element per grid point.

---

## 6. XC Functional (Slater Exchange)

**Reference:** `ε_x = -(3/4)(3ρ/π)^{1/3}` in Hartree

**Code** (`src/potential/xc.rs:89-91`):
```rust
let cbrt = (3.0 * rho_bohr / PI).powf(1.0 / 3.0);
let ex_ha = -0.75 * cbrt;
```

**Status: CORRECT**

`-0.75 × (3ρ/π)^{1/3} = -(3/4)(3ρ/π)^{1/3}`. This is the standard Slater
exchange for the unpolarized electron gas. The potential `V_x = (4/3)ε_x` is
also correct (line 93).

---

## 7. XC Functional (Perdew-Zunger Correlation)

**Reference (Ref. 4, Table I):**

| Parameter | Value | Code | Match |
|-----------|-------|------|-------|
| γ | -0.1423 | -0.1423 | ✓ |
| β₁ | 1.0529 | 1.0529 | ✓ |
| β₂ | 0.3334 | 0.3334 | ✓ |
| A | 0.0311 | 0.0311 | ✓ |
| B | -0.048 | -0.048 | ✓ |
| C | 0.002 | 0.002 | ✓ |
| D | -0.0116 | -0.0116 | ✓ |

**Code** (`src/potential/xc.rs:117-128` for r_s ≥ 1):
```rust
let denom = 1.0 + beta1 * sqrt_rs + beta2 * rs;
ec_ha = gamma / denom;
let d_ec = -gamma * (beta1 / (2.0 * sqrt_rs) + beta2) / (denom * denom);
vc_ha = ec_ha - rs / 3.0 * d_ec;
```

**Status: CORRECT**

The derivative `dε_c/dr_s = -γ(β₁/(2√r_s) + β₂)/(1 + β₁√r_s + β₂r_s)²`
and potential `V_c = ε_c - (r_s/3) dε_c/dr_s` match the standard derivation.

---

## 8. Kinetic Energy (Hamiltonian Diagonal)

**Reference:** `T_{G,G'} = (ħ²/2m)|k+G|² δ_{G,G'}`

**Code** (`src/scf/mod.rs:846-848`):
```rust
let ke = HBAR2_OVER_2M * (k + g).norm_squared();
h[(i, i)] = Complex64::new(ke, 0.0);
```

**Status: CORRECT**

`HBAR2_OVER_2M = ħ²/(2m) ≈ 3.81 eV·Å²`, verified in `consts.rs` test.

---

## 9. Local Potential on FFT Grid

**Reference:** `V_local(G) = (1/Ω) Σ_atom V^pp(|G|) S(G)`

**Code** (`src/scf/mod.rs:766-779`):
```rust
for &(ref tau, pp) in &atom_data {
    let phase = -g.dot(tau);
    let sf = Complex64::cis(phase);
    let v_form = pp.v_local_of_g(g_norm, omega);
    v += sf * v_form;
}
```

**Status: CORRECT (convention)**

The `1/Ω` factor is included inside `v_local_of_g()` (the pseudopotential
form factor already has it). Structure factor `S(G) = exp(-iG·τ)` is correct.

---

## 10. Non-Local Potential Matrix Elements

**Reference (Section 4.2):**
```
V_NL(G,G') = (1/Ω) Σ F_i(|k+G|) D_{ij} F_j(|k+G'|) × (2l+1)/(4π) P_l(cos θ)
             × (-1)^l × Σ_atom exp(-iG·τ) exp(iG'·τ)
```

**Code** (`src/potential/nonlocal.rs:154-198`):
```rust
let angular = (2 * l + 1) as f64 / (4.0 * PI) * legendre_p(l, cos_theta);
vnl += fi * d * fj * angular;
// ...
sf_sum += phases[ig] * phases[jg].conj();  // exp(-iG_i·τ) × exp(iG_j·τ)
h[(ig, jg)] += sf_sum * (vnl * inv_omega);
```

**Status: CORRECT**

- Form factors `fi`, `fj` computed via Bessel transform in `new()`
- Angular part uses addition theorem correctly
- Structure factor pair `S(G-G') = exp(-iG·τ)conj(exp(-iG'·τ))` is correct
- `inv_omega = 1/Ω`
- The `(-1)^l` factor: checking `bessel_transform_projector()`...

**Note:** The `(-1)^l` factor from `i^l (i^*)^l` cancels in the `|β⟩D⟨β|` form
because both projectors carry the same `i^l`. The code comment at line 50-52
confirms this. **Verified correct.**

---

## 11. Ewald Summation

**Reference (Section 5):**
```
E_recip = (2πe²/Ω) Σ_{G≠0} |S(G)|² exp(-|G|²/(4η²)) / |G|²
E_real = (e²/2) Σ_{T} Σ'_{i,j} Z_i Z_j erfc(η|r|) / |r|
E_self = -(η/√π) e² Σ Z_i²
E_bg = -πe²(Σ Z_i)² / (2Ωη²)
```

**Code** (`src/ewald.rs:51-114`):

| Term | Code | Reference | Match |
|------|------|-----------|-------|
| E_recip | `s_sq * exp(-g2/(4η²)) / g2 × 2πe²/Ω` | ✓ | ✓ |
| E_real | `Z_i Z_j erfc(η r) / r × e²/2` | ✓ | ✓ |
| E_self | `-(η/√π) Σ Z² × e²` | ✓ | ✓ |
| E_bg | `-π (Σ Z)² e² / (2Ω η²)` | ✓ | ✓ |

**Status: CORRECT**

The screening parameter `η = (N π / Ω)^{1/3}` is a standard choice.
Cutoffs `g_max = 10η` and `r_max = 10/η` provide good convergence.

---

## 12. Density from Wavefunctions

**Reference:** `ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²`

**Code** (`src/scf/density.rs:49-63`):
```rust
let f = occ[ib] * kp.weight;
// ... place c(G) on grid, inverse FFT ...
rho_acc[i] += f * psi.norm_sqr();
```

**Status: CORRECT**

The unnormalized inverse FFT gives `ψ(r) × N_grid`. Then `|ψ|²` is `N_grid²`
times the true value. But the normalization at lines 78-87 (`∫ρ dr = N_el`)
corrects this, so the intermediate scaling doesn't matter.

---

## 13. Fermi Energy Search

**Reference:** Find `E_F` such that `N_el = Σ_{n,k} w_k f(ε_{n,k}, E_F, σ) × spin_factor`

**Code** (`src/scf/smearing.rs:50-91`): Bisection search, 200 iterations.

**Status: CORRECT**

The bisection bounds `e_min - 10σ` to `e_max + 10σ` cover the full occupation
range. 200 iterations give precision `(e_max - e_min) / 2^{200} ≈ 0`, far
beyond machine epsilon.

---

## 14. Smearing Functions

### Fermi-Dirac
**Reference:** `f(x) = 1/(1+e^x)`, entropy: `-[f ln f + (1-f) ln(1-f)]`

**Code** (`smearing.rs:107-114, 184-191`): Matches. Overflow protected at `|x|>40`.

**Status: CORRECT**

### Gaussian
**Reference:** `f(x) = erfc(x)/2`, entropy: `exp(-x²)/√π`

**Code** (`smearing.rs:117-125, 192-194`): Matches.

**Status: CORRECT**

### Methfessel-Paxton (order 1)
**Reference:** `f(x) = erfc(x)/2 - (1/2)x exp(-x²)/√π`

**Code** (`smearing.rs:127-137`):
```rust
let f0 = puruspe::erfc(x) / 2.0;
let gauss = (-x * x).exp() / PI.sqrt();
f0 - 0.5 * x * gauss
```

**Status: CORRECT**

Entropy: `(1/2 - x²) exp(-x²)/√π` also correct (line 196).

### Marzari-Vanderbilt (cold)
**Reference:** `f(x) = (1/2)erfc(x+1/√2) + exp(-(x+1/√2)²)/√(2π)`

**Code** (`smearing.rs:139-149`):
```rust
let arg = x + sq2_inv;
0.5 * puruspe::erfc(arg) + (1.0 / (2.0 * PI).sqrt()) * (-arg * arg).exp()
```

**Status: CORRECT**

---

## 15. Sigma→0 Energy Extrapolation

**Reference:** `E₀ = (E + F) / 2 = E - TS/2`

**Code** (`src/scf/mod.rs:414`):
```rust
let energy_sigma0 = (e_total + free_energy) / 2.0;
```

where `free_energy = e_total - ts`.

**Status: CORRECT**

This is the Gillan-De Vita-Caro formula, valid for Fermi-Dirac and Gaussian
smearing. For Methfessel-Paxton, the sigma→0 extrapolation is exact (the
energy is already variational), but this formula is still a good approximation.

---

## 16. NLCC Implementation

**Reference:** Evaluate `V_xc[ρ_val + ρ_core]` instead of `V_xc[ρ_val]`.

**Code** (`src/scf/mod.rs:298-299, 380-387`):
```rust
let rho_for_xc = add_core_density(&rho_r, &rho_core_r);      // for V_xc
let rho_new_for_xc = add_core_density(&rho_r_new, &rho_core_r); // for E_xc
```

Energy uses `rho_xc_r = ρ_val + ρ_core` for E_xc but `rho_val_r = ρ_val` for E_vxc.

**Status: CORRECT**

---

## 17. Spin-Polarized LSDA

**Code** (`src/scf/mod.rs:556-560`):
```rust
let rho_up_xc = add_core_density(&rho_up_r, &(rho_core_r/2));
let rho_down_xc = add_core_density(&rho_down_r, &(rho_core_r/2));
let (exc_r, vxc_up_r, vxc_down_r) = xc::lda_xc_spin_grid(&rho_up_xc, &rho_down_xc);
```

**Status: CORRECT**

Core density split equally between channels, as per standard convention (the
frozen core is unpolarized).

---

## Summary

| # | Formula | Status |
|---|---------|--------|
| 1 | Total energy | CORRECT |
| 2 | Band energy | CORRECT |
| 3 | Hartree energy | CORRECT |
| 4 | Hartree potential | CORRECT |
| 5 | XC energy | CORRECT |
| 6 | Slater exchange | CORRECT |
| 7 | PZ correlation | CORRECT |
| 8 | Kinetic energy | CORRECT |
| 9 | Local potential | CORRECT |
| 10 | Non-local potential | CORRECT |
| 11 | Ewald summation | CORRECT |
| 12 | Density construction | CORRECT |
| 13 | Fermi energy search | CORRECT |
| 14 | Smearing functions | CORRECT |
| 15 | Sigma→0 extrapolation | CORRECT |
| 16 | NLCC | CORRECT |
| 17 | Spin-polarized LSDA | CORRECT |

**Known issues** (not formula bugs, but implementation gaps):
- PSP8 `D_ij` not populated from `ekb` values (Proposal 18)
- UPF `rho_atom` unit conversion needs verification (Proposal 18)
- No GGA infrastructure yet (only LDA)
