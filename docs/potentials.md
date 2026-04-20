# Potentials: Local, Hartree, Exchange-Correlation

## Local Pseudopotential

```text
V_local(G) = Σ_atom S(G) × v_form(|G|)
```

where `S(G) = exp(-iG·τ)` is the structure factor and `v_form` is the
radial form factor from the pseudopotential.

### Form factor (Bessel transform)

For G ≠ 0:

```text
v_form(G) = (4π/Ω) ∫₀^∞ r² [V_loc(r) + Ze²/r] sin(Gr)/(Gr) dr − 4πZe²/(ΩG²)
```

For G = 0:

```text
v_form(0) = (4π/Ω) ∫₀^∞ r² [V_loc(r) + Ze²/r] dr
```

The Coulomb singularity is handled by adding `Ze²/r` to make the integrand
short-ranged, then subtracting the Fourier transform of `Ze²/r` analytically.

**Code:** `src/pseudopotential/mod.rs:98-133`

**Known issue:** Our integrand `V_loc(r) + Ze²/r` diverges at r=0 because the
pseudopotential is smooth at the origin (V_loc(0) is finite). QE uses
`erf(r)/r` subtraction instead, which stays bounded. See
[Radial Integration](radial-integration.md) and Proposals 38-39.

### G=0 convention

V_local(G=0) is excluded from the Hamiltonian and added to the total energy
as `V_local(G=0) × N_electrons`. This matches QE's convention
(`setlocal.f90:92`).

**Code:** `src/scf/context.rs:80-83`

### Assembly on FFT grid

**Code:** `src/scf/potentials.rs:19-52`

---

## Hartree Potential

In reciprocal space (Poisson's equation):

```text
V_H(G) = 4πe² ρ̃(G) / |G|²     for G ≠ 0
V_H(G=0) = 0                    (neutralizing background)
```

Hartree energy:

```text
E_H = (Ω/2) Σ_{G≠0} |ρ̃(G)|² × 4πe² / |G|²
```

The factor Ω appears because ρ̃(G) is the 1/N-normalized FFT (spatial
average), and |ρ̃(G)|² must be integrated over the cell volume.

**Code:** `src/potential/hartree.rs` (potential), `src/scf/energy.rs:34-49` (energy)

**Constant:** `E2 = 14.3996 eV·Å` (Coulomb constant, `e²/(4πε₀)`)

---

## Exchange-Correlation (LDA)

Computed in real space on the FFT grid. Unit conversion chain:

1. Input density ρ in e/Å³
2. Convert to e/Bohr³: `ρ_Bohr = ρ × BOHR3_TO_ANG3` (≈ 0.148)
3. Compute ε_xc, V_xc in Hartree
4. Convert to eV: `× HA_TO_EV` (27.2114)

### Slater exchange

```text
ε_x = -(3/4)(3ρ/π)^{1/3}     [Hartree, per electron]
V_x = (4/3) ε_x
```

Reference: Slater, Phys. Rev. 81, 385 (1951)

**Code:** `src/potential/xc.rs:80-96`

### Perdew-Zunger correlation

Two regimes based on Wigner-Seitz radius r_s = (3/(4πρ))^{1/3} in Bohr:

```text
r_s ≥ 1:  ε_c = γ / (1 + β₁√r_s + β₂ r_s)
r_s < 1:  ε_c = A ln(r_s) + B + C r_s ln(r_s) + D r_s
```

Parameters (unpolarized, Hartree):

| γ | β₁ | β₂ | A | B | C | D |
|---|----|----|---|---|---|---|
| -0.1423 | 1.0529 | 0.3334 | 0.0311 | -0.048 | 0.002 | -0.0116 |

Potential: `V_c = ε_c - (r_s/3) dε_c/dr_s`

Reference: Perdew & Zunger, Phys. Rev. B 23, 5048 (1981), Table I

**Code:** `src/potential/xc.rs:98-147`

### Spin-polarized LSDA

Exchange: `ε_x(ρ↑,ρ↓) = (ρ↑ε_x(2ρ↑) + ρ↓ε_x(2ρ↓)) / ρ`

Correlation: interpolated with von Barth-Hedin function:

```text
f(ζ) = [(1+ζ)^{4/3} + (1-ζ)^{4/3} - 2] / [2^{4/3} - 2]
ε_c(r_s, ζ) = ε_c^unpol + f(ζ) [ε_c^pol - ε_c^unpol]
```

Fully polarized PZ parameters: γ=-0.0843, β₁=1.3981, β₂=0.2611

**Code:** `src/potential/xc.rs:175-345`

### Audit status

| Item | Status | Verified against |
|------|--------|-----------------|
| V_local formula | CORRECT | QE `vloc_mod.f90` |
| V_local(G=0) convention | CORRECT | QE `setlocal.f90` |
| Hartree potential | CORRECT | Standard Poisson |
| Hartree energy | CORRECT | Positive-definite |
| Slater exchange | CORRECT | Slater 1951 |
| PZ correlation (all params) | CORRECT | PZ 1981 Table I |
| PZ derivative and V_c | CORRECT | Standard derivation |
| LSDA exchange | CORRECT | Passes ζ=0 limit |
| LSDA correlation | CORRECT | Von Barth-Hedin f(ζ) |
| Unit conversion chain | CORRECT | eV/Å throughout |
| **V_local quadrature** | **ISSUE** | See [Radial Integration](radial-integration.md) |
