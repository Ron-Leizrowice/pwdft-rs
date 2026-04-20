# Non-Local Pseudopotential (Kleinman-Bylander)

Reference: Kleinman & Bylander, Phys. Rev. Lett. 48, 1425 (1982)

## The KB Separable Form

```text
V̂_NL = Σ_atom Σ_{i,j} |β_i⟩ D_{ij} ⟨β_j| × S_atom
```

where β_i are projector functions with angular momentum l_i, D_{ij} is the
coupling matrix (eV), and S_atom is the structure factor.

## Matrix Elements in PW Basis

```text
V_NL(G,G') = (1/Ω) Σ_atom S(G-G') × Σ_{i,j, l_i=l_j} F_i(|k+G|) D_{ij} F_j(|k+G'|)
             × (2l+1)/(4π) P_l(cos θ)
```

where:

- `F_i(q) = 4π ∫ [r·β_i(r)] j_l(qr) r dr` — form factor (Bessel transform)
- `P_l(cos θ)` — Legendre polynomial, `cos θ = q̂·q̂'`
- `S(G-G') = Σ_atom exp(-i(G-G')·τ)`

The sum over m is eliminated by the spherical harmonic addition theorem:
`Σ_m Y_lm(q̂) Y_lm*(q̂') = (2l+1)/(4π) P_l(cos θ)`

The phase factors `i^l` from bra and `(i*)^l` from ket give `|i|^{2l} = 1` (cancel).

**Code:** `src/potential/nonlocal.rs:105-202` (matrix element), lines 53-103 (precomputation)

## Form Factor (Bessel Transform)

```text
F_i(q) = 4π ∫₀^∞ [r·β_i(r)] × j_l(qr) × r × dr
```

UPF files store `χ(r) = r·β(r)` in Bohr^{-1/2}, converted to Å^{-1/2} at
parse time. The extra `× r` in the integrand comes from the 3D Fourier
transform of a function with angular structure `f(r) Y_lm(r̂)`:

```text
f̃(q) = 4π (-i)^l Y_lm(q̂) ∫₀^∞ r² f(r) j_l(qr) dr
```

With `f(r) = β(r)` and UPF storing `r·β(r)`, the integral becomes
`∫ r × [r·β(r)] × j_l(qr) dr`.

**Code:** `src/potential/nonlocal.rs:217-238`

**QE match:** Confirmed identical to `qe-7.5/upflib/beta_mod.f90:112-113`:

```fortran
aux(ir) = upf(nt)%beta(ir,nb) * besr(ir) * rgrid(nt)%r(ir)
```

## Spherical Bessel Functions

Explicit formulas for l=0,1, upward recurrence for l≥2:

```text
j_0(x) = sin(x)/x
j_1(x) = sin(x)/x² - cos(x)/x
j_{l+1}(x) = (2l+1)/x × j_l(x) - j_{l-1}(x)
```

Small-x limit: `j_0(0) = 1`, `j_l(0) = 0` for l>0.

Upward recurrence is stable for l < x. For DFT pseudopotentials, l ≤ 6 and
qr values are moderate, so this is safe.

**Code:** `src/potential/nonlocal.rs:248-268`

## Legendre Polynomials

Bonnet's recurrence: `(n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)`

Tested: P(1)=1, P(-1)=(-1)^l, orthogonality ∫P_l P_m = 2δ/(2l+1).

**Code:** `src/potential/nonlocal.rs:274-290`

## D_ij Matrix

For norm-conserving pseudopotentials, D_ij is stored in PP_DIJ block of UPF
(in Ry, converted to eV at parse time).

**Code:** `src/pseudopotential/upf/convert.rs` (`PP_DIJ` block — Ry→eV scalar multiply).

## Dimensional Analysis

```text
F_i: [Å^{-1/2}] × [Å] × [Å] = [Å^{3/2}]
D_ij: [eV]
F_i × D × F_j: [Å³] × [eV] = [eV·Å³]
(1/Ω): [Å^{-3}]
V_NL: [eV] ✓
```

## Audit Status

| Item | Status | Verified against |
|------|--------|-----------------|
| KB matrix element formula | CORRECT | Standard, QE init_us_2 |
| Form factor integrand (r·χ·j_l·r) | CORRECT | QE beta_mod.f90:112-113 |
| Phase i^l cancellation | CORRECT |  |
| Angular factor (2l+1)/(4π) P_l | CORRECT | Addition theorem |
| Structure factor S(G-G') | CORRECT | exp(-iG·τ)×conj(exp(-iG'·τ)) |
| Bessel function j_l(x) | CORRECT | Tested to 1e-10 |
| Legendre polynomial P_l(x) | CORRECT | Orthogonality verified |
| Unit conversions (β, D_ij) | CORRECT | UPF Ry→eV, Bohr→Å |
| **Form factor quadrature** | **ISSUE** | See [Radial Integration](radial-integration.md) |
