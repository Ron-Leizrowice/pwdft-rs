# Plane-Wave DFT: Theory and Implementation Reference

This document collects the core formulas underlying pwdft-rs, verified against
the canonical references below. Each section states the formula, the convention
we use, and any pitfalls. It is intended as a living reference to prevent
re-deriving solved problems.

## Key References

1. **Payne, Teter, Allan, Arias, Joannopoulos** — *Iterative minimization
   techniques for ab initio total-energy calculations: molecular dynamics and
   conjugate gradients*, Rev. Mod. Phys. **64**, 1045 (1992).
   The canonical review of plane-wave pseudopotential DFT.
   DOI: [10.1103/RevModPhys.64.1045](https://doi.org/10.1103/RevModPhys.64.1045)

2. **Martin, R. M.** — *Electronic Structure: Basic Theory and Practical
   Methods*, 2nd ed., Cambridge University Press (2020).
   Textbook covering all aspects from Hohenberg-Kohn to practical PW-PP
   implementation. Chapters 11-13 are especially relevant.

3. **Kleinman & Bylander** — *Efficacious Form for Model Pseudopotentials*,
   Phys. Rev. Lett. **48**, 1425 (1982).
   DOI: [10.1103/PhysRevLett.48.1425](https://doi.org/10.1103/PhysRevLett.48.1425)

4. **Perdew & Zunger** — *Self-interaction correction to density-functional
   approximations for many-electron systems*, Phys. Rev. B **23**, 5048 (1981).
   Table I gives the LDA correlation parameters used in this code.
   DOI: [10.1103/PhysRevB.23.5048](https://doi.org/10.1103/PhysRevB.23.5048)

5. **Louie, Froyen, Cohen** — *Nonlinear ionic pseudopotentials in
   spin-density-functional calculations*, Phys. Rev. B **26**, 1738 (1982).
   The NLCC paper: core charge must be added for XC when core/valence overlap.
   DOI: [10.1103/PhysRevB.26.1738](https://doi.org/10.1103/PhysRevB.26.1738)

6. **ABINIT theory documentation** — Pseudopotentials:
   https://docs.abinit.org/theory/pseudopotentials/

7. **Theoretical Physics Reference** — DFT chapter:
   https://www.theoretical-physics.com/dev/quantum/dft.html

8. **MAVENs DFT Notes** — SCF convergence:
   https://mavens-group.github.io/dft-notes/08-SCF.html

---

## 1. Kohn-Sham Total Energy

The central quantity. The Kohn-Sham total energy functional is:

```
E[ρ] = T_s[ρ] + E_H[ρ] + E_xc[ρ] + E_ext[ρ] + E_ion-ion
```

where:
- `T_s[ρ]` = kinetic energy of non-interacting electrons
- `E_H[ρ]` = Hartree (classical Coulomb) electron-electron energy
- `E_xc[ρ]` = exchange-correlation energy
- `E_ext[ρ]` = electron-ion interaction (local + non-local pseudopotential)
- `E_ion-ion` = ion-ion electrostatic energy (Ewald)

### 1.1 Band Energy and Double-Counting Correction

The Kohn-Sham eigenvalues already contain T_s, E_ext, V_H, and V_xc through
the effective potential. Using the identity (Ref. 7):

```
T_s = Σ_{n,k} f_{n,k} w_k ε_{n,k} - ∫ V_eff(r) ρ(r) dr
```

where `V_eff = V_ext + V_H + V_xc`, we can express total energy as:

```
E_total = E_band - E_H + E_xc - E_vxc + E_ion-ion
```

where:
- `E_band = Σ_{n,k} f_{n,k} w_k ε_{n,k}` (sum of occupied eigenvalues)
- `E_H = (1/2) ∫ V_H(r) ρ(r) dr` (Hartree double-counting subtraction)
- `E_xc = ∫ ε_xc(r) ρ(r) dr` (XC energy density integrated)
- `E_vxc = ∫ V_xc(r) ρ(r) dr` (XC potential double-counting subtraction)
- `E_ion-ion` = Ewald energy

**Why this works:** E_band counts E_H once (through V_H in V_eff) and E_vxc
once (through V_xc). We subtract E_H (which over-counts by 1/2) and replace
E_vxc with E_xc (the true XC energy, not the potential integral).

**In code:** `src/scf/mod.rs`, function `compute_total_energy_from_components()`:

```rust
e_band - e_hartree + e_xc - e_vxc + e_ewald
```

Plus the `v_local_g0 * n_electrons` constant shift (see Section 3.1).

### 1.2 NLCC Modification

With nonlinear core correction (Ref. 5), the XC functional is evaluated on
the augmented density `ρ_total = ρ_valence + ρ_core`:

```
E_xc = ∫ ε_xc[ρ_val + ρ_core] × (ρ_val + ρ_core) dr    (XC energy)
E_vxc = ∫ V_xc[ρ_val + ρ_core] × ρ_val dr                (double-counting)
```

The key subtlety: `E_xc` uses the total density in both the functional and the
integration measure, but `E_vxc` uses only the valence density in the
integration measure. This is because the band eigenvalues contain `V_xc × ρ_val`
(only valence electrons feel V_eff), so the double-counting correction must
also use `ρ_val`.

**In code:** `src/scf/mod.rs`, lines 380-387:
```rust
let rho_new_for_xc = add_core_density(&rho_r_new, &rho_core_r);
let e_xc = xc::lda_xc_energy(rho_xc_r, exc_r, omega);       // ε_xc × ρ_total
let e_vxc = Σ V_xc × ρ_val × dvol;                           // V_xc × ρ_val
```

---

## 2. Reciprocal-Space Formulation

### 2.1 Bloch Waves and Plane-Wave Expansion

In a periodic crystal, the wavefunction satisfies Bloch's theorem:

```
ψ_{n,k}(r) = e^{ik·r} u_{n,k}(r)
```

where `u_{n,k}(r)` has the periodicity of the lattice. Expanding in plane waves:

```
ψ_{n,k}(r) = (1/√Ω) Σ_G c_{n,k}(G) e^{i(k+G)·r}
```

where `G` are reciprocal lattice vectors and `Ω` is the cell volume.

### 2.2 Kohn-Sham Matrix Equation

Substituting into the KS equation gives a matrix eigenvalue problem:

```
Σ_{G'} H_{G,G'}(k) c_{n,k}(G') = ε_{n,k} c_{n,k}(G)
```

with matrix elements:

```
H_{G,G'}(k) = T_{G,G'}(k) + V_eff(G-G')

T_{G,G'}(k) = (ħ²/2m) |k+G|² δ_{G,G'}     (kinetic, diagonal)
V_eff(G-G') = V_local(G-G') + V_H(G-G') + V_xc(G-G')    (local potentials)
```

Plus the non-local pseudopotential contribution (Section 4).

**In code:** `build_hamiltonian_with_v_eff()` constructs this matrix.
- Kinetic: `HBAR2_OVER_2M * |k+G|²` on the diagonal
- Potential: `v_eff_fft[miller_to_idx(G-G')]` for each (G, G') pair

### 2.3 FFT Convention

We use the convention:

```
Forward:   f̃(G) = Σ_r f(r) e^{-iG·r}           (unnormalized)
Inverse:   f(r) = Σ_G f̃(G) e^{+iG·r}           (unnormalized)
Normalized: f̃(G) = (1/N) × Forward[f(r)]        (divide by N_grid)
```

The 1/N normalization on the forward FFT means:
- `ρ̃(G=0) = (1/N) Σ_r ρ(r) = ⟨ρ⟩` (spatial average of density)
- `V(G-G') = (1/N) Σ_r V(r) e^{-i(G-G')·r}` (Fourier coefficient)

**Why:** With this convention, `V_eff(G-G')` in the Hamiltonian has the correct
units (eV) without extra volume factors. The density `ρ̃(G=0)` equals the
average density `n_electrons/Ω`.

---

## 3. Potentials

### 3.1 Local Pseudopotential

```
V_local(G) = (1/Ω) Σ_atom V^pp(|G|) × S_atom(G)
```

where `V^pp(|G|)` is the Bessel transform of the radial local potential and
`S_atom(G) = exp(-iG·τ_atom)` is the structure factor.

**G=0 convention:** `V_local(G=0)` is an arbitrary constant that depends on
the pseudopotential construction. QE excludes it from the Hamiltonian and adds
`V_local(G=0) × N_electrons` to the total energy. We follow the same convention.

**In code:** `src/scf/mod.rs`, lines 214-216:
```rust
let v_local_g0 = v_local_fft[0].re;
v_local_fft[0] = Complex64::new(0.0, 0.0);  // zero in Hamiltonian
// ... later: e_total += v_local_g0 * n_electrons
```

### 3.2 Hartree Potential

In reciprocal space, Poisson's equation gives:

```
V_H(G) = 4πe² ρ̃(G) / |G|²    for G ≠ 0
V_H(G=0) = 0                   (absorbed into G=0 convention)
```

The Hartree energy:

```
E_H = (Ω/2) Σ_{G≠0} |ρ̃(G)|² × 4πe² / |G|²
```

The factor of Ω appears because our `ρ̃(G)` is the spatial average (divided
by N_grid during FFT normalization), and `|ρ̃(G)|²` needs to be integrated
over the cell volume.

**In code:** `hartree_on_fft_grid()` computes `V_H(G) = 4πe² ρ̃(G)/|G|²`.
Energy uses `0.5 * omega * Σ |ρ̃(G)|² × 4πe²/|G|²`.

### 3.3 Exchange-Correlation (LDA)

Slater exchange (Ref. 4, but following standard LDA):

```
ε_x(ρ) = -(3/4)(3ρ/π)^{1/3}    [Hartree]
V_x(ρ) = (4/3) ε_x(ρ)          [Hartree]
```

Note: `ε_x` is energy per electron, so `E_x = ∫ ρ(r) ε_x(ρ(r)) dr`.

Perdew-Zunger correlation (Ref. 4, Table I, unpolarized):

```
r_s ≥ 1:  ε_c = γ / (1 + β₁√r_s + β₂ r_s)
r_s < 1:  ε_c = A ln(r_s) + B + C r_s ln(r_s) + D r_s
```

with `γ = -0.1423, β₁ = 1.0529, β₂ = 0.3334, A = 0.0311, B = -0.048,
C = 0.002, D = -0.0116` (all in Hartree).

The XC potential: `V_xc = d(ρ ε_xc)/dρ = ε_xc - (r_s/3) dε_xc/dr_s`.

**Unit conversion:** Code works in eV and Å. Density `ρ` in e/ų is converted
to e/Bohr³ via `ρ_Bohr = ρ_Å × Bohr³/ų = ρ_Å × 0.529177³`. Results converted
from Hartree to eV via `HA_TO_EV = 27.211386`.

### 3.4 Nonlinear Core Correction (NLCC)

When the core and valence charge distributions overlap significantly (transition
metals, alkali atoms), the XC potential computed from `ρ_valence` alone is
inaccurate because XC is nonlinear in ρ (Ref. 5).

The fix: evaluate XC on the augmented density `ρ_val + ρ_core`, where `ρ_core`
is the frozen core charge density from the pseudopotential file (PP_NLCC in UPF
format, or the core charge block in PSP8).

```
V_xc[ρ_val + ρ_core]  replaces  V_xc[ρ_val]
```

This applies to both the potential (in the Hamiltonian) and the energy.
Without NLCC, transition metal eigenvalues are shifted by a constant
(~15 eV for Fe), producing totally wrong energies and band structures.

**Common pitfall:** The core density is split equally between spin channels
in spin-polarized calculations: `ρ_core_up = ρ_core_down = ρ_core/2`.

---

## 4. Non-Local Pseudopotential (Kleinman-Bylander)

### 4.1 The KB Separable Form

The full non-local pseudopotential operator (Ref. 3):

```
V̂_NL = Σ_atom Σ_{i,j} |β_i⟩ D_{ij} ⟨β_j| × S_atom
```

where `β_i` are projector functions with angular momentum `l_i`, `D_{ij}` is the
coupling matrix, and `S_atom` is the structure factor for atomic positions.

### 4.2 Matrix Elements in PW Basis

The matrix element between plane waves `|k+G⟩` and `|k+G'⟩` (Ref. 6):

```
V_NL(G,G') = (1/Ω) Σ_type Σ_{i,j with l_i=l_j} F_i(|k+G|) D_{ij} F_j(|k+G'|)
             × (2l+1)/(4π) P_l(cos θ) × (-1)^l
             × Σ_atom S(G_i - G_j)
```

where:
- `F_i(q) = 4π ∫ [r β_i(r)] j_l(qr) r dr` is the form factor (Bessel transform)
- `P_l(cos θ)` is the Legendre polynomial, with `cos θ = q̂·q̂'`
- The sum over `m` is eliminated by the spherical harmonic addition theorem:
  `Σ_m Y_lm(q̂) Y_lm*(q̂') = (2l+1)/(4π) P_l(cos θ)`
- The `(-1)^l` factor comes from `i^l × (i*)^l`
- `S(G-G') = Σ_atom exp(-i(G-G')·τ) = Σ_atom exp(-iG·τ) exp(iG'·τ)`

**Important:** UPF files store `r × β(r)`, so the integral is
`∫ [rβ(r)] j_l(qr) r dr` (note the extra `r`).

**In code:** `NonlocalPotential::new()` precomputes form factors;
`add_to_hamiltonian()` accumulates the matrix elements with structure factors.

### 4.3 D_ij Matrix

For norm-conserving pseudopotentials:
- **UPF format:** D_ij is stored in the PP_DIJ block (in Ry, converted to eV)
- **PSP8 format:** D_ij is diagonal with `D_ii = ekb_i` (the KB energies from
  each projector block header). **Note:** Off-diagonal D_ij requires ONCVPSP
  formalism with multiple projectors per angular momentum channel.

---

## 5. Ewald Summation (Ion-Ion Energy)

The Coulomb energy between periodic point charges is decomposed (standard Ewald):

```
E_ewald = E_real + E_recip + E_self + E_background
```

### 5.1 Reciprocal Space

```
E_recip = (2πe²/Ω) Σ_{G≠0} |S(G)|² exp(-|G|²/(4η²)) / |G|²
```

where `S(G) = Σ_i Z_i exp(iG·r_i)` is the charge-weighted structure factor
and `η` is the Ewald screening parameter.

### 5.2 Real Space

```
E_real = (e²/2) Σ_{T} Σ'_{i,j} Z_i Z_j erfc(η|r_ij + T|) / |r_ij + T|
```

where the prime means exclude `i=j` when `T=0` (self-interaction).

### 5.3 Self-Energy Correction

```
E_self = -(η/√π) e² Σ_i Z_i²
```

### 5.4 Background Charge

```
E_bg = -π e² (Σ_i Z_i)² / (2Ω η²)
```

This neutralizes the divergent G=0 term. For charge-neutral cells
(`Σ Z_i = N_electrons`), this contribution is small but non-zero.

### 5.5 Screening Parameter

```
η = (N_atoms × π / Ω)^{1/3}
```

This balances computational cost between real and reciprocal sums.

---

## 6. Electron Density

### 6.1 From Wavefunctions

```
ρ(r) = Σ_{n,k} f_{n,k} w_k |ψ_{n,k}(r)|²
```

where `f_{n,k}` is the occupation (including spin factor: 0 to 2 for nspin=1,
0 to 1 for nspin=2) and `w_k` is the k-point weight (sum to 1).

In practice:
1. Place PW coefficients `c_{n,k}(G)` onto FFT grid at positions `G`
2. Inverse FFT → `ψ_{n,k}(r)` (unnormalized)
3. Accumulate `ρ(r) += f_{n,k} × w_k × |ψ(r)|²`
4. Normalize so `∫ ρ(r) dr = (Ω/N_grid) Σ_r ρ(r) = N_electrons`

### 6.2 Initial Density (SAD)

Superposition of Atomic Densities: sum isolated-atom charge densities.
For each atom with pseudopotential data `4πr²ρ_atom(r)`:

```
ρ_atom(G) = (1/Ω) ∫ [4πr²ρ(r)] j₀(|G|r) dr × S(G)
```

Note: UPF `PP_RHOATOM` stores `4πr²ρ(r)`, so no extra `4πr²` factor needed.

---

## 7. SCF Convergence

### 7.1 The Fixed-Point Problem

The SCF cycle is: `ρ_out = F[ρ_in]` where F maps an input density through
one complete KS cycle. Convergence requires finding `ρ* = F[ρ*]`.

### 7.2 Linear Mixing

```
ρ_in^{n+1} = ρ_in^n + α R^n
```

where `R^n = ρ_out^n - ρ_in^n` and `α` is the mixing parameter (`mixing_beta`).

### 7.3 Anderson/Pulay (DIIS) Mixing

Construct optimal linear combination of history:

```
ρ̄_in = Σ_i c_i ρ_in^{n-m+i},   Σ c_i = 1
```

Minimize `|R̄|² = Σ_{ij} c_i c_j ⟨R^i|R^j⟩` subject to constraint.
Solve the DIIS matrix equation:

```
| B₁₁  B₁₂  ...  1 | | c₁ |   | 0 |
| B₂₁  B₂₂  ...  1 | | c₂ | = | 0 |
| ...              1 | | .. |   | . |
|  1    1    ...  0 | | -μ |   | 1 |
```

where `B_{ij} = ⟨R^i|R^j⟩`.

Update: `ρ_in^{n+1} = Σ_i c_i [ρ_in^i + α R^i]`

### 7.4 Kerker Preconditioning

Damp long-wavelength charge sloshing (Ref. 8):

```
R̃(G) = [|G|² / (|G|² + q_TF²)] R(G)
```

where `q_TF = √(4πe² N(ε_F))` is the Thomas-Fermi wavevector. This suppresses
the G→0 component of the residual where the dielectric response is largest.

---

## 8. Smearing and Entropy

### 8.1 Mermin Free Energy

At finite electronic temperature, the relevant functional is the Mermin free
energy `F = E - TS`, where S is the electronic entropy from fractional
occupations. The sigma→0 extrapolated energy:

```
E₀ = (E + F) / 2 = E - TS/2
```

This is the best estimate of the T=0 total energy from a finite-σ calculation.

### 8.2 Occupation Functions

All give occupation in [0, 1] for a single state. Multiply by `spin_factor`
(2/nspin) for the full occupation.

**Fermi-Dirac:** `f(x) = 1/(1 + e^x)` where `x = (ε - E_F)/σ`

**Gaussian:** `f(x) = erfc(x)/2`

**Methfessel-Paxton (order 1):** `f(x) = erfc(x)/2 - (1/2) x exp(-x²)/√π`

**Marzari-Vanderbilt (cold):** `f(x) = (1/2) erfc(x + 1/√2) + exp(-(x+1/√2)²)/√(2π)`

### 8.3 Entropy Weights

The entropy contribution per state `S_n = σ × s(x)` where:

- **Fermi-Dirac:** `s = -(f ln f + (1-f) ln(1-f))`
- **Gaussian:** `s = exp(-x²)/√π`
- **Methfessel-Paxton:** `s = (1/2 - x²) exp(-x²)/√π`
- **Cold:** `s = (x+1/√2) exp(-(x+1/√2)²)/√π`

---

## 9. Units Convention

This code works in **eV and Ångströms** throughout (not Hartree/Bohr or Rydberg):

| Quantity | Unit |
|----------|------|
| Energy | eV |
| Length | Å |
| Density | e/ų |
| Potential | eV |
| Wavevector | 1/Å |
| ħ²/2m | eV·Å² (= 3.81 eV·Å²) |
| e² | 14.4 eV·Å (Coulomb constant) |

Conversions:
- `1 Ry = 13.606 eV`
- `1 Ha = 27.211 eV`
- `1 Bohr = 0.529177 Å`
- `1 Bohr³ = 0.148185 ų`

Pseudopotential files (UPF: Ry/Bohr, PSP8: Ha/Bohr) are converted at parse time.

---

## 10. Common Pitfalls (Lessons Learned)

### 10.1 NLCC Omission (Fe bug)
**Symptom:** Constant eigenvalue shift (~15 eV for Fe) relative to QE.
**Cause:** `PP_NLCC` data ignored; XC evaluated on `ρ_val` instead of
`ρ_val + ρ_core`. Si works because it has no NLCC.
**Fix:** Parse `PP_NLCC`, Bessel-transform to FFT grid, add to density before
XC evaluation. See Proposal 25.

### 10.2 V_local(G=0) Convention
**Symptom:** All eigenvalues shifted by a large constant.
**Cause:** `V_local(G=0)` is pseudopotential-dependent and arbitrary.
**Fix:** Exclude from Hamiltonian, add `V_local(G=0) × N_electrons` to total
energy. This matches QE's convention.

### 10.3 FFT Normalization
**Symptom:** Energies off by factors of N_grid or Ω.
**Cause:** Inconsistent 1/N vs 1/Ω normalization between density, potential,
and energy formulas.
**Convention:** Our forward FFT includes 1/N normalization, so `ρ̃(G=0)` is the
spatial average density. Hartree energy includes explicit `× Ω` factor.

### 10.4 PSP8 D_ij
**Symptom:** Zero non-local energy with PSP8 pseudopotentials.
**Cause:** D_ij matrix not populated from `ekb` values in PSP8 format.
**Fix:** See Proposal 18.

### 10.5 Spin Factor in Occupations
**Symptom:** Wrong electron count or doubled energies in spin-polarized.
**Cause:** Occupation function returns `f ∈ [0, spin_factor]` where
`spin_factor = 2/nspin`. Must be consistent in density, energy, and Fermi
energy search.
**Convention:** `occupation()` includes spin_factor. The Fermi energy search
uses the same function. Density accumulates `f × w_k × |ψ|²`.

### 10.6 E_xc vs E_vxc with NLCC
**Symptom:** Small but systematic energy error with NLCC pseudopotentials.
**Cause:** Using `ρ_total` in both `E_xc` and `E_vxc`.
**Fix:** `E_xc = ∫ ε_xc × ρ_total dr` but `E_vxc = ∫ V_xc × ρ_val dr`.
The double-counting correction uses valence-only because that's what the
eigenvalues contain.
