# Proposal 26: Math Audit, Reference Documentation, and Comment Fix

## Problem

The pwdft-rs codebase implements many core DFT formulas — total energy, Hartree, XC, non-local pseudopotential, Ewald, smearing, NLCC — but there is no single document verifying these against canonical references, and no reference to consult when debugging discrepancies (like the Fe NLCC bug in Proposal 25). Stumbling over pre-existing solved problems wastes time.

Additionally, one comment in the codebase is mathematically incorrect and could mislead future contributors.

## What Was Done

A full audit of every formula in the codebase was performed against canonical DFT theory references:

1. **Payne, Teter, Allan, Arias, Joannopoulos** — Rev. Mod. Phys. **64**, 1045 (1992)
2. **Martin, R. M.** — *Electronic Structure*, 2nd ed., Cambridge (2020)
3. **Kleinman & Bylander** — Phys. Rev. Lett. **48**, 1425 (1982)
4. **Perdew & Zunger** — Phys. Rev. B **23**, 5048 (1981)
5. **Louie, Froyen, Cohen** — Phys. Rev. B **26**, 1738 (1982) [NLCC]
6. **ABINIT documentation** — Pseudopotential theory
7. **MAVENs DFT Notes** — SCF convergence algorithms

### Audit Results

All 17 formulas verified **correct**:

| # | Formula | File | Status |
|---|---------|------|--------|
| 1 | Total energy (band - Hartree + XC - VXC + Ewald) | `scf/mod.rs:1020` | CORRECT |
| 2 | Band energy Σ f w ε | `scf/mod.rs:989` | CORRECT |
| 3 | Hartree energy (Ω/2) Σ \|ρ(G)\|² 4πe²/G² | `scf/mod.rs:998` | CORRECT |
| 4 | Hartree potential V_H(G) = 4πe² ρ(G)/G² | `scf/mod.rs:786` | CORRECT |
| 5 | XC energy ∫ ε_xc ρ dr | `potential/xc.rs:64` | CORRECT |
| 6 | Slater exchange -0.75(3ρ/π)^{1/3} | `potential/xc.rs:91` | CORRECT |
| 7 | PZ correlation (all 7 params match Table I) | `potential/xc.rs:117-135` | CORRECT |
| 8 | Kinetic energy (ħ²/2m)\|k+G\|² | `scf/mod.rs:847` | CORRECT |
| 9 | Local pseudopotential on FFT grid | `scf/mod.rs:766` | CORRECT |
| 10 | Non-local KB matrix elements | `potential/nonlocal.rs:154-198` | CORRECT |
| 11 | Ewald (reciprocal + real + self + background) | `ewald.rs:51-114` | CORRECT |
| 12 | Density ρ(r) = Σ f w \|ψ\|² | `scf/density.rs:50-62` | CORRECT |
| 13 | Fermi energy bisection | `scf/smearing.rs:50-91` | CORRECT |
| 14 | Smearing (FD, Gaussian, MP, Cold) | `scf/smearing.rs:107-149` | CORRECT |
| 15 | Sigma→0 extrapolation E₀ = (E+F)/2 | `scf/mod.rs:414` | CORRECT |
| 16 | NLCC (ρ_xc for E_xc, ρ_val for E_vxc) | `scf/mod.rs:380-387` | CORRECT |
| 17 | Spin-polarized LSDA with VBH interpolation | `potential/xc.rs:212-301` | CORRECT |

### Comment Bug Found

**`src/potential/nonlocal.rs`, lines 50-52:**

```rust
///     F_i(|k+G|) D_{ij} F_j(|k+G'|) × (2l+1)/(4π) P_l(cos θ) × (-1)^l
///
/// (The i^l × (i*)^l = (-1)^l factor)
```

The comment claims `i^l × (i*)^l = (-1)^l`. This is wrong.

**Correct identity:** `i^l × (i*)^l = |i|^{2l} = 1` for all `l`.

Proof by cases:
- l=0: 1 × 1 = 1
- l=1: (-i) × (i) = -i² = 1
- l=2: (-1) × (-1) = 1
- l=3: (i) × (-i) = -i² = 1

The full derivation: the projector FT is `⟨k+G|β_lm⟩ = (-i)^l Y_lm(q̂) F_l(q)` (Rayleigh expansion), and its conjugate is `⟨β_lm|k+G'⟩ = (i)^l Y_lm*(q'̂) F_l(q')`. The product `(-i)^l (i)^l = 1`, so **no `(-1)^l` factor appears** in the matrix element.

The **code is correct** — it does NOT apply the `(-1)^l` factor. If it did, the l=1 (p-orbital) channel would have the wrong sign, and Si (which has l=0 and l=1 projectors) would give incorrect band structures. Si matches QE, confirming the code is right and the comment is wrong.

## Implementation

### Step 1: Create `docs/theory.md`

A comprehensive reference document (~18 KB) covering:

1. **Total energy functional** — KS decomposition, double-counting correction, NLCC modification
2. **Reciprocal-space formulation** — Bloch waves, PW expansion, KS matrix equation, FFT convention
3. **Potentials** — Local PP (with G=0 convention), Hartree in G-space, LDA XC (Slater + PZ with all parameters), NLCC
4. **Non-local pseudopotential** — KB separable form, PW matrix elements, form factor Bessel transform, angular sum via addition theorem, D_ij for UPF and PSP8
5. **Ewald summation** — All four terms with formulas, screening parameter choice
6. **Electron density** — From wavefunctions, SAD initial density
7. **SCF convergence** — Linear mixing, Anderson/Pulay DIIS with matrix equation, Kerker preconditioning
8. **Smearing and entropy** — All four schemes (FD, Gaussian, MP, Cold) with occupation and entropy formulas, Mermin free energy, sigma→0 extrapolation
9. **Units convention** — eV/Å table with all conversion factors
10. **Common pitfalls** — NLCC omission, V_local(G=0), FFT normalization, PSP8 D_ij, spin factor

Each section cites specific references (Payne 1992, Martin textbook, PZ 1981, etc.) and links to the relevant code locations.

### Step 2: Create `docs/math-audit.md`

A line-by-line verification (~11 KB) of every formula against the code:
- Each entry shows the canonical formula, the exact code with line numbers, and a CORRECT/BUG status
- Summary table of all 17 verified formulas
- Notes on known implementation gaps (PSP8 D_ij, rho_atom units) with proposal references

### Step 3: Fix the comment in `src/potential/nonlocal.rs`

Lines 49-52, change from:

```rust
///     F_i(|k+G|) D_{ij} F_j(|k+G'|) × (2l+1)/(4π) P_l(cos θ) × (-1)^l
///
/// (The i^l × (i*)^l = (-1)^l factor)
```

To:

```rust
///     F_i(|k+G|) D_{ij} F_j(|k+G'|) × (2l+1)/(4π) P_l(cos θ)
///
/// (The phase factors i^l from bra and (i*)^l from ket give i^l(i*)^l = |i|^{2l} = 1)
```

## Files Modified

- `docs/theory.md` — **new** (reference document)
- `docs/math-audit.md` — **new** (verification document)
- `src/potential/nonlocal.rs` — **comment fix** (lines 50-52)

## Acceptance Criteria

1. **`docs/theory.md` exists** with formulas for all 10 sections, citing canonical references.
2. **`docs/math-audit.md` exists** with line-by-line verification of all 17 formulas.
3. **Comment fixed** in `nonlocal.rs` — no longer claims `(-1)^l`.
4. **No code changes** beyond the comment fix — this is documentation only.
5. **`cargo test` passes** unchanged (comment fix has no effect on compilation or behavior).
6. **Future debugging aided** — the pitfalls section covers the Fe NLCC bug, V_local(G=0), FFT normalization, PSP8 D_ij, and spin factor conventions, so these are not re-derived.
