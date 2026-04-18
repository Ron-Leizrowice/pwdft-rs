---
id: MAUD
title: Mathematical accuracy audit of core physics modules (post-MADOC-A cold read)
priority: medium
complexity: small
risk: low
depends_on: [MADOC]
blocks: []
status: active
author: Researcher
date: 2026-04-18
---

# MAUD — Mathematical accuracy audit of core physics modules

## §1 Scope and methodology

### §1.1 Motivation

MADOC Phase A (PR #87, commit `045b68d`) replaced the terse one-line
docstrings in `src/scf/energy.rs`, `src/scf/driver.rs`, and
`src/scf/driver_spin.rs` with full mathematical specifications:
defining equation, variable definitions with units, and primary-
literature citations. Those docstrings are now the *authoritative*
spec — the engineers' rewrite target if implementation and doc ever
diverge. Before GGAP Phase B starts extending the XC surface into GGA
territory, this audit does a cold line-by-line comparison of every
post-MADOC docstring against the code that purports to implement it,
plus spot-checks against QE 7.5 source and the primary-literature
citations MADOC named.

This is a **scoping + reporting pass**. No code is fixed here.

### §1.2 Modules in scope

| Module | File | What was checked |
|--------|------|------------------|
| Energy assembly | `src/scf/energy.rs` (post-MADOC) | Every public helper's docstring vs. code; `EnergyComponents` identity |
| Non-spin driver | `src/scf/driver.rs` (post-MADOC) | Module header KS equations, per-iteration pipeline |
| Spin driver | `src/scf/driver_spin.rs` (post-MADOC) | LSDA + CCMX derivation |
| LDA / LSDA XC | `src/potential/xc.rs` | Perdew-Zunger formulas; GGAP Phase A `XcEvaluator` |
| Non-local PP (KB) | `src/potential/nonlocal.rs` | `V_NL = Σ\|β⟩D⟨β\|` separable assembly; VNMT test |
| Local PP | `src/potential/local.rs`, `src/pseudopotential/mod.rs::v_local_of_g` | `V_local(G)` spherical Bessel transform |
| Ewald sum | `src/ewald.rs` | Real + recip + self + background; η choice |
| UPF boundary | `src/pseudopotential/upf/convert.rs` | Every Ry/Bohr → eV/Å conversion |
| PCFX G-space sym | `src/symmetry/density/g_space.rs` | Phase convention, rotation transpose |

### §1.3 Reference ground truth

- **Docstring (post-MADOC):** treated as authoritative. If code and
  docstring disagree, one of them is wrong — this proposal flags which.
- **Primary literature:** Perdew-Zunger 1981 PRB 23 5048; Kleinman-
  Bylander 1982 PRL 48 1425; Louie-Froyen-Cohen 1982 PRB 26 1738;
  Harris 1985 PRB 31 1770 + Foulkes-Haydock 1989 PRB 39 12520.
- **QE 7.5 source:** cited by filename and line number whenever a
  specific convention is invoked (e.g., PZ constants as Ha vs Ry,
  V_local G=0 formula, density symmetrization phase).

### §1.4 Finding categories

- **A (bug):** code diverges from docstring or literature in a way
  that changes numerical results.
- **B (ambiguity):** code and docstring consistent, but literature has
  multiple conventions; annotate which one pwdft-rs chose and why it
  still matches QE.
- **C (doc gap):** code correct but docstring is silent or imprecise
  on a subtle point.
- **D (unit-boundary concern):** possible missed or miscomputed
  conversion at the UPF / SI / Ry→eV interface.

---

## §2 Findings

Ordered by module. Within each module, findings are grouped by
category. Code snippets quote the file on `origin/main` @
`1f568be` (MAUD branch point).

### §2.1 `src/scf/energy.rs`

Read the post-MADOC docstrings carefully; the full equation for every
energy helper is now present. No **Category A** findings in this file:
every formula in the code matches the MADOC docstring's stated
equation.

#### Category C (doc gap)

**C-ENG-1. `band_energy` docstring does not distinguish between
nspin=1 and nspin=2 occupation conventions.**

The docstring says `f_{n,k} dimensionless occupation in [0, spin_factor]`
with `spin_factor = 2` for nspin=1 and `1` for nspin=2. In nspin=2 the
`eigenvalues_all` array is `[up_k0, up_k1, …, down_k0, down_k1, …]`
(driver_spin.rs:186), and the matching `weights_all` repeats the
k-weights. That works because the band sum `Σ f·w·ε` doesn't care
about the index structure. Fine as coded — but the docstring does not
mention that the spin driver lays out eigenvalues this way, and a
future change that factors `Vec<Vec<Vec<f64>>>` with an explicit
spin axis would silently change the normalization. *Doc-only.*

**C-ENG-2. `hartree_energy` docstring omits the "G=0 excluded"
qualifier in the displayed formula.**

The code at `scf/energy.rs:94-109` correctly tests
`g2 > G2_ZERO_THRESHOLD` to skip the DC mode, and the surrounding
prose explains the compensating background. But the displayed
equation reads

```text
E_H = (Ω/2) · Σ_{G ≠ 0} |ρ(G)|² · 4πe² / |G|²
```

which does carry the `G ≠ 0` subscript — actually fine. Re-read: not
a gap. **Dismissed.**

**C-ENG-3. `total_energy` / `harris_foulkes_energy` docstrings do not
explicitly state that `e_ewald` is already in eV and already includes
the screened-charge background.**

`ewald.rs::ewald_energy` returns `E_real + E_recip + E_self + E_bg`
with the `E_bg = -πe²(Σ Z_i)²/(2Ωη²)` background-charge term baked in.
The total-energy caller in both drivers simply passes `ctx.e_ewald`
through. If someone reads the `total_energy` docstring and assumes
they should manually add a background term, they'll double-count. A
one-line reference to `crate::ewald::ewald_energy` from inside the
`total_energy` docstring would close this. *Doc-only.*

**C-ENG-4. `add_core_density` docstring doesn't justify the `.max(0.0)`
clamp's effect on Harris-Foulkes stationarity.**

At `energy.rs:355`: `(v + c).max(0.0)`. The clamp is there to prevent
`ρ^(1/3)` NaN in exchange. But negative-density regions are physical
in well-converged SCF only at round-off amplitude; at poorly
converged points (`Δρ_in → ρ_out` large) they can be non-trivial.
Clamping introduces an O(|negative tail|) non-differentiability
which in principle breaks the stationarity argument that makes E_HF
quadratic in Δρ. In practice the clamp threshold is `0`, the tail
amplitude at Si convergence is ~1e-16, and nothing goes wrong; but
the MADOC docstring asserts O(Δρ²) stationarity without flagging
that the clamp is a mild nonlinearity. *Doc-only.*

### §2.2 `src/scf/driver.rs`

Module header now states the full Kohn-Sham equations and 8-step
per-iteration pipeline (see `/tmp/maud_driver.rs:1-40` against
`src/scf/driver.rs:1-10` on disk — the file on disk at MAUD branch
point pre-dates the MADOC commit but HEAD has the updated version;
reviewers should check HEAD:src/scf/driver.rs).

No **Category A** findings. Code matches the pipeline the module
header declares, step for step.

#### Category C (doc gap)

**C-DRV-1. Step 6 symmetrization happens BEFORE the output density is
used in any energy term.**

The module header (post-MADOC):

> 6. Reconstruct the output density and symmetrize in G-space via
>    analytic fractional-translation phase factors (PCFX).

In the code at `driver.rs:232-255`:

```rust
let mut rho_r_new = density::compute_density(...);   // raw from bands
crate::symmetry::density::symmetrize_density_g(&mut rho_r_new, ...);  // in place
```

then `rho_r_new` is fed to `hartree_energy`, `xc_energy_corrected`,
and the VGC5 `e_kinetic` / `e_local` / `e_nonlocal` evaluations. The
docstring should state explicitly that **all energy components at
convergence are evaluated on the symmetrized output density**, not on
the raw band-reconstructed density. This matters for the direct-sum
identity in `EnergyComponents`: if someone ever swapped the order,
the identity would only hold modulo symmetry residual (which, pre-
PCFX, was as large as 1.2 eV on Si). *Doc-only.*

**C-DRV-2. Harris-Foulkes input-density identification is implicit.**

The code evaluates `hartree_energy(&rho_g, …)` for E_HF, where
`rho_g` was last updated at the end of the previous iteration's
mix step (line 392). This is `ρ_in` for the current iteration,
which matches the HF definition — but the docstring on
`harris_foulkes_energy` calls this "input density" without tying it
to the variable name in the driver. A reader navigating the driver
may not immediately map `rho_g` → ρ_in. *Doc-only.*

### §2.3 `src/scf/driver_spin.rs`

Module header (post-MADOC) derives the LSDA + CCMX basis change in
full. See `/tmp/maud_driver_spin.rs:1-57`.

No **Category A** findings; the docstring-vs-code match is clean for
the CCMX basis inverse, NLCC ρ_core/2 per-channel split, and SPNC
per-spin Δρ convergence scalar.

#### Category C (doc gap)

**C-DRS-1. Initial split uses average magnetization, not per-atom.**

At `driver_spin.rs:80`:

```rust
let avg_mag: f64 = per_atom_mag.iter().sum::<f64>() / per_atom_mag.len().max(1) as f64;
let mut rho_up_r = rho_total.iter().map(|&r| r * (1.0 + avg_mag) / 2.0).collect();
```

The initial split is uniform: every point in real space gets the same
ρ_up/ρ_down ratio `(1+avg_mag)/2`. For an alloy with different
`starting_magnetization` per element (e.g., Fe3Al), this loses the
atom-resolved magnetization at iteration 0. The module header
references SPNC and CCMX but does not flag that the SAD split is
coarser than per-atom. In practice: after a few iterations the SCF
finds the right per-atom magnetization and the initial-condition loss
is forgotten. But for delicate magnetic cases (antiferromagnetic FeO)
a proper per-atom initial split would seed the right ordering much
faster. *Doc-only for now — potential follow-up for Researcher if a
real AFM case fails to converge.*

**C-DRS-2. VGC5 spin-polarized `e_xc` uses the coarse spin-averaged
ρ_core.**

At `driver_spin.rs:433`:

```rust
let e_xc_term = xc_energy_bare(&rho_xc_total, &exc_r_out, ctx.omega);
```

`rho_xc_total = ρ_total_new + ρ_core` (full ρ_core, not split). This
is algebraically correct since the two spin channels' `ρ_σ + ρ_core/2`
each carry half of ρ_core, and `Σ_σ (ρ_σ + ρ_core/2) = ρ_total +
ρ_core`. But the VGC5 docstring on `EnergyComponents::e_xc` says
"ρ here is ρ_val + ρ_core", which is ambiguous about whether `ρ_val`
refers to total valence (yes) or per-spin valence (no). *Doc-only.*

### §2.4 `src/potential/xc.rs`

The core PZ/Slater formulas are unchanged by MADOC or GGAP Phase A;
GGAP Phase A only added the `XcEvaluator` dispatch enum without
touching the per-point math. Constants (γ = −0.1423, β₁ = 1.0529,
β₂ = 0.3334 for paramagnetic; γ = −0.0843, β₁ = 1.3981, β₂ = 0.2611
for ferromagnetic; A = 0.0311, B = −0.048, C = 0.0020, D = −0.0116
for high-density paramagnetic; A = 0.01555, B = −0.0269, C = 0.0007,
D = −0.0048 for high-density ferromagnetic) match PZ 1981 Table XII.

#### Category A (bug — docstring only, but flagged because it will
mislead future maintainers)

**A-XC-1. Spin exchange potential docstring has a spurious `2^{1/3}`
factor.**

`xc.rs:302-306` (function `slater_exchange_spin`) docstring reads:

```text
/// V_x_σ = (4/3)·ε_x(2ρ_σ)·2^{1/3}  [derivative of the spin-scaled exchange]
```

The code at line 334 computes:

```rust
let vx_up_ha = (4.0 / 3.0) * ex_up_ha;
```

where `ex_up_ha = -0.75 · (6ρ_up/π)^(1/3) = ε_x(2ρ_up)` (the
unpolarized exchange evaluated at `2ρ_up`). So the code computes
`V_x_σ = (4/3)·ε_x(2ρ_σ)` with **no** `2^{1/3}` factor.

Derivation that code is right, docstring is wrong:

```text
  E_x^σ[ρ_σ] = ρ_σ · ε_x(2ρ_σ)
  V_x_σ = δE_x/δρ_σ = ε_x(2ρ_σ) + ρ_σ · 2 · dε_x/du|_{u=2ρ_σ}
        = ε_x(2ρ_σ) + u · dε_x/du
        = ε_x(2ρ_σ) + ε_x(u)/3           since ε_x = -C u^{1/3} ⇒ u·dε_x/du = ε_x/3
        = (4/3) · ε_x(2ρ_σ)
```

No `2^{1/3}`. This is a MADOC follow-up: the docstring text in
`slater_exchange_spin` was written for the `ε_x(ρ_σ)` parametrization
(not `ε_x(2ρ_σ)`) and the extra `2^{1/3}` would have been the chain-rule
factor if the scaled-density form had been used. Confusing to anyone
trying to reproduce the derivation.

**Fix category:** doc-only; the actual SCF math is correct. The test
`test_lsda_unpolarized_limit` passes because it only checks that ζ=0
reproduces the unpolarized result — it would pass even if the code
were wrong, as long as both ran the same (wrong) formula.

**Cross-reference:** MADOC follow-up. Flagged in the MADOC
completion note as a possible docstring-cleanup task.

#### Category C (doc gap)

**C-XC-1. Spin correlation potential docstring uses ambiguous `(±1 - ζ)`.**

`xc.rs:349-351` reads:

```text
/// V_c_σ = ε_c - (rs/3)·dε_c/drs + (±1 - ζ)·dε_c/dζ
```

The correct formula is:

```text
V_c_up   = V_c_rs + (1 - ζ)·dε_c/dζ
V_c_down = V_c_rs − (1 + ζ)·dε_c/dζ
```

The docstring's `±1 - ζ` form yields `1 - ζ` for up and `−1 - ζ` for
down, which matches `-(1 + ζ)` — so mathematically equivalent. But
the `±1 - ζ` notation is unusual and the reader has to work out the
sign. The code (line 385-386) writes the two cases out explicitly:

```rust
let vc_up_ha   = (1.0 - zeta).mul_add( dec_dzeta, vc_rs_ha);
let vc_down_ha = (1.0 + zeta).mul_add(-dec_dzeta, vc_rs_ha);
```

*Doc-only.*

**C-XC-2. PZ constants stated without the "Hartree vs Rydberg"
resolution.**

PZ 1981 Table XII lists constants that are routinely reported in the
literature as Rydberg-valued (Payne et al. RMP 64 1045 explicitly
writes `ε_c = -0.1423/(1 + β₁√rs + β₂ rs) Ry`). QE's PZ implementation
(`qe-7.5/XClib/qe_funct_corr_lda_lsda.f90:16-73`) returns the bare
`ec = gc/ox` and the caller multiplies by `e2 = 2` at
`qe-7.5/PW/src/v_of_rho.f90:311` — i.e., QE treats these constants as
**Hartree**. pwdft-rs follows the same Hartree convention (`ec_ha · HA_TO_EV`).

A Si-equilibrium sanity check (rs ≈ 3.32 Bohr at ρ = 0.044 e/Å³):

| Interpretation | ε_c (QE, pwdft-rs) |
|----------------|--------------------|
| Hartree        | −0.961 eV          |
| Rydberg        | −0.481 eV          |

Both QE and pwdft-rs return −0.961 eV, matching the
`test_lda_known_values` pin at `xc.rs:429-448`. So the convention is
**Hartree, not Rydberg**, despite what many papers claim. This is
non-obvious and worth a one-line docstring pin. *Doc-only.*

#### Category B (convention pin, no bug)

**B-XC-1. LSDA interpolation is von Barth-Hedin, not
Vosko-Wilk-Nusair.**

`xc.rs:358-385` uses the classic VBH interpolation

```text
f(ζ) = [(1+ζ)^(4/3) + (1-ζ)^(4/3) - 2] / [2^(4/3) - 2]
```

from von Barth & Hedin J. Phys. C 5 (1972) 1629, not the VWN5
interpolation in Vosko, Wilk, Nusair CanJP 58 (1980) 1200. Both QE's
PZ path (`xc_lsda` → `pz_polarized`) and pwdft-rs share the VBH
choice, so no discrepancy. PZ 1981 itself used VBH. Consistent.
*Docstring should state this explicitly — the reader has no way to
tell from the existing text that VWN was not used.*

### §2.5 `src/potential/nonlocal.rs`

KB assembly via expanded real-Y_lm channel basis, then single GEMM
(VNLM, 2026-04-18). D_ij application is block-diagonal in
(atom, l, m) (`nonlocal.rs:283-306`). No **Category A** findings.

#### Category C (doc gap)

**C-NL-1. Docstring silent on phase factor `i^l` cancellation.**

The module docstring (line 41) says:

> Phase factor `i^{l_α}` cancels exactly: it appears as
> `i^{l_α} · (i^{l_β})*` in H_NL, and the block-diagonal structure
> of D forces `l_α = l_β`, so the combined phase is 1. We omit it
> from B entirely.

This is explained in the module header but repeated nowhere in the
per-projector F_i(q) assembly `bessel_transform_projector`. Anyone
auditing `bessel_transform_projector` in isolation would not see
why the `i^l` factor from the textbook `β_i(k+G) = F_i·i^l·Y_lm` is
missing. *Doc-only; no numerical consequence.*

**C-NL-2. Addition-theorem test pins the sum, not individual
m-channels — defense already exists.**

MADOC notes (via VNMT) that `test_ylm_addition_theorem` is only a
sum-level pin; the complementary per-m-channel test
`test_single_channel_l2_m_isolation` (lines 784-894) pins the
Y_{l,m} normalization individually for `(l=2, m ∈ {0, +1, +2})`. This
is the specific defense against the "trace-equivalent-but-projector-
wrong" failure mode called out in Researcher's agent definition. **Good
as-is.** Noting here for audit completeness.

### §2.6 `src/potential/local.rs` and `src/pseudopotential/mod.rs::v_local_of_g`

No **Category A** findings. VGCMP Phases 1-4 already pinned
`V_local(G)` bit-for-bit vs. QE for Si to 1e-9 Ry.

#### Category B (convention pin)

**B-VL-1. erf-subtraction Gaussian width is 1 Å in pwdft-rs vs.
1 Bohr in QE.**

`pseudopotential/mod.rs:167` computes `Z·e²·erf(r)/r` with `r` in Å,
i.e., rc = 1 Å. QE (`qe-7.5/upflib/vloc_mod.f90:138`) uses `erf(r)`
with `r` in Bohr, i.e., rc = 1 Bohr.

The decomposition `V_loc(r) = short(r) − Z·e²·erf(r/rc)/r` is exact
for **any** rc, with the compensating tail
`FT[Z·e²·erf(r/rc)/r] = 4π·Z·e²·exp(-G²·rc²/4)/G²`. pwdft-rs carries
`exp(-G²/4)` (G in 1/Å, rc = 1 Å); QE carries `exp(-G²·tpiba2/4)`
(G in lattice units, tpiba2 translates to Bohr so rc = 1 Bohr). Both
agree to machine precision after unit conversion — see
`tests/vloc_erf_consistency.rs` and VGCMP.

**Pin:** docstring already says "different convention from QE's 1
Bohr, but yields identical V_local(G) values because the decomposition
is exact for any Gaussian width" (`pseudopotential/mod.rs:128-130`).
Already well-documented.

#### Category D (potential unit-boundary concern)

**D-VL-1. `v_local_of_g` at heavy atoms (Z > 14) — VGCMP assumed
weak Z-dependence.**

The erf-subtracted form is numerically stable *for our log mesh* at
Si (Z = 14) and Fe (Z = 26), where the bracketed integrand
`r·V_loc(r) + Z·e²·erf(r)` is bounded everywhere by O(Z) on the PP
log mesh. For Z → 92 (U) the first grid points see integrand values
of order Z·e² ≈ 1e3 eV·Å, and the `sin(Gr)/Gr` beating against a
log-mesh Simpson with dr ≈ 1e-3 Å at small r may lose O(1e-4) relative
precision. VGCMP Phase 1 only verified Si (Z = 14). Heavy-element
pseudopotentials (Cs, W, Pt, U) are **not** in the validation set yet.

**This cross-links to the existing VGCH proposal** (High — Validation
in INDEX.md): VGCH already owns "Heavy-atom V_local(G) residual —
post-VGCMP continuation (closes 5 Z>14 `#[ignore]`s)" and is the
right home for a numerical QE pin at Z ≥ 50. No new proposal needed;
MAUD just annotates that the Z-dependence concern is real and
motivates VGCH. Flagged as **cross-link** in §4 below.

### §2.7 `src/ewald.rs`

#### Category A (borderline — behavior diverges under pathological
cells, unused in production)

**A-EW-1. Cutoff heuristic is fixed at `10η` regardless of truncation
error; QE uses an iterative α-optimization to target 1e-7.**

`ewald.rs:60-65`:

```rust
let eta = (n_atoms as f64 * PI / omega).cbrt();
let g_max = 10.0 * eta;
```

vs. QE's `qe-7.5/PW/src/ewald.f90:93-101`:

```fortran
alpha = 2.9d0
100 alpha = alpha - 0.1d0
...
upperbound = 2.d0 * charge**2 * SQRT(2.d0 * alpha / tpi) * &
      erfc( SQRT(tpiba2 * gcutm / 4.d0 / alpha) )
IF (upperbound > 1.0d-7) goto 100
```

QE iteratively decreases α until the upper bound on truncation error
falls below 1e-7. pwdft-rs uses a fixed heuristic. At `10η`,
`exp(-(10η)²/(4η²)) = exp(-25) ≈ 1.4e-11` — well below 1e-7 for
normal cells — but for highly anisotropic cells (slab geometry, 10×10×2
unit cell, `test_ewald_anisotropic_cell`) the η = (Nπ/Ω)^(1/3) = 1.16
Å⁻¹ is not optimal for both real and reciprocal directions simultaneously
and `g_max = 11.6 Å⁻¹` may leave a larger tail than the 1e-7 target.

**Classification:** this is **not a Category A bug in practice** —
the test suite already exercises anisotropic cells and the test
`test_ewald_nacl` pins to 1% of the Madelung result, which neither
formulation misses. But it could be flagged at future validation
thresholds (phonon calculations, stress tensor for slabs).

**Fix category:** follow-up proposal candidate. Not MAUD-priority.

#### Category C (doc gap)

**C-EW-1. No reference cited for the η = (Nπ/Ω)^{1/3} heuristic.**

The docstring (lines 30-31) gives the formula but no citation. This
is the "optimally balanced" choice for isotropic cubic cells (roots
of d/dη of the combined real+recip cost function), but standard
texts (Martin ch. 13, Payne et al. §IV.B) give it a name and a
derivation. *Doc-only.*

**C-EW-2. Background term `E_bg` sign and zero-net-charge note
not in the equation block.**

Docstring shows `E_bg = -πe²(Σ Z_i)²/(2Ωη²)`. For neutral cells
`Σ Z_i = 0` and `E_bg = 0`. This is worth stating — the
`test_ewald_zero_charges` test pins the trivial case but future
charged-cell (defect calculation) users will need to know
exactly what E_bg does. *Doc-only.*

### §2.8 `src/pseudopotential/upf/convert.rs`

#### Category D (potential unit-boundary concern — no actual bug found,
but audit completeness worth a pin)

**D-UPF-1. Every UPF block's conversion is covered; audit
enumeration.**

| UPF block | Storage unit | pwdft-rs conversion | Line |
|-----------|--------------|---------------------|------|
| `PP_R`    | Bohr         | × BOHR_TO_ANG       | 48   |
| `PP_RAB`  | Bohr         | × BOHR_TO_ANG       | 52   |
| `PP_LOCAL`| Ry           | × RY_TO_EV          | 56   |
| `PP_BETA` | Bohr^{-1/2} (stores χ = r·β) | ÷ √BOHR_TO_ANG | 72-75 |
| `PP_DIJ`  | Ry           | × RY_TO_EV          | 82   |
| `PP_RHOATOM` | e/Bohr (stores 4πr²·ρ_at) | ÷ BOHR_TO_ANG | 92 |
| `PP_NLCC` | e/Bohr³ (bare ρ_core) | ÷ BOHR3_TO_ANG3 | 121 |
| `PP_Q`    | USPP — not supported | — | N/A |
| `PP_PSWFC`| Atomic wavefunctions — not used in LDA SCF path | — | N/A |

Every block that enters the SCF path is converted. `PP_Q` (augmentation
charges) is not parsed because USPP is not supported; if USPP is ever
enabled, this audit must be redone. *Audit clean. No fix needed.*

**D-UPF-2. `PP_NLCC` conversion is the **only** place in the codebase
where `BOHR3_TO_ANG3` is applied to a density.**

Cross-checked: `PP_RHOATOM` uses `BOHR_TO_ANG` (linear, because it
stores `4πr²·ρ`, which has dimensions `1/length` when ρ is in
`e/Bohr³` and `4πr²` in `Bohr²`). `PP_NLCC` uses `BOHR3_TO_ANG3`
(volumetric, because it stores bare `ρ_core` in `e/Bohr³`). This
asymmetry is the NCFX bug in reverse; NCFX caught it going the
wrong way. The current code is correct and the NLCC ρ_core(G)
pins at `convert.rs:303-411` (Si + Fe) regression-guard it. Good.

### §2.9 `src/symmetry/density/g_space.rs`

No **Category A** findings. Convention matches QE after the `S → S⁻¹`
relabelling the docstring explains.

#### Category B (convention pin — pwdft-rs direct vs. QE transpose)

**B-SYM-1. Rotation sign convention differs between pwdft-rs and QE.**

`SpaceGroupOp.rotation` in pwdft-rs is the **direct-space fractional**
rotation `R` acting as `r' = R·r + τ`. QE stores its `s(:,:,ns)` as
**R^T** (see `qe-7.5/PW/src/symm_base.f90:533` where atoms rotate as
`rau = s^T · xau`). The G-space action in both codes is the same —
`ρ(G) → e^{-iG·τ} · ρ(R^T · G)` — but that happens to be
`e^{-iG·τ} · ρ(s·G)` in QE notation and `e^{-iG·τ} · ρ(R_direct^T · G)`
in pwdft-rs notation, i.e., the *implementations* of "apply the
rotation" differ by a transpose that cancels exactly.

This is pinned by the test `test_symmetrize_g_matches_real_space_on_compatible_grid`
(lines 427-461) which cross-validates the G-space form against the
deprecated real-space form on a symmorphic grid. Already well-documented.

**Pin:** see docstring at `g_space.rs:112-117`. Already annotated.

#### Category C (doc gap)

**C-SYM-1. Band-limitation requirement only enforced informally.**

`g_space.rs:124-146` states that for `P² = P` to hold exactly the
input density must be band-limited to `|G|² < |G|²_Nyquist / |R|²_max`
to prevent aliasing. The default `ecutrho_ratio = 4` ensures this.
But there's no runtime `debug_assert!` checking that the input
density actually satisfies the cutoff — a user who passes a
hand-crafted density through `symmetrize_density_g` could silently
alias. Since this helper is only called by `run_scf` it's not
exposed, but the docstring "To avoid this, the input density must be
band-limited …" suggests a contract the code does not enforce.
*Doc-only — not reachable from user input today.*

---

## §3 Recommended fix priority

### §3.1 Category A findings (1)

Only one true A-level finding:

- **A-XC-1** (spin exchange potential docstring has spurious `2^{1/3}`):
  this is a *docstring* bug, not a math bug. The code is right; the
  docstring misleads anyone trying to reproduce the derivation or who
  generalizes the formula for LSDA + gradient correction. **Fix in a
  MADOC follow-up PR** (or as part of GGAP Phase B, where the LSDA
  potential formulas get re-read end-to-end).

- **A-EW-1** (Ewald fixed-cutoff heuristic): classified borderline
  above. In practice passes all tests. **Not a MAUD fix; candidate for
  a separate follow-up proposal if/when validation tightens (phonon
  dispersions, slab-geometry total energies).**

### §3.2 Category B findings (3): pin as doc-only cleanups

- B-XC-1 (VBH vs VWN): add one line in `pz_correlation_spin`.
- B-VL-1 (1 Å vs 1 Bohr Gaussian width): already well-documented.
- B-SYM-1 (rotation direct vs QE transpose): already well-documented.

### §3.3 Category C findings (9): roll up into a MADOC-B cleanup

Most C findings are minor clarifications — 1 to 3 sentences per
docstring. Recommend a single PR that sweeps all of them as
"MADOC-B: follow-up docstring clarifications from MAUD audit."

Priority subset if MADOC-B is deferred:

1. **C-DRV-1** — state that VGC5 energy terms are post-symmetrization.
   Highest-value doc clarification because a future maintainer who
   rearranges the symmetrize-vs-diagonalize order would silently break
   the per-component identity.
2. **C-ENG-4** — note the `max(0.0)` clamp's interaction with HF
   stationarity.
3. **C-XC-2** — pin "Hartree, not Rydberg" for PZ constants.

### §3.4 Category D findings (2): VGCH and heavy-atom validation

- **D-VL-1**: the `v_local_of_g` erf-subtracted form is pinned at
  Z = 14 (Si) and Z = 26 (Fe) but not above. The VGCMP Phase 1 comment
  in `.claude/logbooks/researcher.md` notes this as a known gap. The
  existing **VGCH** proposal (High — Validation) owns this continuation;
  MAUD's D-VL-1 is a motivating annotation, not a new action.
- **D-UPF-1 / D-UPF-2**: audit clean; document the per-block
  conversion table (above) in a module docstring for
  `src/pseudopotential/upf/convert.rs`. One-time PR.

### §3.5 Suggested landing order

1. **No action required** on any Category A code path — SCF numerics
   are correct.
2. **MADOC-B doc sweep** (highest value / lowest risk): sweep all C
   findings plus A-XC-1 (which is a doc fix only). One PR, one
   reviewer.
3. **VGCH heavy-atom proposal** (D-VL-1): to be filed by Researcher as
   separate proposal when time permits.
4. **A-EW-1 follow-up** (Ewald iterative α): defer until a validation
   failure surfaces.

---

## §4 Cross-links

### §4.1 MADOC follow-up (MADOC-B)

- **A-XC-1** (spin exchange docstring `2^{1/3}` spurious factor):
  docstring fix only; the code is correct. Add to MADOC's "Phase B
  follow-up" list.
- **C-ENG-3, C-ENG-4, C-DRV-1, C-DRV-2, C-DRS-1, C-DRS-2, C-NL-1,
  C-XC-1, C-XC-2, B-XC-1, C-EW-1, C-EW-2, C-SYM-1**: all doc-only,
  roll up into a single "MADOC-B MAUD cleanup sweep" PR.

### §4.2 VGCH cross-link (existing proposal)

- **D-VL-1**: the existing **VGCH** proposal already owns the heavy-
  atom V_local validation. MAUD's D-VL-1 adds "Z-dependence of the
  erf-subtracted form factor precision" as one concrete motivation for
  VGCH to include at least one Z ≥ 50 PP (Cs, W, or Pt) in its test
  matrix.

### §4.3 Ewald follow-up (low priority)

- **A-EW-1**: if anyone wants tighter Ewald accuracy for slab or
  charged-defect calculations, port QE's iterative α-optimization
  (`ewald.f90:93-101`). Not MAUD-blocking.

### §4.4 GGAP Phase B

- **A-XC-1** and **C-XC-2** should be revisited when GGAP Phase B
  ports the PBE potential formulas — those derivations reuse the
  spin-scaled exchange bookkeeping and the Hartree convention for the
  correlation constants. A clean docstring at this branch point saves
  GGAP Phase B one round of confused reverification.

---

## §5 Flagged for follow-up (out of scope for MAUD)

- **src/scf/initial_density.rs:80** — uniform initial spin split from
  *average* `starting_magnetization` loses per-atom resolution for
  mixed magnetic systems. Researcher follow-up if a real AFM case
  fails to converge from SAD. (See C-DRS-1.)

- **src/ewald.rs:60** — fixed `η = (Nπ/Ω)^{1/3}`, cutoff `10η` is
  non-adaptive. Researcher follow-up after phonon or slab validation
  surfaces the limitation. (See A-EW-1.)

- **src/pseudopotential/mod.rs (v_local_of_g)** — erf-subtracted form
  not validated for Z ≥ 50. Researcher follow-up as VGCH proposal.
  (See D-VL-1.)

No Rust-idiom / dead-code / perf / production-code issues surfaced
during this audit; the modules in scope have already been cleaned
through the TACC, QLNT, QLN2, and CAST passes.

---

## §6 What MAUD explicitly does NOT cover

- GGAP Phase B (PBE exchange-correlation) — out of scope; MAUD is a
  cold read of current state only.
- HYBR (hybrid functionals) — out of scope.
- USPP / PAW — not implemented in pwdft-rs; no code to audit.
- Stress / forces — separate helper layer, not touched by MADOC Phase A.
- Spin-orbit coupling — not implemented.

---

## §7 Acceptance

This proposal is a *report*, not a code change. Acceptance = EM reads
§2 findings, agrees on the §3 priority, and decides which (if any) of
the §4 follow-ups to file as standalone proposals. No implementation
work falls directly out of MAUD itself.

Expected downstream PRs:
1. **MADOC-B cleanup sweep** (one PR, MADOC owner, bounds A-XC-1 + C-*)
2. **VGCH heavy-atom validation** (separate proposal, Researcher)
3. **(deferred)** Ewald iterative α (separate proposal if needed)
