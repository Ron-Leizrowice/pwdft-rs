---
id: GGAP
title: GGA/PBE exchange-correlation functional (umbrella, phased A-F)
priority: high
complexity: large
risk: medium
depends_on: []
blocks: []
status: draft
author: Researcher
date: 2026-04-18
---

# GGAP — GGA/PBE exchange-correlation functional

## Problem

pwdft-rs currently ships a **single XC functional**: Perdew-Zunger 81 LDA
(`src/potential/xc.rs`, `lda_xc_grid` / `lda_xc_spin_grid`). Every modern
pseudopotential library (PseudoDojo, SSSP, GBRV, ONCVPSP releases since
~2016) ships its *primary* set at the GGA/PBE level. Our own tree
already contains a complete PBE PP set at `pseudopotentials/nc/pbe/`
(71 elements, Ag–Zr) that the solver cannot correctly consume.

Missing PBE support means:

1. **No comparison to published DFT data.** Virtually every materials
   database (Materials Project, AFLOW, OQMD) reports PBE total energies,
   band gaps, and formation enthalpies. LDA total energies are not
   directly comparable — they systematically over-bind by ~0.5–1 eV per
   bond.
2. **Wrong physics for magnetic systems.** LDA famously predicts Fe
   ground state as FCC non-magnetic (vs. the experimentally observed BCC
   ferromagnetic). PBE corrects this. Our Fe BCC test
   (`qe_validation.rs::test_fe_bcc_fm_vs_qe`, currently `#[ignore]`) may
   benefit from PBE in ways unrelated to the CCMX mixer investigation.
3. **Footgun already armed.** `src/settings.rs:335-345` defines
   `XcFunctional::{Pz, Pbe, Pbe0, Hse06}` — a YAML-parseable enum with a
   `pbe_functional_parses()` test at `src/settings.rs:1098`. The driver
   ignores the enum entirely and always calls LDA. A user who writes
   `functional: pbe` in their YAML gets LDA results silently. This is
   exactly the class of footgun NCFX/PCFX/CCMX trained us to look for.

This proposal phases in PBE support across six focused sub-PRs, each
independently reviewable but coordinated under one umbrella so the
reviewer can see the whole shape before Phase A starts.

## Scope

**In scope:** PBE (Perdew-Burke-Ernzerhof 1996, `iflag=1` in QE's family
of GGA variants) for non-spin-polarized (`nspin=1`) and spin-polarized
(`nspin=2`, collinear) SCF. CPU path only; GPU deferred to Phase E.
Pseudopotential functional-tag consistency check at load time.

**Out of scope** (explicitly):

- PBEsol, revPBE, RPBE, BLYP, PW91 and other GGA siblings. Trivial to
  add *once* PBE's scaffolding is in place (QE's `pbex`/`pbec` switch on
  `iflag`), but each needs its own validation pass. **Add on demand.**
- Meta-GGA (SCAN, r²SCAN, TPSS). Requires kinetic-energy-density τ(r)
  from |∇ψ|² — architectural change, separate proposal.
- Hybrid functionals (PBE0, HSE06) — require exact-exchange integrals,
  very large separate proposal.
- libxc C library integration. See §"Build vs buy" below.

## Physics

All equations below cross-reference the **paper** (Perdew-Burke-Ernzerhof,
*Phys. Rev. Lett.* **77**, 3865 (1996); erratum *PRL* **78**, 1396
(1997)) AND the **QE 7.5 source of truth** (per `.claude/agents/researcher.md`:
QE's Fortran is canonical, the paper is a reference). Both point to the
same numbers — matching both is our definition of "correct".

### PBE exchange (non-spin)

**Enhancement-factor form (PBE §III, eq. 14):**

    ε_x^PBE(ρ, ∇ρ)  =  ε_x^LDA(ρ) · F_x(s)

where

    ε_x^LDA(ρ)      =  −(3/4π)·k_F                    (eV per electron)
    k_F             =  (3π²ρ)^{1/3}                   (Å⁻¹, Bohr units)
    s               =  |∇ρ| / (2 k_F ρ)               (dimensionless)
    F_x(s)          =  1 + κ − κ/(1 + μs²/κ)          (eq. 14)

**Constants (PBE eq. 13, Table I of `Appendix` footnote of PRL 77:3865):**

    κ  =  0.804                                       (the "LO" bound)
    μ  =  β · π² / 3  =  0.21951492776...             (PBE choice)
    β  =  0.06672455060314922                         (PBE eq. 11)

**QE ground-truth:** `qe-7.5/XClib/qe_funct_exch_gga.f90:111-331`,
subroutine `pbex` with `iflag=1`. The constant `k(1) = 0.804_DP` at
line 160 and `mu(1) = 0.2195149727645171_DP` at line 162 match our
paper numbers to machine precision. QE returns three outputs per
grid point:

- `sx`  = ρ · (ε_x^PBE − ε_x^LDA), i.e., the *gradient correction* to
  the exchange energy density (the LDA slater part is added elsewhere).
- `v1x` = ∂(ρ·ε_x^PBE) / ∂ρ        (Ry/Å³ → eV/Å³ after our unit conv.)
- `v2x` = ∂(ρ·ε_x^PBE) / ∂(|∇ρ|²)  (eV·Å² per density unit)

This two-potential shape (`v1x`, `v2x`) is the standard GGA interface
and matches libxc's `vrho` / `vsigma` outputs. The total semilocal
potential is assembled by eq. (XX) below.

### PBE correlation (non-spin)

**Form (PBE §IV, eq. 7):**

    ε_c^PBE(ρ, ζ=0, ∇ρ) = ε_c^LDA(ρ) + H(r_s, ζ=0, t)

where the gradient correction is (PBE eq. 7):

    H(r_s, ζ, t) = γ · φ(ζ)³ · ln[ 1 + (β/γ) t² · (1 + A t²) / (1 + A t² + A² t⁴) ]
    A            = (β/γ) / [ exp(−ε_c^LDA / (γ φ(ζ)³ )) − 1 ]      (eq. 8)
    t            = |∇ρ| / (2 φ(ζ) k_s ρ)                           (eq. 9)
    k_s          = sqrt(4 k_F / π)                                 (TF wavenumber)
    φ(ζ)         = ½[(1+ζ)^{2/3} + (1−ζ)^{2/3}]                    (eq. 3)

**Constants (PBE eq. 9/10, footnote):**

    β  = 0.06672455060314922       (same as PBE exchange)
    γ  = (1 − ln 2) / π²           ≈ 0.031090690869654895

**QE ground-truth:** `qe-7.5/XClib/qe_funct_corr_gga.f90:195-259`,
subroutine `pbec` with `iflag=1`. Note QE's internal name for γ is
`ga` (line 214, `ga = 0.0310906908696548950_DP`) and for β is `be(1)`
(line 217, `0.06672455060314922_DP`). These match the PBE paper
values and our formulas to all digits. QE's `pbec` calls `pw(rs, 1,
ec, vc)` (line 229) for the LDA correlation — **we must use the same
PW92 correlation, not PZ**, for PBE consistency (see "Hidden PZ-vs-PW92
issue" below).

### Spin-polarized PBE (`nspin=2`, collinear)

**Exchange (spin scaling — exactly Perdew-Wang 1988):**

    ε_x^PBE(ρ↑, ρ↓, ∇ρ↑, ∇ρ↓) = (1/ρ) [ρ↑ · ε_x^PBE(2ρ↑, 2∇ρ↑)
                                         + ρ↓ · ε_x^PBE(2ρ↓, 2∇ρ↓)]

— i.e., each spin channel evaluates the unpolarized formula at 2ρ_σ
with its own ∇ρ_σ, weighted by ρ_σ/ρ. This is the Oliver-Perdew spin
scaling relation. QE calls `pbex(2*rho_up, 4*grho_up_sq, iflag, ...)`
and `pbex(2*rho_down, 4*grho_down_sq, iflag, ...)` separately (see the
wrappers in `qe-7.5/XClib/qe_drivers_gga.f90`, `gcxc_spin`). The factor
of 4 on |∇ρ_σ|² is because QE's input is σ=|∇ρ|², and ∇(2ρ_σ)² =
4|∇ρ_σ|².

**Correlation (ζ-dependent H):**

    ε_c^PBE(ρ, ζ, ∇ρ) = ε_c^LDA(ρ, ζ) + H(r_s, ζ, t)

with the *same* H(r_s, ζ, t) formula, now using the actual ζ and the
scaled t (t = |∇ρ_total| / (2 φ(ζ) k_s ρ_total)). **Note** that only
∇ρ_total enters t — the gradient of the total density, not ∇ρ↑ and
∇ρ↓ separately. This is a subtle point:

- Exchange: needs ∇ρ↑ and ∇ρ↓ individually (per-channel).
- Correlation: needs ∇ρ_total = ∇ρ↑ + ∇ρ↓ (global).

QE ground-truth for the spin path: `qe-7.5/XClib/qe_funct_corr_gga.f90:446-541`,
subroutine `pbec_spin`. Input is `grho = |∇ρ_total|²` (line 462,
rho_tot = ρ↑ + ρ↓). Output `v1c_up` and `v1c_dw` differ because the
φ(ζ) and dε_c^LDA/dρ_σ factors are spin-dependent.

**σ definition (common source of sign/factor bugs):** QE's
`gcxc_spin` driver uses the four gradient-squared quantities
(`qe-7.5/XClib/qe_drivers_gga.f90`):

    sigma_up_up   = |∇ρ↑|²
    sigma_up_down = ∇ρ↑ · ∇ρ↓        (signed, not absolute-squared)
    sigma_down_down = |∇ρ↓|²
    sigma_total   = |∇ρ_total|² = σ_up_up + 2·σ_up_down + σ_down_down

Our implementation should mirror this: four buffers during the spin
path, pass only the needed ones into each sub-function.

### The semilocal V_xc (the footgun)

For a GGA functional ε_xc(ρ, ∇ρ), the XC potential is:

    V_xc(r) = δE_xc/δρ(r) = ∂(ρ·ε_xc)/∂ρ  −  ∇ · [∂(ρ·ε_xc) / ∂(∇ρ)]
            = v1xc(r)                       −  ∇ · [ 2 v2xc(r) · ∇ρ(r) ]   (*)

where `v1xc ≡ v1x + v1c` and `v2xc ≡ v2x + v2c` are the two partial
derivatives returned by `pbex` / `pbec`. The **sign of the divergence
term (*) is the most common GGA implementation bug.** See QE
`qe-7.5/XClib/qe_drivers_gga.f90` where `v2c` is assembled into `h(r)`
then `∇·h` is formed; the convention is:

    V_xc = v1xc(r)  −  ∇ · [ h(r) ]
    h(r) = 2 · v2xc(r) · ∇ρ(r)          (note the factor of 2; v2 is
                                          ∂/∂σ, σ = |∇ρ|², chain rule
                                          gives 2∇ρ)

We need two FFTs (or three — see below) per SCF iteration for this:

1. **∇ρ(r) step:** FFT ρ_r → ρ_G; multiply by (iG_x, iG_y, iG_z);
   IFFT three times → (∇ρ)_x(r), (∇ρ)_y(r), (∇ρ)_z(r).
2. **Assemble `h(r) = 2 v2xc(r) · ∇ρ(r)`** as a 3-component vector
   field in real space.
3. **∇·h step:** FFT each component → h_x(G), h_y(G), h_z(G);
   sum `iG_x · h_x(G) + iG_y · h_y(G) + iG_z · h_z(G)`; IFFT to real
   space → (∇·h)(r).

Total: 3 forward + 4 inverse FFTs per SCF iter (one ρ→∇ρ set + one
h→∇·h set). Compare LDA: 1 forward FFT (V_xc^LDA into G-space for
V_eff assembly). Cost increase: ~7× FFT work for the XC step alone.
On production grids (48³ = 110k points) this is ~O(few ms) per iter
— negligible vs. eigensolve.

### NLCC + GGA

Non-linear core correction (NCFX, landed): the XC functional sees
ρ_total + ρ_core, not just valence. For PBE, the gradient must also
include the core: **∇(ρ + ρ_core) = ∇ρ + ∇ρ_core**. Since ρ_core
is time-independent across the SCF loop, its gradient is precomputed
once (same FFT trick: ρ_core(G) → multiply by iG → IFFT). Stored
alongside `ScfContext.rho_core_r` as `rho_core_grad_r: [Vec<f64>; 3]`.

For `nspin=2`: ρ_core is split equally (half to each spin channel)
per NCFX. Gradient inherits the same split.

### ε_c^LDA inside PBE (hidden PZ-vs-PW92 issue)

**This is a subtle correctness trap.** PBE correlation uses the
**Perdew-Wang 1992** (PW92) parametrization of the Ceperley-Alder
correlation energy, **not** Perdew-Zunger 81 (PZ). The two LDA
parametrizations agree to ~0.1 meV/electron on typical densities but
they are *not identical*, and the PBE gradient correction was fitted
assuming PW92. QE is strict about this — `pbec` at line 229 calls
`pw(rs, 1, ec, vc)`, which is the PW92 subroutine, never `pz(rs, ...)`.

We currently have only PZ (`perdew_zunger_correlation` and
`pz_correlation_rs` in `src/potential/xc.rs`). **Phase C of this
proposal adds PW92 as a private helper inside the GGA module**, used
only by PBE. The existing LDA path continues to call PZ (unchanged,
regression-safe).

The effect on PBE total energies of swapping PZ for PW92 is
<10 meV/electron but **it is the difference between "matches QE" and
"doesn't"**, and we will spend days hunting it if we don't get it right
the first time. The VGC5-style per-component validation is unforgiving
on this scale.

## Numerics

### Density gradients via FFT

`∇ρ(r) = iG · ρ(G)` is five lines of code reusing `src/fft.rs`:

```rust
// Pseudocode; Phase A lands the real version.
pub fn density_gradient(
    rho_r: &[f64],
    fft: &mut FftPlan,
    g_vectors: &[[f64; 3]],      // G in inv-Å, native FFT ordering
) -> [Vec<f64>; 3] {
    let mut rho_g = vec![Complex64::zero(); rho_r.len()];
    fft.forward_real(rho_r, &mut rho_g);
    // (∇ρ)_α(G) = i G_α ρ(G), α ∈ {x,y,z}
    [0, 1, 2].map(|alpha| {
        let mut grad_alpha_g: Vec<Complex64> = g_vectors.iter().zip(&rho_g)
            .map(|(&g, &rho_g)| Complex64::i() * g[alpha] * rho_g)
            .collect();
        let mut grad_r = vec![0.0f64; rho_r.len()];
        fft.inverse_real(&grad_alpha_g, &mut grad_r);
        grad_r
    })
}
```

Not a library. Do not propose one. (Explicit scope caveat from the
researcher brief.)

### Divergence (the back-leg)

Symmetric: `∇·h(r) = Σ_α ∂h_α/∂r_α`, and `∂h_α/∂r_α (G) = iG_α · h_α(G)`,
so one forward FFT per component + a sum in G-space + one IFFT. Same
code path as the gradient, other direction. **Sign check**:
computing ∇·(∇ρ) of a Gaussian test density should recover −|G|²·ρ(G)
in G-space, i.e., the Laplacian; unit test on Phase A.

### Density floor

PBE numerics are singular as ρ → 0 (k_F → 0, r_s → ∞, F_x(s) → κ+1 with
s → ∞ when |∇ρ|/ρ^{4/3} is finite but ρ is tiny). The LDA floor
(`RHO_FLOOR = 1e-30` e/Å³ in `src/consts.rs`) is fine for PZ but PBE
needs a second check: if `|∇ρ|/ρ > SOMETHING`, clamp s. QE's
convention (see `qe_funct_exch_gga.f90` line ~183, the `small=1.E-10`
check) is:

    if (rho < small) return zeros
    if (sqrt(grho)/rho^{4/3} < very_small) fall back to LDA (s=0 branch)

Mirror this exactly. Phase B includes the tests (`s=0`, `s=10`,
`ρ=1e-12`).

### Simpson vs trapezoidal (reminder)

SIMP already replaced trapezoidal with Simpson for the radial
quadrature (pseudopotential Bessel transforms). PBE does not add any
radial integrals — all our work is grid FFTs. No quadrature concerns.

## API shape

### HYBR compatibility note (EM-applied 2026-04-18)

Phase A's `XcFunctional` dispatch must remain a **data** enum (variants
hold parameters, not closures or trait objects). The hybrid-functional
proposal HYBR (`proposals/HYBR-hybrid-functional-support.md`, landed
2026-04-18) extends the enum with variants that need access to the
wavefunction basis during Hamiltonian construction — if Phase A embeds
the semilocal computation inside a closure, that path isn't reachable
from HYBR's Fock integrator (the closure only sees ρ, never ψ). Keep
`match xc_functional { ... }` dispatch strictly data-driven. `eval` /
`eval_spin` should take `&XcFunctional` parameters, not `Box<dyn ...>`
trait objects, and must not move the functional into `Fn`-typed fields.
See HYBR §3 for the three specific architectural traps to avoid.

### The dispatcher

We follow the **enum-dispatch pattern** (matches `MixingMode`/`Mixer`
from MODR-A and `SmearingScheme`). Single public module method per
call site, internal match on functional type:

```rust
// src/potential/xc.rs (Phase A re-shape)

pub enum XcEvaluator {
    Pz,
    Pbe,
    // future: Pw92, PbeSol, RevPbe, ...
}

impl XcEvaluator {
    pub fn from_settings(s: XcFunctional) -> Result<Self> {
        match s {
            XcFunctional::Pz  => Ok(Self::Pz),
            XcFunctional::Pbe => Ok(Self::Pbe),
            XcFunctional::Pbe0 | XcFunctional::Hse06 =>
                Err(PwdftError::NotImplemented(
                    "hybrid functionals (PBE0, HSE06) not yet supported".into()
                )),
        }
    }

    /// Non-spin: returns (exc_r, v1_r, v2_r) where v2_r is None for LDA.
    pub fn eval(
        &self,
        rho_r: &[f64],
        rho_grad_r: Option<&[[f64; 3]]>,     // only needed for PBE
    ) -> XcGridResult { ... }

    /// Spin-polarized version.
    pub fn eval_spin(
        &self,
        rho_up_r: &[f64],
        rho_down_r: &[f64],
        rho_grad_up_r: Option<&[[f64; 3]]>,
        rho_grad_down_r: Option<&[[f64; 3]]>,
    ) -> XcSpinGridResult { ... }

    /// Whether this functional needs the density gradient.
    /// Drivers use this to skip the ∇ρ FFT for LDA (no regression).
    pub fn needs_gradient(&self) -> bool {
        matches!(self, Self::Pbe)
    }
}

pub struct XcGridResult {
    pub exc_r: Vec<f64>,
    pub v1_r: Vec<f64>,                      // ∂(ρ·ε_xc)/∂ρ
    pub v2_r: Option<Vec<[f64; 3]>>,         // 2·(∂/∂σ) · ∇ρ  (the h vector)
}
```

**Driver call site** (non-spin, phase D):

```rust
let rho_for_xc = add_core_density(&rho_r, &ctx.rho_core_r);
let grad_for_xc = if xc.needs_gradient() {
    let grad = compute_density_gradient(&rho_for_xc, &mut ctx.grid.fft, &ctx.g_vectors);
    Some(grad)
} else {
    None
};
let xc_result = xc.eval(&rho_for_xc, grad_for_xc.as_deref());

// Assemble V_xc(r) = v1_r − ∇·h_r
let vxc_r = assemble_semilocal_vxc(&xc_result.v1_r, xc_result.v2_r.as_deref(),
                                    &mut ctx.grid.fft, &ctx.g_vectors);
```

For LDA, `needs_gradient()` is false, `v2_r` is `None`, and
`assemble_semilocal_vxc` short-circuits to `v1_r` with zero FFT work —
**exact regression safety** for the LDA path.

### Alternatives considered and rejected

- **Two sibling functions** (`lda_xc_grid` + `pbe_xc_grid`, driver
  matches): duplicates the driver dispatch in every call site. The
  CCMX mixer abstraction is our template and it's an enum, not a pair.
  Enum wins.
- **Generics over functional type**: dispatch resolves at compile time;
  would explode the SCF driver into a monomorphized forest. Also
  precludes runtime functional choice (user YAML). Rejected.
- **Trait-object dispatch**: functionally equivalent to enum, adds
  virtual call overhead on every grid point. Enum wins.

## Phasing

| Phase | Scope | Est. time | PR size | Dep | Risk |
|-------|-------|-----------|---------|-----|------|
| **A** | Gradient infrastructure: `compute_density_gradient`, `assemble_semilocal_vxc`, unit tests (Laplacian of Gaussian, linearity). Refactor `lda_xc_grid` call sites to thread through the `XcEvaluator` enum (LDA-only dispatch; no PBE yet). | 1-2 days | ~400 LoC | — | Low — all LDA regressions are test-pinned |
| **B** | PBE exchange (non-spin): port `pbex` from QE with `iflag=1`. Unit tests vs F_x(s) at s ∈ {0, 0.1, 1, 5}. QE cross-check at ρ = 0.1 e/Å³, ∇ρ = 0.05 e/Å⁴. | 1 day | ~200 LoC | A | Low — analytic function, easy to pin |
| **C** | PBE correlation (non-spin) + PW92 LDA helper. Unit tests + Si SCF non-spin regression vs QE PBE. First end-to-end result. | 2-3 days | ~350 LoC | B | **Medium** — PW92 is a separate parametrization; sign/unit checks |
| **D** | Spin-polarized PBE (`nspin=2`). Adds `pbex` spin wrapper + `pbec_spin`. Fe BCC FM SCF regression vs QE PBE. | 2 days | ~300 LoC | C | Medium — correlation uses ∇ρ_total + ζ |
| **E** | GPU PBE shader (`src/gpu/shaders/pbe_xc.wgsl`). Takes `rho_r` and `rho_grad_r` as f32 inputs; returns `v1_r`, `v2_r`. Gradient FFTs run on CPU (wgpu FFT is a separate can of worms, see CUCL). | 2-3 days | ~250 LoC shader + ~200 LoC driver | D | **Medium-high** — f32 precision on F_x(s) at small s; shader branchiness |
| **F** | QE validation suite expansion. Add `tests/qe_validation.rs::test_si_pbe_vs_qe`, `test_al_pbe_vs_qe`, `test_fe_bcc_pbe_vs_qe`. Reference QE runs via `qe-runner` skill. New CSVs under `scripts/validate/gga_pbe_reference.csv`. | 1-2 days | ~200 LoC tests + CSV | D (E optional) | Low — validation-only |

**Total:** ~9-13 days of Core Engineer work (Phases A-D) + 2-3 days
GPU (Phase E) + 1-2 days Researcher validation (Phase F). Single
engineer can do Phases A-D sequentially in ~2 weeks; Phase E parallel
to Phase F once Phase D lands.

**Recommended PR cadence:** one PR per phase, stacked. Phase A stands
alone (pure infrastructure). Phases B+C can merge as one PR if the
reviewer prefers (both are non-spin, small, tightly coupled); I'd split
them because C drags in PW92. Phase D is its own PR. Phase E is
standalone. Phase F is standalone.

## Tests

### Unit tests (Phase B, C, D)

Each added to `src/potential/xc.rs`'s `mod tests`:

1. **F_x(s) pinned values** (exchange only):

   | s   | F_x(s)     | Source |
   |-----|------------|--------|
   | 0   | 1.0        | PBE eq. 14 limit |
   | 0.1 | 1.00273... | Our code == QE `pbex(rho=0.1, grho=sigma, iflag=1)` |
   | 1.0 | 1.2145...  | Same |
   | 5.0 | 1.7951...  | High-s, near LO bound κ=0.804 asymptote |

   Reference values computed from a Python port in
   `scripts/validate/pbe_reference.py` (Phase F).

2. **H(r_s, ζ=0, t) at (ρ = 0.02 e/Å³, |∇ρ| = 0.01 e/Å⁴):** compare
   our H vs QE's `pbec` output. Tolerance 1e-9 eV.

3. **Spin limits:** `pbe_spin(ρ/2, ρ/2, ∇ρ/2, ∇ρ/2) == pbe(ρ, ∇ρ)`
   to machine precision (ζ=0 limit, same shape as our existing
   `test_lsda_unpolarized_limit`).

4. **Zero-gradient limit:** `pbe(ρ, 0)` reduces exactly to PW92
   (not PZ — deliberate). Pin PW92 vs `pw` from QE's `qe_funct_corr_lda.f90`.

5. **Laplacian-of-Gaussian unit test** (Phase A): for
   `ρ(r) = exp(-r²/2σ²)`, check ∇·(∇ρ) matches the analytic
   `(r²/σ⁴ − 3/σ²)·ρ(r)` inside the grid, off-edge to avoid BC artifacts.

### Integration tests (Phase C/D/F)

1. **Si diamond PBE** (non-spin, Phase C): ecut = 30 Ry, 4×4×4
   shifted MP, ONCVPSP PBE Si.upf. Pin total energy vs QE PBE
   reference to ±50 meV (PBE PPs typically need higher ecut than LDA
   for same tolerance; tighten later). Per-component (VGC5-style)
   tolerances: E_kinetic ±10 meV, E_hartree ±5 meV, E_xc ±50 meV
   (most of the slack budget absorbed here), E_local ±10 meV,
   E_nonlocal ±10 meV, E_ewald bit-exact (no PP change).

2. **Al FCC PBE** (non-spin, Phase C): metal with fractional
   occupations. Tests that the smearing/Fermi energy path still works
   with GGA (it should — unaffected).

3. **Fe BCC FM PBE** (spin, Phase D): 8×8×8 MP, currently
   `#[ignore]`'d in our LDA suite. Hope is that PBE converges
   here where LDA fights the mixer (well-known: LDA predicts FCC-NM
   Fe, PBE predicts BCC-FM). Accept convergence at 1e-6 density
   threshold; pin magnetization M = 2.22 μB (PBE reference, vs
   experiment 2.22 μB).

4. **LDA regression sweep**: run *all* existing QEVL LDA tests
   through the new dispatcher, confirm zero change in energies
   (bit-for-bit on Si, sub-meV on Fe).

### GPU consistency (Phase E)

Add `tests/gpu_consistency.rs::test_pbe_xc_gpu_vs_cpu` — same shape
as the existing LDA version. f32 tolerance: 1e-4 eV/electron on
E_xc (vs LDA 1e-5; GGA has one more multiply-with-possible-cancellation
inside F_x(s)).

## Validation (Phase F)

Runs via the `qe-runner` skill. Inputs at `qe_validation/pbe/`:

- `si_pbe.in` — Si diamond, 2 atoms, ecut=30 Ry, 4×4×4 shifted MP,
  `input_dft = 'PBE'`. Reference output: `total_energy`, eigenvalues
  at Γ, X, L, per-component decomposition via `vgc5_per_component.py`.
- `al_pbe.in` — Al FCC, 1 atom, ecut=30 Ry, 8×8×8 shifted MP,
  Gaussian smearing σ=0.01 Ry.
- `fe_pbe.in` — Fe BCC, 1 atom, ecut=50 Ry (Fe wants high ecut),
  8×8×8 shifted MP, `nspin=2`, `starting_magnetization(1)=0.7`.

Reference CSV at `scripts/validate/gga_pbe_reference.csv`, same
schema as `vgc5_qe_si_components.csv`. Python reference port at
`scripts/validate/pbe_reference.py` (ports QE's `pbex`/`pbec` to
Python + NumPy for per-point analytic cross-checks).

## GPU path (Phase E)

CPU PBE lands in Phases A-D. Phase E is optional — CPU PBE is not
gated on it.

New shader: `src/gpu/shaders/pbe_xc.wgsl`. Takes `rho_r`,
`rho_grad_r_x`, `rho_grad_r_y`, `rho_grad_r_z` (four f32 arrays) and
returns `exc_r`, `v1_r`, `v2h_x`, `v2h_y`, `v2h_z` (five f32 arrays —
the assembled `h = 2 v2 ∇ρ` in the output, so the divergence step
runs in G-space on CPU / second GPU pass).

Gradient FFTs: stay on CPU for Phase E (wgpu FFT is a separate
decision, see CUCL proposal). The ∇ρ → h → ∇·h pipeline crosses
CPU↔GPU three times per SCF iter in the mixed CPU-FFT / GPU-xc path.
On M2 8-core, expected wall-time gain over CPU: 1.2-1.5× on the xc
step alone (smaller than LDA's 2-3× because the xc step is no longer
the per-grid-point bottleneck — the FFTs are). Worth it on large
grids (96³+) only. Flag for benchmark-driven re-evaluation.

## Pseudopotentials

### Functional-tag consistency check

**Current state (audited):** `src/pseudopotential/upf/convert.rs` parses
the UPF header but **does not read `functional=` at all.** The field
is silently dropped. A PBE UPF and an LDA UPF for the same element
go into the solver indistinguishably.

**Proposed fix (Phase A):** add `functional_tag: String` field to
`PseudopotentialData`, populated from the UPF `functional=` attribute.
At SCF load time (`ScfContext::new`), check that the user's
`xc.functional` matches every PP's tag. If not, return
`PwdftError::FunctionalMismatch { pp_element, pp_tag, requested }`.

The matching table:

| YAML `functional:` | Accepted UPF `functional=` tags |
|--------------------|-----------------------------------|
| `pz`               | `SLA  PW   NOGX NOGC`, `PZ`, `LDA-PZ` |
| `pbe`              | `PBE`, `GGA-PBE`, `SLA  PW   PBX  PBC` |
| `pbe0`, `hse06`    | (error: not yet supported) |

The 4-token QE form (`SLA PW NOGX NOGC`) is QE's internal shorthand
for Slater exchange + PW92 correlation + no gradient exchange + no
gradient correlation. For PBE it's `SLA PW PBX PBC` = Slater +
PW92 + PBE-exchange-gradient + PBE-correlation-gradient. See
`qe-7.5/Modules/funct.f90` for the full decoder ring.

Accept both the compact form (`PBE`) that ONCVPSP generates and the
QE 4-token form (`SLA PW PBX PBC`) that older GBRV/SSSP PPs use.

### PP files

Already present in the tree at `pseudopotentials/nc/pbe/` (71
elements, checked Phase F's three test systems — Si, Al, Fe — are all
there). No external downloads needed. **This is a pleasant surprise
vs. the researcher brief, which flagged "test suite needs new PP files
added".**

## Build vs buy — libxc

**libxc** (https://tddft.org/programs/libxc/) is the C library that
QE, VASP, ABINIT, and most every other DFT code use for functionals.
~600 functionals, single integration. FFI from Rust is mature
(libxc-sys on crates.io).

**Cost of libxc integration:**
- Native dependency (C library) — conflicts with "zero system deps"
  goal in CLAUDE.md.
- FFI boundary for every grid point (or chunk). Inversion of control
  pattern (libxc drives, we provide callbacks).
- Builds must link libxc on macOS; wgpu/Metal adds Apple-specific
  linker pain.
- Licensing: libxc is MPL-2.0. Compatible with our code but a new
  third-party license to manage.

**Cost of from-scratch:**
- PBE implementation is ~400 LoC total (pbex + pbec + pbec_spin).
  QE port is mechanical transcription. We have the ground-truth source
  on disk.
- Each new GGA variant (PBEsol, RPBE, BLYP, ...) is 50-150 LoC each.
- No new build-time dependencies; pure-Rust stack preserved.

**Recommendation: from-scratch for PBE-only scope.** Revisit libxc
when we want ≥3 GGA functionals or our first meta-GGA (SCAN). The
crossover is: "when the effort to port the next functional exceeds
the effort to integrate libxc". PBE alone is well below that bar.

**If libxc becomes attractive later:** it's a drop-in replacement for
the enum-dispatch layer. `XcEvaluator::Pbe` today, `XcEvaluator::Libxc(id)`
in the future. No call-site changes.

## Risks and mitigations

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| Sign error on ∇·h term (the semilocal potential footgun) | E_total right at convergence, but O(Δρ) error spoils quadratic |E_HF−E_KS| convergence; false "density converged, energy not" warnings | **High** — this is *the* historical GGA pitfall | Phase B test: pin V_xc ≡ ∂(ρ·ε)/∂ρ − ∇·h against QE's output at 3 test points before any SCF run. Phase C: Si |E_HF−E_KS| < 1e-4 eV at convergence, same tolerance as LDA |
| PZ-vs-PW92 mix-up inside PBE correlation | E_total wrong by ~1-10 meV/electron; all our validation thresholds blown | **Medium** — easy to miss if Phase C reuses `perdew_zunger_correlation` | Phase C adds PW92 as a *new* helper; unit test pins PW92 against QE's `pw` subroutine output directly (not via PBE). PBE path calls PW92 only |
| NLCC gradient integration mistake (forget ∇ρ_core, or double-count) | E_xc wrong on every Si/Fe/any-PP-with-NLCC run; Si-PBE validation fails | Medium | Phase A unit test: compute ∇(ρ + ρ_core) two ways (direct vs. ∇ρ + ∇ρ_core), assert equal. Phase C Si test fails hard if this regresses (NLCC contributes ~0.5 eV to Si E_xc) |
| Spin correlation uses ∇ρ_σ instead of ∇ρ_total (common misread of PBE paper) | Fe BCC E_xc wrong by ~10-100 meV | Medium | Phase D test: unpolarized limit (ζ=0) must match Phase C non-spin to 1e-10 eV — any slip from ∇ρ_σ vs ∇ρ_total shows up as a ζ=0 discrepancy because the two forms are only identical at ζ=0 |
| f32 precision on GPU PBE at small-s (Phase E) | GPU vs CPU disagreement > 1e-4 eV/electron; forces us to drop back to CPU | Low-medium | Phase E: use f64 on the GPU path for the `1/(1 + μs²/κ)` denominator (the only non-benign operation). f32 elsewhere. Still net f32 throughput since the expensive ops (cbrt, log, exp) remain f32 |

**Highest-risk item: the ∇·h sign.** It has the same shape as every
GGA bug in every textbook appendix ever written. Phase B's very first
test pins V_xc at three grid points with known (ρ, ∇ρ) against QE —
this catches the sign on day one, before any SCF runs.

## Out-of-scope follow-ups

Capture these as FLUP-family proposals:

- **PBEsol / revPBE / RPBE** (Phase G, FLUP): `iflag=2,3,8` in QE.
  ~30 LoC each on top of Phase C's `pbex`/`pbec` scaffolding.
- **PW91** (FLUP): QE's `ggax`/`ggac` — same structure, different
  enhancement factor.
- **libxc integration** (separate proposal, probably `LIBX`): revisit
  when we want ≥3 GGAs or our first meta-GGA.
- **SCAN / r²SCAN meta-GGA** (separate proposal): requires τ(r) —
  architectural work.
- **MADOC extension** (flag in MADOC's scope file): PBE's own
  mathematical-documentation needs — docstrings for `pbex`, `pbec`,
  `pbec_spin`, `semilocal_vxc_assembly`. Do not modify MADOC in this
  proposal; MADOC owner folds in when PBE lands.

## Acceptance criteria

- Phase A: LDA regression tests pass (bit-exact); Laplacian-of-Gaussian
  unit test passes; `XcFunctional::Pbe` returns `NotImplemented` from
  `from_settings` (placeholder).
- Phase B: F_x(s) unit tests pass at all four s values; three-point
  V_xc pin against QE passes.
- Phase C: Si PBE non-spin SCF converges; |E_total − QE| < 50 meV
  (likely < 5 meV based on LDA comparison — the 50 meV is slack);
  |E_HF − E_KS| < 1e-4 eV at convergence.
- Phase D: Fe BCC PBE FM SCF converges; M = 2.22 μB ± 0.05;
  |E_total − QE| < 100 meV.
- Phase E: GPU PBE matches CPU PBE within 1e-4 eV/electron on E_xc
  over the existing GPU consistency test suite.
- Phase F: three QE reference calcs cached at `scripts/validate/`;
  Python reference port matches QE to machine precision.

## References

- Perdew, Burke, Ernzerhof, "Generalized Gradient Approximation Made
  Simple", *Phys. Rev. Lett.* **77**, 3865-3868 (1996). Erratum
  *PRL* **78**, 1396 (1997).
- Perdew, Wang, "Accurate and simple analytic representation of the
  electron-gas correlation energy", *Phys. Rev. B* **45**, 13244
  (1992). (PW92 LDA correlation used inside PBE.)
- Oliver, Perdew, "Spin-density gradient expansion for the kinetic
  energy", *Phys. Rev. A* **20**, 397 (1979). (Spin-scaling for
  exchange.)
- Martin, *Electronic Structure: Basic Theory and Practical Methods*
  (Cambridge, 2004), Ch. 8.3 "Generalized Gradient Approximation".
- Kresse, Furthmüller, *Phys. Rev. B* **54**, 11169 (1996). §III.C
  on GGA implementation in plane-wave codes — references for the
  FFT-based ∇ρ / ∇·h approach.
- QE 7.5 source (canonical convention):
  - `qe-7.5/XClib/qe_funct_exch_gga.f90` (`pbex`, lines 111-331)
  - `qe-7.5/XClib/qe_funct_corr_gga.f90` (`pbec`, `pbec_spin`)
  - `qe-7.5/XClib/qe_drivers_gga.f90` (`gcxc`, `gcxc_spin` — the ∇·h
    assembly path)
  - `qe-7.5/Modules/funct.f90` (functional-tag decoder ring)
