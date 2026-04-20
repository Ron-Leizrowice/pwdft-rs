---
id: HYBR
title: Hybrid exchange-correlation functionals (PBE0, HSE06) — scoping proposal
priority: medium
complexity: large
risk: high
depends_on: [GGAP]
blocks: []
status: draft
author: Researcher
date: 2026-04-18
---

## HYBR — Hybrid functional support (PBE0, HSE06)

### Problem

Hybrid functionals systematically close the DFT band-gap problem.
For Si, LDA ~0.5 eV, PBE ~0.6 eV, experiment 1.17 eV, HSE06 ~1.17 eV.
Every serious materials application (semiconductor band alignments,
defect levels, surface workfunctions, molecular HOMO-LUMO) wants
hybrids. Once GGAP lands (PBE infrastructure), PBE0 and HSE06 are
the obvious next physics proposal.

XCNI already parses `XcFunctional::{Pbe0, Hse06}` from YAML and returns
`PwdftError::NotImplemented` from the driver. HYBR scopes turning
those two enum variants into working SCF operators.

This is scoping only. No code. No GGAP edits. The EM decides
activation after GGAP Phase A lands.

### Why hybrids are architecturally different

LDA and GGA are *density* functionals: `V_xc(r) = f(ρ(r), ∇ρ(r))`. The
XC contribution to H is a grid quantity folded into V_eff once per
SCF iter. Dispatch is "ρ-in → V_xc-out."

**Hybrid functionals are *wavefunction* functionals.** The Fock term

    V_x^Fock ψ_i(r) = − Σ_j f_j ψ_j(r) ∫ ψ_j*(r') ψ_i(r') w(|r−r'|) dr'  (1)

with `w(r) = 1/r` (PBE0) or `w(r) = erfc(ωr)/r` (HSE06, ω=0.106 Bohr⁻¹)
acts *on* ψ and cannot be precomputed as a grid potential. The
Hamiltonian becomes

    H |ψ⟩ = T |ψ⟩ + V_eff |ψ⟩ + V_NL |ψ⟩ + V_x^Fock |ψ⟩                  (2)

with V_eff still on the grid, V_x^Fock applied band-by-band. Naive
cost O(n_occ² · n_pw² · n_kpts²), ~100× the rest of H combined.
Production codes use **ACE** (Lin Lin, *J. Chem. Theory Comput.* **12**,
2242 (2016)) to reduce wall-time by O(100). QE's implementation is
`qe-7.5/PW/src/exx.f90` (~5000 lines; `use_ace` flag at line 63).

### 1. Physics scope

**PBE0** (Perdew-Burke-Ernzerhof, *J. Chem. Phys.* **105**, 9982 (1996)):

    E_xc^{PBE0} = 0.25 E_x^{Fock} + 0.75 E_x^{PBE} + E_c^{PBE}            (3)

Unscreened Coulomb for Fock. Adamo-Barone (*JCP* **110**, 6158 (1999))
for the 0.25 mixing justification.

**HSE06** (Heyd-Scuseria-Ernzerhof, *JCP* **118**, 8207 (2003); erratum
*JCP* **124**, 219906 (2006)):

    E_xc^{HSE06} = 0.25 E_x^{Fock,SR}(ω) + 0.75 E_x^{PBE,SR}(ω)
                 + E_x^{PBE,LR}(ω) + E_c^{PBE}                            (4)

SR/LR split via `1/r = erfc(ωr)/r + erf(ωr)/r`. The screening
parameter ω = 0.106 Bohr⁻¹ is the 2006-erratum value (the 2003 paper
had 0.15). QE's G-space kernel (`exx_base.f90:809`) is

    w_SR(G) = (4π/G²) · [1 − exp(−G²/(4ω²))]     (Ry, G in Bohr⁻¹)       (5)

with finite G→0 limit `π/ω²` — **screening kills the divergence**,
which is why HSE06 avoids Gygi-Baldereschi / Martyna-Tuckerman
machinery and ends up ~O(10) cheaper than PBE0 overall.

Both functionals add `V_x^Fock` to H; they differ only in `w(G)`. The
architecture must treat PBE0 as "HSE06 with ω=0" (unscreened) — one
`ExxOperator` with a kernel parameter, not two.

### 2. Architectural decision — extending SCF dispatch

GGAP Phase A threads `params.xc_functional` through the driver and
dispatches semilocal XC (LDA | PBE). HYBR needs a second axis: "does
this functional add a wavefunction-level Fock term?"

**Option A — nested enum (recommended):**

```rust
enum XcSpec {
    Pz, Pbe, Hybrid(HybridSpec),
}
struct HybridSpec {
    fraction: f64,               // 0.25 for PBE0/HSE06
    base: GgaVariant,            // PBE for both
    screening: Option<f64>,      // Bohr⁻¹: None=PBE0, Some(0.106)=HSE06
}
```

Driver matches `xc_functional`; on `Hybrid(h)` it (1) computes V_eff
from the scaled semilocal part, (2) applies `ExxOperator` to each ψ
before eigensolving. The operator owns the kernel `w(G)` and current
occupied wavefunctions; exposes `apply(&self, psi: &ArrayView2<_>) ->
Array2<_>`.

**Option B — flat enum + capability methods.** Keep `XcFunctional::{Pz,
Pbe, Pbe0, Hse06}`, add `needs_fock_exchange()` + `fock_kernel()`. Less
type-safe (invalid "Pbe0 with ω>0" representable) but zero churn to
GGAP Phase A.

**Recommendation: Option A, conditional on GGAP Phase A's shape.** The
nested enum makes HSE06 = PBE0 + screening explicit and forbids
ill-formed states. Fall back to B if A's enum resists extension.

The crucial requirement either way: GGAP's enum stays a *data* enum
(match-dispatched), not a trait object or closure. See §3.

### 3. GGAP amendment (propose; do not apply here)

Amendment text for the EM to add to GGAP post-approval:

> **GGAP → HYBR compatibility note.** Phase A's `XcFunctional` dispatch
> must remain a *data* enum (variants hold parameters, not closures or
> trait objects). HYBR extends the enum with variants that need access
> to the wavefunction basis during Hamiltonian construction — if Phase
> A embeds the semilocal computation inside a closure, that path isn't
> reachable from HYBR's Fock integrator (the closure only sees ρ). Keep
> `match xc_functional { ... }` dispatch strictly data-driven. `eval`
> / `eval_spin` should take `&XcFunctional` parameters, not Box or
> move it into Fn-typed fields.

Do not edit GGAP in this PR — advisory text only.

### 4. Phasing (rough; hybrids are large)

| Phase | Scope | CE-days | Risk |
|-------|-------|---------|------|
| **0** | GGAP Phase A shape check; confirm hybrid-extensible. | 0.5 | Low |
| **1** | `ExxOperator` infrastructure. CPU, single k-point, non-spin, no ACE. Bare O(n_occ² · n_pw²) G-space convolution. Si 2-atom HSE06 end-to-end. | 7-10 | High — G=0, periodic-image cutoff, gamma extrapolation pitfalls |
| **2** | Range-separated kernel `w(G)` (eq. 5) unit-tested vs QE `exx_base.f90:803-824`. Adds PBE0 with Gygi-Baldereschi G=0 extrapolation. | 3-4 | Medium — G=0 special case is a footgun |
| **3** | **ACE** (Lin Lin 2016). Rebuild every N_ace iters; apply cheaply between. ~10-100× speedup. | 10-15 | **Highest** — algorithmic density; ratio experiments against QE `use_ace=.TRUE.` |
| **4** | Multiple k-points. O(n_kpts²) scaling; natural distribution over (k,q) pairs. Γ-only first (phases 1-3); MP grids here. | 7-10 | High — q-sum convergence is slow |
| **5** | Spin-polarized (`nspin=2`). Each channel gets its own Fock operator. Fe BCC HSE06 target. | 3-5 | Low-medium |
| **6** | Validation. Si HSE06 gap (1.17 ± 0.05 eV expected); C diamond (~5.4 eV); Al PBE0 total. QE refs via `qe-runner`. | 3-5 | Low |

**Total: ~35-55 CE-days (≈7-11 CE-weeks).** Phase 3 (ACE) is both the
dominant line item and the dominant risk.

### 5. Anti-scope

Each requires its own proposal:

- **libxc.** Same reasoning as GGAP: pure-Rust. Port QE directly.
- **Double-hybrids** (B2PLYP) — add MP2 correlation machinery.
- **Meta-GGAs** (SCAN, TPSS, r²SCAN) — need τ(r); different axis.
- **RPA / GW / BSE** — different theoretical framework.
- **GPU Fock exchange.** CPU-only all phases. Separate proposal.
- **USPP / PAW hybrids.** Norm-conserving only. See `us_exx.f90`,
  `paw_exx.f90`.

### 6. Open questions

1. **Is ACE the right compression?** Alternatives: linear-scaling
   exchange (Wu/Selloni/Car 2009), tablewise EXX. ACE is the
   QE/VASP/CP2K default. Recommend ACE; revisit if memory blocks.
2. **Memory budget.** Bare Fock at n_pw=500, n_occ=4, n_kpts=8, no ACE
   ~O(n_occ² · n_pw · n_kpts²) Complex64 — GB-scale intermediates.
   Production targets (n_pw~5000, n_occ~40, n_kpts~64) approach ~100 GB
   without ACE. ACE compresses to n_occ Complex64 vectors length n_pw
   per k (~10-100 MB). **ACE is not optional for production.**
3. **k-point EXX parallelization.** Natural is rayon over (k,q)
   pairs; duplicates wavefunctions. MPI-style decomp is the long-term
   answer; single-node rayon fine for Phase 4.
4. **Gamma extrapolation for PBE0.** Unscreened Coulomb diverges at
   G=0 in periodic systems; QE uses Gygi-Baldereschi (1986) or
   Martyna-Tuckerman (1999) (`exx_base.f90:820-830`). Pick one. HSE06
   dodges by screening.

### Risks

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| GGAP Phase A paints enum into closure/trait corner | Re-open Phase A; 2-week delay | Medium | §3 amendment; EM applies on GGAP merge |
| ACE bugs produce silently-wrong totals | Hybrid results untrustworthy; gap off 0.6-1.5 eV | Medium-high | Phase 3 bit-compare QE `use_ace=.T.`; cross-check Phase 1 bare-Fock |
| G=0 divergence (PBE0) mishandled | Total energy off ~eV | Medium | Ship HSE06 first (Phase 1, no divergence); PBE0 in Phase 2 with explicit ω→0 limit test |
| k-point EXX memory explodes at n_kpts>32 | Production blocked | High (ACE mitigates) | Phase 4 benchmark-gated; Γ-only flag for applicable systems |

**Highest-risk technical item: ACE integration at production memory
scale.** Phase 3 gates whether HYBR becomes useful or stays a demo.
k-point EXX parallelization (Phase 4) is a close second.

### Acceptance

- **0:** GGAP Phase A landed, data-dispatched.
- **1:** Si 2-atom HSE06 SCF converges; |E − QE (use_ace=.F.)| < 50 meV.
- **2:** `w(G)` matches QE at three G-points to 1e-12 Ry.
- **3:** ACE ≥ 10× faster than Phase 1 bare Fock; |ΔE| < 1 meV.
- **4:** Si 4×4×4 HSE06 gap within 0.05 eV of QE.
- **5:** Fe BCC HSE06 FM converges; M within 0.05 μB.
- **6:** Si HSE06 1.17 ± 0.05 eV; C diamond 5.40 ± 0.1 eV.

### Surprises worth noting

- **HSE06 is architecturally simpler than PBE0.** Screening kills the
  G=0 divergence (eq. 5 has finite limit π/ω²); no Gygi-Baldereschi /
  Martyna-Tuckerman needed. Cost is also lower. **Implement HSE06
  first** — the XCNI enum ordering (PBE0 before HSE06) is wrong from
  an implementation sequencing standpoint.
- **QE's `use_ace=.FALSE.`** path is a genuine bare-Fock reference for
  Phase 1 validation. Slow (minutes/iter) but exists and is correct.
  No need to invent a Phase 1 reference.

### References

- Perdew, Burke, Ernzerhof, *J. Chem. Phys.* **105**, 9982 (1996) —
  PBE0.
- Adamo, Barone, *JCP* **110**, 6158 (1999) — PBE0 mixing parameter.
- Heyd, Scuseria, Ernzerhof, *JCP* **118**, 8207 (2003); erratum *JCP*
  **124**, 219906 (2006) — HSE06.
- Lin Lin, *J. Chem. Theory Comput.* **12**, 2242 (2016) — ACE.
- Paier et al., *JCP* **124**, 154709 (2006) — HSE for solids,
  k-convergence.
- Martin, *Electronic Structure* (Cambridge 2004) Ch. 7.6, 8.4.
- QE 7.5 source:
  - `qe-7.5/PW/src/exx.f90` (main driver; `use_ace` at line 63)
  - `qe-7.5/PW/src/exx_base.f90` (`fac(ig)` kernel at 803-824; G=0
    divergence at 900-953)
  - `qe-7.5/PW/src/exx_band.f90` (band-parallel operator application)
