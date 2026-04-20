---
id: VNLM-CUD
title: VNLM-CUD — Cu Γ-point per-m pin for two-radial-d-projector sum
status: active
priority: high
complexity: small
risk: low
depends_on: [VGCH-2F, VNMT]
blocks: [VQEF]
owner: core-engineer
author: Researcher (2026-04-20)
---

## VNLM-CUD — Cu two-radial-d-projector sum diagnostic

### Problem

VGCH-2F Part C session-2 (PR #179,
`proposals/VGCH-2F-part-c-session-2-findings.md`) narrowed the Cu
shared-density Kleinman-Bylander hypothesis **H-C5** (`D_ij · Σ_lm β_i · β_j`
contraction on Cu d-projectors) to a single structural question:

> **`PP_DIJ` is strictly diagonal 6×6 on Cu (same as Si), but Cu has
> two l=2 radial projectors summed into the same angular channel.
> Does `D_00·β_{l=2,ν=0}·β_{l=2,ν=0} + D_11·β_{l=2,ν=1}·β_{l=2,ν=1}`
> assemble correctly when both radial channels are active?**

Verified facts (from VGCH-2F § H-C5 evidence):

- Cu.upf (`pseudopotentials/nc/lda/Cu.upf`) carries 6 β projectors
  in the layout `(l=0,ν=0), (l=0,ν=1), (l=1,ν=0), (l=1,ν=1),
  (l=2,ν=0), (l=2,ν=1)` (same (l, ν) shape as Si ONCVPSP).
- Cu.upf `PP_DIJ` is strictly diagonal: every off-diagonal entry
  is zero (`grep -cE 'angular_momentum="2"' Cu.upf` → 2, and the
  PP_DIJ block has six nonzero diagonals with 30 zero off-diagonals).
- Si and Cu therefore have identical D_ij structural layouts. The
  difference is physical: Cu's two l=2 radial projectors are both
  materially significant (d-electron valence manifold) while Si's
  l=2 projectors are smaller-magnitude energy-gap corrections.

The VNMT test `test_single_channel_l2_m_isolation` (`potential/nonlocal.rs:914`)
pins Si's **single** l=2 radial projector against an explicit per-m
hand-computed reference, using the addition-theorem identity on Y_{2,m}.
It explicitly zeros five of six D_ij diagonal entries to isolate one
radial channel. This **correctly pins the angular surface for any
single-radial-ν l=2 projector**, but by construction does not exercise
the sum-over-ν that Cu requires.

VGCH-2F identified this as a **trace-equivalent-but-projector-wrong**
risk per the Researcher's agent-prompt warning: a per-ν asymmetry in
the radial magnitude that cancels correctly in the sum Σ_ν over m
(via addition theorem) would pass VNMT's m-isolation test while still
producing an error on individual m-channels. Cu is the canonical
witness for this class of bug: its d-manifold sits on the Fermi
surface and any per-m asymmetry in V_NL matrix elements would
manifest as the +0.26 eV/band Γ-point eigenvalue offset observed
in the Cu transplant (PR #167, `Δε/band = +0.26 eV` across all 5
d-states at Γ).

**Magnitude plausibility.** The Cu shared-density gap of +8.87 eV
in E_xc + partial cancellation yielding +16.34 eV in E_HF (VGCH-2F
Note 1) is an order of magnitude larger than a pure angular-surface
bug (≤ keV atomic-scale shifts). A two-radial sum bug on the d-
channel that flipped the relative weight of the two projectors
would shift the V_NL matrix element by ~D_{00} (in eV) — for Cu,
`D_{00,l=2} ≈ O(100–300 eV)`, so a 10% sign-weight flip on one ν
gives O(10 eV) matrix-element shift, consistent with the observed
gap magnitude.

### Research

#### Why this is physics, not just testing

Implementations of `Σ_i D_ii |β_i⟩⟨β_i|` with multiple radial
projectors per l channel (ν > 1) are rare in semiconductor PPs —
Si ONCVPSP has them but at low amplitude. On d-electron metals
with multi-radial PPs they are load-bearing. Cu is the canonical
case and the Cu transplant fingerprint has exactly this shape.

If a bug exists in pwdft-rs's `NonlocalPotential::new`
(`potential/nonlocal.rs:195–XXX`) — e.g. a missed inner-loop index
when n_proj_l > 1, or a missed ν-index in the pre-computed
`D·B^H` contraction at `nonlocal.rs:364-365` — the bug would
manifest only on d-electron metals with multi-radial PPs. Si
would silently pass (small magnitude), Fe LDA would silently
pass (M=0 and valence-at-d limited), but Cu would display
exactly the +8.87 eV E_xc + one-e opposite-sign fingerprint
observed. **This is not idle hypothesis-testing; it's the
last structural feature unique to Cu that VGCH-2 has not pinned.**

#### Current V_NL assembly — what VNLM-CUD needs to cover

Reading `potential/nonlocal.rs:195–400`:

- Per-k-point, the constructor builds `B[G, channel]` where
  `channel = (atom α, radial projector ν, m)`. `B[G, ch]` =
  `(1/√Ω) · exp(−iG·τ_α) · F_ν(|k+G|) · Y_{l_ν,m}(q̂_{k+G})`.

- `D·B^H` is pre-applied in the constructor (nonlocal.rs:364–365
  loop). The comment at nonlocal.rs:364 says:
  `DB_H[channel(a,i,m), :] += (D_ij / Ω) · conj(B[:, channel(a,j,m)])`
  which implicitly indexes per (a, i, m) and (a, j, m) over i,j
  ∈ [0, n_proj_type). For Cu this is a 6×6 sum per atom per m,
  and because D is diagonal reduces to 6 diagonal terms per atom
  per m — but the 2 terms at l=2 both contribute to the **same m
  slot** (since `channel(a, i=4, m) ≠ channel(a, i=5, m)` but both
  get the same Y_{2,m} angular factor in B).

- `H_NL = B · (D·B^H)` in `add_to_hamiltonian` (nonlocal.rs:412–452)
  via a single GEMM. This is VNLM (`proposals/completed/VNLM-*.md`),
  the 5.2× speedup over the scalar form; it has been regression-
  pinned on Si but Cu's two-ν case is not explicitly pinned on
  the sum-over-ν axis.

The primitive identity the test needs to pin on Cu's d-channel at
Γ, for a single atom at τ = 0:

```text
⟨G₁ | V_NL | G₂ ⟩
  = Σ_ν Σ_ν' D_{νν'} · (1/Ω) · F_ν(|G₁|) F_ν'(|G₂|)
    · Σ_m Y_{2,m}(Ĝ₁) Y_{2,m}(Ĝ₂)                           [diagonal D: ν=ν']
  = (1/Ω) · (Σ_m Y_{2,m}(Ĝ₁) Y_{2,m}(Ĝ₂))
    · (D_{44} · F_4(|G₁|) · F_4(|G₂|) + D_{55} · F_5(|G₁|) · F_5(|G₂|))
```

where `F_4, F_5` are the Bessel transforms of Cu's two l=2 radial
projectors (`pp.beta_projectors[4]`, `pp.beta_projectors[5]`).
The angular sum is the addition theorem `(2l+1)/(4π) · P_2(Ĝ₁·Ĝ₂)`.

**The key observation:** VNMT's test zeros `D[55]` (leaving only
one radial projector active). The structural question Cu raises is
what happens when both are active — i.e. does the Rust assembly
return `D_44·F_4·F_4 + D_55·F_5·F_5` or something else?

### Implementation

One phase, one new test, one manual-break sanity check.

#### Phase 1 — add `test_two_radial_l2_cu_per_m` (1 CE-day)

Land a new test in `pwdft/pwdft-core/src/potential/nonlocal.rs::tests`
as a sibling of `test_single_channel_l2_m_isolation` (line 914):

```rust
/// VNLM-CUD — Cu l=2 two-radial-ν sum per-m pin.
///
/// Complements test_single_channel_l2_m_isolation (which zeroes
/// D[55] to isolate a single radial projector). Here we pin the
/// full-D_ij assembly on Cu, where the l=2 channel has two
/// radial projectors (ν = 0, 1) and both carry the same Y_{2,m}
/// angular factor. A per-ν magnitude error (e.g. a missed sum over
/// ν in `D·B^H`, or a shifted array index at the (a, i=5, m)
/// channel layout) would cancel in the VNMT single-channel test
/// via Y_{2,m} addition-theorem but be detectable here.
///
/// Setup identical to VNMT except:
///   - Load Cu.upf (not Si.upf).
///   - Zero D except the two l=2 diagonals D[4,4] and D[5,5].
///   - Per-m reference is the sum of two radial form factors:
///     F(ν=4; |G|) · F(ν=4; |G₂|) · D_{44}
///     + F(ν=5; |G|) · F(ν=5; |G₂|) · D_{55}
///     (both multiplied by the VNMT angular product for each m).
///   - Pin H_NL[G₁, G₂] per m against the explicit two-ν sum.
#[test]
fn test_cu_l2_two_radial_per_m() {
    // ... see VNMT's structure above ...
    // After loading Cu.upf:
    //   assert_eq!(pp.n_projectors(), 6);
    //   assert_eq!(pp.beta_projectors[4].l, 2);
    //   assert_eq!(pp.beta_projectors[5].l, 2);

    // Zero off-(l=2) D elements; keep D[4,4] and D[5,5].
    let d_44 = pp.dij[4 * 6 + 4];
    let d_55 = pp.dij[5 * 6 + 5];
    assert!(d_44.abs() > 1e-6);
    assert!(d_55.abs() > 1e-6);
    pp.dij.iter_mut().for_each(|d| *d = 0.0);
    pp.dij[4 * 6 + 4] = d_44;
    pp.dij[5 * 6 + 5] = d_55;

    // Use the same G_1 = (1,0,1), G_2 = (2,0,1) pair VNMT pins, so
    // the addition-theorem angular sum matches VNMT's 17/(16π).

    // Expected H_NL[G_1, G_2]:
    //   = (17 / (16π)) · (1/Ω) · (D_44·F_4(|G_1|)·F_4(|G_2|)
    //                              + D_55·F_5(|G_1|)·F_5(|G_2|))
    //
    // F_ν computed via bessel_transform_projector (same production
    // helper VNMT uses).

    // Pin at 1e-10 (same ULP headroom as VNMT).
}
```

**Geometry choice — why Cu FCC instead of the VNMT artificial
cubic.** VNMT uses a = 2π Å (cubic, Miller → Cartesian identity)
to make the hand-computed angular sum clean. For Cu, preserving
that a = 2π Å simple-cubic artificial-geometry is the cleanest
path: it lets us directly reuse VNMT's angular-sum constants
(17/(16π) at G₁ = (1,0,1), G₂ = (2,0,1)), and the only new pin
is the two-ν radial sum. Single Cu atom at origin; geometry is
diagnostic-only.

**Why not load Cu's real FCC geometry.** Cu FCC a = 6.8219 Bohr
= 3.61 Å; reciprocal lattice is FCC with b_i ≠ ê_i in Cartesian.
Miller indices no longer coincide with Cartesian G — any per-m
angular sum has to be recomputed in rotated coordinates. This
adds no value to the two-ν radial sum check and doubles the
hand-computation surface area.

#### Phase 2 — manual-break sanity check (~15 min)

As in VNMT, verify the test is sensitive to a specific bug class.
Two mutations, each must trigger the test at >1e-5 eV residual:

**Mutation A: shift Cu's second radial projector.**
In the test (not production), before calling `NonlocalPotential::new`,
scale `pp.beta_projectors[5].r_beta` by √2. If the test
passes under this mutation, it is not detecting per-ν amplitude errors.

**Mutation B: flip sign on `D_55`.**
Replace `pp.dij[5*6+5]` with `-d_55`. The expected H_NL value also
flips sign on the second-ν term; the test must detect the sign
flip.

If mutation A triggers but B does not, the test is too sensitive
to F₅'s magnitude and not sensitive to the D₅₅ sign — file a
tolerance adjustment and re-pin. If either mutation does not
trigger, the test is defective; abandon VNLM-CUD Phase 1 and
escalate.

Mutations are test-only and are removed before landing. Document
the two mutation invariants in the test docstring as a **"must be
detected by"** contract — future edits that weaken the test break
the contract.

#### Phase 3 — VQEF Cu cell re-check (after Phase 1 passes)

Before PR submission, re-run VQEF Cu LDA end-to-end under `/test
--tier2` and include pre/post numbers (E_total, shared-density
transplant ΔE_xc at iter-1) in the PR body. No behavior change
expected — Phase 1 is a pure regression pin, not a fix.

**Outcome gate (physics interpretation):**

- **Phase 1 test passes + mutations trigger + VQEF Cu numbers
  unchanged.** Cu's two-radial-d assembly is correct. H-C5 is
  fully refuted. Escalate to H-C2 (`n_bands` margin audit on Cu's
  3d DOS tail) as the next Class A Cu candidate. Close VNLM-CUD.

- **Phase 1 test fails from the start (before any mutation) on
  vanilla Cu.** We have a localized KB assembly bug. File a
  VNLM-CUD-FIX follow-up proposal with the exact per-m residual
  table as the diagnosis surface. The Cu +8.87 eV ΔE_xc gap is
  now explained; expect the fix to close Cu in VQEF.

- **Phase 1 test passes but mutations don't trigger.** The test
  is insensitive; its contract is violated. Iterate on the per-m
  reference construction until the sanity checks bite. Do not
  treat this as evidence for "Cu two-radial-d is correct" — the
  evidence is conditional on the mutations firing.

### Verification

- `test_cu_l2_two_radial_per_m` passes at 1e-10 tolerance (Tier-1).
- Mutation A (β₅ × √2) causes the test to fail at >1e-5 eV
  residual. Documented in the PR as a preserved-as-comment log.
- Mutation B (D₅₅ sign flip) causes the test to fail at >1e-5 eV
  residual. Documented in the PR as a preserved-as-comment log.
- VNMT's existing Si test (`test_single_channel_l2_m_isolation`)
  still passes untouched.
- `/test --tier2` on the qe_validation suite: VQEF Cu LDA pre/post
  numbers are identical (pin is a test, not a fix). If they
  differ, the pin's geometry seeded a different production code
  path — regression.

### Non-goals

- **Does NOT fix any Cu residual.** Strict diagnostic.
- **Does NOT touch `NonlocalPotential::new` or
  `add_to_hamiltonian` production code.** Only adds a test.
- **Does NOT audit multi-radial l=0 or l=1 projectors on Cu.** If
  Phase 1 passes, those sum-over-ν surfaces are analogous but
  less load-bearing for Cu's d-state physics; file VNLM-CUP/CUS
  follow-ups only if the d-channel closes and the total residual
  doesn't.
- **Does NOT generalize to other multi-ν-per-l elements** (Fe,
  Ga, As, Mn all have 2-ν-per-l layouts). The test covers Cu
  specifically because Cu is the canonical VGCH-2 witness. If
  Cu's assembly is correct, we have strong indirect evidence
  that other multi-ν elements are correct too.
- **Does NOT touch the angular surface** — VNMT already pins
  Σ_m Y_{l,m} Y_{l,m} via the addition theorem. This test rides
  on VNMT's angular correctness and isolates the ν-summation
  correctness.

### Risks

1. **Test is trace-equivalent but projector-wrong.** The reason
   the Researcher prompt flags this whole class — a test that pins
   the sum might still pass with an asymmetric per-term bug.
   Mitigated by Phase 2's two mutation sanity checks: mutation A
   forces a per-ν magnitude asymmetry; mutation B forces a per-ν
   sign asymmetry. Both must detect. **The test is not valid
   without the mutation-firing evidence in the PR body.**
2. **VNMT's cubic a = 2π Å geometry plus Cu.upf might load-exercise
   code paths that aren't structurally representative of production
   Cu.** Specifically: `BasisSet::new(&lattice, 50.0 eV)` on a = 2π
   Å has `n_pw ≈ 20`, vs production Cu (a = 3.61 Å, ecut = 408 eV)
   at `n_pw ≈ 1500`. This is intentional — we're testing the
   per-channel contraction logic, which is n_pw-independent. But
   if `NonlocalPotential::new` or `add_to_hamiltonian` have
   n_pw-dependent code paths (bound checks, GEMM size thresholds,
   etc.), the low-n_pw test might miss them. Mitigation: Phase 3
   VQEF Cu cell re-check catches any production-size-dependent
   divergence.
3. **Cu.upf's `D_44` and `D_55` magnitudes.** If either is very
   small (< 1 eV), the test's per-m residual at 1e-10 may sit
   at floating-point noise. Check magnitudes in the test setup
   before running hand-computation. If D_55 is < 1 eV, replace
   the test's `assert!(d_55.abs() > 1e-6)` with a magnitude
   guard and bump the tolerance to 1e-8.
4. **Complex64 vs real handling of τ = 0 phase.** VNMT keeps
   τ = 0 so `exp(−iG·τ) = 1`; Cu test must preserve this.
   Already codified — single-atom-at-origin geometry.

### Provenance

- VGCH-2F: `proposals/VGCH-2F-part-c-session-2-findings.md` §
  H-C5 evidence + § Scope decision — H-C5 narrowed to the
  Cu-specific two-ν sum.
- VGCH-2 Part B: `proposals/completed/VGCH-*-part-b-*.md` (Cu
  transplant); +8.87 eV ΔE_xc at shared density.
- VNMT: `proposals/completed/VNLM-*.md` — single-GEMM V_NL
  assembly (5.2×); VNMT test at `potential/nonlocal.rs:914`
  (`test_single_channel_l2_m_isolation`).
- Cu.upf structural verification:
  `pseudopotentials/nc/lda/Cu.upf` — 6 β projectors,
  (l=0,0,1,1,2,2), diagonal PP_DIJ (verified via VGCH-2F inspection).
- Si.upf structural parity:
  `pseudopotentials/nc/lda/Si.upf` — 6 β projectors, diagonal
  PP_DIJ; identical layout to Cu.upf.
- Rust V_NL assembly sites:
  - `potential/nonlocal.rs:111` (`NonlocalPotential` struct).
  - `potential/nonlocal.rs:195` (`NonlocalPotential::new`).
  - `potential/nonlocal.rs:364–365` (`D·B^H` per-channel loop
    comment).
  - `potential/nonlocal.rs:412` (`add_to_hamiltonian`, GEMM).
- Paper citations:
  - Kleinman & Bylander, *Phys. Rev. Lett.* **48**, 1425 (1982).
  - Blöchl, *Phys. Rev. B* **41**, 5414 (1990), Eq. 12 (explicit
    plane-wave KB matrix element with Σ_{i,j} D_{ij}).
  - Gonze *et al.*, *Comput. Mater. Sci.* **25**, 478 (2002) §II.C
    (GEMM-friendly B·D·B^H layout used by VNLM).

### Cost

- Phase 1: 1 CE-day (test authoring + hand-computed per-m reference).
- Phase 2: ~15 min (two mutation sanity checks, in-test only, no
  production impact).
- Phase 3: Tier-2 re-run (automated; ~10 min wall).

Total: **1 CE-day**.
