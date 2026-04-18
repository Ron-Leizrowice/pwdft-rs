---
id: TACC
status: proposed
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# TACC: Test Accuracy + Relevance Audit (post-TAUD refresh)

## Motivation

TAUD (~6 weeks ago) removed silent-pass arms, zip-code tolerances, and
stale assertion-comments pre-NCFX. Since then, 10 PRs landed ~15 tests
(NCFX ρ_core pin, NLCC A/B/C, CCMX Fe free-mag, PCFX identity guard, ITEV
8 eigensolver tests, PRPL siblings). CCMX shipped with a stale
`final_delta` comment caught in review.

Three prompts: (a) ITEV carries one upstream-blocked `#[ignore]`; (b) PRPL's
`(Err,Err)=>{}` was fixed in review — have siblings survived? (c) several
`#[ignore]` reasons cite VERF, which landed "cosmetic-only" and did *not*
close the Si gap.

## Findings

### 1. Silent-pass arms that slipped past TAUD

Two sites. PRPL's fix (`src/scf/mixing/anderson.rs:913`) wasn't applied to
its siblings:

| File:line | Test | Current arm |
|-----------|------|-------------|
| `src/scf/mixing/anderson.rs:548-550` | `test_kerker_vs_plain_scf_convergence` | `(Err(_), Err(_)) => { /* acceptable */ }` |
| `src/scf/mixing/broyden.rs:348-350` | `test_broyden_vs_plain_scf_convergence` | same pattern |

If both mixers fail (kinetic/initial-density/V_local regression), both tests
pass silently. No other silent-pass arms found in `tests/` or `src/**`.

### 2. Stale `#[ignore]` reasons (VERF archived)

VERF is in `proposals/completed/` (2026-04-17) and its completion note
says it did *not* close the Si gap. Six QE-validation tests still cite it:

| File:line | Current reason | Accurate blocker |
|-----------|----------------|------------------|
| `tests/qe_validation.rs:272` (C) | "SCF stalls…tracked with VERF" | SYKP/MPSH + VGCMP |
| `:307` (Al) | "depends on VERF/Si root-cause fix" | SYKP/MPSH |
| `:394` (GaAs) | same | VGCMP heavy-atom |
| `:432` (Cu) | same | VGCMP heavy-atom |
| `:465` (NaCl) | same | VGCMP heavy-atom |
| `:503` (MgO) | same | VGCMP heavy-atom |

Si (line 235) and Fe (line 348) were already refreshed to cite
SYKP/MPSH/VGCMP; propagate that style to the six above.

`tests/itev_iterative_eigensolver.rs:100` ignore reason (faer 0.24
Lanczos reorthogonalization) is accurate — no change.

### 3. Zip-code tolerance

`tests/lapack_smoke.rs:60-69` (`test_faer_complex_hermitian_eigen`): asserts
only `trace == 7.0` and `ev > -1e-10`. Passes for `[0, 0, 7]`. Must assert
eigenvalues match `[1, 2, 4]` with `relative_eq!(…, epsilon=1e-10)`,
mirroring the sibling real-symmetric test.

Other eV-scale tolerances (`gpu_consistency.rs` 0.1 eV) are defensible: f32
noise ~1e-3 eV/iter × ~10 iters ≈ 1e-2 eV, so 0.1 eV is ~10× expected.

### 4. Weak/dead file: `tests/fe_debug.rs`

223 LOC, residue of the Fe 210 eV investigation. Of 6 tests: 1 zip-code
(`|basis-79|<10`), 2 trivial (band-0≈0, `is_finite()`), 1 duplicate
(`Hermiticity` — already covered at `nonlocal_symmetry.rs:30`), 1 dead
(`test_fe_full_hamiltonian_eigenvalues` prints QE comparison with no
`assert!`), 1 keeper (`test_fe_ewald_energy`, <0.01 eV vs QE). VGCMP
Phases 1–4 supersede the meaningful checks. Delete file; migrate the
Ewald pin into `qe_validation.rs`.

### 5. Regression-guard coverage matrix

Every 2026-04 physics PR has a guard that trips on revert:

| PR | Guard | File:line |
|----|-------|-----------|
| NCFX | Si Q_core ≈ 0.74 (pre: 0.207) | `src/pseudopotential/upf/convert.rs:242` |
| NCFX | Fe E_xc vs QE < 1 eV (pre: ~49) | `tests/qe_validation.rs:612` |
| NLCC | 4 ρ_core(G) pins Si+Fe | `convert.rs:304–413` |
| CCMX | Fe free-mag \|HF-KS\| < 1e-3 eV (pre: ~13) | `tests/spin_polarization.rs:295,313` |
| PCFX | Si per-component Σ < 1e-5 eV (pre: 1.204) | `tests/vgc5_per_component_si.rs:302` |
| ITEV | Γ-degeneracy preserved + real-Si dense match | `src/eigensolver/iterative.rs:488,585` |

TACC does **not** need to propose new integration tests.

### 6. Harness duplication (flagged, out of scope)

9+ tests re-declare Si Γ-only ecut=100, 16³, 4-band fixtures: parallel/
spin/gpu × 2–3 each plus 3 in-module mixer tests. Shared `tests/common/`
helper would cut ~200 LOC. Flag to Core Engineer.

### 7. Stale comment spots

- `tests/spin_polarization.rs:244-246` refers to `ScfResult.final_delta` as
  "not currently exposed". Verify against post-MODR state — either add the
  tighter `Δρ` assertion the note anticipates or drop the note.
- `tests/vgc5_per_component_si.rs`: pins correctly labelled
  "pre-PCFX/PRE-NCFX"; no stale numbers found.

### 8. GPU-gate coverage

NCFX ρ_core tests in `src/pseudopotential/upf/convert.rs::tests` run
unconditionally (CPU + `--features gpu`). No `#[cfg]` oversights.

## Recommended fixes (atomic, prioritized)

1. **Silent-pass sibling fix** (`anderson.rs:548`, `broyden.rs:348`).
   Replace `(Err,Err)=>{}` with `panic!` matching PRPL at `anderson.rs:913`.
   Impact: test quality. Risk: low.
2. **Six `#[ignore]` reason strings in `qe_validation.rs`** (lines 272, 307,
   394, 432, 465, 503). Point at SYKP/MPSH (light) or VGCMP heavy-atom
   (Z>14). Impact: documentation. Risk: low.
3. **Tighten `lapack_smoke::test_faer_complex_hermitian_eigen`** to match
   `[1, 2, 4]` eigenvalues. Impact: quality. Risk: low.
4. **Delete `tests/fe_debug.rs`**, migrate only `test_fe_ewald_energy` into
   `qe_validation.rs`. −222 LOC, -5 weak tests. Impact: quality + CI.
   Risk: low.
5. **Resolve ITEV `#[ignore]` integration**: unit tests at
   `src/eigensolver/iterative.rs:538,593` already cover n_pw=89 Si. EM call:
   delete the `#[ignore]`'d integration (Option A) or shrink to 1–2 iters
   below the reorth hang regime (Option B). Impact: coverage. Risk: low/med.
6. **Strip `final_delta` TODO comment** in `tests/spin_polarization.rs:244`
   after verifying MODR state. Impact: hygiene. Risk: low.

Fixes 1–4 are independent ≤30-line diffs. 5 and 6 depend on state
verification. All six fit a single PR if EM prefers, or one PR per fix.

## What this is NOT

- **Not** replacing `relative_eq!` with exact equality.
- **Not** un-ignoring upstream-blocked tests (ITEV line 100 stays ignored
  until faer `iterate_lanczos` is fixed).
- **Not** adding new integration tests — regression-guard matrix (§5) is
  complete for every 2026-04 physics PR.
- **Not** refactoring Si-Γ fixture duplication (flagged to Core Engineer).

## Open questions

- Fix #4: keep `test_fe_full_hamiltonian_eigenvalues` as a commented-out
  debug fixture, or delete wholesale? (Recommend: delete; `git log` brings
  it back if ever needed.)
- Fix #5: Option A (delete `#[ignore]`'d integration) vs B (shrink)?
  Recommend A — zero unique coverage today.
- Worth retrying Al/C at small ecut (10 Ry, 2×2×2) to un-ignore
  *something*? Researcher call, separate proposal if pursued.

## Flagged for follow-up

- `tests/fe_debug.rs:211-223` — prints QE eigenvalues without `assert!`.
  **Code Reviewer** (this proposal, fix #4).
- Si Γ-only SCF harness duplication across 9+ tests. **Core Engineer**
  (`tests/common/` helper module).
- `ScfResult.final_delta` exposure status post-MODR. **Core Engineer**.
