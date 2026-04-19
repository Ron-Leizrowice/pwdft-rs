---
id: MIXA
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# MIXA: Pin Anderson-stall-on-wide-gap-insulators as negative regression

## Problem

Plain Anderson mixing (`MixingMode::Plain`) stalls on wide-gap insulators.
Two independent observations establish the pattern as a real, latent,
mixer-robustness limitation rather than a one-off anomaly:

1. **C diamond** (VQEF / `test_c_diamond_vs_qe`, 2026-04-19 sweep): at
   ecut = 30 Ry on a Γ-centered 4×4×4 grid, plain Anderson
   (β = 0.3, ndim = 8) plateaus at Δρ ≈ 1.75e-8 and does not decay
   within 150 iterations (driver returns `ConvergenceFailure` against a
   `conv_threshold = 1e-8` target). With `β = 0.1` or `ndim = 16`
   Anderson stalls even higher at Δρ ≈ 1e-5. Every alternative
   converges cleanly in 10–15 iterations: Kerker auto (15 iters,
   Δρ 2.2e-10), Kerker `q_tf = 0.5` (12 iters), Broyden (12 iters),
   Broyden + Kerker (12 iters, Δρ 8e-11), PeriodicPulay `p = 4`
   (10 iters).
2. **Si at FFT grid 20/24** (PCRS investigation, 2026-04-17, Researcher
   logbook line 94): "Anderson mixer loses conditioning and oscillates"
   at `conv_threshold = 1e-8`. Same pathology on a different system.

The underlying mechanism is long-wavelength charge sloshing that Anderson
alone — without Kerker's `|G|² / (|G|² + q_TF²)` preconditioner or
Broyden's approximate inverse Jacobian — cannot damp. Wide-gap insulators
are susceptible because their dielectric response is weak at small |G|.

## Recommendation

**Document and pin. Do not fix right now.** The robust-mixer path
(Kerker / Broyden / Periodic Pulay) already exists and converges to the
same physics — so the Anderson limitation is a *performance* issue
(forces the user to pick the right mixer), not a *correctness* issue.
Fixing plain Anderson to match Kerker/Broyden robustness is a nice-to-
have but out of scope for the current VQEF validation push; the VQEF
C-diamond test already pins `Broyden { kerker: true }` as the mixer
under which the QE comparison runs, so the blocker there is the
residual 1.45 eV energy gap (VGCH light-atom), not the mixer choice.

Scope of MIXA:

1. **Regression test** — `tests/mixer_robustness.rs` with one
   `#[test]` that runs C diamond SCF under `MixingMode::Plain` with
   `max_iter = 150` and asserts the driver returns
   `ConvergenceFailure { delta, .. }` rather than `Ok(..)`. The test is
   `#[ignore]`'d (Tier-2) so default `cargo test` stays green. Two
   outcomes flag human attention: `Ok(res)` (someone fixed Plain
   Anderson — rewrite as positive convergence pin), or `Err(other)`
   (failed for an unexpected reason — don't paper over). A secondary
   `delta > 1e-8` assertion guards against the driver reporting
   ConvergenceFailure on an already-converged state.
2. **Docstring update** on `MixingMode::Plain` in
   `src/scf/mixing/mod.rs`: short `///` note documenting the
   wide-gap-insulator limitation and pointing users at Kerker / Broyden
   / PeriodicPulay. No proposal-ID tokens in public rustdoc (DCLN).
3. **This file + INDEX entry** — one-line reference.

## Non-goals

- **Not fixing Anderson.** A real fix would likely mean either
  auto-enabling Kerker-like preconditioning inside Plain when the user
  doesn't supply one, or swapping the default mixer to Broyden + Kerker.
  Both are ELMN-scale design changes that deserve their own proposal.
- **Not tightening the C diamond QE validation.** That arm already
  pins `Broyden { kerker: true }` and the residual is VGCH territory.
- **Not a per-system mixer-selection heuristic.** The user picks the
  mixer via YAML; MIXA only documents the tradeoffs.

## Success criteria

- `cargo test --test mixer_robustness -- --ignored` runs the new test
  and it passes (i.e. the documented stall still reproduces).
- `cargo test` (Tier-1) unchanged — new test is `#[ignore]`'d.
- Public rustdoc on `MixingMode::Plain` mentions the limitation and
  suggests alternatives.

## Follow-up triggers

Reopen MIXA (or spawn a follow-up `MIXB`) when:

- Someone empirically fixes plain Anderson (e.g. a small default change
  to DIIS re-initialization or history trimming) and the test flips —
  then the test becomes a positive convergence pin.
- A third observation of the stall lands on yet another system,
  strengthening the case that plain Anderson is too weak a default.
- The default `MixingMode` changes (currently `Plain` via `#[default]`)
  — the new default will need its own regression-pin decision.

## References

- `proposals/VQEF-qe-validation-full-matrix.md` § 2a mentions the
  C-diamond mixer sweep.
- `proposals/completed/PCRS-per-component-residual.md` — Si at FFT
  grid 20/24, same stall pattern (Researcher logbook line 94).
- `tests/qe_validation.rs::test_c_diamond_vs_qe` — full sweep table
  in the docstring.
- `tests/mxba_adaptive_beta_fe.rs` — prior art for the
  document-a-failure-as-negative-regression pattern.
