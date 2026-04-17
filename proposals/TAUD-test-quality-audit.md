---
id: TAUD
status: active
priority: high
complexity: medium
risk: low
depends_on: []
blocks: []
---

# TAUD: Test Suite Quality Audit (Silent-Pass and Tolerance Gaps)

## Origin

**2026-04-17.** SPNC (PR landed same day) exposed that
`tests/spin_polarization.rs::test_fe_ferromagnetic_fixed_moment` had been
silently passing while Fe BCC fixed-mag=2 was not actually converging: the
per-spin Anderson mixers entered a `±ε` limit cycle that cancelled in the
total density, and the total-only convergence criterion hid it. Pre-SPNC
|HF-KS| was 22.2 eV (pre-SPXC) / 13.0 eV (post-SPXC) — both at a "converged"
state that was anything but.

If *one* test was lying for months, others almost certainly are. This
proposal audits the whole suite looking for the same pattern and related
smells: `match`/`if let` on fallible SCF calls that only `eprintln!` on the
`Err` branch without failing the test, tolerances set to an eV when the
natural scale is µeV, stale `#[ignore]`s gated on problems already fixed,
and assertions that check a tangential property instead of the one named in
the test.

## Executive Summary

Audit scope: 12 integration-test files in `tests/` (63 `#[test]`s) plus 27
`#[cfg(test)]` modules in `src/` (~195 `#[test]`s). Total ~258 tests.

**Findings:**

| Category | Count | Severity of worst offender |
|----------|-------|----------------------------|
| Silent-pass SCF result patterns | **5** | critical — would have hidden SPNC |
| Tolerance gaps (placeholder or off-scale) | **8** | high |
| Stale / reason-no-longer-valid `#[ignore]`s | **1** confirmed + **8** conditional | medium |
| Duplicate coverage | **2** clusters | low |
| Convergence-criterion blind spots | **3** | high |

The 5 silent-pass SCF patterns are the most alarming: they are all copies of
exactly the bug that hid the Fe BCC false convergence. Each one
`eprintln!`s on `Err` and falls off the end of the test function, which
`cargo test` treats as pass.

---

## Per-Finding Table

Severity legend:
- **critical** — would mask a bug that silently changes physics / total energy
- **high** — masks bugs that change assertions by >10× the natural scale
- **medium** — documentation drift; test still runs correctly
- **low** — cleanup / future-proofing

### Category 1 — Silent-pass SCF result patterns

| # | Severity | File:line | Description | Suggested fix |
|---|----------|-----------|-------------|---------------|
| 1.1 | **critical** | `tests/spin_polarization.rs:78-101` | `test_si_nspin2_matches_nspin1`: 4-arm `match (result1, result2)` where only `(Ok, Ok)` runs assertions; `(Err, _)` and `(_, Err)` arms just `eprintln!` then fall through. Test passes even if the whole SCF pipeline returns `Err`. | Replace with `let r1 = result1.expect("nspin=1 must converge"); let r2 = result2.expect("nspin=2 must converge");` and then assert. |
| 1.2 | **critical** | `tests/spin_polarization.rs:239-259` | `test_fe_ferromagnetic_fixed_moment`: `match result { Ok(r) => {...assertions...}, Err(e) => eprintln!(...) }`. This is the exact test SPNC proved was lying. Since SPNC, Fe fixed-mag=2 is *known* to not converge on this PP — the `Err` branch is guaranteed to fire, but the test still reports pass. The `Ok` arm's magnetization check is unreachable. | Either (a) delete the test (SPNC already documents it as unreachable for this PP — see `SPNC-spin-per-density-convergence.md` §Empirical Result), or (b) change to `assert!(result.is_err(), "Fe fixed-mag=2 is known not to converge under nc/lda/Fe.upf; delete this test or change the pseudopotential");` so the day the PP changes, somebody notices. |
| 1.3 | **critical** | `tests/parallel_consistency.rs:178-194` | `test_scf_serial_vs_parallel`: `if let (Ok(s), Ok(p)) = (&result_s, &result_p)` — eigenvalue + energy assertions only run if both converge. If either path fails, test passes trivially. The serial-vs-parallel property is exactly what needs asserting even if only one converges (that would itself be a bug). | `let s = result_s.expect("serial SCF must converge"); let p = result_p.expect("parallel SCF must converge");` then compare. |
| 1.4 | **critical** | `tests/parallel_consistency.rs:238-253` | `test_scf_kerker_serial_vs_parallel`: identical pattern to 1.3. Kerker path is specifically the one that has historically been flakier — silencing its failures is exactly backwards. | Same fix as 1.3. |
| 1.5 | **critical** | `tests/gpu_consistency.rs:366-415` | `test_gpu_vs_cpu_scf_direct_comparison`: 4-arm match. `(Err, Err)` arm just `eprintln!`s "both did not converge" and passes. "GPU and CPU both diverge" is a critical bug, not a pass condition. | Change `(Err, Err)` arm to `panic!("both GPU and CPU SCF diverged: gpu={e1}, cpu={e2}");`. The other three arms are correctly handled (two panic, one passes). |

**Additional, lower-severity silent-pass patterns:**

| # | Severity | File:line | Description | Suggested fix |
|---|----------|-----------|-------------|---------------|
| 1.6 | medium | `tests/gpu_consistency.rs:242-248` | `test_gpu_vs_cpu_scf_eigenvalues`: `match gpu_result { Ok(r) => ..., Err(e) => { eprintln!; return; } }`. Early-return on GPU divergence is OK if the test is truly best-effort, but it is not documented as such. | Add `assert!(false, "GPU SCF did not converge; this likely indicates a real regression")` or a `#[ignore]` with explanatory comment if the test is intentionally best-effort on hosts without GPU (it already has `try_new()` short-circuit for that). |
| 1.7 | medium | `tests/gpu_consistency.rs:456-468` | `test_gpu_scf_kerker_converges`: `match result { Ok(r) => assert_energy_range, Err(e) => eprintln! }`. Test is literally named "…_converges" yet passes when it does not. | Replace with `let r = result.expect("GPU Kerker SCF must converge");`. |

---

### Category 2 — Tolerance gaps

| # | Severity | File:line | Description | Suggested fix |
|---|----------|-----------|-------------|---------------|
| 2.1 | **high** | `tests/spin_polarization.rs:87-90` | `assert!(de < 0.5, ...)` — 0.5 eV tolerance on a Si nspin=1 vs nspin=2 energy comparison where the expected difference is machine-zero (same physics in the unpolarized limit). SPNC demonstrated this system reaches O(µeV). | Tighten to `de < 1e-3` (or, better, measure the actual value and set at 10× empirical). |
| 2.2 | **high** | `tests/spin_polarization.rs:93-96` | `assert!(r2.magnetization < 0.1, ...)` — Si is non-magnetic; M should be exactly 0 at convergence. 0.1 µB is 10% of a full electron spin. | Tighten to `r2.magnetization < 1e-3` or `.abs() < 1e-3`. Note: the test uses `<` not `.abs() < ` which allows arbitrarily-negative magnetization. Fix both. |
| 2.3 | high | `tests/spin_polarization.rs:247-249` | `assert!((r.magnetization - 2.0).abs() < 0.5, ...)` — fixed-mag=2 should give *exactly* 2.0 by construction (it's a constraint, not an observable). 0.5 µB is huge. Becomes moot if 1.2 is fixed (test deleted). | Tighten to `< 1e-6` if the test is kept at all. |
| 2.4 | medium | `tests/gpu_consistency.rs:100-103` | `assert!(rel < 1e-3, ...)` for Hartree GPU-vs-CPU. Hartree is a linear operator — its GPU/CPU comparison should be bounded by f32 relative precision (~1e-6) per op times a few ops, so 1e-5 is the right scale. 1e-3 is 100× looser than needed. | Tighten to `1e-5` and watch for which case (if any) triggers failure — that will tell us where f32 loses precision. |
| 2.5 | medium | `tests/gpu_consistency.rs:255-259` | `gpu_result.total_energy < -100.0 && > -300.0` — a 200 eV-wide window for a Si energy that is known to be around -212 eV. This is an "is it in the zip code" check. | Replace with a specific value ± tolerance, e.g. `(total_energy + 212.0).abs() < 1.0`. Cross-reference against the single-threaded `parallel_consistency` golden value. |
| 2.6 | medium | `tests/gpu_consistency.rs:262-266` | `fermi_energy > -5.0 && < 10.0` — 15 eV-wide window for Si E_F. Same smell as 2.5. | Replace with a specific value ± tolerance. |
| 2.7 | medium | `tests/gpu_consistency.rs:302-308` | `e < gpu_result.fermi_energy + 5.0 * sigma` — 5σ above E_F for "occupied" bands. σ=0.05 eV so this allows 0.25 eV above E_F. Thermodynamically justified, but the comment just says "within smearing width" — the 5× factor is uncommented magic. | Document the 5× explicitly (`5σ ≈ f=exp(-5) ≈ 0.7% Fermi tail`) or tighten to 3σ. |
| 2.8 | medium | `tests/gpu_consistency.rs:460-463` | `r.total_energy < -100.0 && > -300.0` — duplicate of 2.5 in `test_gpu_scf_kerker_converges`. | Same fix. |

**Also surfaced but acceptable with their current comments:**

- `tests/spin_polarization.rs:187-193` — `1e-5 eV` tolerance with a specific documented empirical value (7.19e-7). Properly calibrated. *No action.*
- `tests/gpu_consistency.rs:146-154` — XC `exc_err < 0.01` / `vxc_err < 0.015`. Appropriate for f32 LDA. *No action.*

---

### Category 3 — Stale `#[ignore]`s

| # | Severity | File:line | Description | Suggested fix |
|---|----------|-----------|-------------|---------------|
| 3.1 | **medium** | `tests/kb_projector_validation.rs:934` | `test_vloc_comparison_with_qe`: `#[ignore = "VERF: v_local_of_g uses bare Coulomb subtraction; QE uses erf(r)/r..."]`. VERF landed on 2026-04-17 (see `proposals/INDEX.md` Completed table). The test's stated blocker is gone. | Unignore. Run the test. If it fails, the QE reference numbers in the test (lines 976-978) have a legitimate grid/quadrature mismatch that needs a new tolerance or fix — document whichever. |
| 3.2 | conditional | `tests/qe_validation.rs:229,266,301,334,380,418,451,489` | All 8 QE-validation tests `#[ignore]` pending VGCMP/Si-13.4-eV fix. Still valid for now. | Keep as-is; unignore one-by-one as VGCMP phases close the gap. Add a follow-up TODO comment at the top of the file to re-audit this list when VGCMP Phase 4 lands. |

---

### Category 4 — Duplicate coverage

| # | Severity | File:line | Description | Suggested fix |
|---|----------|-----------|-------------|---------------|
| 4.1 | low | `tests/kb_projector_validation.rs:404-514` (test_05) and `:649-755` (test_08) | test_05 checks diagonal V_NL(G,G) against manual computation; test_08 checks off-diagonal V_NL(G,G') against manual. Together they're good. Not duplicated, keeping for tracking. | **No action** — flagged for awareness only. |
| 4.2 | low | `tests/fe_debug.rs` entire file (6 tests) | This file is a diagnostic-era artifact from Fe VGCMP investigation. Every test only `eprintln!`s diagnostic info plus one or two weak sanity asserts (e.g. `test_fe_v_local_at_g0` only asserts `v_g0.is_finite()`). All the information it carries is now subsumed by `tests/vgcmp_vloc_cross_check.rs` and `tests/vgcmp_beta_q_cross_check.rs`, which compare against Python references. | Delete `tests/fe_debug.rs` entirely after VGCMP Phase 2 lands and the per-shell comparison is trusted. **Do not delete in this audit** — gate on VGCMP. Open as follow-up. |

---

### Category 5 — Convergence-criterion blind spots (tests calling SCF)

| # | Severity | File:line | Description | Suggested fix |
|---|----------|-----------|-------------|---------------|
| 5.1 | high | `tests/spin_polarization.rs:75-77` | `run_scf` called twice; neither result's `.n_iterations` nor `delta` is ever asserted against `max_iter` or `conv_threshold`. A run that took 39/40 iters (one iter away from `ConvergenceFailure`) passes the same as a run that took 8. | After unwrapping, assert `result.n_iterations < params.max_iter` (with margin). |
| 5.2 | high | `tests/gpu_consistency.rs:206-314` | `test_gpu_vs_cpu_scf_eigenvalues`: converges with `conv_threshold = 1e-6` but never asserts delta reached that threshold. Same blind spot as 5.1. | Same fix. |
| 5.3 | medium | `tests/parallel_consistency.rs:155-156` | `assert!(result_serial.is_err()); assert!(result_parallel.is_err());` — this correctly asserts `ConvergenceFailure` (the test intentionally runs only 5 iters). The comment says "Won't converge in 5 iters — that's fine" but the assertions don't distinguish `ConvergenceFailure` from any other `Err` variant. If the SCF started throwing `Gpu` or `Eigensolver` errors for unrelated reasons, the test would still pass. | Assert on the specific error variant: `assert!(matches!(result_serial, Err(PwdftError::ConvergenceFailure { .. })));`. |

---

## Prioritized Fix Plan

Group into PR-sized chunks (each < 1 hour of Core Engineer work):

### PR A — Silent-pass critical fixes (5 tests)
Handles findings **1.1, 1.3, 1.4, 1.5, 1.7**. Straight mechanical edits:
replace `match`/`if let` on `Result<ScfResult, _>` with `.expect("...")`.
For the "Both Err" cases (1.5) swap `eprintln!` for `panic!`.

Add a follow-up assertion for 5.1 and 5.2 at the same time: `assert!(r.n_iterations < params.max_iter)`.

**Files:** `tests/spin_polarization.rs`, `tests/parallel_consistency.rs`, `tests/gpu_consistency.rs`.
**Expected risk:** some tests may now fail — that is the *point* of the
audit. If a test fails, the root cause is a real bug that needs a follow-up
proposal (the bug was there all along; we just started looking).
**Validation:** `cargo test` (CPU) and `cargo test --features gpu` (GPU) must both pass post-fix. If any fail, open a follow-up.

### PR B — Delete or invert the Fe fixed-mag test (1 test)
Finding **1.2**. Per SPNC `proposals/SPNC-spin-per-density-convergence.md`
§Empirical Result, Fe BCC fixed-mag=2 is known to diverge under the current
pseudopotential. The test's `Ok` branch is unreachable; the test exists
only as a form of regression detector, but it regresses silently. Either
delete or invert (`assert!(result.is_err())`).

Recommend **invert** over delete — the day somebody tries a different Fe
PP or a better mixer, the test reverses and flags the improvement for
review.

**Files:** `tests/spin_polarization.rs`.

### PR C — Tolerance tightening (6 call sites)
Findings **2.1, 2.2, 2.3, 2.4, 2.5 (and 2.8), 2.6, 2.7**.
Mechanical constant edits, but each requires running the affected test
once to confirm the tighter tolerance still passes (or to set a
data-driven value).

**Files:** `tests/spin_polarization.rs`, `tests/gpu_consistency.rs`.
**Expected risk:** some tests may now fail → either the current assertion is wrong (code bug) or the empirical value is noisier than we thought. Both are valuable discoveries.

### PR D — Unignore `test_vloc_comparison_with_qe`
Finding **3.1**. Remove the `#[ignore]`, run the test. If it passes, great —
the VERF work transitively fixed KB test 11 as well. If it fails, tighten
the tolerance to something data-driven (current hardcoded `0.5 eV` from
line 1015-1022 is already loose) or open a follow-up proposal.

**Files:** `tests/kb_projector_validation.rs`.

### PR E — Convergence-variant specificity (1 file)
Finding **5.3**. Replace `.is_err()` with `matches!(..., Err(PwdftError::ConvergenceFailure { .. }))` in `test_scf_serial_vs_parallel`.

**Files:** `tests/parallel_consistency.rs`.

### Deferred (not in this audit's output PRs)
- Finding **4.2** (delete `tests/fe_debug.rs`): wait for VGCMP Phase 2 landing.
- Finding **3.2** (QE validation `#[ignore]`s): unignore as VGCMP phases close the gap.

---

## Explicitly NOT Fixed in This Audit

This audit **identifies** problems. Fixes are scheduled as PRs A-E above.
Specifically, this proposal does not:

- Modify any test file, any `src/` file, or any other proposal.
- Attempt to diagnose *why* any assertion is currently loose — that's up
  to the Core Engineer implementing the fix PR.
- Delete `tests/fe_debug.rs` (defer to VGCMP Phase 2 completion).
- Re-examine the 8 QE-validation `#[ignore]`s beyond confirming VGCMP is
  still live (per `proposals/INDEX.md` 2026-04-17 note).

## References

- `proposals/SPNC-spin-per-density-convergence.md` — the SPNC case study
  this audit generalizes from (per-spin cancellation false-convergence).
- `proposals/completed/SPXC-spin-xc-consistency.md` — the XC-input/output
  bug SPNC's convergence criterion was previously masking.
- `proposals/completed/VERF-vloc-erf-subtraction.md` — the landed fix that
  stales `#[ignore]` 3.1.
- `.claude/logbooks/core-engineer.md` 2026-04-17 SPNC/SPXC entries —
  corroborating empirical numbers.
- `proposals/completed/ERRH-error-handling-cleanup.md` — the `PwdftError`
  enum that 5.3 will pattern-match on.

## Status (2026-04-17)

- **PR A done** (PR #20 merged). 5 silent-pass SCF tests now hard-assert. No bugs unmasked.
- **PR B pending** — Fe fixed-mag test: invert to `assert!(result.is_err())` (per SPNC, current PP can't sustain fixed-mag=2; CCMX may fix later).
- **PR C pending** — tighten 6 tolerance call sites (findings 2.1–2.7).
- **PR D pending** — unignore `test_vloc_comparison_with_qe` (post-VERF).
- **PR E pending** — pattern-match `ConvergenceFailure` variant in `test_scf_serial_vs_parallel`.
