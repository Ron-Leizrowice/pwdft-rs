# Code Reviewer Logbook

Entries: date, metrics (actual counts), findings, proposals affected. Track quality trends over time.

## 2026-04-19 — FLP2 FLUP sweep

Reconciled `proposals/FLUP-followup-backlog-seeding.md` against main at commit `1f568be`. PR #89 (FLP2/flup-sweep-2026-04-19).

Counts: 1 struck (DWGT, landed 5456c80 with CLAUDE.md + agent defs propagated), 0 obsolete, 6 still-live (G2ZT, DFLT, MXB2, ITVF, EIGV, EIGW), 0 promoted to standalone proposal, 1 added (TYPE-AX).

TYPE-AX added: 5 `try_from` expect sites from TYPE-A PR #80 at `src/basis.rs:65`, `src/symmetry/operations.rs:{71,120,151}`, `src/symmetry/detect.rs:185`. Not its own proposal — ERR2 P1 decision item. Structural-bound candidates, i.e. add `reason` comment, not Result.

New "Status summary" section at top of FLUP so next sweep gets one-glance state instead of grepping for `~~`.

Nothing promoted: MXB2 is small-medium but already has full RCA + three candidate tunings inline, no scoping gap. EIGV/EIGW may self-resolve on next bench pass. ITVF genuinely blocked on faer 0.25.

Did NOT add: "EnergyComponents double-counting test gap" (prompt flagged conditional on MAUD/TRV2 surfacing it — they haven't landed, so I don't preempt).

## 2026-04-18 — ERR2 panic-free production audit

Wrote `proposals/ERR2-panic-free-production-audit.md`. Post-ERRH census:

- **Category A (unwrap/expect/panic):** 0 unwrap, 15 expect (all `BUG:` prefix), 0 panic — all legitimate invariants, no Result conversions needed. +2 since ERRH (`scf/mixing/broyden.rs:66` BROY, `scf/driver.rs:103` driver-split).
- **Category B (assert/debug_assert):** 16 production sites. All keep-as-assert or keep-as-debug_assert after triage; 0 user-reachable panics.
- **Category C (unreachable/todo):** 1 `unreachable!` in `scf/mixing/anderson.rs:87` (documented). 0 `todo!`.
- **Category D (`InvalidInput` catch-all):** 12 distinct failure kinds stringly-typed into one `InvalidInput(String)` — biggest finding. Proposed split into `InvalidParam`/`InvalidCrystal`/`UnknownElement`.

Enforcement: proposed `unwrap_used`/`expect_used`/`panic`/`unreachable` = `warn`, `todo`/`unimplemented` = `deny`. Test modules get `#[allow(..., reason = "ERR2 Phase 0")]` banners (CAST discipline).

Surprise: UPFV (PR #72) landed mid-audit and obsoleted the only Category B "should-be-Result" candidate (`potential/nonlocal.rs:132`), which is now belt-and-suspenders. Updated proposal accordingly.

Cosmetic nit flagged for Phase 0: `scf/mixing/broyden.rs:66` missed the `BUG:` prefix that all other production expects use.

## 2026-04-17 — TAUD test-suite audit

Wrote `proposals/TAUD-test-quality-audit.md`. Top 3 silent-pass findings:
1. `tests/spin_polarization.rs::test_fe_ferromagnetic_fixed_moment` — `match result { Ok=>assert, Err=>eprintln }` means the test passes vacuously; the `Ok` arm is unreachable post-SPNC.
2. `tests/parallel_consistency.rs:178,238` — serial-vs-parallel SCF tests guard real assertions with `if let (Ok, Ok)`. Any regression silently passes.
3. `tests/gpu_consistency.rs:406` — `(Err, Err)` treated as pass condition. Critical.

Counts: 5 critical silent-pass, 8 tolerance gaps, 1 confirmed stale `#[ignore]`, 8 conditional stale (QE validation, gated on VGCMP).

## 2026-04-17 — MUST / ERRH / SDED

- **MUST:** Applied `#[must_use]` to 70/80 candidate functions; 10 deferred for parallel-agent safety. Remaining 10 should close after SPXC/XCPR/VGCMP merge.
- **ERRH:** 16 production unwrap/panic/expect → 0. 13 remaining `expect()` all documented with `BUG:` prefix.
- **SDED:** `MixingModeType` and `OccupationType` dedup closed; `From` conversions preferred over inline match in YAML settings paths.

## 2026-04-16 — Baseline sweep

Clippy 0; production unwrap 11, panic 2, expect 3; 185 tests, 3 failing (kb_projector_validation — later fixed by KBTF). DDUP 5 issues present.
