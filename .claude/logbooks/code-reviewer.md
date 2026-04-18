# Code Reviewer Logbook

Entries: date, metrics (actual counts), findings, proposals affected. Track quality trends over time.

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
