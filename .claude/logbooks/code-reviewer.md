# Code Reviewer Logbook

Entries: date, metrics (actual counts), findings, proposals affected. Track quality trends over time.

## 2026-04-16 — Baseline sweep

Clippy warnings 0; production unwrap 11, panic 2, expect 3; `#[allow()]` 3; 0 TODOs; 185 tests, 3 failing (kb_projector_validation). Proposals confirmed: CNST/CBRT/DLTB done; CLEN has 5 remaining; ERRH count refresh 25→16; SDED reduced (MixingModeType+OccupationType remain); DDUP 5 issues present.

## 2026-04-16 — SDED (PR #5)

MixingModeType dedup: added `From<MixingModeType> for MixingMode`, replaced inline match in `to_scf_params` with `.into()`, added conversion test. OccupationType confirmed no-op. 160 tests pass, 0 clippy. SmearingType was done prior; this closes SDED.

## 2026-04-16 — ERRH (PR #9)

Replaced all 16 production unwrap/panic/expect with Result/Option propagation. Added `Eigensolver` and `Gpu` variants to `PwdftError`. Signatures updated: `diagonalize_hermitian`, `diagonalize_lowest`, `find_for_atom`, `to_crystal`, `LocalPotential::new`, `NonlocalPotential::new`, `ScfContext::new`, `compute_band_structure`. 13 remaining `expect()` calls documented with `BUG:` prefix. Metrics: unwrap 11→0, panic 2→0, undocumented expect 3→0, BUG-tagged expect 0→13. 183 tests pass.

## 2026-04-17 — MUST implementation

**Done:** Applied `#[must_use]` to 70 functions/methods flagged by `-W clippy::must_use_candidate` across 19 files. Skipped 10 hits in 3 files deferred for parallel-agent safety (`src/potential/xc.rs`: 5, `src/pseudopotential/mod.rs`: 5, `src/scf/mod.rs`: 0). PR opened on branch `MUST/must-use-attributes`.

**Counts:** 80 hits total → 70 applied, 10 deferred. Follow-up PR needed after SPXC/XCPR/VGCMP merge.

**Quality gate:** `cargo test` all pass (same 8 ignored pre-existing), `cargo clippy -q --all-targets` = 0 warnings. `-W clippy::must_use_candidate` goes from 80 → 10 (the 10 deferred).

**No restructurings needed.** Every hit was a straight `#[must_use]` add; no borderline cases required `#[allow]` or signature changes. Diff is purely additive.

## 2026-04-17 — TAUD test-suite audit

**Wrote:** `proposals/TAUD-test-quality-audit.md`. Investigation only — zero `src/`/`tests/` edits.

**Top 3 most alarming findings:**
1. `tests/spin_polarization.rs:239` — `test_fe_ferromagnetic_fixed_moment` is a direct clone of the bug SPNC just fixed. SCF is known (post-SPNC) to diverge on this PP, but `match result { Ok=>assert, Err=>eprintln }` means the test passes anyway. The `Ok` arm is unreachable. Either delete or invert to `assert!(result.is_err())`.
2. `tests/parallel_consistency.rs:178,238` — both serial-vs-parallel SCF tests guard their real assertions with `if let (Ok, Ok)`. A regression in either path makes the test pass silently. Kerker variant is particularly exposed.
3. `tests/gpu_consistency.rs:406-407` — `(Err(e1), Err(e2)) => eprintln!("both did not converge")` treats simultaneous GPU+CPU divergence as a pass condition. That is a critical bug, not a pass condition.

**Counts:** 5 critical silent-pass, 8 tolerance gaps, 1 confirmed stale `#[ignore]` (kb_projector test_vloc — VERF landed), 8 conditional stale (QE validation, gated on VGCMP), 3 convergence-criterion blind spots.

**Fix plan:** 5 PR-sized chunks (A-E) grouped in the proposal. Each < 1 hr. Total ~5 follow-up PRs to address everything actionable; 2 items deferred (fe_debug.rs deletion waits for VGCMP Phase 2; 8 QE-validation `#[ignore]`s unlock as VGCMP phases close).

**Quality gate:** N/A (investigation only, no code change). Branch `TAUD/test-quality-audit` off `origin/main`.
