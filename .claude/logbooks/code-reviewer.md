# Code Reviewer Logbook

Entries: date, metrics (actual counts), findings, proposals affected. Track quality trends over time.

## 2026-04-16 — Baseline sweep

**Quality metrics:**
| Metric | Count |
|--------|-------|
| Clippy warnings | 0 |
| Production unwrap() | 11 |
| Production panic! | 2 |
| Production expect() | 3 |
| #[allow()] | 3 |
| TODOs/FIXMEs | 0 |
| Tests total | 185 |
| Tests failing | 3 (kb_projector_validation) |

**Proposal status:** CNST/CBRT/DLTB confirmed done and archived. CLEN has 5 remaining items (item #2 done). ERRH needs count refresh (25→16). SDED reduced scope (SmearingType done, MixingModeType+OccupationType remain). DDUP all 5 issues confirmed present.

## 2026-04-16 — SDED implementation

**Done:** Implemented SDED (MixingModeType dedup). Added `From<MixingModeType> for MixingMode` impl, replaced inline match in `to_scf_params` with `.into()`, added conversion test. OccupationType confirmed no-op (no SCF-internal duplicate exists). PR #5 created on branch `SDED/settings-enum-dedup`.

**Quality gate:** 160 unit tests pass, 0 clippy warnings. 2 pre-existing QE validation failures (Si energy 13.33 eV off, C diamond convergence) unrelated.

**SDED status:** Complete pending merge. SmearingType was done prior session, MixingModeType done now, OccupationType needs no action.

## 2026-04-16 — ERRH implementation

**Done:** Implemented ERRH (error handling cleanup). Replaced all 16 production unwrap/panic/expect calls with proper error propagation. Added `Eigensolver` and `Gpu` variants to `PwdftError`. Changed signatures of `diagonalize_hermitian`, `diagonalize_lowest`, `find_for_atom`, `to_crystal`, `LocalPotential::new`, `NonlocalPotential::new`, `ScfContext::new`, `compute_band_structure` to return Result/Option. Updated all callers in src/, tests/, and benches/. Documented 13 remaining `expect()` calls with `BUG:` prefix. PR #9 created on branch `ERRH/error-handling-cleanup`.

**Quality gate:** 183 tests pass, 0 clippy warnings. 2 pre-existing QE validation failures (Si energy, C diamond convergence) unrelated.

**Metrics:** Production unwrap: 11->0, panic: 2->0, undocumented expect: 3->0, documented BUG expect: 0->13.

**ERRH status:** Complete pending merge.

## 2026-04-17 — MUST implementation

**Done:** Applied `#[must_use]` to 70 functions/methods flagged by `-W clippy::must_use_candidate` across 19 files. Skipped 10 hits in 3 files deferred for parallel-agent safety (`src/potential/xc.rs`: 5, `src/pseudopotential/mod.rs`: 5, `src/scf/mod.rs`: 0). PR opened on branch `MUST/must-use-attributes`.

**Counts:** 80 hits total → 70 applied, 10 deferred. Follow-up PR needed after SPXC/XCPR/VGCMP merge.

**Quality gate:** `cargo test` all pass (same 8 ignored pre-existing), `cargo clippy -q --all-targets` = 0 warnings. `-W clippy::must_use_candidate` goes from 80 → 10 (the 10 deferred).

**No restructurings needed.** Every hit was a straight `#[must_use]` add; no borderline cases required `#[allow]` or signature changes. Diff is purely additive.
