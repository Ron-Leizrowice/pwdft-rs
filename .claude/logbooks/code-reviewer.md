# Code Reviewer Logbook

Entries: date, metrics (actual counts), findings, proposals affected. Track quality trends over time.

## 2026-04-19 — TSPL landed (PR #150)

18 new `#[ignore = "TSPL Tier-2: ..."]` tags. Tier-1 wall: 95 s → **12 s warm** (8× on top of TPRF's 7×; net 55× vs. opt-level=0 baseline). Tier-2 wall: 58 s warm (passing tests only). Gate: 271 unit + 49 integration pass; clippy 18/24 baseline unchanged; rustdoc clean.

Breakdown of new Tier-2:
- `parallel_consistency.rs`: 2 SCF tests (FFT round-trips stay Tier-1)
- `spin_polarization.rs`: 4 (all SCF-running tests)
- `vgc5_per_component_si.rs`: 4 (vgc5_si/fe_per_component + madoc_band_sum_si/fe)
- `wfrx_subspace_consistency.rs`: 3 (all — each runs 2 SCFs)
- `itev_iterative_eigensolver.rs`: 1 end-to-end (defect-1 single-shots stay Tier-1)
- `qe_validation.rs`: 1 (Si-LDA only; Fe/Cu/GaAs/NaCl/MgO/C/Al arms reserved for other agents; Si-Fermi already has VGCH ignore)
- `gpu_consistency.rs`: 3 SCF tests (kernel-unit tests stay Tier-1)

Fe NLCC regression guard (`test_fe_bcc_xc_nlcc_regression_guard`) **not touched** — uses Fe fixture but is a defensive ecut=15/4×4×4 NLCC guard, not a Fe-vs-QE comparison arm. Left Tier-1 even though it runs SCF — it's fast enough at nspin=1 and catches the NCFX 49 eV pathology on every commit, worth the ~10 s cost.

Limitation: `cargo test -- --ignored` runs the union of TSPL Tier-2 + physics-blockers (VGCH heavy-atom, Al ecut, C mixer, MXBA adaptive-β failure). The 8 physics-blocker tests still fail by design when run with `--ignored`. PR body documents the skip list; reviewer discipline required.

Unexpected rebase: `qe_validation.rs` was modified by a concurrent PR (GGAP Phase B #145) mid-session. TSPL Si ignore absorbed cleanly via `rebase origin/main`; no manual fixup needed.

## 2026-04-18 — DEAD landed all three sub-items (PR #127)

DRSD + SMRT + DHPC all shipped in one bundle after rebasing around HKIN/XCTH/ELMN concurrent merges.

Net LOC delta: **−239** (98 ins / 337 del). `real_space.rs` gone (−269). Gate: 306 tests pass, clippy 17/23 unchanged baseline, doc clean.

Notes for next sweep:
- HKIN's landing mid-session freed `tests/free_electron_bands.rs` so DHPC could land fully. If blocking constraints ever force DHPC-only deferral, keep the bench `diagonalize_lowest(h, h.nrows())` change — it's equivalent but survives on its own.
- `symmetrize_real_ref` inline in g_space.rs tests = ~75 LOC (proposal said 30). The extra LOC is per-operation loop body; non-negotiable — can't shrink without losing readability.
- `src/main.rs:88` comment still says "symmetrize_density" (identifier gone). Flagged for Technical Writer as single-line cleanup.

## 2026-04-19 — TRV2 fresh test-suite review

Wrote `proposals/TRV2-test-suite-review-2026-04-19.md`. PR #91. Targeted gaps TACC/TAUD left. 15 findings, 5 categories.

Counts: Cat1=5 (physics gap), Cat2=1, Cat3=3 (brittle), Cat4=2 (org), Cat5=4 (bench). 6-PR fix sequence.

Top-3 Cat1 gaps (true physics coverage holes):
1. **MADOC band-sum identity** `E_band = e_kin + e_loc + e_nl + 2·e_H + e_vxc` documented `src/scf/energy.rs:596-608` but unpinned. VGC5 self-check is trivial (`Σ = total` by construction). Factor-2 Hartree bug would pass.
2. **CCMX basis-change round-trip** inlined `src/scf/driver_spin.rs:537-577`, zero `#[test]`. Sign-swap on ρ↓ branch that converges wrong evades all assertions.
3. **NLCC Cu/Mn blind** — NCFX fix is universal; pins are Si/Fe only (`convert.rs:305/332/369/395`). 7 other LDA PPs have `core_correction=T`.

Fix PRs sequenced biggest-win-first: PR 1 MADOC (~80 LOC) → PR 2 CCMX helper extract (~40) → PR 3 NLCC Cu/Mn (~40) → PR 4 GPU nspin=2 (~70) → PR 5 fixture dedup (~700 net deletion) → PR 6 bench expansion (~100).

Incidental: `benches/scf_benchmarks.rs:69` has stale Accelerate TMO comment (no Accelerate in stack post-faer). Flagged for §5.1.

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
