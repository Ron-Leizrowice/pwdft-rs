# Code Reviewer — 2026-04-21 — URES unused_results lint

Implemented URES: replaced 87 `#[must_use]` annotations across 22 files with a
single `#![warn(unused_results)]` in `lib.rs`. Surfaced 7 new warnings; all triaged
and fixed without `let _ = ...`.

## Key numbers

- `#[must_use]` annotations removed: 87 (proposal said 82 — 5 added since proposal was written)
- New `unused_results` warnings surfaced: 7
- Warnings by category:
  - binding-with-comment (category 3): 7 — all were value-returning mutations where the return was intentionally discarded
  - Result-propagation (category 1): 0
  - Iterator-fix (category 2): 0
  - test-allow (category 4): 0 — handled via the existing `cfg_attr(test, allow(...))` block

## Warning breakdown

| File | Lines | Type | Fix applied |
|------|-------|------|-------------|
| `potential/nonlocal.rs:292` | `VacantEntry::insert()` returns `&mut V` | binding | `let _entry = e.insert(...)` |
| `potential/xc.rs:1518` | `AtomicUsize::fetch_add()` debug counter | binding | `let _prev = ...` |
| `potential/xc.rs:1640` | same pattern for spin counter | binding | `let _prev = ...` |
| `scf/mixing/anderson.rs:201,202` | `Vec::remove(0)` evicts oldest history | binding | `let _evicted_in/res = ...` |
| `scf/mixing/broyden.rs:178,179` | same pattern | binding | `let _evicted_dv/df = ...` |

## EM recovery note (2026-04-21)

Original agent stalled at quality gate (stream-idle watchdog). EM found 6 bare `self.queue.submit()` calls in `gpu/mod.rs` that triggered `unused_results` under `--features gpu`. Fixed by adding `let _submission = ...` bindings at lines 436, 471, 552, 593, 655, 695. Commit amended. Quality gate re-run: all green. PR #186.

## Blocked / unfinished

Nothing. Quality gate passes clean.

## Tangential

- `atoms.rs` was listed in the proposal as having 2 annotations but had 0 at implementation time — likely removed by an intermediate PR.
- `must_use_candidate` clippy lint is confirmed not enabled in any `Cargo.toml` or `.clippy.toml`.
- Adding `unused_results` to the existing `cfg_attr(test, allow(...))` block is cleaner than per-module `#[allow]` in test code.
