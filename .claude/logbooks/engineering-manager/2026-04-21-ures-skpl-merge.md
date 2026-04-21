# 2026-04-21 — Merge: URES (#186) + SKPL (#187)

## URES PR #186

- Squash-merged: `aeed6c1`
- Worktree `agent-a5be0f91` removed; branch `URES/unused-results-lint` deleted
- Main rebased onto `origin/main` — clean
- Stash complications: `cargo clippy --fix` side-effects on `basis.rs` and `ewald.rs`
  in the main checkout collided with URES changes during stash pop; both stash entries
  dropped (changes already subsumed by URES merge or minor style only)
- Proposal archived: `proposals/completed/URES-unused-results-lint.md`

## SKPL PR #187

- Squash-merged: `e4e4216`
- Worktree `agent-ab7bd2ea` removed; branch `SKPL/tier2-skip-list` deleted
- Main rebased onto `origin/main` — clean
- Proposal archived: `proposals/completed/SKPL-tier2-skip-list-automation.md`

## INDEX + shipping-log

- Removed URES and SKPL rows from `proposals/INDEX.md`
- Shipping-log counter: 75 → 77 PRs
- Added URES + SKPL one-liners under Infra + quality

## Open items

| Item | Status |
|------|--------|
| PZPW | Worktree `agent-ae5acc28` has Phase 0 uncommitted (xc.rs, settings.rs); no PR yet; needs fresh dispatch |
| MLRW | Two stalled worktrees (a0147660 wrong branch, a639cb1b partial Python library); deferred |
| ERR2 P1.e + Phase 2 + 2.5 | Ready to dispatch after PZPW clears |
