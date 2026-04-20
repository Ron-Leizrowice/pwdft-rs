# ROTI: revert SpaceGroupOp rotation i8 → i32 — recovery

**Date:** 2026-04-21
**Proposal:** `proposals/ROTI-revert-symm-rotation-i8.md`
**Branch:** `ROTI/revert-symm-rotation-i8`
**PR:** #184

## Scope inherited

Previous core-engineer agent (`a7fe9ca547748e383`) implemented the proposal's Steps 1-7 and pushed checkpoint `8b6fbbc` titled "ROTI: WIP checkpoint from aborted agent run (2026-04-20)". CI passed (SUCCESS, mergeable), but no `/quality-gate` or `/test --tier2` run was reported, the PR body was a recovery note rather than a proper Summary/Test Plan, and the title still carried "WIP". My job was to validate + promote out of WIP.

Diff (confirmed unchanged on inspection): +48 / −151 LOC across

- `pwdft/pwdft-core/src/symmetry/operations.rs`
- `pwdft/pwdft-core/src/symmetry/detect.rs`
- `pwdft/pwdft-core/src/symmetry/density/g_space.rs`
- `pwdft/pwdft-core/src/symmetry/density/mod.rs`

Scope is exactly the i8→i32 revert + removal of the `#[expect]`/`try_from` noise. No creep.

## Rebase

`git rebase origin/main` was a clean fast-forward equivalent — no conflicts (the recent main advance — STYS proposal, gitignore cleanup, no-backcompat policy, TDBG CI — doesn't touch symmetry or the vendored faer files that ROTI is concerned with).

## Validation

### Quality gate (194 s, M3 Max)

All five steps green:

- `cargo clippy -q --fix --allow-dirty --allow-staged --all-targets` — no fixes applied
- `cargo clippy -q --all-targets` — 0 warnings
- `cargo clippy -q --all-targets --features gpu` — 0 warnings
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` — green **after** a one-line fix to vendored `pwdft/faer/faer/src/mat/mod.rs` (see FLUP #1)
- `cargo test` (Tier-1) — pass

### Tier-2 (symmetry is a trigger per CLAUDE.md)

Ran `cargo test --no-fail-fast -- --ignored`. Outcome recorded in the PR body. The authoritative designed-to-fail skip list in the `#[ignore]` reason strings (C-diamond Plain Anderson, Al ecut, VGCH heavy-atom cells, MXBA) behaved exactly as documented; no *new* failures attributable to this diff.

## Vendored faer rustdoc fix

The `/quality-gate` `cargo doc --no-deps` step was failing on `pwdft/faer/faer/src/mat/mod.rs:146` — a pre-existing malformed doc block where `[\`reborrow::Reborrow\`] \`\`\`` sat on one line, leaving rustdoc to close the previous fence and open an empty one. The fix is a single prose edit (split the inline fence onto its own line). The error is present on `origin/main` at `2b2487f` (and current `3a5ca46`) — not caused by ROTI. CI doesn't catch it because `.github/workflows/rust.yml` runs `cargo clippy -p pwdft-core` and never invokes `cargo doc`. Shipping the faer fix in this PR so the local gate is green; CLAUDE.md § Vendored dependencies explicitly sanctions editing vendored faer on a feature branch.

## Flagged for follow-up

1. **CI doesn't run `cargo doc -D warnings`.** The vendored faer rustdoc break was on main for at least one merge cycle because CI doesn't exercise the rustdoc step that `/quality-gate` runs. A cheap workflow addition (`RUSTDOCFLAGS="-D warnings" cargo doc --no-deps -p pwdft-core` — note the `-p pwdft-core` scope to match local intent) would catch future regressions without re-triggering this same vendored-faer issue. Worth a small proposal.
2. **`/quality-gate` scope is wider than CI.** `cargo doc --no-deps` (unscoped) documents every workspace member including vendored faer, while CI only clippies `-p pwdft-core`. This divergence means an innocent PR can get blocked by vendored-faer prose. Either narrow `/quality-gate` to `-p pwdft-core` (mirror CI) or widen CI to match local gate. The vendored-faer `mat/mod.rs` prose was not the only sloppy rustdoc in that tree — a full rustdoc sweep of vendored faer likely turns up more.
3. **CLAUDE.md cache-line math never fully audited.** The ROTI proposal found one comment off by 15× in `operations.rs`. Worth a one-shot sweep by the technical-writer for other layout/alignment claims in symmetry/ and basis/.

None of the above are in ROTI's scope.
