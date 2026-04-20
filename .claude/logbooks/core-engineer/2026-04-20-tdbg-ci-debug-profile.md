# TDBG: CI Tier-1 cargo test --profile=dev

**Date:** 2026-04-20
**Proposal:** `proposals/TDBG-tier1-debug-build-in-ci.md`
**Branch:** `TDBG/tier1-debug-build-in-ci`

## Scope landed

- Phase A: `.github/workflows/rust.yml` Tier-1 step switched from `cargo test -p pwdft-core` to `cargo test --profile=dev -p pwdft-core`. Comment block added explaining why (Tier-1 is structurally compile-bound on cold runners; `#[ignore]`-gated SCF loops live in Tier-2, which stays on O3).
- Phase B: `Swatinem/rust-cache@v2` step gets `shared-key: pwdft-ci-tier1-dev` so the debug-profile `target/` can't collide with any future release-profile cache (PYQE Phase D nightly, PMTL workspace split).
- Phase C: CLAUDE.md § "Tests & tiers" — one-paragraph "Profile split: local vs. CI" note pinning the divergence. Baseline run `24663826105` (5 m 43 s = 343 s for the test step) and 4-minute (240 s) budget both cited.
- Phase D: Guardrail wall-timer inline in the single test step; fails the step if wall > 240 s.

Explicit non-goal: **TPRF's completed proposal was not prepended with a "Superseded-in-CI-context" note.** Per the task instructions, that sub-step was conditional on having measured CI numbers from a live PR run. Deferring until the PR run lands. CLAUDE.md text is enough auditability for now.

## Measurements

### Baseline (main, run 24663826105, 2026-04-20 push of #180)

- `clippy + test` job total: **10 m 43 s** (11:24:17 → 11:35:00).
- `cargo test (tier 1)` step alone: **5 m 43 s** (343 s), 11:29:02 → 11:34:45.
- Clippy default + GPU combined: 4 m 26 s.

### Post-change

Filled in after the first CI run on this branch lands. See PR body.

## Surprises / friction

1. **Pre-existing rustdoc failure on origin/main.** `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` fails on vendored faer `pwdft/faer/faer/src/mat/mod.rs:146` — an empty triple-backtick block that rustdoc treats as an invalid Rust code block. Reproduced on a clean `git stash`. Not introduced by TDBG. Flagging it for a separate tiny proposal (annotate as ```` ```text ```` per rustdoc's own hint).
2. Quality gate ran successfully for clippy (default + gpu + auto-fix) and Tier-1 `cargo test` (all green). The failure above only hits the rustdoc step. Since TDBG's diff is workflow YAML + CLAUDE.md prose only, no new doctests were introduced — the rustdoc failure is orthogonal.
3. **Worktree isolation hook fired once** when I first tried to edit `.github/workflows/rust.yml` at the main-checkout path instead of the worktree copy. Fixed by using the worktree path. Noting here so the next Core Engineer remembers that path resolution must be inside `.claude/worktrees/agent-*/...`.

## FLUP

- **Rustdoc empty-codeblock in vendored faer** (not a TDBG blocker; tiny proposal). One-line fix: change ` /// ``` ` to ` /// ```text ` at `pwdft/faer/faer/src/mat/mod.rs:146`. Would unblock the EM's standard `/quality-gate` on any PR that didn't hit the preexisting cache.
- **If TDBG's measured win is <30%**, escalate to the proposal's option C (scope `profile.test` to the workspace member only). To be decided from the PR's CI run.
