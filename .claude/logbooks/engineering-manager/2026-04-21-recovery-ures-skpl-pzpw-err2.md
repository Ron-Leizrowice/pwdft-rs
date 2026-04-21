# 2026-04-21 — Recovery wave: URES/SKPL/PZPW + ERR2 refresh

## Context

Resumed from a prior context that dispatched PZPW, SKPL, URES, and MLRW agents in
parallel. All four stalled from the stream-idle watchdog before completing. Partial
work was preserved in their worktrees.

## Worktree audit

| Agent | Branch | State on arrival |
|-------|--------|-----------------|
| ae5acc28 | PZPW/slater-pw92-lda-correlation | Phase 0 uncommitted (pw92_correlation un-gated + LdaPw92 in settings.rs + spin test) |
| ab7bd2ea | SKPL/tier2-skip-list | All 4 parts done, uncommitted |
| a5be0f91 | URES/unused-results-lint | Full commit `de9fd75` — quality gate not run |
| a0147660 | worktree-agent-a0147660 | On WRONG branch; rewrote machine-lock Bash (Part D of MLRW) without Python library |
| a639cb1b | MLRW/reader-writer-machine-lock | Partial Python library (lock/__init__.py, paths.py, pid.py, state.py) untracked |

## ERR2 proposal refresh

Fresh production census on `e113379`: 18 `.expect()` sites (+3 vs prior 15 — three
Mutex poison sites in gpu/mod.rs, legitimately implausible). Added 2026-04-21 refresh
section to `proposals/ERR2-panic-free-production-audit.md`:
- User's standard: "only expect for extremely implausible failure modes, Results for
  genuinely fallible paths"
- P1.e: transplant.rs:121,130 → InvalidParam (2 missed sites, ~5 LOC)
- Phase 2.5: gpu/mod.rs:755 read_staging_buffer → Result<Vec<f32>, PwdftError>
- Next move: bundle P1.e + Phase 2 + Phase 2.5 into one small PR after current wave clears
- Frontmatter updated: status draft → active, complexity medium → small

## SKPL recovery

Prior agent had all changes uncommitted. EM:
1. Wrote logbook `.claude/logbooks/core-engineer/2026-04-21-skpl-tier2-skip-list.md`
2. Committed all 8 files as commit `cac6510` (pre-commit hook reformatted CLAUDE.md; staged + recommitted)
3. Verified consistency script: 9 tests enumerated correctly
4. Pushed branch + dispatched quality gate (background)
5. PR pending quality gate completion

## URES recovery

Prior agent had full commit `de9fd75` but quality gate not run. Recovery agent found
5 `gpu/mod.rs` warnings under `--features gpu` (bare `queue.submit()` calls).
EM applied the fix directly (6 occurrences, `let _submission =` binding pattern),
amended commit → `500dd0c`. Quality gate re-run: all green.
PR #186 created: https://github.com/Ron-Leizrowice/pwdft-rs/pull/186

## PZPW recovery

Prior agent had Phase 0 uncommitted (xc.rs + settings.rs, ~57 LOC). Phase 1+2 not
started. Dispatched new core-engineer agent with precise instructions:
- Phase 1: add lda_xc_pw92, lda_xc_spin_pw92, grid lifts, XcEvaluator::LdaPw92
- Phase 2: simplified Fe BCC diagnostic test (no vgch_transplant_fe.rs dependency)
- Tier-2 required (xc.rs, settings.rs on trigger list)

## MLRW status

- a0147660 rewrote machine-lock Bash shim (Part D) but on wrong branch
  (`worktree-agent-a0147660`). Changes not usable until Python library is complete.
- a639cb1b has partial Python library (4 files, all untracked) on the correct branch.
- Decision: MLRW deferred. The a0147660 Bash shim design is useful but dependent on
  the Python library completing. Will start fresh dispatch after current wave clears.

## CWD drift incident

After committing SKPL files with `cd ...agent-ab7bd2ea`, EM's CWD was set to the
SKPL worktree. The check-worktree.sh hook blocked subsequent Edit calls targeting
URES worktree. Fixed by `cd /main/checkout` in a Bash call. This is a known risk
(CWDL) — EM must use `cd` back to main after any worktree operation.

## Open items at this handoff

| Item | Status |
|------|--------|
| URES PR #186 | Open, quality gate green, awaiting review |
| SKPL PR | Quality gate running; PR to be created on completion |
| PZPW PR | Agent dispatched (Phase 1+2), awaiting completion |
| ERR2 mop-up (P1.e + Ph2 + Ph2.5) | Proposal updated; dispatch after wave clears |
| MLRW | Deferred; assess after URES/SKPL/PZPW clear |
