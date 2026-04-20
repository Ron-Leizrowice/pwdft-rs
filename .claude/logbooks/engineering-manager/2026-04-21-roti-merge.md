# 2026-04-21 — ROTI merge + FLUP seeds + MLRW/MLDX proposal drafts

## Shipped

- **PR #184 ROTI merged.** Squash SHA `43ef482`. `+48 / −151` LOC on the symmetry diff plus one incidental faer rustdoc fix (`pwdft/faer/faer/src/mat/mod.rs:135` — malformed inline code fence; pre-existing on `origin/main`, caught by local `/quality-gate` but not by CI). Tier-2 wall 19 min cold cache on the sub-agent's M3 Max, 14 pass / 6 fail, all 6 failures on the authoritative designed-to-fail list in the `#[ignore]` reason strings (measured ΔE matched predictions to 3–4 sig figs across Fe LDA/PBE, GaAs, Cu, MgO, C diamond).
- Proposal archived to `proposals/completed/ROTI-revert-symm-rotation-i8.md`; frontmatter flipped to `status: completed`; INDEX row removed; shipping-log Type-audit bucket extended (`TYPE, TYPB, CAST → TYPE, TYPB, CAST, ROTI`); cumulative PR count 73 → 74.

## In flight at session end

- **PR #185 ESPL.** Resumed via second core-engineer agent (`add8fe19`) after the first stalled on worktree-isolation confusion (`aae36465`, stalled at orient step). Second agent rebased, worked in its assigned worktree, still running at session end — last observed with uncommitted `pwdft/pwdft-core/src/settings.rs` edits, lock free, no active cargo processes. Tier-2 run still pending on that branch.
- **User pushed local main to origin mid-session.** Local main had been 6 commits ahead (STYS proposal, gitignore/pr-workflow tightening, no-backcompat policy + ESPL `spin_polarized` amendment, TDBG archive, agent-config updates). After push, `origin/main` ≡ local main ≡ `3a5ca46`. The amended ESPL proposal + `shared/no-backcompat.md` are now reachable on the PR base.

## Recovery pattern — sub-agent stalls on worktree orientation

Both the first ESPL agent (`aae36465`) and first ROTI agent (`a2dd3c91`) stalled at the 600 s watchdog threshold while still trying to figure out their worktree assignment. Both started in a placeholder worktree (`worktree-agent-<id>`) rather than on the PR branch, and the stall happened between "orient" and "cd/reset to branch." The retry-pair (`add8fe19`, `a8f56f3c`) got tighter prompts explicitly forbidding `git worktree add`/`switch branches`, instructing `git reset --hard origin/<branch>` + `git rebase origin/main` from within the assigned worktree. Both retries proceeded past orient on the first try. Pattern worth capturing if it recurs: sub-agents fight their own worktree isolation when told to "work on PR branch X."

## Worktree cleanup

Removed: `agent-a2dd3c91`, `agent-aae36465` (stalled placeholders, cleaned mid-session), `agent-a7fe9ca5` (ROTI WIP orphan, held the branch after merge — blocked local `branch -D` until removed), `agent-a8f56f3c` (successful ROTI recovery, merge trilogy).

Still on disk, all orphan ESPL attempts from the original aborted run: `agent-a8d588f7` (ESPL-v2 branch), `agent-ad8d4aaf` (ESPL-final), `agent-af025721` (ESPL original), `agent-espl-recovery` (stalled-agent's partial rebase branch). These are safe to remove once the active ESPL agent (`add8fe19`) returns with its PR green. Defer cleanup until then to avoid racing the live agent.

## FLUP seeds this session

- **CIDO** — CI rustdoc gate + `/quality-gate` scope harmonization. CI doesn't run rustdoc, and `/quality-gate` runs rustdoc workspace-wide so benign pwdft-core PRs can be blocked on vendored-faer prose. The ROTI agent had to fix faer's `mat/mod.rs:135` to make the gate pass — exactly the drift this seed captures. Recommended fix: add `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps -p pwdft-core` to `rust.yml`.
- **SYMC** — `symmetry/` / `basis.rs` / mixer layout-claim sweep. ROTI diff found one layout claim in `symmetry/operations.rs` off by ~15×. One-shot technical-writer sweep.
- **CWDL** — EM cwd sanity check against worktree drift. Session hit a `check-worktree.sh` denial on an Edit to `proposals/FLUP...` because the EM shell's cwd had silently drifted into `agent-a8f56f3c`. Cousin of the 2026-04-20 memory-path drift incident. Recommended: both a `pwd` assertion at EM session-start and explicit `cd "$MAIN_REPO"` preambles in `/merge`/`/pr-review`.

## Drafted

- **MLRW** (`proposals/MLRW-machine-lock-reader-writer.md`) — reader-writer machine lock in `pwdft_validation.lock` (Python library primitive + bash shim). Exclusive for bench/profile/QE; shared for test/clippy/build/doc; `exclusive_pending` flag for priority-promotion anti-starvation. State file moves from per-field directory to single JSON so the bash hook reads it via `python3 -c`. Extends hook-matched commands to include `cargo doc` and all QE binaries (pre-existing drift). Motivated directly by the current session — two concurrent recovery agents spent ~22 min serialized on the exclusive lock when both cargo commands (test + doc) would have been shared-class. Medium complexity, medium risk. Blocks MLDX.
- **MLDX** (`proposals/MLDX-machine-lock-daemon-dispatch.md`) — FastAPI daemon layer (Phase 1) + heavy-command dispatch (Phase 2). Depends on MLRW. Phase 1 motivations: sub-ms acquire (vs. ~500 ms–1 s `uv run`), event-driven waits, heartbeat reclaim. Phase 2 motivations: bench env reproducibility, wall-time history persistence. Explicit library-fallback contract: daemon down ⇒ CLI reverts to MLRW direct, no functional difference. CI uses `PWDFT_LOCK_DIRECT=1` to skip daemon entirely. Large complexity, medium-high risk — gated on MLRW stability.

## Main state at handoff

`origin/main = 43ef482b ROTI: revert SpaceGroupOp rotation i8 → i32 (#184)`. Local main ≡ origin/main post-rebase. Uncommitted WIP popped from merge stash: modified INDEX.md + FLUP; untracked MLRW/MLDX proposal files + `.claude/scheduled_tasks.lock` (harness state, not for commit). Next user action likely: review MLRW/MLDX, approve or request changes. ESPL PR #185 to land next — waiting for the active agent to finish.
