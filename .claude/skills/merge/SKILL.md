---
name: merge
description: Engineering Manager only. Merge an approved PR, delete its branch, remove the worktree, and fast-forward the main checkout. Wraps the three-step merge trilogy. Trigger on "/merge <PR#>".
user_invocable: true
---

# /merge — Engineering Manager merge trilogy

Run the squash-merge + branch delete + worktree cleanup + main fast-forward sequence for PR number `$ARGUMENTS`. EM-only — other roles should not invoke this.

## Prereqs

- You are acting as the Engineering Manager.
- The PR has been reviewed and approved against the § PR review checklist in `.claude/agents/engineering-manager.md`.
- You know the absolute path of the main checkout (call it `MAIN_CHECKOUT`). This is typically `/Users/<user>/…/pwdft-rs` without any `.claude/worktrees/agent-*` segment.

## Steps

1. **Verify the PR state.**

   ```bash
   gh pr view "$ARGUMENTS" --json state,mergeable,title,headRefName,author
   ```

   Confirm `state=OPEN` and `mergeable=MERGEABLE`. If not mergeable, stop and report — do not force-merge.

2. **Squash-merge and delete the remote branch** (the repo has `deleteBranchOnMerge: true`, so this is the canonical shape):

   ```bash
   gh pr merge "$ARGUMENTS" --squash --delete-branch
   ```

3. **Remove the sub-agent's worktree.** Find it by matching the PR's head ref to `git worktree list` output:

   ```bash
   git worktree list
   git worktree remove -f -f .claude/worktrees/agent-<...>
   ```

   Double `-f` is needed for worktrees that were never pushed or have dirty state left behind.

4. **Fast-forward the main checkout.** Run from the main checkout's absolute path, **not** from `$(pwd)` — if you `cd`-ed into a now-removed worktree, `$(pwd)` is a dangling directory.

   ```bash
   git -C "$MAIN_CHECKOUT" pull --ff-only origin main
   ```

5. **Update INDEX.md** if the proposal is now complete. Move the proposal row from Active to Completed; archive the proposal file to `proposals/completed/`. The `/proposal complete <ID>` skill wraps this.

6. **Append a logbook entry** to `.claude/logbooks/engineering-manager.md` noting the merge, any INDEX updates, and any FLUP seeds captured from the PR.

7. **Report** the merged SHA, INDEX changes, and any new stub proposals to the user.

## Errors

- `gh pr merge` fails with "not mergeable" → check for conflicts, failing CI, or missing approvals. Fix at source; do not `--admin`-override without explicit user direction.
- `git worktree remove` fails with "dirty worktree" → the sub-agent may have uncommitted scratch. Investigate before forcing (extra `-f`); don't discard unknown work silently.
- `pull --ff-only` fails with "non-fast-forward" → someone pushed to main concurrently. `git fetch`, inspect, and either accept the divergence or retry.

## See also

- `.claude/agents/engineering-manager.md` § The merge trilogy
- `.claude/skills/proposal/SKILL.md` § `complete <id>`
