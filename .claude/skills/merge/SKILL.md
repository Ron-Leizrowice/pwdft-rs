---
name: merge
description: Engineering Manager only. Merge an approved PR, remove the sub-agent worktree, delete the local branch, and bring the main checkout up to date. Wraps the full merge trilogy. Trigger on "/merge <PR#>".
user_invocable: true
---

# /merge — Engineering Manager merge trilogy

Run the squash-merge + worktree cleanup + local branch delete + main rebase sequence for PR number `$ARGUMENTS`. EM-only — other roles should not invoke this.

## Prereqs

- You are acting as the Engineering Manager.
- The PR has been reviewed and approved against the § PR review checklist in `.claude/agents/engineering-manager.md`.
- You know the absolute path of the main checkout (call it `MAIN_CHECKOUT`). This is typically `/Users/<user>/…/pwdft-rs` — the path without any `.claude/worktrees/agent-*` segment. If you're running inside a worktree, it's the worktree's parent-of-parent-of-parent.

## Why this order matters

GitHub's `gh pr merge --delete-branch` tries to delete the local branch right after merging. If the sub-agent's worktree still has that branch checked out, the local delete fails with "cannot delete branch used by worktree". So the correct order is: merge-remote first, then remove the worktree, then delete the local branch. The repo has `deleteBranchOnMerge: true`, so the **remote** branch auto-deletes regardless of whether we pass `--delete-branch` — we drop the flag and handle the local side ourselves.

Similarly, `pull --ff-only` assumes local main is behind origin. In practice the user makes direct-to-main commits (CI renames, agent-def tweaks) that leave local ahead. Use `fetch + rebase origin/main` instead — it handles both ahead and behind cases. Stash/pop around the rebase if the main checkout has uncommitted WIP.

## Steps

1. **Verify the PR state.**

   ```bash
   gh pr view "$ARGUMENTS" --json state,mergeable,title,headRefName,author
   ```

   Confirm `state=OPEN` and `mergeable=MERGEABLE`. If not mergeable, stop and report — do not force-merge.

2. **Identify the sub-agent worktree.**

   ```bash
   git worktree list
   ```

   Find the row whose branch matches the PR's `headRefName`. Capture that worktree path; you'll remove it in step 4. If no matching worktree shows up, the branch may have been created outside the worktree protocol — skip step 4 and go straight to the local delete in step 5.

3. **Squash-merge.** Drop `--delete-branch` — the remote deletes on its own via `deleteBranchOnMerge: true`, and we handle the local side after the worktree is gone.

   ```bash
   gh pr merge "$ARGUMENTS" --squash
   ```

4. **Remove the sub-agent's worktree.** Use the path captured in step 2.

   ```bash
   git worktree remove -f -f .claude/worktrees/agent-<...>
   ```

   Double `-f` is needed for worktrees that were never pushed or have dirty state left behind. If the agent reported leftover work, investigate before forcing — see Errors below.

5. **Delete the local branch.** Safe now that nothing holds it.

   ```bash
   git -C "$MAIN_CHECKOUT" branch -D <head-ref-name>
   ```

   Use the exact branch name from step 1's `headRefName`.

6. **Bring the main checkout up to date.**

   Run from the main checkout's absolute path, not `$(pwd)` — if you merged from a worktree, `$(pwd)` may now be a dangling directory. Use rebase, not ff-only, because local main may be ahead of origin (direct-to-main commits from the user).

   If the main checkout has uncommitted WIP, stash it first:

   ```bash
   git -C "$MAIN_CHECKOUT" status --short
   # If anything shows:
   git -C "$MAIN_CHECKOUT" stash push -m "EM merge-trilogy stash" -- <modified paths>
   ```

   Then:

   ```bash
   git -C "$MAIN_CHECKOUT" fetch origin
   git -C "$MAIN_CHECKOUT" rebase origin/main
   ```

   If you stashed, pop it:

   ```bash
   git -C "$MAIN_CHECKOUT" stash pop
   ```

   Do **not** push. Direct-to-main commits on local main are the user's territory; the EM never pushes main.

7. **Archive the proposal if it's now complete.**

   Move `proposals/<ID>-*.md` to `proposals/completed/<ID>-*.md`, flip the frontmatter `status` to `completed`, remove the row from `proposals/INDEX.md`, and append a one-line entry to `.claude/logbooks/engineering-manager/shipping-log.md` under the appropriate theme bucket. The `/proposal complete <ID>` skill wraps this. INDEX stays concise — it only lists ongoing work.

   Skip this step for proposal-drafting PRs (the landing creates new proposals; it doesn't complete existing ones) and for PRs that only partially address a multi-phase proposal.

8. **Write a merge logbook entry** at `.claude/logbooks/engineering-manager/YYYY-MM-DD-<pr-slug>.md` noting:

   - PR number + merged SHA.
   - Any stash/pop or rebase conflict handling you did in step 6.
   - Any archival done in step 7.
   - Any FLUP items from the sub-agent's return message — whether you actioned them inline, deferred them, or filed them as new stub proposals.

9. **Report** the merged SHA, INDEX/shipping-log changes, and any new stub proposals to the user.

## Errors

- **`gh pr merge` fails with "not mergeable"** → check for conflicts, failing CI, or missing approvals. Fix at source; do not `--admin`-override without explicit user direction.
- **`git worktree remove` fails with "dirty worktree"** → the sub-agent may have uncommitted scratch. Check `git -C <worktree-path> status` before forcing a second `-f`. If there's real work (a test file, a script, a partial commit), ask the user what to do. Do not discard unknown work silently.
- **`git branch -D` fails with "checked out at ..."** → another worktree still holds the branch (maybe a cousin agent). Remove that worktree first, or leave the branch for the user to clean up.
- **`git rebase origin/main` fails with conflicts** → the user's direct-to-main commits collided with a merged PR. Abort with `git rebase --abort`, pop any stash, and report to the user. Don't improvise a conflict resolution — direct-to-main commits are the user's, not yours to rewrite.
- **`git stash pop` fails with conflicts** → the rebase moved the base in a way that collides with the user's WIP. Report and stop; let the user resolve.

## See also

- `.claude/agents/engineering-manager.md` § The merge trilogy
- `.claude/skills/proposal/SKILL.md` § `complete <id>`
