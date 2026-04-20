---
name: pr-draft
description: Push the current branch to origin and open (or update) a draft PR. Used as a mid-implementation checkpoint — runs only `cargo check`, not the full quality gate. Trigger on "/pr-draft".
user_invocable: true
---

# /pr-draft — Early checkpoint PR

Commit, push, and open (or update) a **draft** PR as soon as the implementation compiles. This is the safety checkpoint against worktree loss and the first handoff surface the user can eyeball. It is **not** the end-of-work submission — that's `/pr-submit`.

Use this as soon as your first file change compiles. Re-run it after each logical increment. Convert the draft to ready with `/pr-submit` when implementation is done + quality gate green.

## Steps

1. **Verify worktree + branch.** `pwd` must be under `.claude/worktrees/agent-*`; current branch must be `<PROPOSAL-ID>/<slug>`. Abort if on `main` or a stale branch.

2. **Cheap sanity check — Rust-only.** If the diff touches any `.rs` file, run `cargo check -q` under the machine lock. Skip the check entirely if the diff is proposal-only, doc-only, fixture-only, or workflow-only — `cargo check` is wasted wall-time on a non-Rust diff, and this skill is meant to be cheap.

   ```bash
   if git -C "$(pwd)" diff --name-only HEAD 2>/dev/null | grep -q '\.rs$' || \
      git -C "$(pwd)" ls-files --others --exclude-standard 2>/dev/null | grep -q '\.rs$'; then
     .claude/bin/machine-lock run "<role>" "pr-draft cargo check" -- cargo check -q
   fi
   ```

   If `cargo check` fails, stop and fix; a broken-compile commit is not useful as a checkpoint.

3. **Stage everything in the worktree and commit.** `git add .` — the worktree is isolated, so whatever is in it belongs to this PR: primary code changes, logbook entries, agent-memory writes, test fixtures. No need to enumerate paths; worktree isolation already scopes the sweep. Commit with a WIP message — either a summary from the caller argument, or `<ID>: checkpoint` if none:

   ```bash
   git -C "$(pwd)" add .
   git -C "$(pwd)" commit -m "<ID>: checkpoint — <summary or auto-generated>"
   ```

   Multiple small checkpoints are fine — the final `/pr-submit` does not squash locally; the EM's merge squashes at merge time anyway.

4. **Push with `--force-with-lease`.** Drafts get rewritten often; the lease protects against overwriting a concurrent push you didn't see.

   ```bash
   git -C "$(pwd)" push --force-with-lease -u origin "$(git -C "$(pwd)" branch --show-current)"
   ```

5. **Open or update the draft PR.**

   ```bash
   existing=$(gh pr list --head "$(git -C "$(pwd)" branch --show-current)" --json number --jq '.[0].number')
   ```

   - **If `$existing` is empty** (first `/pr-draft` on this branch): open a draft.

     ```bash
     gh pr create --draft --title "<ID>: <description> (draft)" --body "$(cat <<'EOF'
     ## Draft — do not merge

     Checkpoint PR opened by `/pr-draft`. Implementation in flight; `/pr-submit` will convert to ready when the quality gate + Tier-2 (if applicable) pass and the session logbook is written.

     ## Proposal

     <ID>: <title>

     ## Progress

     - <what compiles so far>
     EOF
     )"
     ```

   - **If `$existing` is non-empty**: the push in step 4 already updated the PR server-side. No further action needed unless you want to edit the body with new progress:

     ```bash
     gh pr edit "$existing" --body "<updated body with new Progress bullets>"
     ```

6. **Report** the PR URL, number, and "Draft — checkpoint only, not ready for review" to the user. The EM should not review a draft PR against the merge checklist — that happens after `/pr-submit`.

## Errors

- **`cargo check` fails** → stop; fix the compile error. Do not push a broken branch even as a draft — it defeats the point of the checkpoint being a safe recovery point.
- **Nothing to commit** → print "no checkpointable changes" and exit cleanly. Re-running `/pr-draft` after a prior checkpoint with no new changes is a no-op, not an error.
- **`gh pr create` fails because a PR already exists** → the `gh pr list --head` check in step 5 should catch this; if it slipped through, run `gh pr list --state all --head <branch>` to inspect and route to the `gh pr edit` path.

## See also

- `.claude/skills/pr-submit/SKILL.md` — end-of-work skill; promotes the draft to ready after full quality gate.
- `.claude/agents/core-engineer.md` § Implementing an approved proposal — when to call `/pr-draft` vs `/pr-submit`.
