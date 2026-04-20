---
name: worktree-start
description: Enter a worktree and start a new feature branch from origin/main. Use at the start of every implementation session to verify isolation and create the branch. Trigger on "/worktree-start PROP/slug".
user_invocable: true
---

# /worktree-start — Enter worktree and branch from origin/main

Run the standard worktree-entry check + branch creation for `$ARGUMENTS` (a branch name in `<PROPOSAL-ID>/<slug>` form).

## Steps

1. **Verify isolation.** Run:

   ```bash
   pwd
   git worktree list
   ```

   `pwd` must resolve to something under `.claude/worktrees/agent-*`. If it is the main checkout (`/Users/.../pwdft-rs` without a `.claude/worktrees/agent-*` segment), STOP and report a harness failure — the user needs to spawn you with `isolation: "worktree"` or enter a worktree first.

2. **Parse the argument.** `$ARGUMENTS` should look like `PROP/some-slug`. Reject empty or whitespace-only input. Do not silently default the branch name.

3. **Fetch and branch from `origin/main`** (not local `main`, which can lag):

   ```bash
   git -C "$(pwd)" fetch origin
   git -C "$(pwd)" checkout -b <$ARGUMENTS> origin/main
   ```

4. **Report** the branch name, current `HEAD` SHA, and the worktree path back to the user.

## Errors

- Branch already exists → report and ask the user whether to switch to it (`git checkout <branch>`) or pick a new name.
- `origin/main` missing → run `git remote -v` and tell the user; don't guess the remote name.
- Not in a worktree → hand back the harness-failure message; do not try to `EnterWorktree` on the user's behalf.

## See also

- `.claude/agents/shared/worktree.md` — full worktree-isolation protocol.
