---
name: User occasionally pushes to main; agents must sync before work
description: The user makes direct-to-main commits (CI tweaks, agent-def edits) that may be pushed between agent sessions. Sub-agents must fetch + rebase before starting, not assume their local view of origin/main is fresh.
type: feedback
---

The user occasionally pushes directly to `main` outside the PR flow — CI workflow tweaks, agent-def edits, small config changes. These appear as unexpected commits on `origin/main` that sub-agents started behind on.

**Why:** Confirmed by user 2026-04-20 after an EM session hit a "Not possible to fast-forward" during the merge trilogy — local main had user's direct-to-main commits and `origin/main` had both those commits (pushed) plus a freshly-merged PR. Combined divergence.

**How to apply:** When spawning any sub-agent that will branch from `origin/main` (which is all of them), the briefing must instruct the agent to `git fetch origin && git checkout origin/main` (or `git pull --rebase origin main` on the worktree's main tracking branch) before creating their feature branch. Don't assume the worktree's view of `origin/main` matches the true remote.

At merge time, the EM's `/merge` trilogy already uses `fetch + rebase origin/main` (not `pull --ff-only`) to handle the same divergence on the main checkout. See `.claude/skills/merge/SKILL.md` § Step 6.

For session-start (EM): run `git -C <MAIN> fetch origin && git -C <MAIN> status -sb` as part of orientation so the EM knows whether there's an unseen push to reconcile.
