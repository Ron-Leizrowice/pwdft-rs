# PR workflow — draft early, promote when ready

Applies to every sub-agent role that branches off `origin/main` and opens a PR: core-engineer, performance-engineer, researcher, technical-writer, code-reviewer.

## The rule

**Open a draft PR as soon as you have your first file on disk.** Don't wait until the whole implementation is done. The draft PR is:

1. **A visible checkpoint** — the user (and EM) can eyeball your WIP without interrupting you.
2. **A backup** — if the worktree gets corrupted, force-removed, or the agent session crashes, your work is on `origin` and recoverable. Without the draft, a wiped worktree = lost work.
3. **The place CI runs early** — long-running CI (clippy + test, eventually PYQE's QE reference regen) starts burning down before the final `/pr-submit`, so end-of-work latency is lower.

## The cadence

1. **First file change that compiles** → `/pr-draft`. The skill runs `cargo check` (auto-skipped if no `.rs` touched), commits, pushes, opens a draft PR.
2. **After each logical increment** (a test added, a migration finished, a phase of a multi-phase proposal done) → `/pr-draft` again. Same command; the skill detects the existing PR and just pushes.
3. **When implementation is done** — quality-gate green, Tier-2 outcome in hand (if triggered), session logbook written — run `/pr-submit`. This promotes the draft to ready-for-review, updates the PR body with the final Test Plan, and hands off to the EM.

## What `/pr-draft` runs

- `cargo check -q` — only when the diff touches `.rs` files. For proposal-only, doc-only, or fixture-only PRs the check is skipped.
- `git add .` on your worktree root, commit with a `<ID>: checkpoint — <msg>` message. The worktree is isolated, so every file in it belongs to this PR — memories, logbooks, and scratch included.
- `git push --force-with-lease` — safe for drafts; subsequent checkpoints rewrite freely.
- `gh pr create --draft` (first time) or no-op (subsequent calls — the push updates the existing PR server-side).

It does **not** run clippy, rustdoc, tests, or Tier-2. Those belong to `/pr-submit`.

## What `/pr-submit` runs

The full quality gate + rebase on `origin/main` + promote-or-create PR. See `.claude/skills/pr-submit/SKILL.md` for the step list. If a draft PR exists on your branch, `/pr-submit` promotes it via `gh pr ready` instead of opening a duplicate.

## When NOT to use `/pr-draft`

- **Nothing committed yet** — `/pr-draft` without changes is a no-op; you don't need to pre-open an empty PR.
- **Reviewer comments** (code-reviewer role) — comments go via `gh pr review` on someone else's PR. If the code-reviewer role is *filing* a PR (audit finding, doc cleanup), the draft-early rule applies like any other role.
- **Experiments you intend to abandon** — if you're scoping a proposal and may throw away the exploratory branch, don't open a PR. Write the proposal file first, open the PR when you commit to implementing.

## Recovery when a worktree dies

If your worktree is gone and you have a draft PR:

```bash
gh pr checkout <N>                                  # from the main checkout
cp -r <files you need> <new worktree>               # or just keep working in the main checkout; the branch is local now
# continue; next /pr-draft or /pr-submit uses the new worktree
```

Without a draft PR, your only recovery is local branches / reflog, which may have been pruned when the worktree was removed. Draft PRs are the durable backup.
