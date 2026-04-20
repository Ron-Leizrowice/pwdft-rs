# PR workflow — draft early, promote when ready

Applies to every sub-agent role that branches off `origin/main` and opens a PR: core-engineer, performance-engineer, researcher, technical-writer, code-reviewer.

## The rule

**Every atomized change that compiles gets committed, and the
first commit opens a draft PR.** There's no timer. "Atomized"
just means the smallest self-contained increment you'd be
comfortable leaving the branch at — the code compiles, the
worktree isn't mid-rewrite. It's the unit you'd stage + commit
anyway.

The draft PR gives you:

1. **A visible checkpoint** — the user (and EM) can eyeball
   your WIP without interrupting you.
2. **A backup** — if the worktree is removed, the agent
   crashes, or a hook destroys local state, every committed
   increment is on `origin` and recoverable.
3. **Early CI** — clippy + test (and eventually PYQE's QE
   reference regen) start burning down before the final
   `/pr-submit`.

## The cadence

1. **First atomized commit** → `/pr-draft`. The skill runs
   `cargo check` (auto-skipped if no `.rs` touched), commits
   what's unstaged, pushes, opens a draft PR.
2. **Each subsequent atomized commit** → `/pr-draft` again.
   Same command; the skill detects the existing PR and just
   pushes.
3. **When implementation is done** — quality-gate green,
   Tier-2 outcome in hand (if triggered), session logbook
   written — run `/pr-submit`. Promotes the draft to
   ready-for-review, updates the body, hands off to the EM.

What counts as "atomized" is your call. A file rename + its
import updates = one commit. A test added for the helper you
just wrote = one commit. A failed debugging detour you
immediately backed out = zero commits. Don't commit code that
doesn't build — if you're not sure whether it builds, run
`cargo check` first.

## Always `git add .` before committing

Never compose a partial-path `git add <path1> <path2>`
sequence. prek's pre-commit stash/restore dance corrupts
overlapping edits between staged and unstaged changes — the
fix is to leave nothing unstaged. Worktrees are isolated;
everything in yours belongs to this PR anyway, including agent
memory writes and logbook entries. The `/pr-draft` and
`/pr-submit` skills already do `git add .`; the rule is for
hand-crafted commits.

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
