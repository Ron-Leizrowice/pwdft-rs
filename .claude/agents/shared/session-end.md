# Session-end logbook

Every session produces exactly **one new logbook file** in `.claude/logbooks/<role>/` before ending.

## File naming and location

- Path: `.claude/logbooks/<your-role>/YYYY-MM-DD-<short-slug>.md` — e.g. `.claude/logbooks/core-engineer/2026-04-20-pzpw-rewire-correlation.md`.
- **Write it inside your worktree.** The file is a fresh creation, so the worktree-isolation hook permits it; the file lands as part of your PR. No paste-into-PR-body dance is needed.
- **EM sessions** running directly in the main checkout write to `.claude/logbooks/engineering-manager/` on the main checkout in place.

## File structure

```markdown
# <Role> — YYYY-MM-DD — <short title>

<what was done or decided, 1–3 sentences>

## Key numbers

<metrics, measurements, discrepancies — not prose>

## Blocked / unfinished

<anything the next session needs to pick up>

## Tangential

<one-liners worth keeping, if any>
```

Keep entries under ~30 lines. If a section grows longer, extract it into a proposal — logbooks are handoffs, not spec docs.

## Search before asking

At session start, search your role's directory for prior work on the same area:

```bash
rg -l <keyword> .claude/logbooks/<role>/
rg <keyword>    .claude/logbooks/<role>/
```

Cross-read other roles' directories when the topic spans them (physics question in a perf session → `rg <keyword> .claude/logbooks/researcher/`).

Also skim `history.md` in your role's dir — that's the consolidated pre-refactor flat logbook. Useful context, read-only.

## What to record

- What was done or decided
- What's blocked or unfinished
- Key numbers (metrics, measurements, discrepancies — not prose)
- Tangential ideas worth capturing (one line each)

## What not to record

- Detail that belongs in the PR body or a proposal.
- Prose summaries of the diff — `git log` does that.
- Duplicates of information that's already in an existing logbook entry (link instead).
