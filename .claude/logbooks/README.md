# Logbooks

Per-role logbook directories. Each directory holds the session-end notes for one role, one file per session.

## Layout

```text
.claude/logbooks/
  README.md                         ← you are here
  <role>/
    history.md                      ← pre-refactor consolidated entries (read-only archive)
    YYYY-MM-DD-<slug>.md            ← one file per session; chronological via filename
  archive/
    <role>-YYYY-MM.md               ← older entries rotated out
```

Roles: `engineering-manager`, `core-engineer`, `performance-engineer`, `researcher`, `code-reviewer`, `technical-writer`.

## Session-end convention

**Every session produces one new file** named `YYYY-MM-DD-<short-slug>.md` in the role's directory. Write it **inside your worktree** (not to the main checkout); the file lands as part of your PR. No paste-into-PR-body dance, no EM-copies-later handoff.

File structure:

```markdown
# <Role> — YYYY-MM-DD — <short title>

<what was done or decided>

## Key numbers

<metrics, measurements, discrepancies — not prose>

## Blocked / unfinished

<anything the next session needs to pick up>

## Tangential

<one-liners worth keeping, if any>
```

Keep entries under ~30 lines. If an entry grows past that, extract the long material into a proposal — logbooks are handoffs, not spec docs.

## Search before asking

At session start, **search your role's directory** for prior work on the same area:

```bash
rg -l <keyword> .claude/logbooks/<role>/         # filename match
rg <keyword> .claude/logbooks/<role>/            # content match
```

Grep the other roles' directories too when the topic crosses boundaries (e.g. a performance question that involves physics → search both `performance-engineer/` and `researcher/`).

## Why per-session files

- **No more rebase conflicts.** Each PR writes a unique file; concurrent PRs don't collide on a shared `<role>.md`.
- **Logbooks land with PRs.** Sub-agents in a worktree can write to their worktree copy of the logbook tree without tripping the main-checkout hook, because they're creating a new file, not modifying a shared one.
- **Clean attribution.** Filename encodes date and topic; git log per file gives session provenance.
- **Cheap archiving.** Old months can be moved into `archive/` without rewriting a big flat file.

## History before the refactor

`<role>/history.md` holds the pre-refactor flat logbook for each role. Read-only — new entries always get their own file.

Older monthly rollups live in `archive/`. When a role's directory accumulates more than ~20 dated files, the EM rotates the oldest into a monthly file there.
