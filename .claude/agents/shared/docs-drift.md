# Documentation drift protocol

If during your session you encounter a reference to a function, module, type, file path, or repository layout that no longer matches the code, **do not leave it uncorrected.** A wrong docstring or a stale path reference is worse than no documentation — the next reader will act on it.

## Two paths, depending on scope

- **Within your PR's scope:** fix the reference inline as part of your current PR. A one-line docstring correction or a stale path edit does not count as scope creep.
- **Outside your PR's scope** (the fix would touch files unrelated to your proposal, or would widen the diff meaningfully): add a line to your `## Flagged for follow-up` block with target role **Technical Writer** and a specific file:line citation. The EM will triage it.

## Where drift shows up most

- Function / type names that were renamed but still appear in docstrings, `README.md`, `CLAUDE.md`, module-level `//!` headers, or proposal text.
- Repository paths that moved during a refactor but are still cited in comments or in the skills (`.claude/skills/*/SKILL.md`).
- Equation citations where the implementation was updated but the paper / section / equation reference was not.
- URLs to external resources that 404 or redirect.

## What counts as "a reference"

Anything a reader could plausibly follow: a backticked identifier, a path, a citation, a link, a line number in a comment. If you're not sure whether it still matches reality, `rg` / `git log` / open the target file — don't leave a guess in the docs.
