# Session-end logbook

Before ending your session, append a dated entry to `.claude/logbooks/<role>.md`.

Logbooks are handoff documents, not diaries. The next session in this role should be able to read it in 30 seconds and know:

- What was done or decided
- What's blocked or unfinished
- Key numbers — metrics, measurements, discrepancies (not prose)
- Tangential ideas worth capturing (one line each)

If an entry exceeds ~30 lines, you're writing too much. Compress or move to a proposal.

## Worktree sessions can't write to the main logbook

Sub-agents in a worktree cannot Write to `.claude/logbooks/<role>.md` — the hook blocks cross-worktree writes. Options:

- **Preferred:** paste the handoff text into your PR body or your final return message. The EM appends to the logbook when merging.
- **Acceptable:** write the note to `/tmp/<role>-handoff-YYYYMMDD.md` and reference its path in your return message.
- **Never:** try to `cp` or `mv` the note into the main checkout. The hook and your lack of write access both stop this.

Sessions running directly in the main checkout (no worktree) can append to logbooks in place.
