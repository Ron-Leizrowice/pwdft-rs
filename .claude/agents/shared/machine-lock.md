# Machine lock

Multiple agents share this machine. The lock serializes CPU-bound work so `cargo bench` wall-time numbers stay clean. It is a benchmark-contamination guard, not a test-isolation guard — two heavy jobs running concurrently are usually correct; they just invalidate each other's timings.

## Acquire before

- Every `cargo` subcommand that compiles or runs code: `build`, `test`, `bench`, `clippy`, `run`, `doc`.
- Every Quantum ESPRESSO invocation: `pw.x`, `ph.x`, `pp.x`, `bands.x`, `projwfc.x`, `q2r.x`, `matdyn.x`, `dos.x`.
- Any long-running CPU-bound helper (e.g. `pwdft-validate` sub-commands that spin up BLAS).

File reads, edits, git operations, and proposal writing do **not** need the lock.

## Usage

Prefer the `run` one-liner — it acquires, runs, and releases atomically:

```bash
.claude/bin/machine-lock run "<your-role>" "<what you're doing>" -- <command>
```

For multi-command sessions, acquire once at session start and release when done:

```bash
.claude/bin/machine-lock acquire --wait --timeout=900 "<role>" "<desc>"
trap '.claude/bin/machine-lock release' EXIT
# ... work ...
```

Check status without acquiring: `.claude/bin/machine-lock status`.

## Rules

- Acquire from your worktree; the lock is worktree-scoped (MLFX). Running cargo from a different worktree while the lock is held is denied by the Bash PreToolUse hook.
- Never force-remove another agent's lock. If `status` shows a stale lock (dead PID), the next `acquire` will reclaim it automatically.
- Release promptly. Don't hold the lock while reading code or drafting proposals.
