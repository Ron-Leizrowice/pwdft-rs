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

## Cargo and QE: synchronous vs. background

The right choice depends on expected wall-time.

**Short commands (< 5 min expected): synchronous Bash.** `run_in_background: false` (the default). The harness returns stdout when cargo finishes. Includes: Tier-1 `/test`, single-file `/lint`, `/pr-draft` (`cargo check`), warm-cache clippy.

**Long commands (> 5 min expected): background Bash + do other work.** `run_in_background: true`. The harness emits a stream-idle watchdog if your conversation shows no assistant activity for 600 seconds, and a synchronous cargo that compiles silently for 10+ minutes will kill your session. The 2026-04-21 ESPL and MLRW recovery agents both died this way — cold-cache Tier-2 / pytest compilations ran ~10–15 min before emitting any output.

Includes: `/test --tier2` or `--all` on a cold cache, `/bench`, `/profile`, `/quality-gate` on a cold cache, full `cargo doc`, pytest suites on first run after `uv sync`, QE SCF for heavy-atom cells.

**Rules for background jobs:**

- Do *real* parallel work between launch and completion: draft the PR body, write the session logbook, re-read the proposal's Acceptance section, review the diff. This keeps your conversation stream alive and the wall-time isn't wasted.
- Don't use tight `sleep 10 && tail` / `until ! ps -p X; do sleep 5; done` / `while ! grep -qE` polling loops. Each iteration has hook + shell-startup overhead, and `sleep` burns real wall-time. The 2026-04-21 ROTI agent lost ~5 min of its session to this pattern.
- If you need to check completion, do it after a sizable interval (60–120 s) and only when you've run out of other work. The harness will also emit a completion notification when the bg bash exits; prefer waiting for that.
- A well-structured background workflow: launch cargo bg → write PR body → write logbook → draft FLUP section → check for completion notification or (if nothing came) one `TaskOutput`-style read.

**Warm the cache before long runs.** If you've just `git reset --hard` + `git rebase`, the worktree's `target/` is invalidated and the next cargo call pays the full compile. Run one Tier-1 test first (synchronous, ~30 s) to warm the cache, then run Tier-2 or bench commands after — they'll start emitting test output in seconds instead of silently compiling for 10+ minutes.
