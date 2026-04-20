---
id: MLRW
status: active
priority: medium
complexity: medium
risk: medium
depends_on: []
blocks: [MLDX]
---

# MLRW: Reader-writer machine lock — Python library in `pwdft_validation.lock`; exclusive (bench/profile/QE) vs shared (test/clippy/build/doc)

## Problem

The machine lock today is a single mutex implemented as a bash script
(`.claude/bin/machine-lock`) with atomic `mkdir` as the primitive: one
CPU-bound job at a time, full stop. That was the right first-cut because
benchmark contamination is the worst-case failure mode — a `cargo bench`
run sharing CPU with `cargo test` produces garbage wall-time numbers,
which is the thing we're most trying to prevent (see WFRX #99
post-mortem referenced in `check-cargo-lock.sh`'s MLFX fix comment).

But the policy is stricter than the invariant it protects. Two
`cargo clippy` runs from different worktrees don't invalidate each
other. Neither do two `cargo test` runs, nor a test + a clippy. The
only combinations that actually require exclusion are:

- **anything + benchmark** (benchmark wall-time is the protected quantity)
- **anything + samply profile** (profile traces conflate with background noise)
- **anything + QE** (QE is typically a timing-sensitive reference run)

With 2–4 concurrent recovery agents (as in the 2026-04-21 ESPL + ROTI
session), a single `/test --tier2` blocks every other agent's `/lint`
or `/test` for the full Tier-2 wall-time. That's a real throughput tax
for runs that don't care about each other's CPU share.

**Secondary problem.** The current shell implementation is hard to test.
`.claude/bin/tests/machine-lock.test.sh` exercises the happy path but
cannot cleanly simulate concurrent acquires, stale-PID races, or the
priority-promotion machinery a reader-writer lock needs. The
reader-writer algorithm itself is moderately tricky — the kind of code
that wants pytest + threading + a clock monkeypatch, not nested bash
`case` blocks.

**Tertiary problem.** `pwdft_validation/qe/lock.py` (planned in PYQE
Phase A) shells out to the bash `machine-lock` script to coordinate QE
runs. If the primitive lives in Python instead, PYQE gets a direct
import and loses a subprocess hop.

## Proposal

Two changes bundled because they land together:

1. **Move the lock primitive to Python**, living at
   `pwdft/pwdft-validation/pwdft_validation/lock/`.
   `.claude/bin/machine-lock` becomes a thin bash shim that execs the
   Python CLI. The `check-cargo-lock.sh` hook stays pure bash (it fires
   on every Bash tool call; must remain ~50 ms) and reads the JSON
   state file the Python layer maintains.
2. **Add reader-writer classes.** Exclusive (bench/profile/QE) stays
   as-is. Shared (build/test/clippy/doc) becomes unbounded-concurrent
   as long as no exclusive holder exists. Exclusive is priority-
   promoting: once an exclusive acquire posts intent, new shared
   acquires must wait. Prevents starvation under a steady stream of
   cheap shared acquires.

### Truth table

| Caller          | Shared holders present | Exclusive holder present | `exclusive.pending` set | Action     |
|-----------------|:----------------------:|:------------------------:|:-----------------------:|------------|
| **Shared**      | yes                    | no                       | no                      | join       |
| **Shared**      | no                     | no                       | no                      | join       |
| **Shared**      | any                    | yes                      | any                     | wait       |
| **Shared**      | any                    | no                       | yes                     | wait       |
| **Exclusive**   | no                     | no                       | no                      | acquire    |
| **Exclusive**   | yes                    | no                       | no                      | wait+drain |
| **Exclusive**   | any                    | yes                      | any                     | wait       |

### File layout — state

```text
.claude/locks/
├── state.json              single-file state, written atomically via tempfile + rename
└── history.jsonl           append-only audit log (acquire/release/timeout events)
```

`state.json` schema (sketch):

```json
{
  "schema_version": 1,
  "exclusive": {
    "agent": "performance-engineer",
    "desc": "cargo bench scf_iter",
    "worktree": "/Users/.../pwdft-rs/.claude/worktrees/agent-abc",
    "pid": 12345,
    "ts_acquired": 1745320000
  },
  "exclusive_pending": {
    "agent": "core-engineer",
    "worktree": "/Users/.../pwdft-rs/.claude/worktrees/agent-def",
    "pid": 23456,
    "ts_pending": 1745320050
  },
  "shared": [
    {"agent": "core-engineer", "desc": "cargo test", "worktree": "...", "pid": 34567, "ts_acquired": 1745319990},
    {"agent": "code-reviewer", "desc": "cargo clippy", "worktree": "...", "pid": 45678, "ts_acquired": 1745319995}
  ]
}
```

Any of `exclusive`, `exclusive_pending`, `shared` may be absent/empty.
Atomicity: all writes go through `write-tempfile-fsync-rename` so a
concurrent reader never sees a torn file.

### Python package layout

```text
pwdft/pwdft-validation/pwdft_validation/lock/
├── __init__.py          public API: acquire, release, status, run, with-statement helper
├── state.py             LockState dataclass, JSON ↔ dataclass round-trip, atomic write
├── acquire.py           class-aware acquire algorithm (shared + exclusive + pending)
├── pid.py               PID-liveness (kill -0) + staleness (3 h hard cap)
├── paths.py             LOCK_DIR anchored at `git rev-parse --git-common-dir`'s parent
└── cli.py               cyclopts subcommands, wired into existing pwdft-validate app
```

Tests at `pwdft/pwdft-validation/tests/test_lock_*.py` using pytest +
`threading.Barrier` for concurrent-acquire coverage.

### Acquire algorithm — shared

```python
def acquire_shared(agent, desc, worktree, pid, timeout):
    deadline = monotonic() + timeout
    while True:
        with state_lock():             # filesystem-level critical section (mkdir-atomic)
            state = read_state()
            if not state.exclusive and not state.exclusive_pending:
                state.shared.append(SharedHolder(agent, desc, worktree, pid, now()))
                write_state(state)
                return SharedToken(pid)
        if monotonic() > deadline:
            raise AcquireTimeout("shared")
        sleep(POLL_SECONDS)
```

### Acquire algorithm — exclusive

```python
def acquire_exclusive(agent, desc, worktree, pid, timeout, drain_timeout):
    # Post intent under the filesystem lock so shared acquires back off.
    with state_lock():
        state = read_state()
        if state.exclusive_pending and state.exclusive_pending.pid != pid:
            ...  # wait for competing exclusive; see below
        state.exclusive_pending = Pending(agent, worktree, pid, now())
        write_state(state)

    # Drain in-flight shared holders with a bounded wait.
    drain_deadline = monotonic() + drain_timeout
    while True:
        with state_lock():
            state = read_state()
            prune_dead(state)              # kill -0 check on shared holders
            if not state.exclusive and not state.shared:
                state.exclusive = ExclusiveHolder(agent, desc, worktree, pid, now())
                state.exclusive_pending = None
                write_state(state)
                return ExclusiveToken(pid)
        if monotonic() > drain_deadline:
            with state_lock():
                state = read_state()
                if state.exclusive_pending and state.exclusive_pending.pid == pid:
                    state.exclusive_pending = None
                    write_state(state)
            raise DrainTimeout()
        sleep(POLL_SECONDS)
```

`state_lock()` is a filesystem-level mutex serializing the read-modify-
write cycle: `mkdir .claude/locks/state.lock.d` on entry, `rmdir` on
exit, with PID-based staleness for the lock itself (identical to the
current MLFX primitive — unchanged hardening).

### Defaults

| Constant                | Value      | Rationale                                          |
|-------------------------|-----------:|----------------------------------------------------|
| `POLL_SECONDS`          | 2          | Unchanged from current shell impl.                 |
| `DEFAULT_TIMEOUT`       | 600        | Unchanged; overall acquire wait.                   |
| `DRAIN_TIMEOUT`         | 120 (new)  | How long an exclusive waits for shared drain.      |
| `HARD_TIMEOUT_SECONDS`  | 10800      | Unchanged; stale-PID safety net.                   |

### CLI surface

The bash shim `.claude/bin/machine-lock` becomes:

```bash
#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
exec uv run --project "$REPO_ROOT/pwdft/pwdft-validation" \
    pwdft-validate lock "$@"
```

`pwdft-validate lock` subcommands (default class stays **exclusive** so
pre-migration invocations remain safe):

```text
pwdft-validate lock run [--shared|--exclusive] <role> <desc> -- <cmd...>
pwdft-validate lock acquire [--shared|--exclusive] [--wait] [--timeout=N] <role> <desc>
pwdft-validate lock release                                    # class inferred from PID ownership
pwdft-validate lock status [--json]                            # dumps state.json; --json for machine-readable
```

PYQE Phase A.lock becomes a two-line context manager importing
`pwdft_validation.lock.acquire` directly — no subprocess hop.

### Skill updates

| Skill            | Mode       | Rationale                                        |
|------------------|------------|--------------------------------------------------|
| `/bench`         | exclusive  | Benchmark measurements are protected.            |
| `/profile`       | exclusive  | samply trace fidelity.                           |
| `/qe-runner`     | exclusive  | QE reference runs.                               |
| `/test`          | shared     | Tier-1/2 correctness; timing not tracked.        |
| `/lint`          | shared     | clippy.                                          |
| `/quality-gate`  | shared     | Aggregates clippy + tests.                       |
| `/pr-submit`     | shared     | Wraps `/quality-gate`.                           |
| `/pr-draft`      | shared     | Runs `cargo check` only.                         |
| `/cargo`         | exclusive (default) + `--shared` flag | Safe-default conservative.  |

### Hook changes (`check-cargo-lock.sh`)

Stays bash. Three changes, all small:

1. **Read the new JSON state file** via a `python3 -c "import json; ..."`
   one-liner (same pattern the hook already uses for parsing
   `tool_input`). Inspect `state.exclusive` and `state.shared[]`.
2. **Rule reorder**:
   - If `state.exclusive` exists and caller's cwd ≠ `exclusive.worktree` → deny.
   - Else if caller's cwd matches any `shared[].worktree` → allow (caller already holds a shared slot).
   - Else → deny (caller must acquire, shared or exclusive).
3. **Extend matched command set.** Current regex misses `cargo doc` and
   every QE binary (`pw.x`, `ph.x`, `pp.x`, `bands.x`, `projwfc.x`,
   `q2r.x`, `matdyn.x`, `dos.x`) — both listed in `machine-lock.md`'s
   "acquire before" section but un-enforced in the hook. Pre-existing
   drift, closed here to avoid a separate tiny PR.

### Bootstrap / failure modes

- **uv environment broken.** The bash shim's `uv run` fails → cargo
  command fails with a clear error. Mitigation: add a `uv sync` check
  to the repo's setup script; document in CLAUDE.md that `uv sync` is a
  one-time prerequisite after a fresh clone or a dep bump.
- **Python raises mid-acquire.** Releases any state.lock.d held, leaves
  no partial `state.json` (tempfile+rename). Next acquire recovers.
- **Process crash with the exclusive lock held.** PID-liveness check on
  the next acquire reclaims it (same as today).
- **CI environment.** `uv sync` already runs in CI; `uv run pwdft-validate
  lock run ...` works unchanged.

## Risk

**Medium.** Concurrency primitives are easy to get subtly wrong. The
test plan is structured to catch the usual failure modes.

Specific cases in the pytest suite:

- Two parallel `acquire_shared` from different worktrees → both
  succeed, both land in `shared[]`.
- `shared` held, `acquire_exclusive` arrives → posts
  `exclusive_pending`, blocks; concurrent new `acquire_shared`
  attempts block on the pending flag. Shared releases → exclusive
  drains + acquires; pending flag cleared.
- `exclusive` held, new `acquire_shared` → blocks until release.
- Shared holder process dies → next exclusive's drain step prunes via
  `kill -0`.
- Two exclusives race to post `exclusive_pending` → one wins via the
  filesystem `state_lock()`; the loser sees the won flag and waits.
- Exclusive's drain times out → releases `exclusive_pending`, raises;
  shared acquires proceed.
- Legacy state-file format detected → treated as stale, cleared with
  a warning (migration from MLFX shell-era).

Shell-level integration tests retained as
`pwdft/pwdft-validation/tests/integration/test_lock_cli.sh` — exercise
the bash shim + hook interaction end-to-end.

## Non-goals

- No daemon. See MLDX for the FastAPI + dispatch follow-up.
- No numeric concurrency cap on shared. If oversubscription becomes a
  real problem, file follow-up with a measured `SHARED_MAX`.
- No core-affinity work.
- Not rewriting `check-cargo-lock.sh` in Python — hook latency budget
  (~50 ms) rules out uv-startup per Bash tool call.

## Implementation plan

One PR, organized by file:

1. `pwdft/pwdft-validation/pwdft_validation/lock/` — new sub-package
   (`__init__.py`, `state.py`, `acquire.py`, `pid.py`, `paths.py`,
   `cli.py`).
2. `pwdft/pwdft-validation/pwdft_validation/cli.py` — register
   `lock` subcommand group on the cyclopts app.
3. `pwdft/pwdft-validation/tests/test_lock_*.py` — pytest concurrency
   suite.
4. `pwdft/pwdft-validation/tests/integration/test_lock_cli.sh` —
   shell integration test.
5. `.claude/bin/machine-lock` — replace 260-line bash script with
   ~8-line uv-run shim.
6. `.claude/bin/check-cargo-lock.sh` — read `state.json`; extend
   matched command regex to include `doc` + QE binaries; dual-class
   ownership check.
7. `.claude/bin/tests/machine-lock.test.sh` — keep as the end-to-end
   shell smoke test; update expectations for the new state-file
   format.
8. `.claude/agents/shared/machine-lock.md` — document the
   shared/exclusive split; update usage examples.
9. `.claude/skills/*/SKILL.md` — per the table above, thread
   `--shared` or `--exclusive` into each skill's
   `machine-lock run` invocation.

Migration is a hard break per the no-backcompat policy: the state-file
format changes from per-field directory to single JSON file. Any
existing `machine.lock.d/` on disk is detected by the new code and
cleared with a "migrated from MLFX shell-era lock" log line.

## Acceptance

- `pwdft_validation.lock` sub-package exists and is imported by
  `pwdft_validation.qe.lock` (forward PYQE hookup; PYQE may land
  after MLRW and consume the new API).
- Two `machine-lock run --shared` invocations from different
  worktrees run concurrently on warm cache; wall-time ≈ max, not sum.
- `machine-lock run --exclusive` with a shared holder present waits
  (not deny) until the shared holder releases.
- `machine-lock run --shared` while `exclusive_pending` is set waits,
  not joins.
- All existing bench/profile/QE paths remain exclusive (regression
  guard: `/bench` still excludes `/test`).
- `cargo doc` and QE binaries are hook-enforced (new).
- `pytest pwdft/pwdft-validation/tests/test_lock_*.py -q` green.
- `.claude/bin/tests/machine-lock.test.sh` green.
- CLAUDE.md § "Machine coordination" updated to reference the new CLI
  surface.

## Measurement

Before merging, record wall-time of two scenarios on warm cache:

1. **Serial baseline** (today's semantics): two agents sequentially
   run `/lint` + `/test` from separate worktrees. T₁ = sum.
2. **Parallel-shared** (proposal): same two agents, run concurrently.
   T₂ ≤ max + ε.

Expect T₂ ≈ max(T_lint, T_test) + O(5 s) sync noise on M3 Max,
versus T₁ = T_lint + T_test. On current Tier-1 numbers (~12 s warm
test + ~20 s warm clippy) that's a ~50% wall saving for 2-way
concurrency, and the speedup grows with N agents.

## Previous attempt (2026-04-21, stalled mid-pytest)

A first core-engineer agent (ID `a639cb1b4a3237eba`) made partial
progress before stalling on the stream-idle watchdog during pytest
compilation. Work preserved as uncommitted files in the worktree
**`.claude/worktrees/agent-a639cb1b/`**; the branch
`MLRW/reader-writer-machine-lock` on origin/local is a no-op branch
pointing at origin/main. The resume agent should **copy the 4
Python files verbatim** from that worktree as their starting point —
they are well-architected and only need ty-diagnostics cleanup (19
reported, mostly unused imports + `from __future__` placement), not
redesign.

### Files preserved (copy from `agent-a639cb1b/pwdft/pwdft-validation/pwdft_validation/lock/`)

| File           | Lines | Content                                                                                              |
|----------------|------:|------------------------------------------------------------------------------------------------------|
| `__init__.py`  | 44    | Public API re-exports (`AcquireResult`, `Holder`, `LockClass`, `LockState`, `acquire`, `release`, `status`). |
| `paths.py`     | 85    | `LockPaths` dataclass; anchors lock dir at main-repo (parent of `git rev-parse --git-common-dir`); `PWDFT_LOCK_DIR` env override for test redirection to scratch. |
| `pid.py`       | 60    | `pid_alive(pid)` via `os.kill(pid, 0)` (EPERM → alive since pid exists but isn't ours); `is_stale(ts, pid)` = dead PID OR `age ≥ HARD_TIMEOUT_SECONDS` (3h); `self_pid()` = `ML_OWNER_PID` env OR `os.getppid()`. |
| `state.py`     | 251   | `LockState`/`Holder` dataclasses; `STATE_VERSION = 2` single-JSON blob; `_atomic_write` (tempfile + fsync + `os.replace`); `state_lock` ctx mgr (mkdir-atomic, 50 ms poll, 30 s timeout); `clear_legacy` for MLFX migration; `reap_stale(state) → (new_state, messages)` as pure function. |

### Design decisions already made (preserve these)

- **Single JSON blob at `state.json`**, versioned `STATE_VERSION = 2`
  (bumped from MLFX's per-field-directory layout). Versioning in the
  file lets future migrations detect format.
- **`state.lock.d/` mkdir-atomic** is the read-modify-write
  serializer. Same primitive MLFX used successfully. No `fcntl.flock`.
- **Defaults:** `_STATE_LOCK_POLL_S = 0.05` (50 ms), `_STATE_LOCK_TIMEOUT_S = 30`. `HARD_TIMEOUT_SECONDS = 10800` (3 h) — matches MLFX.
- **PID-liveness treats EPERM as "alive"** — pid exists, not ours. Safety-first: never clear a lock owned by another user.
- **`ML_OWNER_PID` env var** for the owning-PID override (shell callers pass `$PPID`). `PWDFT_LOCK_DIR` env var for test redirection to scratch. Both are test ergonomics surface the agent picked, keep them.
- **`reap_stale` as a pure function** returning `(new_state, messages)` — clean test surface, caller decides whether to log the messages.
- **`clear_legacy(paths) → bool`** detects pre-MLRW `machine.lock.d/` or flat-file lock and removes them on first run; returns whether anything was cleared. Callers emit a warning ("migrated from MLFX shell-era lock") when `True`.

### What's NOT done yet

- `acquire.py` — the reader-writer state machine. Follow the
  pseudocode in the proposal's § Proposal § Acquire algorithm
  exactly. Call into `state.state_lock`, `state.reap_stale`,
  `pid.self_pid` for the building blocks.
- `cli.py` — cyclopts subcommands `run`, `acquire`, `release`,
  `status`. Wire into `pwdft_validation.cli` as a sub-app named
  `lock`.
- `pytest pwdft/pwdft-validation/tests/test_lock_*.py` — the truth
  table cases + concurrency tests (use `threading.Barrier`).
- `.claude/bin/machine-lock` rewrite to `uv run` shim.
- `.claude/bin/check-cargo-lock.sh` rewrite to read `state.json` +
  extended command match (`cargo doc` + QE binaries).
- `.claude/agents/shared/machine-lock.md` — document the split.
- 8× skill SKILL.md updates to thread `--shared` / `--exclusive`.
- Session logbook at `.claude/logbooks/core-engineer/...`.
- Fix the 19 `ty` diagnostics in the preserved files (mostly unused
  imports / `from __future__` placement; see `uv run ty check
  pwdft/pwdft-validation/pwdft_validation/lock/` for the list).

### Why the agent stalled

The 10-minute stream-idle watchdog fires when the agent's conversation
shows no activity for 600 s. The first MLRW agent ran a synchronous
pytest whose first-run compilation went silent past that window.

**Resume-agent directive:** when running pytest or cargo on a cold
cache, launch with `run_in_background: true` and use the time to
write the next file (acquire.py, cli.py, etc.) so the stream stays
alive. See `.claude/agents/shared/machine-lock.md` § "Cargo and QE:
synchronous vs. background" for the full protocol (updated 2026-04-21
after the MLRW + ESPL stalls).

## Flagged for follow-up

- `cargo doc` + QE hook coverage drift (items under Hook changes item
  3). Bundled here; if a cleaner split is preferred, file as its own
  tiny PR and land first.
- PYQE Phase A.lock rework: once MLRW lands, amend PYQE to consume
  `pwdft_validation.lock` directly instead of shelling out to
  `machine-lock`. Pure subtraction — drops a subprocess hop.
- MLDX (daemon + dispatch) builds on this library. See that proposal
  for the phased roadmap.
