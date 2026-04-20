---
id: MLDX
status: active
priority: low
complexity: large
risk: medium-high
depends_on: [MLRW]
blocks: []
---

# MLDX: FastAPI daemon for machine-lock + heavy-command dispatch

## Problem

MLRW moves the lock primitive into `pwdft_validation.lock` as a Python
library. Every acquire still costs a `uv run` spinup — ~500 ms–1 s per
invocation. For cargo commands that take ≥ 12 s that's acceptable noise;
for workflows that acquire-release many times (CI-style per-file clippy
sweeps, or a future job-dispatch model) it's a meaningful tax. Poll-
based waits are also coarse: 2 s polling means a shared holder releasing
imposes up to 2 s latency on a waiting exclusive.

A persistent local daemon (FastAPI over a unix socket) erases both
costs, and opens a set of capabilities the filesystem-primitive can't
offer cleanly:

- Event-driven waits (SSE / long-poll) instead of poll loops.
- Structured in-memory contention metrics and acquire history.
- Job dispatch — agents POST `cargo bench …` and the daemon runs it in
  a known environment, streams output, persists wall-time history.
  Benchmark reproducibility jumps because every bench runs under the
  same process parent, same env vars, same working directory
  discipline.

This proposal is explicitly a follow-up to MLRW. It is NOT a
prerequisite for anything on the critical path today, and should only
proceed once MLRW's library API has been stable for a while.

## Proposal

Two phases, each a separate PR.

### Phase 1 — daemon layer (`MLDX-1`)

Wrap `pwdft_validation.lock` in a FastAPI server listening on a unix
socket. CLI / library clients prefer the daemon when it's running;
fall back to direct library use when it isn't. Hook stays filesystem-
only (the daemon maintains the same `state.json` the hook reads).

```text
┌─────────────────────────────────────────────────────────────┐
│ MLDX-1:  FastAPI server (uvicorn, unix socket)              │
│          POST /locks/acquire   {class, role, desc, worktree}│
│          DELETE /locks/{token}                              │
│          POST /locks/{token}/heartbeat                      │
│          GET /status           (JSON)                       │
│          GET /events           (SSE: acquires + releases)   │
├─────────────────────────────────────────────────────────────┤
│ MLRW:    pwdft_validation.lock library                      │
│          acquire_shared / acquire_exclusive / release       │
│          state.json maintained as a shadow for the hook     │
└─────────────────────────────────────────────────────────────┘
```

#### Client behavior

`pwdft_validation.lock.acquire(...)` gains a daemon-first path:

```python
def acquire(class_, agent, desc, ...):
    sock = LOCKD_SOCKET
    if sock.exists() and _ping_daemon(sock):
        return _acquire_via_daemon(sock, class_, agent, desc, ...)
    return _acquire_via_library(class_, agent, desc, ...)
```

Both code paths converge on the same `state.json` and the same MLRW
semantics — just different transport and different wait primitives
(event-driven vs. poll loop).

#### Daemon lifecycle

`launchd` user agent on macOS (+ `systemd --user` on Linux, for
symmetry; production is macOS-only today but the Linux runner in CI
benefits too). Install script under
`pwdft/pwdft-validation/scripts/install-lockd.sh`:

- Writes `~/Library/LaunchAgents/com.pwdft.lockd.plist` with
  `KeepAlive=true`, `RunAtLoad=true`, `WorkingDirectory=<repo root>`,
  `StandardOutPath=<repo>/.claude/locks/lockd.log`.
- `launchctl load` on first run.
- Idempotent re-install on repo update.

Fallback when the daemon isn't installed: CLI logs one-time "running
without lockd — consider `pwdft-validate lock install-daemon` for
event-driven waits" and uses the library path. No functional
difference to the caller.

#### Heartbeats

Acquire returns a `token` + `heartbeat_interval_s`. Client sends
`POST /locks/{token}/heartbeat` every `interval/3` seconds. Daemon
reclaims a lock with no heartbeat for `3 × interval` — 30 s reclaim
window by default. Complements the current `kill -0` PID check and
catches the case where the PID is reused by an unrelated process.

The Python context-manager helper (`with lock(...)`) handles
heartbeats on a background thread automatically; direct CLI users
don't need to manage them.

#### Crash recovery

Daemon writes every state transition to an append-only journal
(`.claude/locks/lockd.journal`) before responding to the client. On
startup:

1. Replay journal to reconstruct in-memory state.
2. Verify every token against PID liveness.
3. Rewrite `state.json` to match.
4. Truncate journal.

Children of dispatched jobs (see Phase 2) survive daemon death; the
journal records their PIDs, and on restart the daemon reattaches via
`/proc`-style enumeration (on macOS: `ps -p <pid>`; on Linux:
`/proc/<pid>`).

#### Hook behavior (unchanged from MLRW)

`check-cargo-lock.sh` reads `state.json` as in MLRW. No HTTP calls
from the hook — we explicitly do not make the hot path depend on
daemon availability.

#### CI

Ephemeral runners skip the daemon. Set `PWDFT_LOCK_DIRECT=1` in the
workflow env; the client-side check honors it and goes straight to
library mode. No daemon install overhead on CI.

#### Observability

- `pwdft-validate lock status --json` proxies to `GET /status` when
  the daemon is running; shows live in-memory data not just shadow
  state.
- `pwdft-validate lock events` tails `GET /events` (SSE) — useful
  during multi-agent sessions to see who's contending on what.
- Acquire history persisted to `.claude/locks/history.jsonl` — already
  planned in MLRW — gains `waited_ms` and `drained_holders` fields
  once the daemon is dispatching.

### Phase 2 — dispatch layer (`MLDX-2`)

Once Phase 1 is stable, extend the daemon with a job-runner API.
Agents submit heavy commands; the daemon executes them in the lock's
protection and streams output back. Benchmark reproducibility and
historical wall-time tracking are the primary wins.

#### API

```text
POST /jobs
  body: {class: exclusive|shared, cmd: [...], cwd: "...", env: {...}, timeout_s: 1800}
  → 201 Created  {job_id, token}

GET /jobs/{id}                   # status
GET /jobs/{id}/output            # SSE stream of stdout + stderr
DELETE /jobs/{id}                # cancel (SIGTERM → SIGKILL after 10 s)

GET /jobs                        # history + current
```

#### Executor whitelist

`cmd[0]` must be in a compile-time whitelist:
`cargo`, `mpirun`, `pw.x`, `ph.x`, `pp.x`, `bands.x`, `projwfc.x`,
`q2r.x`, `matdyn.x`, `dos.x`, `samply`. Anything else → `403 Forbidden`.

Further argument constraints:

- `cargo` — first positional must be one of
  `test|build|bench|clippy|run|doc|check`.
- `samply` — must be `samply record …`.
- QE binaries — any args allowed.

Rejects attempts to smuggle shell metacharacters; `cmd` is a list,
never a string.

#### Bench reproducibility

Dispatched jobs run with:

- `PWD` set to the submitting worktree root.
- `CARGO_TARGET_DIR` preserved.
- `RUSTFLAGS` cleared unless explicitly in the submission's env.
- `LANG=C`, no `LC_*` overrides.
- Process group created via `setsid()` so `DELETE` can SIGTERM the
  whole tree reliably.

Wall-time measured from `fork()` to `waitpid()` in the daemon, not
from `cargo`'s self-report. Persisted to `history.jsonl` keyed by
`(cmd_signature, git_sha, host)`. A future `pwdft-validate bench
trend …` CLI reads this for regression detection.

#### Agent UX

New skill `/run-bench <bench-name>`:

```bash
# Under the hood:
curl --unix-socket .claude/locks/lockd.sock \
    http://lockd/jobs \
    -X POST -H 'Content-Type: application/json' \
    -d @- << EOF
{"class":"exclusive","cmd":["cargo","bench","-p","pwdft-core","--","$BENCH"]}
EOF
```

Returns a job URL; agent follows the SSE stream for output. On
completion, a summary is emitted with wall-time + bench-local results.

The agent never holds a lock directly — it asks the daemon to run
something on its behalf, and the daemon coordinates exclusion and
measurement.

## Risk

**Medium-high.** Three distinct failure classes:

1. **Daemon reliability.** Crash recovery is only as good as the
   journal replay. Test plan: fault-injection harness that SIGKILLs
   the daemon mid-acquire, mid-release, and mid-dispatch; asserts
   state.json + in-memory state converge correctly on restart.
2. **Hook divergence from reality.** If the daemon forgets to fsync
   `state.json` before responding, a caller can believe they hold the
   lock while the hook sees the old state. Mitigation: hook write-
   before-respond contract; test via deliberate race.
3. **Dispatch scope creep.** The executor whitelist is the
   blast-radius bound. Adding entries requires a PR touching the
   hard-coded tuple. Keep it explicit.

The fallback path (library-direct when daemon is down) is the primary
safety net. If MLDX-1 turns out to be flaky, disabling the daemon
restores MLRW behavior with no data loss.

## Non-goals

- Not a multi-machine coordination layer. One daemon per host.
- Not authenticated. Unix-socket filesystem permissions gate access;
  anyone who can read the repo can use the daemon.
- Not a replacement for cron / CI scheduling. Heavy jobs dispatched
  here are interactive — submitted by an agent during a session,
  reaped when the session ends.
- Not a general-purpose task queue. Executor whitelist is explicit.

## Implementation plan

### Phase 1 — MLDX-1 (one PR)

1. `pwdft/pwdft-validation/pwdft_validation/lockd/` — FastAPI app
   (`app.py`, `routes.py`, `journal.py`, `heartbeat.py`, `shadow.py`).
2. `pwdft/pwdft-validation/pwdft_validation/lock/client.py` —
   unix-socket client that prefers the daemon and falls back to
   `pwdft_validation.lock.acquire` (library path, unchanged from
   MLRW).
3. `pwdft/pwdft-validation/pwdft_validation/lock/cli.py` — add
   `install-daemon` / `uninstall-daemon` / `events` subcommands.
4. `pwdft/pwdft-validation/scripts/install-lockd.sh` +
   `com.pwdft.lockd.plist` template.
5. `pwdft/pwdft-validation/tests/test_lockd_*.py` — FastAPI testclient
   - fault-injection tests.
6. `.claude/agents/shared/machine-lock.md` — document the daemon,
   when it matters, when to ignore it.
7. CLAUDE.md — one-paragraph note under § Machine coordination.

### Phase 2 — MLDX-2 (one PR, gated on MLDX-1 stability)

1. `pwdft/pwdft-validation/pwdft_validation/lockd/jobs.py` — executor
   whitelist, job registry, SSE streaming.
2. `pwdft/pwdft-validation/pwdft_validation/lockd/history.py` —
   wall-time + signature persistence.
3. `.claude/skills/run-bench/SKILL.md` — new skill.
4. Tests: dispatch fault injection, whitelist bypass attempts,
   concurrent-job coordination with MLRW lock classes.

## Acceptance — MLDX-1

- Daemon runs under `launchd`; survives reboots; restarts cleanly
  after SIGKILL with no state loss.
- Acquire latency with daemon ≤ 10 ms p50 (vs. ~500 ms–1 s library
  uv-run). Measured via a microbench in
  `pwdft/pwdft-validation/benches/lock_acquire.py`.
- Shared→exclusive handoff latency ≤ 50 ms p50 (vs. up to
  `POLL_SECONDS = 2` in MLRW).
- Hook continues to work when daemon is stopped (filesystem-only
  path); no regressions in `.claude/bin/tests/machine-lock.test.sh`.
- Heartbeat reclaim: a client that stops heartbeating loses its lock
  within 30 s.
- Crash: SIGKILL the daemon, restart, state reconstructs identically.

## Acceptance — MLDX-2

- `POST /jobs` with a non-whitelisted executor returns 403.
- Dispatched `cargo bench` completes; output captured; wall-time
  persisted to `history.jsonl`.
- `DELETE /jobs/{id}` terminates the entire process tree.
- Two concurrent `class: shared` jobs run in parallel; two `class:
  exclusive` jobs serialize.
- `pwdft-validate bench trend <bench>` reads history.jsonl and shows
  the last N runs with their git shas.

## Measurement

Phase 1 microbenchmark:

```text
Scenario                          MLRW (no daemon)   MLDX-1 (daemon)
Single acquire+release            ~1000 ms           ~5 ms
10 back-to-back acquires          ~10 s              ~50 ms
Shared-held, exclusive waits      up to 2 s          ~10 ms
```

Phase 2 reproducibility check: run `scf_iter` bench 10 times each via
`/bench` (MLRW) and `/run-bench` (MLDX-2); compare coefficient of
variation. Expect MLDX-2 CoV ≤ 0.5× MLRW CoV because of the
process-environment control.

## Flagged for follow-up

- If the launchd install proves flaky or user-hostile, consider an
  auto-spawn-on-first-use model with self-exit-after-idle. Document
  experience from MLDX-1's first two weeks before deciding.
- Authenticated multi-agent support (tokens, scoped permissions) is
  explicit non-scope; re-open only if we ever run the daemon outside
  a single-user dev machine.
- `pwdft-validate bench trend` UX is out of scope for MLDX-2 and
  should land as its own small proposal when the history JSONL has
  enough data to justify the CLI.
