---
id: PROF
status: active
priority: medium
complexity: trivial
risk: low
depends_on: []
blocks: []
---

# PROF: Adopt `samply` as the canonical profiler; codify the observability vs. profiling split

## Problem

The project's observability + profiling tooling is implicit. The
performance engineer's logbook has ad-hoc references to `cargo bench`
and one-off `Instant::now()` instrumentation; CLAUDE.md mentions
`cargo bench` for benchmarks but never names a profiler. When a
performance question comes up ("where does Si SCF spend its 75 ms per
iter?"), each agent picks a different tool — `cargo flamegraph`,
hand-rolled timers, manually editing the SCF loop with `eprintln!`s
and timestamps.

Two consequences:

1. **Inconsistent profiling answers.** Different tools, different
   sampling rates, different presentation. Hard to compare
   measurements across agents or across weeks.
2. **Easy to reach for the wrong abstraction.** A `tracing`-based stack
   would offer "span-based timing without surgery", which sounds
   appealing until you realize a sampling profiler does the same thing
   with zero code overhead and sees inside `faer` / `ndrustfft` /
   Accelerate-BLAS where `#[instrument]` annotations cannot reach. We
   considered and dropped the `tracing` ecosystem — see § "When to
   reconsider tracing" below.

We need to (a) name the canonical profiler and (b) document when to
use which tool — profiler vs. benchmark vs. log — so future
performance work converges on a single answer.

## Research

### Why samply (not `cargo flamegraph`, not `tracing-flame`, not Instruments.app)

| Tool                 | Platform | Code overhead | Sees libraries | Notes |
|----------------------|----------|---------------|-----------------|-------|
| **samply**           | macOS, Linux | none      | yes              | Sampling profiler, native to macOS. `cargo install samply`. Output is a self-contained HTML profile (Firefox Profiler UI) — shareable, viewable in any browser, no setup on the receiving end. |
| `cargo flamegraph`   | macOS (with sudo / `dtrace`), Linux (`perf`) | none | yes | Works but rougher on macOS — needs sudo for `dtrace`, output is SVG only. |
| `tracing-flame`      | any       | ~30–50 ns per `#[instrument]`, even when unsubscribed | no — only annotated functions | Biases hot-path measurements. Cannot see into `faer`, BLAS, FFT internals. |
| Instruments.app      | macOS    | none          | yes              | Excellent UI, but Apple-only; requires Xcode; opens its own GUI; profiles are not portable. Good for deep dives, bad for "drop into the agent's logbook." |
| `perf`               | Linux    | none          | yes              | Not available on darwin. Skip. |

Samply wins on darwin because: (a) zero overhead, (b) sees library
internals, (c) profile-as-HTML (the perf engineer can paste a profile
URL into a logbook entry; the EM can open it in any browser without
installing tools), (d) one-line install via `cargo install samply`.

The HTML profile is the killer feature for our agent workflow — agent
output needs to be inspectable by the human reviewer, and a Firefox
Profiler URL is more reviewable than a 4 MB SVG.

### Cross-checking against the "stay on `log`" decision

The observability vs. profiling vs. benchmarking split is:

- **Observability** (what is the SCF doing?): `log` + `env_logger` +
  `indicatif`. Answers questions like "did it converge?", "what was the
  final β?", "which mixer ran?". Always-on, structured for human
  reading, low overhead.
- **Profiling** (where does wall-time go?): `samply`. Answers questions
  like "is V_NL apply or eigensolve dominating now?", "did
  GOPT-PR-B regress XC kernel scheduling?". Run on demand, not in
  production, sees everything.
- **Benchmarks** (did this PR regress function X?): `criterion` via
  `cargo bench`. Answers regression questions with statistical
  confidence. Already in place; no change.

Three tools, three jobs, no overlap. Each is the simplest tool that
serves its job; together they cover everything an alternative
unified-observability stack (e.g., `tracing` + `tracing-flame` +
`tracing-indicatif`) would have tried to do under one roof.

### When to reconsider `tracing`

The `tracing` ecosystem is not on this project's roadmap. Reopen the
question only if one of the following becomes true:

1. An automated regression watcher or validation dashboard needs to
   parse per-iteration SCF events as structured records (not regex over
   `log` lines).
2. A dependency starts emitting `tracing` events we want to filter or
   re-route through our own subscriber. (`faer`, `ndrustfft`, `wgpu`
   are all `log`-only today.)
3. Profiling needs evolve beyond on-demand sampling — e.g., we want
   continuous in-process timing aggregates surfaced live to the user
   rather than `samply record` snapshots.

Until then, MIXL covers the observability content gap, LOGH consolidates
the transport on `log`, and this proposal pins down the profiling story.
Each is a focused tool for one job, which is cheaper than one framework
trying to serve all three.

## Implementation

This is a documentation + tooling-discoverability proposal, not a code
change. Three small landings.

### Step 1 — Add samply install and recipe to CLAUDE.md

In `CLAUDE.md` § Tests & Benchmarks (or a new § Profiling), add:

```markdown
## Profiling

Use `samply` for wall-time profiling of `cargo run` or `cargo bench`
invocations. Install once:

    cargo install samply

Run a profile:

    .claude/bin/machine-lock run "Performance Engineer" "samply Si SCF" -- \
      samply record cargo run --release -- --input examples/si_scf.yaml

Samply opens a browser tab with the Firefox Profiler view. Paste the
local URL (or the exported `.json.gz`) into your logbook entry so
reviewers can replay the profile.

For benchmarks specifically:

    .claude/bin/machine-lock run "Performance Engineer" "samply scf_bench" -- \
      samply record cargo bench --bench scf_benchmarks -- --profile-time 10

Do **not** use `tracing-flame` or hand-rolled `Instant::now()` timers
for new profiling work — samply sees inside `faer` / `ndrustfft` /
BLAS and adds zero overhead. Use criterion (`cargo bench`) for
regression checks of specific functions, samply for end-to-end
wall-time questions, and `log` / `info!` for runtime observability.
```

### Step 2 — Update `.claude/agents/performance-engineer.md`

Add a "Profiling stack" section pointing at samply. Replace any existing
references to "cargo flamegraph" or hand-rolled timers with the samply
recipe. Note the lock-acquire pattern (samply runs are CPU-saturating
and need the machine lock just like `cargo bench`).

### Step 3 — Capture the decision in the Performance Engineer logbook

Append one entry to `.claude/logbooks/performance-engineer.md`
recording the decision (samply adopted 2026-04-19, supersedes ad-hoc
profiling) and a link to this proposal so the choice is discoverable
when the next perf engineer reads back.

### Step 4 (optional) — Smoke-test with one real profile

Run samply on a Si SCF (`examples/si_scf.yaml`) and capture the resulting
profile in the logbook entry as a baseline. Future profiles can compare
against this snapshot to validate the toolchain is working consistently.
This is a "verification by use" step — if the recipe in CLAUDE.md
doesn't actually work, this catches it.

## Verification

Acceptance criteria are documentation-existence checks, not code
behaviour:

- `which samply` returns a path (`cargo install samply` succeeded for
  the agent setting up the box).
- CLAUDE.md § Profiling exists and contains the recipe.
- `.claude/agents/performance-engineer.md` references samply by name.
- One Si SCF profile is captured in
  `.claude/logbooks/performance-engineer.md` as a baseline.

## Out of scope

- Integration with CI / automated regression detection.
  `cargo bench`'s criterion JSON output already supports that and is
  the right place to add it; PROF only documents the manual profiling
  story.
- Alternative profilers (Instruments.app, perf, hyperfine) — they may
  be useful for one-off deep dives but should not be the default.
- Documenting `criterion` itself; it's already documented by use in
  `benches/`.
- Migrating any existing benchmark to a new format. Criterion stays.

## Engineering stack — codified

For future-proofing, here's the full observability/measurement stack
this proposal locks in:

| Use case                                         | Tool                          | Status                |
|--------------------------------------------------|-------------------------------|-----------------------|
| Runtime user-facing status (SCF progress, errors) | `log` + `env_logger`          | In place              |
| Per-iteration progress display                    | `indicatif`                   | In place              |
| Wall-time profiling (where does time go?)         | `samply`                      | **PROF adopts**       |
| Per-function regression detection                 | `criterion` via `cargo bench` | In place              |
| Cross-agent narrative + decision history          | `.claude/logbooks/`           | In place              |
| Numerical accuracy validation                     | QE 7.5 reference calculations | In place              |
| Error reporting                                   | `thiserror` + `Result`        | In place              |

Anything outside this set is a tooling decision worth its own proposal.
