---
name: profile
description: Record a samply profile of a command under the machine lock. Trigger on "/profile". Samply is the canonical profiler — cross-platform, unprivileged, zero code overhead, emits a shareable Firefox Profiler HTML.
user_invocable: true
---

# /profile — Samply under the machine lock

Wrap `$ARGUMENTS` in `samply record` under the machine lock. Samply saturates the CPU like `cargo bench`, so the lock is mandatory.

## Prereq

`samply` is installed (`cargo install samply`). If missing, report and stop; do not auto-install.

## Invocation

```bash
.claude/bin/machine-lock run "Performance Engineer" "samply $ARGUMENTS" -- samply record $ARGUMENTS
```

Typical uses:

- `/profile cargo run --release -- --input inputs/si_scf.yaml` — SCF end-to-end profile.
- `/profile cargo bench --bench scf_benchmarks -- --profile-time 10` — bench profile, criterion gated.
- `/profile cargo test --release --test free_electron_bands` — targeted test profile.

## When to reach for something else

- `cargo flamegraph`, `tracing-flame`, hand-rolled `Instant::now()` timers — **don't.** Samply sees inside `faer` / `ndrustfft` / BLAS where those can't, adds zero overhead, and produces a portable artifact.
- Instruments.app — fallback only, when the question is specifically about Metal GPU timelines (Xcode Metal debugger, GPU trace captures).

## Artifact handling

Samply prints the profile URL (or serves it on localhost). Paste the URL (or the exported `.json.gz` path) into your PR body or logbook entry so the reviewer can replay it in Firefox Profiler.

## On failure

Report samply's stderr and stop.
