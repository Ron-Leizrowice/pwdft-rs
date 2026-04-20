---
name: bench
description: Run cargo bench under the machine lock. Accepts a bench name and any trailing criterion args. Trigger on "/bench".
user_invocable: true
---

# /bench — Machine-locked benchmark run

Run `cargo bench` with `$ARGUMENTS` appended, under the machine lock.

## Invocation

```bash
.claude/bin/machine-lock run "Performance Engineer" "cargo bench $ARGUMENTS" -- cargo bench $ARGUMENTS
```

Typical uses:

- `/bench --bench scf_benchmarks` — full SCF macro bench.
- `/bench --bench gpu_benchmarks --features gpu` — GPU-feature bench group.
- `/bench --bench scf_benchmarks foo_n725` — filter to a specific case.
- `/bench --bench scf_benchmarks -- --profile-time 10` — run for samply-style per-case wall time.

## Production-scale reminder

Every new bench group must include at least one production-scale configuration (`n_pw ≥ 500`, `n_atoms ≥ 8`, `n_kpoints ≥ 4×4×4`). Micro-bench-only additions get rejected at review — a 2× speedup on a 50 ms Γ-only Si test is a rounding error; research value lives at the minutes-to-hours scale.

## Measurement hygiene

- Baselines and post-change numbers must come from runs in the **same worktree** (same `Cargo.lock`, same `target/` cache). Don't compare against a number reported from a different commit.
- Criterion variance is real — a single outlier doesn't prove a regression. Diff the code first; if the diff is zero-cost, predict "this re-benches flat" before paying the lock time (see Performance Engineer § Analytical-first investigation).

## On failure

Report the failing bench, the run log, and stop.
