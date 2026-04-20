---
id: LTOB
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# LTOB: Benchmark `lto = "fat"` vs `"thin"` on the SCF hot path

## §1 Scope

Measure the wall-time delta between `lto = "fat"` and the current
`lto = "thin"` release profile on representative SCF workloads, then
decide whether the build-time cost of `fat` is justified by the
runtime win.

Strictly a measurement-and-decision proposal — no code refactors, no
algorithmic changes. The deliverable is a numbers table plus a
one-line `Cargo.toml` flip (or a decision to keep `thin`).

## §2 Why now

Background: during the workspace restructure (2026-04-20) the root
`Cargo.toml` briefly switched to `lto = "fat"` before being reverted
to `"thin"` pending data. `fat` LTO does whole-program optimization
across the entire crate graph (pwdft-core + faer + gemm + ndrustfft
+ nalgebra + …) in a single unit, so hot paths that cross crate
boundaries — especially `faer::…` calls inlined into SCF eigensolve
and V_NL assembly — can see additional vectorization and dead-code
elimination beyond what `thin` achieves.

Published experience elsewhere is mixed:

- **Typical win:** 2-8% on inner-loop-heavy numeric code.
- **Typical cost:** 2-5× longer release link time; worse in a
  workspace where every dep + member participates in the single-unit
  link.

On our codebase the vendored faer is large and the dep graph deep, so
the link-time cost could be substantial. Against that, production
SCF runs take minutes to hours — even a 3% runtime improvement pays
back quickly.

Without measurement this remains a guessing game. This proposal
resolves it.

## §3 Method

1. **Pick the benchmark set.** Reuse the existing criterion harness
   at `pwdft/pwdft-core/benches/scf_benchmarks.rs`. Focus on the
   targets that actually exercise cross-crate hot paths:
   - `eigensolver` (faer-heavy)
   - `vnl_apply` (faer GEMM + pwdft-core assembly)
   - `fft_forward` / `fft_inverse` (ndrustfft)
   - one full-SCF integration bench if present; otherwise add a
     short Si SCF loop bench (20 iters, n_pw≈725, dense eigensolver)
     for this audit only.

2. **Measure baseline (thin).** Acquire the machine lock, run
   `cargo bench -p pwdft-core --bench scf_benchmarks`. Save the
   `target/criterion/` report as the thin baseline.

3. **Flip to fat.** In `Cargo.toml` set `lto = "fat"` under
   `[profile.release]`. Record clean `cargo build --release`
   wall-time (new number; distinct from bench runtime).

4. **Re-measure.** Re-run the same criterion bench set under fat.
   Criterion will surface per-target runtime deltas vs the thin
   baseline and flag statistical significance.

5. **Record production-run deltas.** Run one Si SCF at production
   cutoff (`inputs/si_scf.yaml`, ecutwfc = 200 eV or higher) under
   each setting, three times, record mean ± stdev.

6. **Decision rule.** Adopt `fat` iff:
   - Criterion reports a statistically significant (≥1σ) runtime
     improvement on ≥1 hot-path target, AND
   - Production SCF wall-time improves ≥2% (absolute), AND
   - Clean release link time stays under 5 minutes on an M3 Max.

   Otherwise keep `thin`. Either way, leave a numbers table in the
   logbook.

## §4 Acceptance

- A single logbook entry (`.claude/logbooks/perf.md`) with:
  - Clean release link time, thin vs fat.
  - Criterion output summary for the target set.
  - Si SCF wall-time, thin vs fat (3 runs each).
  - Decision + one-line justification.
- If decision is "adopt fat": `Cargo.toml` flip committed with
  commit message `LTOB: adopt fat LTO (<X>% faster on <target>)` and
  the logbook entry linked from the commit body.
- If decision is "keep thin": proposal archived with a note in the
  logbook; no code change.

## §5 Non-goals

- No other profile tuning (codegen-units, incremental, etc.) — those
  are independent variables and belong in separate audits.
- Not an optimization project. Purely "does the dial deserve
  turning?".
- Does not block on `pwdft-core`'s own benchmark coverage being
  complete — if a hot path isn't yet benched, note it and move on.
