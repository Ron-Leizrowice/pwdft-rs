---
id: CISP
status: active
priority: low
complexity: small
risk: low
depends_on: [TDBG]
blocks: []
---

# CISP: CI speedup follow-ups (mold, nextest, job parallelization)

## Problem

A caching sweep on `.github/workflows/ci.yml` landed alongside the initial CI rollout (PR #180) and took the obvious no-brainer wins:

- `CARGO_INCREMENTAL=0` (drop `.incremental/` fingerprint bloat in the cached `target/`)
- `RUSTFLAGS="-C debuginfo=0"` (strip DWARF from test and clippy artifacts)
- `cache-on-failure: true` on `Swatinem/rust-cache@v2`
- Workflow-level `concurrency` group to cancel superseded runs

TDBG then takes the next layer: switching Tier-1 `cargo test` to `--profile=dev` so LLVM doesn't optimize unit-scale test binaries.

That leaves three measured-but-not-landed levers from the CI-speedup plan, each independent and each worth considering *after TDBG's numbers come back*. They live here so they aren't lost.

## Research

### The three deferred levers

**A. `mold` linker.** `rui314/setup-mold@v1` installs mold in ~5 s, then `RUSTFLAGS="-C debuginfo=0 -C link-arg=-fuse-ld=mold"` routes rustc through it. mold is ~10× faster than the default `ld.bfd` on Linux for Rust workloads. Linking is a meaningful fraction of total wall-clock on this workspace (faer + wgpu + ndrustfft generate many objects per final binary — clippy all-targets alone produces dozens of test/bench/example executables).

Risk: low. mold is production-quality, drop-in, no runtime behavior change. The only failure mode is someone on a non-Linux runner trying to use it — our CI is `ubuntu-latest` only, so a non-issue.

**B. `cargo-nextest` for the test step.** `taiki-e/install-action@nextest` installs a prebuilt nextest binary. Swap `cargo test -p pwdft-core` → `cargo nextest run -p pwdft-core`. nextest parallelizes per-process (not just per-thread), catches hangs with a timeout, and produces cleaner output.

Risk: low-medium. nextest doesn't run doctests — if any doctest provides unique coverage not duplicated by unit tests, we need a second step `cargo test --doc -p pwdft-core`. Audit the doctest population before flipping; add the second step if warranted.

**C. Split the `rust` job into parallel jobs.** Today, clippy (default + gpu), cargo test, and cargo doc run sequentially in one job. Splitting into `clippy` / `test` / `doc` jobs drops wall-clock from sum-of-four to max-of-three (doc shares enough compile surface with clippy that it belongs in the clippy job).

Risk: medium. Each job holds its own `target/` cache, so total GitHub Actions cache storage grows ~3×. GitHub's per-repo cache ceiling is 10 GB; empirical current cache size after CICH tweaks will determine whether this is safe. Also, both clippy invocations (default + gpu) must stay in one job — they share ~90% of the build graph via `target/` and splitting them would double dependency compile time.

### Why measurement-driven, one at a time

The TDBG projection is "≥50% wall reduction on cold-build PR runs." If that number materializes, the Tier-1 test step drops to ~90 s cold and ~30 s warm. At that point:

- **mold** still pays on every build (linking happens regardless of opt-level). Independent of TDBG.
- **nextest** pays less when Tier-1 test compile time dominates wall-clock. Its win is on *run* time, which is already short. Reconsider after TDBG.
- **Job parallelization** pays proportionally to how unbalanced the steps are. Post-TDBG, clippy (dev, but -O3 deps still compile) likely dominates. Parallelization still helps, but may be a smaller win than expected.

Land them in the order of highest-confidence win first (mold), then measure, then decide on the others.

## Implementation

### Phase A — `mold` linker

Add a step before `Swatinem/rust-cache@v2` (mold install is fast enough that it doesn't meaningfully affect cache key timing):

```yaml
- uses: rui314/setup-mold@v1
```

Update the job-level env:

```yaml
env:
  CARGO_INCREMENTAL: "0"
  RUSTFLAGS: "-C debuginfo=0 -C link-arg=-fuse-ld=mold"
```

Keep the `RUSTDOCFLAGS` step-level env on the rustdoc step as-is. The link-arg is ignored by rustdoc (which doesn't link).

**Sanity check.** In a manual run, look for `-fuse-ld=mold` in `cargo build -v` output to confirm mold is actually wired in. If rustc silently ignores an unknown linker argument, we'd get no speedup with no error.

### Phase B — `cargo-nextest`

Audit doctests first:

```bash
cargo test --doc -p pwdft-core 2>&1 | tail -20
```

If the doctest count is zero or all duplicates of unit tests, swap the test step to:

```yaml
- uses: taiki-e/install-action@nextest

- name: cargo nextest (tier 1)
  run: cargo nextest run -p pwdft-core
```

If doctests have unique coverage, add a companion step:

```yaml
- name: cargo test --doc (tier 1)
  run: cargo test --doc -p pwdft-core
```

Interaction with TDBG: nextest honors cargo's `--profile` flag, so if TDBG has landed, the step becomes `cargo nextest run --profile=dev -p pwdft-core`.

### Phase C — Split into parallel jobs

Only pursue if Phases A+B leave a meaningful gap. Sketch:

```yaml
jobs:
  clippy:
    steps:
      - clippy -p pwdft-core --all-targets
      - clippy -p pwdft-core --all-targets --features gpu
      - cargo doc --no-deps -p pwdft-core  # shares compile surface with clippy
  test:
    steps:
      - cargo nextest run -p pwdft-core  # or --profile=dev per TDBG
```

Each job gets its own `Swatinem/rust-cache@v2` with a distinct `shared-key` so caches don't collide. Concurrency group already handles cancellation. The existing `validation` job (Python) stays as-is.

Before landing Phase C, check the post-CICH cache size from a real run's "Post Swatinem" log line. If it's > 3 GB, parallelization would push us into cache-eviction territory on the 10 GB ceiling.

## Verification

Each phase lands in its own PR with a before/after wall-clock comparison pulled from Actions:

- Cold-cache run: compare "Set up job" → "Complete job" total on a cache-invalidation commit (e.g., touching `Cargo.lock`).
- Warm-cache run: same measurement on a second consecutive push that hits the restored cache.

Phase A success: ≥10% wall-clock reduction on the total rust-job timing, with `-fuse-ld=mold` visible in a `cargo build -v` manual run.

Phase B success: either the doctest audit confirms a clean swap (nextest replaces cargo test with no coverage loss) or the two-step variant works. Measure 2-3× test-step runtime reduction on tests that run in practice.

Phase C success: wall-clock of the new slowest job (probably `clippy`) is at least 40% below the pre-split `rust` job total. If the new max is within 20% of the old sum, parallelization is not paying off and the extra cache storage isn't worth it — revert to a single job.

## Non-goals

- **`sccache`** — TDBG mentions it as a fallback if nothing else pans out. Not pursued here unless A+B+C land and still aren't enough. Adds moving parts (remote cache backend, auth, eviction) that `Swatinem/rust-cache@v2` already covers for 90% of cases.
- **Paid larger runners** (`ubuntu-latest-8-cores`) — real cost, marginal win after mold + parallelization on a 4-core runner.
- **Pinning the Rust toolchain** — `@stable` only invalidates on new stable releases (every 6 weeks); not worth the maintenance.
- **Changes to profile definitions in `Cargo.toml`** — TDBG's territory; CISP stays on the workflow side.
