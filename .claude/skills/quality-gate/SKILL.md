---
name: quality-gate
description: Run the full merge-criterion quality gate under the machine lock (fix clippy, default clippy, gpu clippy, rustdoc with -D warnings, and the Tier-1 test suite). Trigger on "/quality-gate". This is the single check every PR must pass before the EM will merge.
user_invocable: true
---

# /quality-gate — Full merge-criterion gate

Run all five checks that the Engineering Manager requires on every PR. The `$ARGUMENTS` slot is ignored (no knobs — this is either green or not).

## What it runs

In order, short-circuiting on the first failure:

1. `cargo clippy -q --fix --allow-dirty --allow-staged --all-targets` — auto-fix and re-lint the default feature set.
2. `cargo clippy -q --all-targets` — default-feature warnings.
3. `cargo clippy -q --all-targets --features gpu` — lint the GPU source tree and GPU-only test binaries (silently missed otherwise).
4. `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` — rustdoc must be warning-clean. The `-D warnings` flag must travel through the env var; current cargo rejects it after `--`.
5. `cargo test` — Tier-1 suite.

## Invocation

Pick the role hint that matches your agent (Core Engineer, Performance Engineer, Researcher, Code Reviewer, Technical Writer, or Agent).

```bash
.claude/bin/machine-lock run "<role>" "quality gate" -- bash -c '
  cargo clippy -q --fix --allow-dirty --allow-staged --all-targets &&
  cargo clippy -q --all-targets &&
  cargo clippy -q --all-targets --features gpu &&
  RUSTDOCFLAGS="-D warnings" cargo doc --no-deps &&
  cargo test'
```

## Rules

- Do **not** suppress clippy warnings (`#[allow(...)]` on a code smell). Refactor instead.
- Do **not** `#[allow]` rustdoc warnings. Fix the prose — escape brackets, remove links at private items.
- Tier-2 tests are separate — see the `/test` skill. The Tier-2 PR policy applies when the diff touches `pwdft/pwdft-core/src/{scf,potential,symmetry,pseudopotential,eigensolver,gpu}/`, `basis.rs`, `fft.rs`, `ewald.rs`, `crystal.rs`, `kpoints.rs`, or bumps a numerics dep.

## On failure

Report the failing step, its stderr, and stop. Do not re-run automatically — the user needs the output.
