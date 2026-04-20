---
name: lint
description: Run both clippy invocations under the machine lock — auto-fix pass plus default-feature and --features gpu warning sweeps. Trigger on "/lint".
user_invocable: true
---

# /lint — Dual clippy sweep

Run the three clippy passes under the machine lock, in order:

1. `cargo clippy -q --fix --allow-dirty --allow-staged --all-targets` — auto-fix first.
2. `cargo clippy -q --all-targets` — default-feature warnings.
3. `cargo clippy -q --all-targets --features gpu` — GPU-feature warnings.

Both default and GPU runs are required: without `--features gpu` the `pwdft/pwdft-core/src/gpu/` source tree and GPU-only test binaries are not linted, and warnings accumulate silently.

## Invocation

```bash
.claude/bin/machine-lock run "<role>" "clippy" -- bash -c '
  cargo clippy -q --fix --allow-dirty --allow-staged --all-targets &&
  cargo clippy -q --all-targets &&
  cargo clippy -q --all-targets --features gpu'
```

## Rules

- Warning count must go down, never up, across both clippy invocations.
- Do not `#[allow]` code-smell warnings like `too_many_arguments` — refactor the code.
- The auto-fix pass can rewrite your working tree. Run `git diff` after to see what changed.

## When to call the full gate instead

If you're about to open a PR, use `/quality-gate` — it runs `/lint` plus rustdoc and the Tier-1 tests. `/lint` alone is for mid-work hygiene.

## On failure

Report the failing pass and its output. Don't re-run automatically.
