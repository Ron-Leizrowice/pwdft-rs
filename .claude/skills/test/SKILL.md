---
name: test
description: Run cargo test under the machine lock with tier awareness. Default runs Tier-1 (fast); pass "--tier2" for Tier-2 heavy SCF suites, "--all" for both, or any other cargo test args. Trigger on "/test".
user_invocable: true
---

# /test — Tier-aware test runner

Parse `$ARGUMENTS` and run `cargo test` under the machine lock.

## Argument handling

- **empty** → `cargo test` (Tier-1).
- `--tier2` → `cargo test -- --ignored` (Tier-2 heavy SCF suites).
- `--all` → `cargo test -- --include-ignored` (both tiers).
- `--features gpu [extra]` → `cargo test --features gpu [extra]`.
- anything else (e.g. `test_name`, `--nocapture`, `--test free_electron_bands`) is passed through verbatim after `cargo test`.

## Invocation

Pick the role hint that matches your agent.

```bash
.claude/bin/machine-lock run "<role>" "cargo test <args>" -- cargo test <resolved-args>
```

## Tier-2 reminders

`cargo test -- --ignored` includes the pre-TSPL physics-blocker ignores (VGCH heavy-atom cells, Al ecut, C mixer, MXBA) that fail by design. The authoritative skip list is in the `#[ignore]` reason strings in `pwdft/pwdft-core/tests/qe_validation.rs` and `pwdft/pwdft-core/tests/mxba_adaptive_beta_fe.rs`. Report outcomes honestly in the PR body — do not hide expected failures.

**Tier-2 PR policy:** any PR that touches `pwdft/pwdft-core/src/{scf,potential,symmetry,pseudopotential,eigensolver,gpu}/`, `basis.rs`, `fft.rs`, `ewald.rs`, `crystal.rs`, `kpoints.rs`, or bumps a numerics dep in Cargo.toml must run `/test --tier2` and include the outcome in the PR body's Test Plan.

## Runtime expectations

Warm-cache, M3 Max: Tier-1 ≈ 12 s wall, Tier-2 ≈ 58 s wall (excluding the designed-to-fail ignores above).

## On failure

Report the failing test(s) and their output. Do not re-run automatically.
