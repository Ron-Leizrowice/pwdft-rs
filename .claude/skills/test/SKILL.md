---
name: test
description: Run cargo test under the machine lock with tier awareness. Default runs Tier-1 (fast); pass "--tier2" for Tier-2 heavy SCF suites, "--all" for both, or any other cargo test args. Trigger on "/test".
user_invocable: true
---

# /test — Tier-aware test runner

Parse `$ARGUMENTS` and run `cargo test` under the machine lock.

## Argument handling

- **empty** → `cargo test` (Tier-1).
- `--tier2` → Tier-2 with auto-skip (see § Invocation). Physics-blocker tests tagged `SKIP-TIER2` are excluded automatically.
- `--tier2 --no-skip-list` → raw `cargo test -- --ignored` with no auto-skipping. Use this only when you explicitly want physics-blockers to surface (e.g. debugging VGCH-MECH). Expect 4–6 min wall on warm cache.
- `--all` → `cargo test -- --include-ignored` (both tiers).
- `--features gpu [extra]` → `cargo test --features gpu [extra]`.
- anything else (e.g. `test_name`, `--nocapture`, `--test free_electron_bands`) is passed through verbatim after `cargo test`.

## Invocation

Pick the role hint that matches your agent.

**Tier-1 (default):**
```bash
.claude/bin/machine-lock run "<role>" "cargo test" -- cargo test
```

**Tier-2 with auto-skip (default since SKPL):**
```bash
SKIP_FLAGS=$(awk '
  /SKIP-TIER2/ { found=1; in_str=(/"\]$/ ? 0 : 1); next }
  found && in_str { in_str=(/"\]$/ ? 0 : 1); next }
  found && /^[[:space:]]*(#\[|pub[[:space:]]+)?fn[[:space:]]+[a-zA-Z_]+/ {
    match($0, /fn [a-zA-Z_]+/)
    print "--skip " substr($0, RSTART+3, RLENGTH-3)
    found=0
  }
  found && /^[[:space:]]*(#\[|\/\/)/ { next }
  found { exit 1 }
' pwdft/pwdft-core/tests/*.rs | sort -u)
.claude/bin/machine-lock run "<role>" "cargo test tier-2" -- \
  cargo test -- --ignored $SKIP_FLAGS
```

**Tier-2 escape hatch (no auto-skip, runs physics-blockers):**
```bash
.claude/bin/machine-lock run "<role>" "cargo test tier-2 --no-skip-list" -- \
  cargo test -- --ignored
```

## Tier-2 reminders

By default, `/test --tier2` auto-excludes every `#[ignore]` test whose reason string begins with `SKIP-TIER2`. These are tests that are *designed to fail or stall* on current `main` (VGCH heavy-atom cells, Fe PBE, MXBA, Plain Anderson stall). They are not bugs in the test runner — they document open physics residuals. The skip is automatic; you never need to hand-craft `--skip` flags.

Tests tagged `TSPL Tier-2: ...` (without the `SKIP-TIER2` prefix) are heavy-but-passing SCF runs. They always execute under `/test --tier2` and must pass.

The authoritative skip list is `scripts/check-tier2-skip-consistency.sh` — run it to enumerate what would be skipped before invoking cargo.

**Tier-2 PR policy:** any PR that touches `pwdft/pwdft-core/src/{scf,potential,symmetry,pseudopotential,eigensolver,gpu}/`, `basis.rs`, `fft.rs`, `ewald.rs`, `crystal.rs`, `kpoints.rs`, or bumps a numerics dep in Cargo.toml must run `/test --tier2` and include the outcome in the PR body's Test Plan.

## Runtime expectations

Warm-cache, M3 Max: Tier-1 ≈ 12 s wall, Tier-2 ≈ 45–60 s wall (with auto-skip, default since SKPL). Without auto-skip (`--no-skip-list`), physics-blockers run SCFs to convergence before asserting failure — expect 4–6 min wall.

**Cold-cache cost, especially after a rebase.** If you just ran `git reset --hard` + `git rebase origin/main`, the worktree's `target/` is invalidated and the first Tier-2 invocation will pay the full workspace re-compile (~5–10 min on M3 Max) *before* the ~58 s of actual test wall lands. Plan for 15–20 min total on a rebased worktree; don't dispatch Tier-2 + wait-in-a-poll-loop as if it's a 1 min task. Run one Tier-1 first to warm the compile cache if you plan to invoke cargo many times in the same session — each cargo call amortizes the previous one's compile work.

## On failure

Report the failing test(s) and their output. Do not re-run automatically.
