---
id: SKPL
status: completed
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# SKPL: `/test --tier2` pre-applies the authoritative skip list

## Problem

Tier-2 cargo tests carry two very different `#[ignore]` reasons that today look identical to the `/test --tier2` skill:

1. **Heavy-but-passing** — `TSPL Tier-2: ...` (runs Si SCF twice, ecut=100, 40 iters, etc.). These are supposed to execute and pass. Typical wall: seconds to tens of seconds each.
2. **Physics-blocker / designed-to-fail** — `VGCH-MECH Class A: ...`, `VGCH-MECH Class B: ...`, MXBA Fe LDA adaptive-β failure, and the `test_plain_anderson_stalls_on_c_diamond` documented-stall case. These start real SCFs and most run to convergence before their reason-string assertion fires; several never converge (mixer stalls, heavy-atom residuals) and consume minutes of wall apiece.

Today `/test --tier2` just runs `cargo test -- --ignored`. Agents are expected to read every `#[ignore]` reason string and hand-craft a `--skip test_name` for each designed-to-fail case. In practice:

- The 2026-04-21 ROTI recovery agent skipped exactly one test (`test_plain_anderson_stalls_on_c_diamond` — the only one that genuinely stalls) and ran all other physics-blockers to completion, costing ~12 extra min of wall beyond the necessary Tier-2 wall.
- The `/test --tier2` skill's SKILL.md points agents at the `#[ignore]` reason strings as the "authoritative skip list," but there is no programmatic link. The skip list drifts every time a new physics-blocker lands.

**Concrete instance.** ROTI PR #184's Tier-2 was 19 min on cold cache. Removing the 5 `VGCH-MECH Class A` cells + Fe-LDA Class B + MXBA would have saved ~8–10 min of that, leaving ~9 min which is all cold compile + genuine Tier-2 wall. The skip decision was structural — ROTI is a type-width revert; every designed-to-fail cell is irrelevant to it — but the agent couldn't express that without hand-coding six `--skip` flags.

## Proposal

Two changes that together automate the skip list:

### Part A — tag physics-blocker `#[ignore]` reasons with a stable prefix

Update the `#[ignore]` reason strings in `pwdft/pwdft-core/tests/qe_validation.rs`, `pwdft/pwdft-core/tests/mxba_adaptive_beta_fe.rs`, and any similar file to prepend `SKIP-TIER2:` to the reason for any test that is *designed to fail or stall* under current `main`. E.g.:

```rust
#[ignore = "SKIP-TIER2 VGCH-MECH Class A: Fe LDA +11.14 eV ..."]
#[ignore = "SKIP-TIER2 mixer-stall: test_plain_anderson_stalls_on_c_diamond ..."]
#[ignore = "SKIP-TIER2 VGCH-MECH Class B: Fe LDA BSUM 0.99× ratio ..."]
```

Regular Tier-2 heavy-but-passing tests keep their existing `TSPL Tier-2: ...` prefix without the `SKIP-TIER2` marker. The marker is the contract.

**Enumerated scope** (as of 2026-04-21): six `VGCH-MECH Class A` cells (Fe/GaAs/Cu/NaCl/MgO/C-diamond LDA), one `VGCH-MECH Class B` cell (Fe LDA BSUM), Fe LDA PBE Class A, `test_plain_anderson_stalls_on_c_diamond`, the MXBA Fe LDA test. Total ≈ 9–10 tests. Specific line numbers enumerated in the PR body at implementation time.

### Part B — `/test --tier2` reads the marker and auto-`--skip`s

Update `.claude/skills/test/SKILL.md` invocation pseudocode:

```bash
# Build skip list from ignore-reason markers.
SKIP_FLAGS=$(rg -n '#\[ignore = "SKIP-TIER2[^"]*"\]' pwdft/pwdft-core/tests/ \
             --type rust --no-heading \
             | python3 -c '
               import sys, re
               names = set()
               path_cache = {}
               for line in sys.stdin:
                   path, lineno, _ = line.split(":", 2)
                   # read the file once, find the fn name under this #[ignore] line
                   if path not in path_cache:
                       path_cache[path] = open(path).read().splitlines()
                   src = path_cache[path]
                   i = int(lineno)
                   # find the next "fn <name>" after this line
                   while i < len(src) and not src[i].lstrip().startswith("fn "):
                       i += 1
                   m = re.match(r"\s*fn\s+(\w+)", src[i])
                   if m: names.add(m.group(1))
               for n in sorted(names):
                   print(f"--skip {n}")
             ')

.claude/bin/machine-lock run "<role>" "cargo test tier-2" -- \
    cargo test -- --ignored $SKIP_FLAGS
```

Or in a cleaner pure-bash variant that does not need the Python post-processor — use `awk` to walk the source from the ignore-line to the next `fn`:

```bash
SKIP_FLAGS=$(awk '/SKIP-TIER2/{found=1; next} found && /^\s*(#\[|pub\s+)?fn\s+[a-zA-Z_]+/{match($0, /fn [a-zA-Z_]+/); print "--skip " substr($0, RSTART+3, RLENGTH-3); found=0}' pwdft/pwdft-core/tests/*.rs | sort -u)
```

Pick whichever is more readable to maintainers; both are ~5 lines. The skill's SKILL.md documents the extraction logic and the `SKIP-TIER2` marker contract.

### Part C — CI consistency check (optional, cheap)

A small `scripts/check-tier2-skip-consistency.sh`:

- Greps `SKIP-TIER2` markers and verifies each tagged test corresponds to a real `#[test] fn` declaration immediately below.
- Prints the enumerated skip list; agents can `diff` this against what `/test --tier2` actually applies.

Run from the existing `.github/workflows/rust.yml` as a fast pre-clippy step (~2 s). Catches drift where someone adds a physics-blocker without the marker.

### CLI surface (no change)

`/test --tier2` remains a single command that Does The Right Thing. `/test --tier2 --no-skip-list` as an escape hatch runs the raw `cargo test -- --ignored` for when an agent explicitly wants to see a physics-blocker surface (e.g. debugging VGCH-MECH itself).

## Risk

**Low.**

- Part A is pure documentation (comment text in `#[ignore = "..."]` strings). Zero runtime impact.
- Part B is a shell wrapper. Worst case: buggy extraction emits wrong `--skip`, tests run as they do today. Smoke test: `scripts/check-tier2-skip-consistency.sh` enumerates what the wrapper would skip before invoking cargo.
- Part C is an additive CI step. If it false-positives it flags a PR; trivial to fix.

Edge cases to verify:

- Test with both `#[cfg(feature = "gpu")]` and `#[ignore = "SKIP-TIER2 ..."]` should still be skipped (the `cfg` gate handles feature exclusion separately; `--skip` is orthogonal).
- Parametrized tests (if any) where the skip list should cover every instantiation. The awk extraction walks to the next `fn`, which handles `#[test_case(...)]` / `rstest` only if they produce a single `fn`. Validate against any existing parametric Tier-2 tests.

## Non-goals

- Not re-categorizing tests. If a current `VGCH-MECH Class A` is reclassified (e.g. PZPW closes Fe LDA Class B), the `SKIP-TIER2` marker is dropped as part of that PR. Normal lifecycle.
- Not changing Tier-1 (default) semantics.
- Not adding a declarative test-matrix DSL. The `#[ignore = "..."]` reason string is already authoritative; we just parse it.

## Implementation plan

One PR:

1. Grep current Tier-2 ignore reasons; tag physics-blockers with `SKIP-TIER2` prefix.
2. Update `.claude/skills/test/SKILL.md` with the auto-skip invocation.
3. Add `scripts/check-tier2-skip-consistency.sh` + CI step in `rust.yml`.
4. Update CLAUDE.md § "Tier-2 PR policy" paragraph to point at the marker contract.

Deletes: the ROTI session's hand-crafted `--skip test_plain_anderson_stalls_on_c_diamond` from any future Tier-2 invocation — the skill handles it automatically.

## Acceptance

- `/test --tier2` with no arguments runs `cargo test -- --ignored` with all `SKIP-TIER2`-tagged tests auto-excluded. Wall-time on warm cache ≈ 30 s (the TSPL Tier-2 suite without the designed-to-fail long-runners).
- Agents invoking `/test --tier2` never need to hand-craft `--skip` flags.
- `scripts/check-tier2-skip-consistency.sh` green in CI.
- CLAUDE.md cross-reference to the marker contract in place.

## Measurement

Pre/post wall-time on warm-cache M3 Max (current `main`):

- Today: `/test --tier2` with no manual skip → ≈ 4–6 min wall (physics-blockers run SCFs to convergence before asserting).
- Today: `/test --tier2` + agent remembers `--skip test_plain_anderson_stalls_on_c_diamond` → ≈ 3–4 min wall (one stalling test excluded, others still run).
- Today: `/test --tier2` + agent hand-crafts the full designed-to-fail skip list → ≈ 45–60 s wall. Few agents do this correctly.
- Post-SKPL: `/test --tier2` auto-skips → ≈ 45–60 s wall, no agent judgment required.

The third row is what CLAUDE.md advertises as "Tier-2 ≈ 58 s warm-cache." SKPL makes that the default, not the aspirational case.

## Flagged for follow-up

- If Part C's consistency check turns up tests that lack a clear category (`SKIP-TIER2` vs `TSPL Tier-2` heavy-but-passing), file as a mini-audit at landing time.
- Once PZPW closes Fe LDA Class B, the corresponding marker comes off. Acts as a progress signal.
