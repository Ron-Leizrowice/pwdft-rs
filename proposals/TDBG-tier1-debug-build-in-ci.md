---
id: TDBG
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# TDBG: Run Tier-1 tests with debug opt-level in CI to shrink compile time

## Problem

CI compile time dominates Tier-1 runtime. The Tier-1 suite is structurally designed around unit checks and single-shot kernel validations (see CLAUDE.md § Test suite tiers: "every test here is either a pure unit check or a single-shot operation on a small matrix — no test runs an SCF loop for more than a handful of iterations at production `n_pw`"). These tests do **not** need `-O3` to finish in seconds — but the current cargo profile forces `opt-level = 3` for every test build because TPRF (`[profile.test] opt-level=3`) optimized for *local* `cargo test` runtime on a warm-cached M3 Max.

The CI machine is neither warm-cached nor M3 Max. It's a GitHub Actions runner doing a cold build every run, and the `opt-level = 3` test profile makes it pay rustc's optimizer cost (~60–70% of cold-build wall-time on a hot-path numerics crate like `pwdft-core`) for unit tests that would run in milliseconds at `opt-level = 0`.

**Concrete state.**

- `Cargo.toml` sets `[profile.test] opt-level = 3` and `[profile.test.package."*"] opt-level = 3`. TPRF's motivation (`proposals/completed/TPRF-test-profile-opt-level.md`) was *local test wall-time* (11 min → 95 s on M3 Max). In CI, the ratio inverts: **build time dwarfs run time**.
- Previous CI run `24651829658` took 3 min 16 s total; of that, the clippy (full workspace, `-O3`) step was the dominant cost. The test step has historically pushed CI past 5 min on cold builds.
- Tier-1 wall-time on M3 Max is 12 s (warm). On a 4-vCPU runner at `opt-level=0` it should be ~20–30 s to *run*, plus ~60–90 s to compile (vs. ~3–4 min at `opt-level=3`).

**Why now.** Two forcing functions:

1. **PYQE gate lands a QE CI job** (~2 min for the fast tier). Total per-PR CI budget is creeping up. Tier-1 must get cheaper or every PR becomes a 5+ min wait.
2. **PMTL splits `pwdft-metal`** into a second crate. Workspace `cargo test` will compile both — a second 40-crate tree through `-O3` — unless we switch CI to debug for Tier-1.

## Research

### What TPRF was actually solving

TPRF set `[profile.test] opt-level=3` because `pwdft-core`'s hot-path unit tests (FFT round-trips, Simpson quadrature smoke tests) hit inner loops that are painfully slow at `opt-level=0`. But "painfully slow" on a warm M3 Max cache means **still sub-millisecond per test** — the 11-minute pre-TPRF number came mostly from `indicatif`/`faer`/`ndrustfft` dependency compile time plus the stray SCF-loop test that pre-TSPL was still in Tier 1.

Post-TSPL (2026-04-19), every SCF-loop test moved to Tier 2 (`#[ignore]`, `cargo test -- --ignored`). The remaining Tier-1 tests are:

- FFT forward/inverse round-trips on 18³ and 24³ grids
- VGCMP Phase 1–4 cross-checks (form-factor arithmetic; no SCF loop)
- Free-electron bands (analytic; no SCF loop)
- KB projector validation (Bessel transform; scipy-parity pins)
- Non-local symmetry checks (structure-factor arithmetic)
- ITEV single-shot eigenvalue check on a defect-1 matrix
- Ewald + LAPACK smoke
- ALOC F-5 cache check

None of these require `-O3` to finish in CI wall-time budget. The hot-path concern TPRF addressed is a Tier-2 problem now.

### Measurement plan

Before committing, measure. On a fresh CI-like machine (or `act pull_request` locally with a cold cache):

| Config | Compile | Tier-1 test wall | Total |
|---|---|---|---|
| `cargo test -p pwdft-core` (current, `O3`) | ~210 s | ~20 s | ~230 s |
| `cargo test -p pwdft-core --profile=dev` (`opt-level=0`, no LTO) | ~70 s | ~35 s | ~105 s |
| `cargo test -p pwdft-core --profile=dev-ci` (`opt-level=1`, no LTO) | ~110 s | ~25 s | ~135 s |

(Numbers are rough projections based on 4-vCPU runner throughput and the TPRF pre-bench ratios.) The middle row — straight debug — is the most likely winner for Tier-1 in CI.

**Tier 2 stays on `-O3`.** A Tier-2 SCF loop at `opt-level=0` would be glacial (minutes per test). Tier-2 runs nightly or on-label, where latency is not user-facing.

### Options

| Option | How | Pros | Cons |
|---|---|---|---|
| A. `--profile=dev` in CI | workflow uses `cargo test --profile=dev` for Tier 1 | zero code change | `dev` profile is also used locally by `cargo build` — not a problem, but couples two contexts |
| B. New `ci-fast` profile | add `[profile.ci-fast] opt-level=0, inherits="dev"` | explicit, documented | adds profile surface; name-bikeshedding |
| C. Scope `profile.test` to workspace only | `profile.test.package.pwdft-core = { opt-level=0 }` then `"*"=3` for deps | targeted | deps still take O3 time in CI; half the win |

Recommendation: **A** for the fastest landing. Tier-1 CI uses `cargo test --profile=dev`; `cargo test` locally remains `-O3` (via `profile.test`) so the M3 Max numbers in CLAUDE.md stay accurate. No code churn, no new profile. If later we want both contexts to share a config, promote to option B.

## Implementation

### Phase A — Switch CI Tier-1 to `--profile=dev`

1. `.github/workflows/ci.yml`:

   ```yaml
   - name: cargo test (tier 1, debug build)
     run: cargo test --profile=dev -p pwdft-core
     # Tier 1 is structurally compile-bound, not runtime-bound.
     # Local `cargo test` stays on the O3 test profile (see Cargo.toml).
   ```

2. Clippy stays on the default `cargo clippy -p pwdft-core --all-targets` (`dev` profile) — that's already the fastest path; no change needed there.

3. Rustdoc step stays as-is.

### Phase B — Cache the debug `target/` separately

`Swatinem/rust-cache@v2` keys on profile by default, so the dev and release caches are isolated. One nit: make sure the cache key string includes the workflow step ID so the Tier-1 dev cache isn't polluted by the Tier-2 release cache when the nightly workflow (PYQE Phase D) runs.

   ```yaml
   - uses: Swatinem/rust-cache@v2
     with:
       shared-key: pwdft-ci-tier1-dev
   ```

### Phase C — Document the policy

1. CLAUDE.md § Tests & Benchmarks: add a note that CI Tier-1 uses `--profile=dev` while local runs use the O3 test profile; explain when each is appropriate. Include the measured wall-time delta from Phase A so the choice is auditable.

2. `proposals/completed/TPRF-test-profile-opt-level.md` — prepend a "Superseded-in-CI-context" note pointing here. The *local* TPRF win stays intact.

### Phase D — Guardrail: runtime regression check

Add a one-line assertion in the CI job that Tier-1 stays under a budget:

   ```yaml
   - name: Tier-1 runtime budget
     run: |
       start=$(date +%s)
       cargo test --profile=dev -p pwdft-core
       elapsed=$(( $(date +%s) - start ))
       if (( elapsed > 240 )); then
         echo "::error::Tier-1 CI wall-time regressed to ${elapsed}s (budget: 240s)"
         exit 1
       fi
   ```

   (Merge with the existing step; don't double-run.) Budget starts at 4 min — comfortable ceiling; tighten on evidence.

### Phase E — Non-goals explicitly called out

- **No change to `cargo test` locally.** CLAUDE.md § Tests & Benchmarks still says `cargo test` for Tier 1. Warm M3 Max stays at 12 s.
- **No change to Tier 2.** `cargo test -- --ignored` continues to use the O3 test profile. The whole point of Tier 2 is production-scale SCF loops, which need optimization.
- **No change to benches.** `cargo bench` already uses the release profile.

**If the number doesn't pan out.** The projections above are ballpark. If Phase A's measured win is <30% of cold-build time, promote to option C (scope `profile.test` to the workspace member only). If still insufficient, keep the O3 Tier-1 and instead invest in `sccache`/`cargo-nextest`/`mold` — these are follow-up options, not blockers for TDBG to land.
