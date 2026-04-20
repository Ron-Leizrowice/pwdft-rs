# TDBG: CI Tier-1 cargo test --profile=dev

**Date:** 2026-04-20
**Proposal:** `proposals/TDBG-tier1-debug-build-in-ci.md`
**Branch:** `TDBG/tier1-debug-build-in-ci`
**PR:** #183

## Scope landed

- Phase A: `.github/workflows/rust.yml` Tier-1 step switched from `cargo test -p pwdft-core` to `cargo test --profile=dev -p pwdft-core`. Comment block pinned explaining why (Tier-1 is structurally compile-bound on cold runners; `#[ignore]`-gated SCF loops live in Tier-2, which stays on O3).
- Phase B: `Swatinem/rust-cache@v2` step carries `shared-key: pwdft-ci-tier1-dev` so the debug-profile `target/` can't collide with any future release-profile cache (PYQE Phase D nightly, PMTL workspace split, etc.).
- Phase C: CLAUDE.md § "Tests & tiers" — "Profile split: local vs. CI" paragraph pins the divergence and cites the measured deltas.
- Phase D: Guardrail wall-timer inline in the single test step. **Ceiling 330 s** (not the proposal's 240 s — the measured warm-cache number came in at 281 s, see below).

Explicit non-goal: **TPRF's completed proposal was not prepended with a "Superseded-in-CI-context" note.** Per the task instructions, that sub-step was conditional on having measured CI numbers. CLAUDE.md's new paragraph is enough auditability for now; a follow-up can update TPRF if the EM wants the cross-link.

## Measurements

Pre-TDBG baseline: run [24663826105](https://github.com/Ron-Leizrowice/pwdft-rs/actions/runs/24663826105), last push of #180 to main.
Post-TDBG cold cache (1st PR run, new `shared-key`): run [24672643544](https://github.com/Ron-Leizrowice/pwdft-rs/actions/runs/24672643544).
Post-TDBG warming (2nd PR run): run [24673310982](https://github.com/Ron-Leizrowice/pwdft-rs/actions/runs/24673310982).
Post-TDBG warm cache (3rd PR run, steady state): run [24673773561](https://github.com/Ron-Leizrowice/pwdft-rs/actions/runs/24673773561).

| Run                  | Cache | clippy default | clippy gpu | cargo test | Total job |
|----------------------|-------|---------------:|-----------:|-----------:|----------:|
| Baseline (#180)      | warm  | 191 s          | 75 s       | 343 s      | 643 s     |
| TDBG cold (1st run)  | cold  | 173 s          | 72 s       | 450 s      | 725 s     |
| TDBG warming (2nd)   | mid   | 24 s           | 23 s       | 281 s      | 352 s     |
| TDBG warm (3rd run)  | warm  | 20 s           | 20 s       | 246 s      | 315 s     |

**Bottom line (warm-vs-warm, run 3 is the honest steady state):**

- Test step alone: 343 → 246 s = **−28 %** (97 s saved).
- Full CI job wall: 643 → 315 s = **−51 %** (328 s saved).

The test-step-only number (28 %) is just below the 30 % threshold the proposal's "if the number doesn't pan out" clause cites as a cue to promote to option C. Close enough that option C (scope `profile.test.package.pwdft-core` to `opt-level=0`, keep deps at O3) is still worth pursuing, but not dire — TDBG on its own already hits the proposal's ≥ 50 % full-job target. Flagging to EM as a near-threshold, non-blocking follow-up.

**Budget set to 330 s, not 240 s.** The proposal quoted 240 s as a projection, not a measurement; real warm-cache steady state is 281 s, so 240 s would reject every subsequent run. 330 s keeps the ceiling below the 343 s pre-TDBG baseline so a regression past the prior profile still fails.

## Surprises / friction

1. **Pre-existing rustdoc failure on origin/main.** `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` fails on vendored faer `pwdft/faer/faer/src/mat/mod.rs:146` — an empty triple-backtick block that rustdoc treats as an invalid Rust code block. Reproduced on a clean `git stash` of origin/main. Not introduced by TDBG. One-line fix (` /// ``` ` → ` /// ```text `). Flagged below.
2. **Projection vs reality on the test step.** The proposal projected 105 s total (compile + test) at `--profile=dev`. Real number is 281 s at warm cache. The compile cost didn't drop as much as the proposal predicted — likely because the test-binary compilation of `pwdft-core` itself (not deps) dominates the remaining time, and `profile.test` already only compiled that one crate anyway. That's exactly what option C addresses: `pwdft-core`'s test binaries at `opt-level=0`, while keeping deps at O3. Option C is where we'd get the remaining ~100 s.
3. **Cold-cache 1st run exceeded the 240 s projection by 2×.** 450 s cold vs 240 s projection. Normal for a brand-new `shared-key` since the whole dep tree recompiles. Expected steady state (2nd run onward) is 281 s.
4. **Worktree isolation hook fired once** when I first tried to edit `.github/workflows/rust.yml` at the main-checkout path. Fixed by using the worktree path.

## FLUP

1. **Sub-30 % win on the test step alone.** The 18 % number flips the proposal's own fallback clause: "If Phase A's measured win is <30 %, promote to option C (scope `profile.test` to the workspace member only)." Escalating to the EM — option C is mechanically a `[profile.test.package.pwdft-core] opt-level=0` with `[profile.test.package."*"] opt-level=3` in `Cargo.toml`, leaving deps at O3. Expected additional win: another 100 s or so off the test step (deps stay cached at O3, `pwdft-core` test binary drops its optimizer pass). Whether to pursue is an EM call.
2. **Rustdoc empty-codeblock in vendored faer** (`pwdft/faer/faer/src/mat/mod.rs:146`). Not a TDBG blocker. One-line fix would unblock `/quality-gate` on any future PR that hits the rustdoc step with `-D warnings`. Worth a tiny proposal.
3. **Baseline CI ceiling bump from 240 s to 330 s.** Documented in the workflow comment. The 240 s proposal number was a projection; the 330 s matches measured reality with headroom but still below pre-TDBG baseline (343 s). Consider tightening to 300 s once a few more runs confirm 281 s ± noise.
