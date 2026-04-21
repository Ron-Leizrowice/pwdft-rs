# 2026-04-21 — SKPL: Tier-2 physics-blocker skip list automation

## Shipped

All 4 parts implemented (session terminated by stream-idle watchdog before commit;
EM completed administrative steps):

- **Part A**: Tagged 9 physics-blocker tests with `SKIP-TIER2` prefix in `#[ignore]` reason strings:
  - `qe_validation.rs`: C diamond, Fe LDA, GaAs, Cu, NaCl, MgO, Fe BCC PBE (7 tests)
  - `mxba_adaptive_beta_fe.rs`: MXBA Fe adaptive-β failure (1 test)
  - `mixer_robustness.rs`: Plain Anderson stall on C diamond (1 test)
- **Part B**: `SKILL.md` updated — `--tier2` invocation now uses awk to extract `SKIP-TIER2`-marked test names and auto-passes `--skip` flags. Escape hatch `--no-skip-list` runs raw `--ignored`.
- **Part C**: `scripts/check-tier2-skip-consistency.sh` — enumerates what would be skipped; exits 1 on orphaned markers. CI step added to `rust.yml` (fast pre-clippy check).
- **Part D**: `CLAUDE.md` Tier-2 section updated with one sentence describing the `SKIP-TIER2` marker convention.

## Quality gate

[Pending — EM dispatched quality gate + PR agent after stream-idle recovery]

## FLUP

None.
