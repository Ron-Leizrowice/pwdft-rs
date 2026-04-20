---
id: ESPL
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# ESPL: Split `ElectronSettings` — system physics vs. convergence knobs; drop default `max_iter` to 50

## Problem

`src/settings.rs:209` defines `ElectronSettings` as a single struct that mixes two orthogonal kinds of parameters:

| Physics (what system are we solving) | Convergence (how we solve it) |
|---|---|
| `nspin` | `mixing_beta` |
| `starting_magnetization` | `mixing_ndim` |
| `tot_magnetization` | `smearing` |
| `occupations` | `smearing_width` |
|  | `mixing_mode` |
|  | `pulay_period` |
|  | `adaptive_beta` |

These have orthogonal rates of change. A caller sweeping `mixing_beta` in a retry loop shouldn't be reaching into the struct that owns magnetization; a caller studying a magnetic system shouldn't need to round-trip through mixing state. Today they do.

## Proposal

Three changes, bundled because they land in the same struct:

### Part A — Split into `ElectronsPhysics` + absorb-rest-into-`ScfSettings`

New `src/settings.rs` layout:

```rust
pub struct ElectronsPhysics {
    pub spin_polarized: bool,                // see Part C — was `nspin: usize` (1|2)
    pub starting_magnetization: HashMap<String, f64>,
    pub tot_magnetization: Option<f64>,
    pub occupations: OccupationType,
}

pub struct ScfSettings {
    // existing
    pub max_iter: usize,
    pub conv_threshold: f64,
    pub energy_threshold: f64,
    pub n_bands: Option<usize>,
    pub eigensolver: EigensolverType,
    pub wfrx_subspace: bool,
    // migrated in from ElectronSettings
    pub mixing_beta: f64,
    pub mixing_ndim: usize,
    pub mixing_mode: MixingModeType,
    pub pulay_period: usize,
    pub adaptive_beta: bool,
    pub smearing: SmearingScheme,
    pub smearing_width: f64,
}
```

`electrons:` block in YAML shrinks to `spin_polarized` (bool) + magnetization. Everything else moves under `scf:`.

### Part B — Drop `scf.max_iter` default from 100 to 50

Current default (`src/settings.rs:171`) is 100. 50 is enough for LDA insulators (typically converge in 10-30 iters under Anderson/Pulay), and when it isn't the `PwdftError::ConvergenceFailure` path fires with a clean message. Condition: the error message must name `scf.max_iter` as the user-facing knob.

### Part C — Rename `nspin: usize` → `spin_polarized: bool`

`nspin` was never really a count — it's a flag that's 1 (spin-unpolarized, no m) or 2 (spin-polarized with up/down channels). The `usize` type allows nonsense values (0, 3, 17) that every consumer then has to reject, and it invites "what about `nspin=4` for non-collinear?" churn that we have no plan to implement. Non-collinear DFT (SOC, 2×2 spinors) is a fundamentally different code path, not a knob on this struct.

Migration at use sites:

| Before | After |
|---|---|
| `nspin == 1` | `!spin_polarized` |
| `nspin == 2` | `spin_polarized` |
| `nspin as f64` (in a few density-normalization spots) | keep internal plumbing that needs `1.0` or `2.0` in a tiny local `let nspin = if spin_polarized { 2 } else { 1 };` — the YAML/`Settings` surface is bool; the numeric constant stays a local derivation where it's actually used. Don't expose `nspin()` as a method on `ElectronsPhysics`; callers who need the integer should derive it at the call site so the type system pushes them toward the bool-shaped logic. |

YAML migration is a hard break. The field is literally named `spin_polarized` (bool) in YAML; any deck still writing `nspin: <anything>` fails to parse with a message naming `spin_polarized`. In-repo decks under `inputs/`, `examples/`, and `tests/` migrate in the same PR.

Downstream callers to update (non-exhaustive; the compiler will find the rest):

- `ScfParams::to_scf_params`
- `scf::driver` vs `scf::driver_spin` dispatch
- `pwdft_core::potential::xc` spin-polarized branch selection
- `initial_density` SAD that scales by `nspin`
- Any test fixture or integration test that constructs `ElectronSettings { nspin: 2, .. }` directly

Non-collinear is explicitly out of scope. If we ever add SOC, it becomes a new enum (`Collinear(bool) | NonCollinear`) or a separate `SpinMode` field — not a resurrection of integer `nspin`.

## Migration

Hard break. pwdft-rs is pre-release with no external users; no deprecation shims, no `#[serde(alias)]`, no one-release warnings.

- All repo-internal YAML fixtures (`inputs/`, `examples/`, `tests/`) migrate in the **same** PR to the new layout and to `spin_polarized`. A deck that still writes `nspin: 2` or puts `mixing_beta` under `electrons:` is a hard parse error pointing at the new location / field.
- `ScfParams::to_scf_params` plumbing: unchanged externally; internally picks from the new location.

## Risk

- **Low.** Settings.rs has 4 existing `full_settings_roundtrip`-style tests that catch schema regressions; extend them.
- **Zero** for callers who use defaults.
- **Hard break** for any stale YAML deck — caught immediately on first parse with a message naming the new field.

## Non-goals

- Not renaming `ScfSettings` or its existing fields. Existing names stay.
- Not touching `ScfParams` (the internal struct that `to_scf_params` produces). Only the YAML-facing split.

## Test plan

- Add roundtrip tests for both new layouts.
- Add a test that the deprecated layout still parses and produces an equivalent `Settings`, emitting a warning.
- `cargo test` green.

## Acceptance

- `ElectronsPhysics` exists, holds the four physics fields.
- `ElectronsPhysics::spin_polarized: bool` (renamed from `nspin: usize`); no `nspin` field or method survives on the new struct.
- `ScfSettings` absorbs mixing/smearing.
- Default `ScfSettings::max_iter == 50`.
- `PwdftError::ConvergenceFailure` display names `scf.max_iter`.
- Repo YAML fixtures migrated (both structural split and `nspin` → `spin_polarized`).
- Any YAML writing `nspin:` at all is a hard parse error naming `spin_polarized`. No `#[serde(alias)]`, no deprecation shim.

## Previous attempt (2026-04-20, aborted)

> **Heads-up for the recovery agent:** Part C (`nspin: usize` → `spin_polarized: bool` rename) was **added to this proposal after the 1ae1416 checkpoint**. The preserved branch keeps `nspin: usize` in the split struct. You inherit the Part A split + Part B max_iter drop from that branch, then layer Part C on top (rename the field, migrate call sites, add the legacy-nspin-integer serde shim, update YAML fixtures a second time). Do not re-derive Parts A/B.

A core-engineer agent (ID `af0257216752d1873`) made substantial progress
but aborted before `/pr-submit`. Work preserved on branch
**`origin/ESPL/electrons-settings-split`** at commit `1ae1416`.

- **Diff**: +461 / −98 LOC across 6 files — `pwdft/pwdft-core/src/settings.rs`,
  `pwdft/pwdft-core/src/error.rs`, and three `inputs/*.yaml` fixtures
  (`si_scf.yaml`, `si_scf_converged.yaml`, `si_scf_qe_match.yaml`).
- **Session logbook**: `.claude/logbooks/core-engineer/2026-04-20-espl-settings-split.md`
  on that branch documents the field-by-field migration map and the
  deprecation-alias approach the agent chose.
- **Not done**: `/quality-gate`, `/test --tier2` (settings.rs cascades
  into SCF driver plumbing), PR creation, verification that every
  non-touched deck under `examples/` and `tests/` still parses.

**Recovery plan for the next agent.** Check out the preserved branch,
rebase on `origin/main`, run `/pr-draft` immediately, then read the
logbook to understand scope. Confirm the deprecation-alias path works
on an old-layout YAML fixture, run `/quality-gate` + `/test --tier2`,
then `/pr-submit`. Don't re-derive the settings split — the 461 lines of
new code are load-bearing.

```bash
git -C <MAIN> fetch origin
git worktree add .claude/worktrees/<new-agent> -b ESPL-resume/electrons-settings-split origin/ESPL/electrons-settings-split
cd .claude/worktrees/<new-agent>
git rebase origin/main
/pr-draft "ESPL recovery pickup — settings split + YAML fixtures inherited from 1ae1416"
# ...read logbook, verify deprecation-alias path, quality gate + tier-2...
/pr-submit
```
