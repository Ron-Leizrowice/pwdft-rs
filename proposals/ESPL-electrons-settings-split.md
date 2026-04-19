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

Two changes, bundled because they land in the same struct:

### Part A — Split into `ElectronsPhysics` + absorb-rest-into-`ScfSettings`

New `src/settings.rs` layout:

```rust
pub struct ElectronsPhysics {
    pub nspin: usize,                        // 1 or 2
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

`electrons:` block in YAML shrinks to just `nspin` + magnetization. Everything else moves under `scf:`.

### Part B — Drop `scf.max_iter` default from 100 to 50

Current default (`src/settings.rs:171`) is 100. 50 is enough for LDA insulators (typically converge in 10-30 iters under Anderson/Pulay), and when it isn't the `PwdftError::ConvergenceFailure` path fires with a clean message. Condition: the error message must name `scf.max_iter` as the user-facing knob.

## Migration

- YAML schema change is strictly breaking (fields move between blocks). Handle it with a one-release serde deprecation warning: accept the old location with a `#[serde(alias)]` that logs a warning and maps to the new location.
- All repo-internal YAML fixtures (`examples/`, `tests/`) migrate in the same PR.
- `ScfParams::to_scf_params` plumbing: unchanged externally; internally picks from the new location.

## Risk

- **Low.** Settings.rs has 4 existing `full_settings_roundtrip`-style tests that catch schema regressions; extend them.
- **Zero** for callers who use defaults.
- **One-release cost** for callers with YAML configs — they see a deprecation warning until they rename.

## Non-goals

- Not renaming `ScfSettings` or its existing fields. Existing names stay.
- Not touching `ScfParams` (the internal struct that `to_scf_params` produces). Only the YAML-facing split.

## Test plan

- Add roundtrip tests for both new layouts.
- Add a test that the deprecated layout still parses and produces an equivalent `Settings`, emitting a warning.
- `cargo test` green.

## Acceptance

- `ElectronsPhysics` exists, holds the four physics fields.
- `ScfSettings` absorbs mixing/smearing.
- Default `ScfSettings::max_iter == 50`.
- `PwdftError::ConvergenceFailure` display names `scf.max_iter`.
- Repo YAML fixtures migrated.
