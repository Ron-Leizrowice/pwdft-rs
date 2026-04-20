# Core Engineer — 2026-04-20 — ESPL: split ElectronSettings into ElectronsPhysics + ScfSettings

Implemented ESPL end-to-end: `ElectronSettings` is gone; system-physics
knobs live in `ElectronsPhysics`, convergence/mixing/smearing knobs live
in `ScfSettings`. `scf.max_iter` default dropped 100 → 50.
`ConvergenceFailure` Display now names `scf.max_iter`.

## Key numbers

**Field migration map (old → new YAML path):**

| Old `electrons.<field>`   | New `scf.<field>`      |
|---------------------------|------------------------|
| `mixing_beta`             | `scf.mixing_beta`      |
| `mixing_ndim`             | `scf.mixing_ndim`      |
| `mixing_mode`             | `scf.mixing_mode`      |
| `pulay_period`            | `scf.pulay_period`     |
| `adaptive_beta`           | `scf.adaptive_beta`    |
| `smearing`                | `scf.smearing`         |
| `smearing_width`          | `scf.smearing_width`   |

Remaining in `electrons:` (ElectronsPhysics): `nspin`,
`starting_magnetization`, `tot_magnetization`, `occupations`.

**Fixtures migrated:** `inputs/si_scf.yaml`, `inputs/si_scf_qe_match.yaml`,
`inputs/si_scf_converged.yaml`. `si_free_electron.yaml` didn't have an
`electrons:` block. No test files carried an embedded YAML `electrons:`
block beyond the unit tests in `src/settings.rs`.

**Default changes:**

- `ScfSettings::max_iter`: 100 → 50.
- `ConvergenceFailure` Display: now reads
  `"SCF did not converge after {n} iterations (delta = {δ}); raise
  scf.max_iter or tighten mixing/smearing in the input YAML"`.

**Deprecation-warning text (one `log::warn!` per migrated field):**

```text
YAML field `electrons.<field>` is deprecated; move to `scf.<field>`
(pre-ESPL layout accepted for one release)
```

**Tier-2 wall-time:** See PR body.

## Tangential

- The pre-ESPL YAML layout is accepted via a `SettingsWire`/`ElectronsWire`
  shim (`#[serde(from = "SettingsWire")]` on `Settings`) because serde's
  `#[serde(alias)]` can't cross struct boundaries between `electrons:` and
  `scf:`. Planned removal: one release after ESPL lands.
- `ScfParams` (the internal SCF-layer struct) is unchanged. The YAML split
  is flattened at `Settings::to_scf_params` boundary, so downstream code
  is unaware of the YAML reorganization.
