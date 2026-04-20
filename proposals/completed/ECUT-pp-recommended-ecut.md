---
id: ECUT
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# ECUT: Per-PP recommended wavefunction cutoff; drop hardcoded 204.09 eV default

## Problem

`src/settings.rs:101` hardcodes `ecutwfc: 204.09` eV (= 15 Ry) as the global default. This is Si's PseudoDojo NC/LDA value and happens to be reasonable for most light elements, but it is:

- Too low for first-row transition metals and heavier (Cu typically wants 40-50 Ry, Fe ~60 Ry for the `.stringent` set).
- Too high for simple cases (H, Li, Na can converge at 8-10 Ry).
- Not derivable from any per-PP metadata the engine currently reads — the number was pinned from "what Si wants".

Users see no diagnostic when the default fires. They get an ecut that is quietly wrong for their system until convergence tests surface the gap.

## Non-problem: the UPF `rho_cutoff` attribute

The PseudoDojo UPF header carries `rho_cutoff=...` (e.g., Si: 15.09 Ry). **This is NOT the recommended ecutwfc.** It is the cutoff used for the atomic PP_RHOATOM radial mesh, unrelated to the plane-wave cutoff. Using it as a default would be a coincidence for Si and wrong for every other element. This must not be the source of truth.

## Proposal

Ship a per-element recommended-ecut table derived from PseudoDojo's publicly documented `.standard` / `.stringent` defaults. Apply max-across-species when the crystal has more than one element.

### Part A — Recommended-ecut table

New `src/pseudopotential/recommended_ecut.rs`:

```rust
/// PseudoDojo NC/LDA `.standard` recommended ecutwfc in eV, per element.
/// Sourced from <http://www.pseudo-dojo.org/> as of 2026-04-19.
/// `.stringent` variant available via `STRINGENT_ECUT[z]`.
pub fn recommended_ecut_ev(z: u32, variant: EcutVariant) -> Option<f64> { ... }

pub enum EcutVariant { Standard, Stringent }
```

Table is ~92 entries, static `const` array. Maintenance cost: re-sync on PseudoDojo revisions (rare).

### Part B — Default uses the table + warns

`BasisSettings::ecutwfc` becomes `Option<f64>` in YAML. When absent:

1. Look up each species' `Standard` ecut from the table.
2. Take the max across all species in the crystal. Optional 1.0× multiplier (configurable later; default 1.0).
3. `log::warn!("ecutwfc not set, using recommended {} eV from PseudoDojo .standard (max over species); set`basis.ecutwfc`to override", ecut)`.

When `basis.ecutwfc` is set, use it unchanged (current behavior). No warning.

### Part C — Remove the 204.09 literal

Delete the hardcoded default. A system with an unknown element in the table falls back to the highest known value with a louder warning — fail loud, not quiet.

## Sanity check against qe_validation

The 8 reference systems in `qe_validation/` use Si/C/Al/Fe/Ga/As/Mg/O/Na/Cl. Compute the table-derived default for each and confirm it matches or exceeds the hand-picked `ecutwfc` in the existing `.in` files. Any mismatch is a table-entry correction, not a physics change.

## Risk

- **Low.** The change is additive: tests/fixtures that set ecutwfc explicitly are unaffected.
- Fixtures that rely on the default (= 204.09) may shift to higher cutoffs for heavy-element systems; re-pin as needed. This is the point.

## Non-goals

- Not shipping a `.stringent` table in the first pass — `.standard` only.
- Not auto-scaling ecutrho from ecutwfc beyond the existing `ecutrho_ratio` knob. That's a separate conversation.
- Not removing `basis.ecutwfc` from YAML. User override must remain.

## Acceptance

- `recommended_ecut_ev(Z, Standard)` returns the PseudoDojo value for every Z the repo's `pseudopotentials/nc/lda/` tree ships (currently up to Ra-88).
- Defaulting path emits `log::warn!` with the derived value and named override knob.
- `src/settings.rs` has no `204.09` literal.
- 8 QE-validation SCFs still produce the same energies when `ecutwfc` is explicit (unchanged path).
