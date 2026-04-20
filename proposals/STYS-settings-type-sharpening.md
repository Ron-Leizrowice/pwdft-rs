---
id: STYS
status: active
priority: medium
complexity: medium
risk: low
depends_on: [ESPL]
blocks: []
---

# STYS: Settings type-sharpening audit — make invariants visible in the type

## Problem

pwdft-rs's `Settings` surface inherits QE-style "integer-as-flag" and
"string-as-enum" habits that Rust's type system can express directly.
Each one is a small pothole, but the cumulative effect is:

- Parse errors come out as validation-rejected values instead of
  "serde says this field must be one of {true, false}."
- Every consumer re-asserts the invariant (`if nspin != 1 && nspin != 2
  { … }`) instead of relying on the compiler.
- "What about nspin=4?" churn gets invited by fields that are
  structurally binary.

ESPL Part C kicked off the pattern (`nspin: usize` → `spin_polarized:
bool`). STYS is the proactive sweep so we don't discover each remaining
case ad-hoc.

## Audit findings

| Field | Current type | Actual domain | Proposed | File |
|-------|--------------|---------------|----------|------|
| `occupations` | `String` (parsed later) | `"fixed" \| "smearing"` | `OccupationType` enum | `settings.rs` |
| `smearing` | `String` (parsed later) | `"gaussian" \| "mp1" \| "fd" \| "cold"` | reuse `SmearingScheme` | `settings.rs` |
| `mixing_mode` | `String` (parsed later) | `"anderson" \| "broyden" \| "pulay" \| "periodic_pulay"` | reuse `MixingModeType` | `settings.rs` |
| `eigensolver` | `String` (parsed later) | `"dense" \| "iterative"` | reuse `EigensolverType` | `settings.rs` |
| `kgrid.shift` | `[usize; 3]` | `{0, 1}^3` | `[bool; 3]` or `Shift` enum | `settings.rs` |
| `SymmetrySettings.use_symmetry` | `bool` | ✓ already good | — | — |
| `SymmetrySettings.tolerance` | `f64` | ≥ 0 physics constraint | stays `f64`; add `NonNegF64` wrapper or validate | `settings.rs` |
| `BasisSettings.ecutrho_ratio` | `f64` | ≥ 1 physics constraint | stays `f64`; document invariant + validate | `settings.rs` |

(Counts are from a grep pass on `settings.rs` @ commit d408eee; refresh
when implementing.)

For the enum-shaped `String` fields: most already have a companion enum
elsewhere in the crate (`OccupationType`, `MixingModeType`, etc.). The
current architecture parses the string at the YAML boundary, then
converts. STYS flattens that — the YAML field is typed directly via
`#[serde(rename_all = "snake_case")]`. Invalid strings become a serde
parse error naming the field and the accepted variants.

## Research

### Why `String` at all?

Reading git history on `settings.rs` — the string-typed fields predate
the companion enums (which were added in MODR/CFGN/TSPL for internal
plumbing). The `Settings` surface was never retrofitted. STYS is the
retrofit.

### What about forward-compat (new variants)?

Adding a new `SmearingScheme` variant right now means: (a) add the
variant to the enum, (b) handle it in the dispatcher. With STYS, the
same change also means: (c) the new variant becomes accepted in YAML
automatically via `#[serde(rename_all = "snake_case")]`. That's one
fewer place to update, not one more.

### Why only settings-surface?

Internal plumbing (`ScfParams.nspin: usize`) stays integer — see
ESPL's non-goals. Settings-surface is where *user* invariants live;
those are the ones worth pinning in the type. Internal plumbing that
happens to use `usize` as a loop bound or array index stays `usize`.

## Implementation

Five independent phases; each lands as its own commit within the
same PR. Compiler-driven — flip the field type, follow the errors.

### Phase A — `ElectronsPhysics.occupations: String` → `OccupationType`

Reuse the existing `OccupationType` enum if one exists; otherwise
add it:

```rust
#[derive(Deserialize, ..)]
#[serde(rename_all = "snake_case")]
pub enum OccupationType { Fixed, Smearing }
