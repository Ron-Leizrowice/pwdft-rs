---
id: SOPT
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# SOPT: Drop `Option<&SymmetryInfo>` from SCF Plumbing

## Problem

`ScfContext`, `run_scf`, and `run_scf_spin` carry symmetry as
`Option<&SymmetryInfo>`. The `Option` exists to support the
`symmetry.enabled = false` YAML setting, which makes
`Settings::to_symmetry_info()` return `None`. Two call sites
(`src/scf/mod.rs:337` and `:631`) then guard `symmetrize_density` with
`if let Some(symm) = ctx.symmetry`.

But every crystal has at least the identity group — the
`disabled` mode is just "user wants only the identity op". There is no
material that genuinely lacks symmetry. The `Option` therefore encodes
a setting (was symmetrization requested?) rather than a structural fact
about the physics, and it leaks into eight type signatures that would
otherwise be uniform.

This is a code smell flagged on 2026-04-17. It is not a correctness
bug — both branches produce the right answer. It is purely a
maintainability and clarity refactor.

## Proposed Change

1. Replace `symmetry: Option<&'a SymmetryInfo>` with
   `symmetry: &'a SymmetryInfo` in:
   - `ScfContext` (`src/scf/context.rs`)
   - `ScfContext::new` parameter list
   - `run_scf` (`src/scf/mod.rs`)
   - `run_scf_spin` (`src/scf/mod.rs`)
   - `compute_band_structure` if it carries symmetry (check
     `src/bandstructure.rs`)
2. When the user sets `symmetry.enabled = false`, construct an
   identity-only `SymmetryInfo` rather than returning `None` from
   `Settings::to_symmetry_info()`. Helper:
   ```rust
   pub fn identity(crystal: &Crystal) -> SymmetryInfo {
       SymmetryInfo {
           ops: vec![SymmOp::identity()],
           has_time_reversal: false,
           ..Default::default()
       }
   }
   ```
   Add this to `src/symmetry/mod.rs` alongside `from_crystal`.
3. In `symmetrize_density` (`src/symmetry/density.rs`), short-circuit
   when `symmetry.ops.len() == 1 && symmetry.ops[0].is_identity()`:
   ```rust
   if symmetry.is_trivial() { return; }
   ```
   `is_trivial()` is the right cheap predicate — add it to
   `SymmetryInfo`. The check costs one comparison per call.
4. Remove the `if let Some(symm) = ctx.symmetry` guards at
   `src/scf/mod.rs:337` and `:631`; call `symmetrize_density` directly
   with `ctx.symmetry`.
5. Update all callers in `src/main.rs`, `tests/`, and `benches/` to
   construct a real `SymmetryInfo` (either via `from_crystal` or the
   new `identity()` helper) instead of passing `None`. Audit:
   ```bash
   grep -rn "symmetry: None\|None,?\s*$" src/ tests/ benches/ | grep -i scf
   ```

## Why this is a good idea

- **Removes a class of bugs:** `Option`-handling forgetting to convert
  the `None` case is a common mistake. Eliminating the `Option`
  eliminates the bug class.
- **Makes signatures uniform:** all SCF entry points then take the
  same `&SymmetryInfo` rather than some taking `Option` and some not.
- **No physics change:** an identity-only `SymmetryInfo` produces
  bit-identical results to the current `None`-skipping path, because
  `symmetrize_density` short-circuits.
- **Follows the principle of "least options":** if the answer is
  always "yes, you have *some* symmetry," express that in the type.

## Why it might not be worth it

- Pure refactor, no perf/correctness win. If we're trying to keep PRs
  focused on user-visible behavior, defer.
- Touches ~8 files for a stylistic improvement. The current `Option`
  pattern works.

## Verification

Mechanical:
- `cargo test` — all 220+ tests pass with bit-identical results
  (modulo runtime variance).
- `cargo clippy -q --all-targets` — clean.
- Specifically, the existing `symmetry.enabled = false` test path (if
  one exists in `tests/`) must still produce the same final energy as
  before.

Numerical equality test: pick one of the QE-validation systems (Si or
Fe), run with old code and new code, assert total energies match to
1e-12 eV.

## Estimated Effort

~1 hour for a Core Engineer who already understands the symmetry
plumbing. Half the time is the mechanical signature changes; the
other half is verifying every caller and writing/updating the
identity-construction helper.

## Origin

Identified by the user during a code review on 2026-04-17, after the
SYKP audit. Logged here as a follow-up rather than coupled to SYKP
(SYKP was an audit-only outcome; this is its own refactor).
