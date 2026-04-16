# Proposal 33: Consolidate Physical Constants

## Problem

The Coulomb constant `e²` is defined in 4 separate places with the same value `14.399_645_351_950_548`. The canonical location `consts::E2_COULOMB` exists but is unused — everyone imports from `hartree::E2` or defines their own copy.

Additionally, `consts::RHO_FLOOR` (1e-20) is defined but never used, while `xc.rs` hardcodes `1e-30` in 4 places. And `consts::PI` re-exports `std::f64::consts::PI` for a single consumer.

| Constant | Location | Status |
|----------|----------|--------|
| `E2_COULOMB` | `src/consts.rs:13` | Defined, never imported |
| `E2` | `src/potential/hartree.rs:11` | `pub`, imported by 6 files |
| `E2` | `src/pseudopotential/mod.rs:101` | Function-local `const` |
| `e2` | `src/pseudopotential/upf.rs:236` | `let` in test |
| `RHO_FLOOR` | `src/consts.rs:19` | Defined (1e-20), never used |
| `1e-30` | `src/potential/xc.rs:27,177,226,269` | Hardcoded density floor |
| `PI` | `src/consts.rs:1` | Re-export, only used by `crystal.rs` |

## Implementation

### Step 1: Consolidate E2

Delete `pub const E2` from `src/potential/hartree.rs:11`. Delete `const E2` from `src/pseudopotential/mod.rs:101`. Replace all uses:

| File | Current | Replacement |
|------|---------|-------------|
| `src/potential/hartree.rs:26,52` | `E2` | `crate::consts::E2_COULOMB` |
| `src/pseudopotential/mod.rs:111,121,130` | `E2` | `crate::consts::E2_COULOMB` |
| `src/ewald.rs:14` | `use potential::hartree::E2` | `use crate::consts::E2_COULOMB as E2` |
| `src/scf/energy.rs:35,156` | `hartree::E2` | `crate::consts::E2_COULOMB` |
| `src/scf/mod.rs:183` | `hartree::E2` | `crate::consts::E2_COULOMB` |
| `src/gpu/mod.rs:489` | `hartree::E2` | `crate::consts::E2_COULOMB` |
| `tests/gpu_consistency.rs:82,172` | `hartree::E2` | `pwdft_rs::consts::E2_COULOMB` |
| `benches/gpu_benchmarks.rs:34` | `hartree::E2` | `pwdft_rs::consts::E2_COULOMB` |

### Step 2: Use RHO_FLOOR in XC

Update `consts::RHO_FLOOR` from `1e-20` to `1e-30` to match actual usage. Then replace hardcoded `1e-30` in `xc.rs` lines 27, 177, 226, 269 with `crate::consts::RHO_FLOOR`.

### Step 3: Remove consts::PI

Delete `pub const PI` from `src/consts.rs:1`. Update `src/crystal.rs:4` to `use std::f64::consts::PI;`.

## Verification

```bash
cargo clippy -q --all-targets
cargo test
cargo test --test qe_validation  # physics unchanged
```

## Estimated Effort

Under an hour. Mechanical search-and-replace with import updates.
