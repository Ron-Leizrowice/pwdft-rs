# Proposal 21: Replace atoms.rs with elements_rs Crate

## Problem

`src/atoms.rs` is a hand-rolled `Element` enum covering Z=1-92 with three capabilities: `atomic_number()`, `symbol()`, and `from_symbol()`/`from_z()` conversions. It uses `unsafe { std::mem::transmute }` in `from_z()` (line 118) and maintains a parallel `SYMBOLS` array that must stay synchronized with the enum variants.

The `elements_rs` crate (v0.2.1) provides all 118 elements with rich metadata that could replace our enum and unlock future features:

| Capability | Current `atoms.rs` | `elements_rs` |
|-----------|-------------------|---------------|
| Elements covered | Z=1-92 | Z=1-118 |
| Symbol/Z conversion | Manual enum + array | Built-in `FromStr`, `From<u8>`, `TryFrom` |
| Atomic mass | Not available | `standard_atomic_weight()` |
| Atomic radii | Not available | Slater, Rahm, Cordero, Bondi, Mantina |
| Element classification | Not available | Metal/nonmetal, category enum |
| Valence electrons | Not available | `valence_electrons()` trait |
| Electron config | Not available | Full `orbitals()` |
| `unsafe` code | `transmute` in `from_z` | None |
| Maintenance | Manual | Community-maintained |

Properties we don't currently use but would benefit from in future proposals:

- **Atomic mass:** Needed for molecular dynamics, phonon calculations, and mass-weighted coordinates.
- **Atomic radii:** Useful for initial guess (overlap radius), bond detection, and visualization. Currently the SAD initial guess (`initial_density.rs`) uses pseudopotential radial grids, but covalent/vdW radii from a reference source would help with grid partitioning and neighbor detection.
- **Element classification:** Could auto-select smearing parameters (metals need smaller sigma) or mixing strategy.

## References

- Crate: [elements_rs on crates.io](https://crates.io/crates/elements_rs) (v0.2.1, published 2026-04-05)
- Repository: <https://github.com/earth-metabolome-initiative/elements-rs>
- License: GPL-3.0 per Cargo.toml (verify compatibility with project license)

## Implementation

### Step 1: Add dependency

```toml
[dependencies]
elements_rs = ">=0.2"
```

**License:** The crate is GPL-3.0. All existing pwdft-rs dependencies are permissive (MIT, Apache-2.0, BSD) and GPL-3 compatible. Adding this dependency means pwdft-rs must be distributed under GPL-3 terms, which is fine for this project.

### Step 2: Replace internal Element type

Remove `src/atoms.rs` entirely. Create a thin re-export wrapper:

```rust
// src/atoms.rs (replacement — re-export + compatibility)
pub use elements_rs::Element;

/// Extension trait for DFT-specific element properties not in elements_rs.
pub trait ElementExt {
    fn atomic_number_u32(self) -> u32;
}

impl ElementExt for Element {
    fn atomic_number_u32(self) -> u32 {
        u8::from(self) as u32
    }
}
```

### Step 3: Update call sites

There are 6 call sites outside `atoms.rs`:

**`src/input.rs:124`** — element lookup from TOML symbol:

```rust
// Before:
let elem = crate::atoms::Element::from_symbol(&ai.symbol)
    .unwrap_or_else(|| panic!("unknown element: {}", ai.symbol));
Atom::new(elem.atomic_number(), ai.position)

// After:
let elem: Element = ai.symbol.parse()
    .map_err(|_| PwdftError::InvalidInput(format!("unknown element: {}", ai.symbol)))?;
Atom::new(u8::from(elem) as u32, ai.position)
```

**`src/settings.rs:343`** — same pattern as input.rs.

**`src/pseudopotential/mod.rs:94-95`** — match element by symbol and Z:

```rust
// Before:
crate::atoms::Element::from_symbol(&pp.element)
    .is_some_and(|e| e.atomic_number() == z)

// After:
pp.element.parse::<Element>().ok()
    .is_some_and(|e| u8::from(e) as u32 == z)
```

**`src/pseudopotential/psp8.rs:148-149`** — derive symbol from Z:

```rust
// Before:
let element = crate::atoms::Element::from_z(z_int)
    .map(|e| e.symbol().to_string())
    .unwrap_or_else(|| format!("Z{z_int}"));

// After:
let element = Element::try_from(z_int as u8)
    .map(|e| e.symbol().to_string())
    .unwrap_or_else(|_| format!("Z{z_int}"));
```

### Step 4: Leverage new properties (optional, future)

Once integrated, new capabilities become available without additional code:

```rust
use elements_rs::{Element, AtomicRadius, Electronegativity, ElementClassification};

// Auto-detect metallic systems for smearing defaults
fn is_metallic_system(atoms: &[Atom]) -> bool {
    atoms.iter().any(|a| {
        Element::try_from(a.z as u8)
            .map(|e| e.is_metal())
            .unwrap_or(false)
    })
}

// Atomic radius for neighbor detection or SAD partitioning
fn atomic_radius_ang(z: u32) -> Option<f64> {
    Element::try_from(z as u8).ok()
        .and_then(|e| e.slater_atomic_radius())
}

// Atomic mass for future MD
fn atomic_mass_amu(z: u32) -> Option<f64> {
    Element::try_from(z as u8).ok()
        .map(|e| e.standard_atomic_weight())
}
```

### Step 5: Remove unsafe code

The current `from_z` uses `unsafe { std::mem::transmute }`. With `elements_rs`, this becomes `Element::try_from(z as u8)` which is safe. The entire `SYMBOLS` array and manual enum can be deleted.

## Acceptance Criteria

1. **All existing tests pass:** Element lookups by symbol and Z produce identical results for Z=1-92.
2. **`unsafe` removed:** No `transmute` in the codebase for element conversion.
3. **`atoms.rs` simplified:** The file is either deleted or reduced to a re-export + extension trait.
4. **Extended range:** Elements Z=93-118 are now recognized (no panic on Np, Pu, etc.).
5. **No behavior change:** SCF results identical — this is a pure refactor of element lookup code.
