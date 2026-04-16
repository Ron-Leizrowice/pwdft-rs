---
id: SDED
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# SDED: Deduplicate Settings Enums

## Problem

The YAML migration introduced serde-friendly enum types in `settings.rs` that duplicate existing types in the SCF modules:

| Settings type | SCF type | Mapping code |
|---------------|----------|-------------|
| `SmearingType` (5 variants) | `SmearingScheme` (4 variants) | `to_scf_params` lines 395-404 |
| `MixingModeType` (2 variants) | `MixingMode` (2 variants) | `to_scf_params` lines 407-412 |
| `OccupationType` (2 variants) | *(not used internally)* | Only checked implicitly |

The `to_scf_params` method has a 20-line match block translating between equivalent types. Every time a new smearing scheme is added, both enums and the match must be updated.

## Implementation

### Option A: Add serde to SCF types directly (preferred)

Add `#[derive(Serialize, Deserialize)]` and `#[serde(rename_all = "snake_case")]` to `SmearingScheme` and `MixingMode` in the SCF modules. Then use them directly in `ElectronSettings`.

**Changes in `src/scf/smearing.rs`:**
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SmearingScheme {
    FermiDirac,
    Gaussian,
    MethfesselPaxton,
    Cold,
}
```

**Changes in `src/scf/mixing.rs`:**
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MixingMode {
    Plain,
    #[serde(deserialize_with = "deserialize_kerker")]
    Kerker { q_tf: Option<f64> },
}
```

The `Kerker` variant has internal state (`q_tf`), so it needs a custom deserializer that maps the YAML string `"kerker"` to `Kerker { q_tf: None }`. Alternatively, keep `MixingModeType` as the serde type and convert only `MixingMode` (2 lines).

**Changes in `src/settings.rs`:**
- Remove `SmearingType` enum entirely
- Use `crate::scf::smearing::SmearingScheme` in `ElectronSettings`
- Either remove `MixingModeType` or keep it as a thin serde adapter for `MixingMode`
- Remove the match blocks from `to_scf_params`

For `SmearingType::Fixed` (which has no counterpart in `SmearingScheme`): this belongs in `OccupationType` logic, not smearing. When `occupations: fixed`, the smearing scheme is irrelevant. Handle this in `to_scf_params` with a default fallback.

### Option B: Keep settings types, add `From` impls

If we want to keep a clean separation between serde types and internal types:

```rust
impl From<SmearingType> for SmearingScheme {
    fn from(st: SmearingType) -> Self {
        match st {
            SmearingType::FermiDirac => Self::FermiDirac,
            SmearingType::Gaussian => Self::Gaussian,
            SmearingType::MethfesselPaxton => Self::MethfesselPaxton,
            SmearingType::Cold => Self::Cold,
            SmearingType::Fixed => Self::FermiDirac, // unused when occupations=fixed
        }
    }
}
```

This centralizes the conversion but still requires maintaining both types.

### Recommendation

Option A for `SmearingScheme` (it's a simple flat enum). Option B for `MixingMode` (the `Kerker { q_tf }` variant makes serde tricky). This eliminates `SmearingType` entirely and keeps `MixingModeType` as a 2-line adapter.

## Verification

```bash
cargo test
```

Existing roundtrip tests (`smearing_types_roundtrip`, `mixing_mode_type_roundtrip`) validate serde behavior. After migration, these should test the unified types.

## Estimated Effort

Under an hour. The `SmearingScheme` change is mechanical. The `MixingMode` serde requires a small custom deserializer or keeping the adapter.
