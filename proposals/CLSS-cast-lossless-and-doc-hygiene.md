---
id: CLSS
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# CLSS: `cast_lossless` + Doc Hygiene (`missing_errors_doc`, `missing_panics_doc`)

## Problem

Three pedantic lints catch small-but-real issues that the CLIP baseline skipped.

### 1. `cast_lossless` — 85 hits

Expressions like `(2 * l + 1) as f64` or `i as f64` (from `i32`/`u32`) silently cast even though the conversion is infallible. `f64::from(x)` is:

- Explicit about losslessness (will fail to compile if the source type becomes lossy, e.g. `i64`).
- Easier to read and grep for.

Representative hits:

```rust
// tests/kb_projector_validation.rs:490
let angular = (2 * l + 1) as f64 / (4.0 * PI);
// → let angular = f64::from(2 * l + 1) / (4.0 * PI);

// tests/kb_projector_validation.rs:595
let q_vals: Vec<f64> = (0..50).map(|i| i as f64 * 0.5).collect();
// → let q_vals: Vec<f64> = (0..50).map(|i| f64::from(i) * 0.5).collect();
```

This is **distinct from** `cast_precision_loss` (intentionally not enabled — see CLIP): lossless-from-i32/u32 is always safe; `usize → f64` is where precision loss can occur, and that's flagged separately.

### 2. `missing_errors_doc` — 7 hits

Public functions that return `Result` without a `# Errors` section. With only 7 items flagged, this is a quick, high-value improvement post-ERRH: callers can see what error conditions exist without grepping the source.

### 3. `missing_panics_doc` — 10 hits

Functions with `panic!` / `unwrap` / `expect` that aren't documented as panicking. After ERRH, most remaining panics are invariant checks, and documenting them makes the invariants visible at the API boundary.

## Implementation

1. Enable the lints:

   ```toml
   [lints.clippy]
   cast_lossless = "warn"
   missing_errors_doc = "warn"
   missing_panics_doc = "warn"
   ```

2. `cargo clippy --fix --allow-dirty --allow-staged --all-targets` handles all 85 `cast_lossless` hits mechanically.

3. Manually add `# Errors` sections to the 7 flagged functions. Format:

   ```rust
   /// ...
   ///
   /// # Errors
   /// Returns `Err(PseudopotentialError::InvalidFormat)` if the UPF header is malformed.
   pub fn parse(...) -> Result<...> { ... }
   ```

4. Manually add `# Panics` sections to the 10 flagged functions. If a panic is truly "this should never happen" (invariant), prefer rewriting to return `Result`; otherwise document the trigger.

## Verification

```bash
cargo clippy -q --all-targets    # zero warnings
cargo test                       # all tests pass
cargo doc --no-deps              # rustdoc builds clean
```

No physics validation needed — these changes are mechanical or pure documentation.
