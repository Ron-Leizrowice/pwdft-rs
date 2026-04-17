---
id: DBGC
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# DBGC: ScfResult Debug derive + GPU-test reference-value const

## Origin

Two nits surfaced in code review of TAUD bundle PR #26 (2026-04-17):

1. **`ScfResult` lacks `#[derive(Debug)]`** (`src/scf/mod.rs:131`). Forced TAUD PR E to use a 3-arm `match` instead of `assert!(matches!(result, Err(...)))` — the natural pattern would have been `result.is_err_and(|e| matches!(...))` with a debug-format diagnostic on failure. Adding `Debug` simplifies all future inverted-test patterns.

2. **`-198.8926` (Si total energy reference) duplicated** in two GPU tests in `tests/gpu_consistency.rs`. Hoist to a module-level `const SI_REFERENCE_TOTAL_EV: f64 = -198.8926;` (and similar for `SI_REFERENCE_FERMI_EV: f64 = 6.969;`). One source of truth for these magic numbers.

## Implementation

1. Add `#[derive(Debug)]` to `ScfResult` in `src/scf/mod.rs`. Verify all transitively-required types also derive `Debug` (likely already do — `Vec<f64>`, `f64`, `usize`, etc.).
2. Hoist the GPU Si reference constants in `tests/gpu_consistency.rs`. Adjust the two tightened-tolerance assertions from PR C to reference the consts.
3. Run the suites as a sanity check — should be no-op behaviorally.

## Verification

- `cargo test` — 222+ pass, no behavior change.
- `cargo test --features gpu` — 231+ pass.
- `cargo clippy -q --all-targets --features gpu` — clean.

## Notes

Tiny PR (~10 lines). Run as a one-shot Code Reviewer task; no need to bundle.
