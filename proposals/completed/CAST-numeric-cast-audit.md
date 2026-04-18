---
id: CAST
status: active
priority: low
complexity: medium
risk: medium
depends_on: []
blocks: []
---

# CAST: Numeric Cast Safety Audit

## Problem

Scientific code uses `as` casts liberally, and CLIP (legacy #27) deliberately skipped the three `cast_*` correctness lints because ~296 warnings were predominantly intentional. That decision was pragmatic at the time but left a real residual risk: with 142 `as` casts in `src/` plus more in tests, a *single* unintentional narrowing cast could silently corrupt a Miller index, grid coordinate, or shape count and produce wrong physics with no error.

Enabling the three lints flags every cast for a **one-time audit**, after which intentional cases carry `#[allow(...)]` with a rationale and future additions get reviewed automatically.

Current hit counts (from `-W clippy::pedantic` on 2026-04-17):

| Lint | Hits | Concern |
|------|------|---------|
| `cast_possible_truncation` | 52+12+1 = 65 | `usize → i32` (52), `f64 → i32` (12), `f64 → usize` (1). Narrowing. |
| `cast_possible_wrap` | 52+14+1 = 67 | `usize → i32` wrap on 32-bit (52), `u32 → i32` (14), `usize → i64` 32-bit target (1). |
| `cast_sign_loss` | 10+1+1 = 12 | `i32 → usize` (10), `i64 → usize` (1), `i32 → u32` (1). |
| `manual_midpoint` (overflow-safe avg) | 5 | Covered by QLNT, listed here for context. |

Total: ~144 sites requiring triage. Most will be `#[allow]`'d with a comment; the audit's value is in the small fraction that reveal actual bugs or surface a cleaner design (e.g. using `i32::try_from(usize)` at a boundary).

## Research

Typical patterns seen in this codebase:

```rust
// Intentional: Miller indices derived from geometric inequality, always fit in i32 for reasonable ecut
let n_max = (g_max / b1.norm()).ceil() as i32;   // basis.rs

// Intentional: grid dims asserted, fits trivially
let nx = dims.0 as usize;                        // fft.rs

// Potentially hidden bug: if atoms.len() > i32::MAX, wraps — unlikely but possible
let n = atoms.len() as i32;
```

The first two get `#[allow(clippy::cast_possible_truncation, reason = "bounded by ecut / asserted grid dim")]`. The third is the kind of case we want to find.

## Implementation

Do this as a **dedicated audit PR**, one reviewer pass. Not incremental across other work — the point is the concentrated review.

1. Enable the lints:

   ```toml
   [lints.clippy]
   cast_possible_truncation = "warn"
   cast_possible_wrap = "warn"
   cast_sign_loss = "warn"
   ```

2. Generate the full list of hits:

   ```bash
   cargo clippy -q --all-targets 2>&1 | \
     grep -E "casting|possible_wrap|sign_loss" > /tmp/cast-audit.txt
   ```

3. For each hit, choose one of:

   - **Rewrite** with a safe conversion (`try_from`, `.min(i32::MAX as usize) as i32`, or a new type-level bound) when the cast can actually fail in realistic inputs.
   - **`#[allow]` with `reason = "..."`** when the cast is bounded by a documented invariant:

     ```rust
     #[allow(clippy::cast_possible_truncation,
             reason = "n_max bounded by ecut; exceeding i32::MAX would need ecut > 10^9 Ry")]
     let n_max = (g_max / b1.norm()).ceil() as i32;
     ```

4. Prefer **module-level** allows over file-level `#![allow(...)]` so the lint keeps working elsewhere.

5. After the audit, the lints remain `warn` so any new cast gets flagged automatically — at which point the author either fixes it or adds a rationale.

## Verification

- `cargo clippy -q --all-targets` — zero warnings after audit.
- All existing tests pass, including QE cross-validation.
- Every `#[allow(clippy::cast_*)]` has a `reason = "..."` attribute explaining the invariant.
- Spot-check: grep for `#[allow(clippy::cast` and ensure ratios look reasonable (if >90% of hits are blanket-allowed, the audit wasn't thorough enough).

## Notes

This is the largest of the five clippy follow-up proposals by effort but the one most likely to expose a real bug. Worth doing after QLNT, FMAD, CLSS, and MUST land so the author isn't fighting multiple lint waves at once.
