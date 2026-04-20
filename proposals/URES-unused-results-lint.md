---
id: URES
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# URES: Replace `#[must_use]` annotations with `unused_results` lint

## Problem

The codebase carries 82 `#[must_use]` annotations across 22 files, all added by the
completed MUST proposal via `must_use_candidate` suggestions. These annotations are
point-in-time — new functions added without them silently opt out. The underlying goal
(warn when a computed return value is discarded) is better served by the rustc
`unused_results` lint, which fires unconditionally on every ignored non-`()` return
value without requiring per-function annotation.

Annotation counts by file:

| File | Count |
|------|-------|
| `src/symmetry/operations.rs` | 20 |
| `src/settings.rs` | 9 |
| `src/basis.rs` | 8 |
| `src/scf/mixing/anderson.rs` | 5 |
| `src/kpoints.rs` | 5 |
| `src/crystal.rs` | 6 |
| `src/scf/mixing/mod.rs` | 3 |
| `src/scf/smearing.rs` | 3 |
| `src/symmetry/mod.rs` | 3 |
| `src/fft.rs` | 4 |
| `src/scf/mixing/broyden.rs` | 2 |
| `src/potential/local.rs` | 2 |
| `src/symmetry/density/mod.rs` | 2 |
| `src/atoms.rs` | 2 |
| remaining 8 files | 1 each |
| **Total** | **82** |

## Implementation

1. **Add `unused_results` to `src/lib.rs`:**

```rust
#![warn(unused_results)]
```

Place it alongside the existing `cfg_attr(test, allow(...))` block at the top of the file.

2. **Compile and collect new warnings.** Run `cargo build 2>&1 | grep unused_results`
   to find all call sites that currently discard a return value without `let _ =`.
   For each warning, either:
   - Add `let _ =` if the discard is intentional (e.g. inserting into a map where
     the old value is irrelevant), or
   - Assign the result to a named binding if it should actually be used.

3. **Strip all 82 `#[must_use]` annotations** from the 22 files listed above.
   They are now redundant — `unused_results` covers every function unconditionally.

4. **Disable `must_use_candidate` if it was enabled.** Grep `src/` for
   `must_use_candidate`; remove any `#![warn(clippy::must_use_candidate)]` found.
   (Current audit shows it is not present, but confirm before shipping.)

## Verification

```bash
cargo build                                     # zero unused_results warnings
cargo clippy -q --all-targets                   # zero must_use_candidate warnings
cargo clippy -q --all-targets --features gpu
cargo test
```

After the change, discarding any non-`()` return value anywhere in `src/` produces a
warning without the author needing to remember to add `#[must_use]`.
