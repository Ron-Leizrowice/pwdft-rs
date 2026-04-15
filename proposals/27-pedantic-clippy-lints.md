# Proposal 27: Enable Pedantic Clippy Lints

**Status:** Not started. Experimental run showed 822 warnings, all auto-fixable.

## Problem

The codebase currently runs default clippy only. Enabling `clippy::pedantic` catches additional code quality issues: missing `#[must_use]`, `From` vs `as` casts, unreadable numeric literals, missing doc sections, redundant closures, and more.

An experimental run with `pedantic = "warn"` in `[lints.clippy]` produced 822 warnings. Running `cargo clippy --fix` resolved all of them automatically.

## Plan

1. Add `pedantic = "warn"` to `[lints.clippy]` in Cargo.toml
2. Run `cargo clippy --fix --allow-dirty --allow-staged --all-targets`
3. Review auto-fixes for correctness (especially cast changes)
4. Selectively `#[allow]` categories that are noise for scientific computing:
   - `cast_possible_truncation` / `cast_sign_loss` on grid index math
   - `cast_precision_loss` on `usize as f64` (grids are never >2^52)
5. Verify zero warnings, all tests pass

## Warning breakdown (from experimental run)

| Count | Warning | Auto-fixable |
|-------|---------|--------------|
| 290 | `missing_docs_in_private_items` → backtick items in docs | Yes |
| 88 | `cast_lossless` → `i32 as f64` should use `From` | Yes |
| 56 | `must_use_candidate` | Yes |
| 52 | `cast_possible_truncation` | Suppress |
| 48 | `cast_precision_loss` | Suppress |
| 35 | `unreadable_literal` → add underscores to constants | Yes |
| 31 | `return_self_not_must_use` | Yes |
| 21 | `uninlined_format_args` | Yes |
| 14 | `float_cmp` | Review case-by-case |
| ~187 | Other (redundant closure, similar names, etc.) | Mostly yes |

## Acceptance Criteria

1. `[lints.clippy] pedantic = "warn"` in Cargo.toml
2. Zero clippy warnings with `--all-targets`
3. All tests pass unchanged
