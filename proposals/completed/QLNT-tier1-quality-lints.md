---
id: QLNT
status: completed
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# QLNT: Tier-1 Quality Lints (Follow-up to CLIP)

## Problem

The `CLIP` proposal (legacy #27) enabled 7 curated clippy lints. Since then the codebase has matured (ERRH, BROY, SIMP, DDUP landed) and a fresh pedantic+nursery audit shows a small set of additional lints with near-zero false-positive rates for this codebase. Enabling them now catches real issues without noise.

Baseline `cargo clippy -q --all-targets` is clean (0 warnings). Running with `-W clippy::pedantic -W clippy::nursery` surfaces the following low-noise lints that are worth promoting to `warn` in the workspace config:

| Lint | Hits | Why it matters here |
|------|------|---------------------|
| `manual_midpoint` | 5 | Overflow-safe midpoint — correctness in numerical code. |
| `redundant_clone` | 6 | Removes wasted allocations (some are in SCF hot path). |
| `match_same_arms` | — | Catches duplicated arms — often a copy-paste bug. |
| `semicolon_if_nothing_returned` | — | Style consistency; flags accidental expression returns. |
| `explicit_iter_loop` | — | Prefer `for x in &vec` over `for x in vec.iter()`. |
| `implicit_clone` | — | Flags `.to_vec()`/`.to_owned()` hiding clones on grid-sized data. |
| `inefficient_to_string` | — | `format!("{x}")` → `.to_string()` where appropriate. |
| `unnecessary_wraps` | — | Functions that always return `Ok(...)`/`Some(...)` — simplifies `Result` plumbing left over from ERRH. |
| `items_after_statements` | — | Makes module/function structure clearer. |
| `doc_link_with_quotes` | — | Fixes broken rustdoc links. |
| `unnecessary_box_returns` | — | Avoids `Box<T>` where `T` suffices. |

These lints were not addressed by CLIP because they are either in the `nursery` group (not part of `pedantic`) or were pre-empted by the "too-noisy-for-domain" filter that legitimately applied to `doc_markdown`/`cast_precision_loss` etc. but didn't apply to these.

## Implementation

Extend the existing `[lints.clippy]` block in `Cargo.toml`:

```toml
[lints.clippy]
# --- Existing (CLIP / legacy #27) ---
float_cmp = "warn"
needless_pass_by_value = "warn"
cloned_instead_of_copied = "warn"
uninlined_format_args = "warn"
redundant_closure = "warn"
unreadable_literal = "warn"
manual_let_else = "warn"

# --- QLNT additions ---
manual_midpoint = "warn"
redundant_clone = "warn"
match_same_arms = "warn"
semicolon_if_nothing_returned = "warn"
explicit_iter_loop = "warn"
implicit_clone = "warn"
inefficient_to_string = "warn"
unnecessary_wraps = "warn"
items_after_statements = "warn"
doc_link_with_quotes = "warn"
unnecessary_box_returns = "warn"
```

Then fix in order:

1. `cargo clippy --fix --allow-dirty --allow-staged --all-targets` — handles most style lints mechanically.
2. Manually audit `manual_midpoint` (5 hits) — confirm each replacement preserves behavior (no underflow in signed cases).
3. Manually audit `redundant_clone` (6 hits) — especially any in SCF iteration loops (`src/scf/mod.rs`, `src/scf/mixing.rs`).
4. `unnecessary_wraps` — review each hit; either drop the `Result`/`Option` or `#[allow]` with a comment explaining why the wrapper is load-bearing.

## Verification

```bash
cargo clippy -q --all-targets    # zero warnings
cargo test                       # full suite (~177 tests)
cargo test --features gpu        # GPU path (~186 tests)
```

No physics re-validation required — these are purely structural / style lints. Tests catch any behavioral regression from `manual_midpoint` or `unnecessary_wraps` edits.
