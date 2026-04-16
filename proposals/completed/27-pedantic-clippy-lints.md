# Proposal 27: Pedantic Clippy Lints

## Problem

The project has no clippy lint configuration beyond the defaults. Enabling `clippy::pedantic` as a blanket group produces ~887 warnings, many of which are noise for a scientific computing codebase (cast lints, similar variable names for physics notation, doc formatting for math terms). A curated subset gives the correctness and quality benefits without the noise.

## Research

Full `clippy::pedantic` warning breakdown (887 total):

| Count | Lint category | Domain relevance |
|-------|---------------|-----------------|
| 363 | `doc_markdown` — backtick formatting | Noise — flags physics terms and math symbols |
| 296 | `cast_*` — numeric casts | Noise — scientific code casts `f64`/`f32`/`usize`/`i32` constantly |
| 85 | `must_use_candidate` — missing `#[must_use]` | Low — useful for libraries, noisy for applications |
| 35 | `unreadable_literal` — unseparated constants | Useful — physical constants like `27.211386245988` |
| 21 | `uninlined_format_args` — `format!("{}", x)` | Useful — modernization, auto-fixable |
| 15 | `float_cmp` — `==` on floats | Critical — almost always a bug in numerical code |
| 11 | `similar_names` — `q_i` vs `q_j` | Noise — standard physics notation |
| 11 | `missing_panics_doc` — undocumented panics | Medium — documents crash conditions |
| 8 | `unnecessary_wraps` — raw string hashes | Useful — auto-fixable |
| 7 | `docs_for_function_returning` — missing docs | Medium |
| 6 | `redundant_closure` — `\|x\| foo(x)` | Useful — auto-fixable |
| 6 | `missing_errors_doc` — undocumented errors | Medium |
| 5 | `wildcard_enum_match_arm` | Useful — catches missing future variants |
| 5 | `manual_*` implementations | Useful — use std library equivalents |
| 2 | `manual_let_else` | Useful — modern Rust pattern |
| 2 | `too_many_lines` | Noise — SCF loop is inherently complex |
| 1 | `single_char_binding_names` | Noise — physics variables (x, k, q, r) |

## Recommended Lints

### Tier 1 — Bug prevention

| Lint | Warnings | Rationale |
|------|----------|-----------|
| `float_cmp` | 15 | Catches `==` / `!=` on floats. Critical for numerical code — these are almost always bugs or should use `approx`. |

### Tier 2 — Performance

| Lint | Warnings | Rationale |
|------|----------|-----------|
| `needless_pass_by_value` | ~55 | Flags functions taking owned `T` where `&T` suffices. Avoids unnecessary clones, especially in hot paths. |
| `cloned_instead_of_copied` | ~5 | `.cloned()` on `Copy` types should be `.copied()` — clearer intent, signals the type is cheap. |

### Tier 3 — Code modernization (auto-fixable)

| Lint | Warnings | Rationale |
|------|----------|-----------|
| `uninlined_format_args` | 21 | `format!("{}", x)` → `format!("{x}")`. Auto-fixable via `--fix`. |
| `redundant_closure` | 6 | `\|x\| foo(x)` → `foo`. Auto-fixable. |
| `unreadable_literal` | 35 | `27.211386245988` → `27.211_386_245_988`. Improves readability of physical constants. Auto-fixable. |
| `manual_let_else` | 2 | Modernizes `match x { Some(v) => v, None => return }` → `let Some(v) = x else { return }`. |

### Skipped — too noisy for this domain

| Lint | Warnings | Why skip |
|------|----------|----------|
| `cast_possible_truncation`, `cast_sign_loss`, `cast_precision_loss` | ~296 | Scientific code casts between `f64`/`f32`/`usize`/`i32` constantly. Nearly all are intentional. |
| `similar_names` | 11 | Flags standard physics notation (`q_i`/`q_j`, `z_sum`/`z2_sum`, `jl`/`jlm1`/`jlp1`). |
| `doc_markdown` | 363 | Flags physics terms and math symbols as needing backticks. Massive noise. |
| `must_use_candidate` | 85 | Suggests `#[must_use]` on many functions. Useful in library crates, noisy for an application. |
| `too_many_lines` | 2 | The flagged functions (`run_scf`, etc.) are inherently complex; splitting them would obscure the algorithm flow. |
| `single_char_binding_names` | 1 | Physics code uses single-letter variables by convention (x, k, q, r). |

## Implementation

Add to `Cargo.toml`:

```toml
[lints.clippy]
# Tier 1 — bug prevention
float_cmp = "warn"

# Tier 2 — performance
needless_pass_by_value = "warn"
cloned_instead_of_copied = "warn"

# Tier 3 — code modernization
uninlined_format_args = "warn"
redundant_closure = "warn"
unreadable_literal = "warn"
manual_let_else = "warn"
```

Then fix in order:

1. `cargo clippy --fix --allow-dirty --allow-staged` — auto-fixes `uninlined_format_args`, `redundant_closure`, `unreadable_literal`
2. Manually review and fix `float_cmp` hits (15) — each one is a potential bug
3. Audit `needless_pass_by_value` hits (~55) — some may need `#[allow]` where ownership is intentional
4. Fix remaining `cloned_instead_of_copied` and `manual_let_else` (trivial)
5. Run `cargo test` to verify nothing broke

## Verification

```bash
cargo clippy -q --all-targets         # zero warnings from enabled lints
cargo test                             # all tests pass
cargo test --features gpu              # GPU tests pass (if applicable)
```

## Estimated Effort

Tier 3 auto-fixes: minutes. Tier 1–2 manual fixes: a couple sessions, mostly the `needless_pass_by_value` audit.
