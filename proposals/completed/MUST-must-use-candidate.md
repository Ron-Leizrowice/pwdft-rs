---
id: MUST
status: completed
priority: low
complexity: medium
risk: low
depends_on: []
blocks: []
---

# MUST: `must_use_candidate` Annotations

## Problem

Clippy's `must_use_candidate` (pedantic) flags pure functions/methods whose return value is easy to accidentally drop. The current codebase has ~80 hits:

- 53 methods (e.g. builders, accessors returning new values)
- 27 free functions

In a DFT solver, dropping a result is a real failure mode: forgetting to assign the output of `assemble_v_eff(...)` or `symmetrize_density(...)` can silently run with stale data. `#[must_use]` makes the compiler catch this.

The CLIP proposal (legacy #27) originally skipped this lint with the rationale "useful in library crates, noisy for applications." That call made sense pre-ERRH when many functions returned `Result` and the `#[must_use]` on `Result` already covered most cases. Post-ERRH, many functions that now return bare `T` (e.g. energy scalars, density arrays) would benefit from explicit `#[must_use]`.

## Research

Hit distribution by kind (from `-W clippy::pedantic`):

| Kind | Count | Example |
|------|-------|---------|
| `this method could have a #[must_use] attribute` | 53 | `ScfParams::with_mixing_beta` style builders |
| `this function could have a #[must_use] attribute` | 27 | Pure numerical helpers in `numerics.rs`, `ewald.rs` |
| `missing #[must_use] attribute on a method returning Self` | 6 | Likely builder methods — these are the strongest case |

The builder-returns-`Self` cases (6) are the highest value: forgetting to chain the builder result is a clear bug. The rest are judgment calls.

## Implementation

1. Enable the lint:

   ```toml
   [lints.clippy]
   must_use_candidate = "warn"
   ```

2. Triage in three passes:

   **Pass A — builders returning `Self` (6 hits):** always add `#[must_use]`. Unambiguously correct.

   **Pass B — pure numerical computations (most of the 27 functions):** add `#[must_use]`. These functions have no side effects; dropping the result is always a bug.

   **Pass C — methods on types (53 hits):** case-by-case. If the method has interior mutability or a side effect, leave it; otherwise add `#[must_use]`.

3. If any hit genuinely should not be `#[must_use]` (e.g. the method is called purely for its effect and the return value is a convenience), apply `#[allow(clippy::must_use_candidate)]` with a one-line comment explaining why.

4. Consider adding a short `#[must_use = "..."]` message on the highest-value items (energy results, convergence checks) for extra friction on accidental drops.

## Verification

```bash
cargo clippy -q --all-targets    # zero warnings
cargo test                       # full suite
```

No physics impact — attributes are advisory only and don't change runtime behavior.
