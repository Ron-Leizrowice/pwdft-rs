---
id: ELMN
status: active
priority: low
complexity: trivial
risk: low
depends_on: []
blocks: []
---

# ELMN: Minimize `src/atoms.rs` — delete `from_symbol` / `from_z`, keep only `pub use`

## Problem

`src/atoms.rs` (72 LOC) is a thin wrapper around `mendeleev::Element`:

```rust
pub use mendeleev::Element;

pub fn from_symbol(s: &str) -> Option<Element> {
    Element::iter().find(|e| e.symbol() == s)
}

pub fn from_z(z: u32) -> Option<Element> {
    Element::list().iter().copied().find(|e| e.atomic_number() == z)
}
```

Both helpers are 6-line wrappers around mendeleev's native iter APIs. They exist because someone preferred `atoms::from_symbol("Si")` to `Element::iter().find(|e| e.symbol() == "Si")`. Not wrong, but the indirection means:

- One more file to grep through when looking for the element data model.
- Test sub-mod (27 LOC) exercises the wrappers, not anything novel.
- Every touch has to decide "does this belong in atoms.rs or use mendeleev directly?"

## Proposal

Shrink `src/atoms.rs` to one line:

```rust
pub use mendeleev::Element;
```

Delete `from_symbol` and `from_z`. Callers migrate to `Element::iter().find(|e| e.symbol() == s)` or `Element::from_atomic_number(z)` (mendeleev's native API — confirm during implementation that this exists with exactly that signature; if not, use `Element::list()`).

## Migration

Grep call sites:

- `src/atoms::from_symbol(...)` — roughly 2-4 sites (PP loader, settings parser).
- `src/atoms::from_z(...)` — roughly 1-2 sites (UPF header parse).

Expected churn: ~10 lines across 3-5 files. No behavioral change.

## Risk

- **Near zero.** Mechanical rewrite. Existing tests that exercised the helpers become tests of mendeleev's native API — redundant, delete.

## Non-goals

- Not removing `pub use mendeleev::Element;`. The `atoms::Element` alias lets the rest of the crate speak one vocabulary; it's cheap and valuable. Keep.
- Not bumping `mendeleev` beyond the current range. Version pin stays.

## Acceptance

- `src/atoms.rs` is ≤ 5 lines (re-export + module doc if any).
- No caller imports `atoms::from_symbol` or `atoms::from_z`.
- `cargo test` green — including any tests that indirectly relied on the wrappers.
