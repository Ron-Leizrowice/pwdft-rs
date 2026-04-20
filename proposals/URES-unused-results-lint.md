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

## Implementation policy

**`let _ = ...` outside `#[cfg(test)]` is a code smell, not the fix.** Every new
`unused_results` warning in production code is telling you one of three things,
none of them routine. Triage each warning with this decision tree:

### 1. `Result<_, _>` discard → real bug

You were silently swallowing an error. Options, in preference order:

- Propagate with `?` if the caller returns a `Result`.
- Handle explicitly with a `match` / `if let Err(e) = ...`.
- `.expect("SAFETY: infallible because <reason>")` when you can prove it can't
  fail (e.g., writing to a `Vec<u8>` or a pre-sized `String`). The `expect`
  string documents why.
- **Never `let _ = ...` a `Result`**. It hides the error path from grep, CI,
  and reviewers.

### 2. Iterator chain without terminal → laziness bug

`x.iter().map(|y| side_effect(y))` doesn't run. Fix by adding the right terminal:

- `.for_each(...)` for side effects only.
- `.collect::<Vec<_>>()` when you want the results.
- `.sum()` / `.count()` / `.find(...)` etc. for aggregations.

### 3. Value-returning mutation whose return you don't want → caller-side API smell

The stdlib returns values from mutation methods to make them composable. If you
don't need the return, the caller is often using the wrong method:

| Current call | Better alternative |
|---|---|
| `vec.remove(i)` (ignoring value) | `vec.swap_remove(i)` if order doesn't matter, or `vec.retain(...)` for filter-style |
| `vec.pop()` (just to shrink) | `vec.truncate(vec.len() - 1)` or `vec.clear()` |
| `option.take()` (just to clear) | `*option = None;` |
| `map.remove(k)` (just to delete) | no alternative; document with a comment why the value is uninteresting |
| `map.insert(k, v)` (discarding prior) | defensible — stdlib has no `set()` — document with a comment if the discard is surprising |
| `write!(s, "...")` on infallible writer | `write!(s, "...").expect("infallible: String writer")` |

For the defensible cases (a handful), write:

```rust
let _prev = map.insert(k, v); // prior value intentionally discarded: first-wins semantics
```

Bind to a `_name` variable, not bare `_`. This keeps the reader oriented on
what's being thrown away and allows the reason to surface in grep.

### 4. Tests get a looser policy

`#[cfg(test)]` code may have legitimate bare `let _ = ...` patterns (e.g.,
consuming a `Result` in a setup helper where the test asserts on later state).
Tests can carry `#[allow(unused_results)]` at the module level if volume is
high; production code cannot.

## Implementation steps

1. **Add `unused_results` to `src/lib.rs`** alongside the existing
   `cfg_attr(test, allow(...))` block:

   ```rust
   #![warn(unused_results)]
   ```

2. **Count new warnings before triaging.** Run `cargo build 2>&1 | grep -c 'unused_results'`
   on a clean build and record the number in the PR body. If it's > 200, stop
   and decide whether to stage the rollout (e.g., enable per-module) or accept
   a larger first-pass burden.

3. **Walk each warning through the decision tree above.** Reject `let _ =` as
   a general response. Expect most hits to resolve via option 3 (caller-side
   API refactor) or option 1 (`?` propagation).

4. **Strip all 82 `#[must_use]` annotations** from the 22 files listed above.
   They are now redundant — `unused_results` covers every function unconditionally.

5. **Disable `must_use_candidate`** if it was enabled (current audit shows it
   isn't; confirm).

## Counter-case: when `let _name = ...` IS the right answer in production

There are a handful of genuinely-fine cases that the decision tree above covers:

- **`HashMap::insert(k, v)` in set-semantics code**: stdlib has no
  `.set(k, v) -> ()` variant, and the prior-value `Option<V>` is usually
  uninteresting. Document with `let _prev = ... // set-semantics: prior value
  intentionally dropped`.
- **`Mutex::lock().unwrap()` / `RwLock::write().unwrap()` held for its
  side-effect of locking**: bind to `let _guard = ...` (NOT bare `_` — that
  drops immediately, defeating the lock).
- **Channel `send` on a flaky remote that can disconnect**: document with a
  reason (`// fire-and-forget: downstream may be gone`).

These are documented exceptions, not the default. Production code review should
flag any new bare `let _ = ...` that isn't in one of these patterns.

## Verification

```bash
cargo build                                     # zero unused_results warnings
cargo clippy -q --all-targets                   # zero must_use_candidate warnings
cargo clippy -q --all-targets --features gpu
cargo test
```

After the change, discarding any non-`()` return value anywhere in `src/` produces a
warning without the author needing to remember to add `#[must_use]`.
