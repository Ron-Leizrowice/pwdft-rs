---
name: No backwards compatibility
description: pwdft-rs is pre-release with no users; proposals must not carry deprecation shims, legacy aliases, or migration warnings
type: feedback
---

Do **not** endorse legacy code paths, backwards compatibility, or deprecation aliases in proposals or implementations. Breaking renames / field moves / schema changes land as hard breaks.

**Why:** pwdft-rs is in active pre-1.0 development with zero external users. The only YAML decks, input fixtures, and API call sites are in-repo; they migrate in the same PR as the breaking change. Carrying `#[serde(alias)]` shims, `#[deprecated]` wrappers, or "one-release warning" paths adds permanent maintenance weight to buy flexibility for users who don't exist.

**How to apply:**

- When reviewing a proposal that adds a shim, `#[serde(alias)]`, or "legacy still parses" clause → strike it. Acceptance criterion becomes: legacy form is a hard parse error with a message naming the new field.
- When a rename / move lands, all in-repo call sites (`inputs/*.yaml`, `examples/`, `tests/`, docs) migrate in the same PR. No follow-up "migrate the last fixture" PRs.
- When a sub-agent proposes "let's add a deprecation warning for one release" — tell them no, hard-break it.
- Move forwards. The current codebase is the only thing that needs to keep working.
