# No backwards compatibility, no legacy code

**pwdft-rs is pre-release with zero external users.** When a change improves
something, the old version is **deleted** in the same change. No parallel
old/new paths, no "keep it alive for one release," no gentle migrations.

## What this forbids

- `#[deprecated]` wrappers that keep the old API callable.
- `#[serde(alias = "...")]`, `#[serde(rename = "...")]` shims, or
  `#[serde(untagged)]` enums that accept both old and new YAML forms.
- "One-release deprecation warning" paths where the legacy form still
  parses but logs a warning.
- Feature flags that preserve old behavior.
- `if use_old_path { … } else { … }` branches carrying behavior from a
  previous design.
- Re-exports or type aliases whose sole purpose is to keep old import
  paths working.

## What this requires

- The old form is a **hard** error: a compile error (renamed / removed
  function, renamed / removed field), or a hard YAML parse error with a
  message naming the new thing.
- Every in-repo call site migrates in the **same** PR — no follow-up
  "migrate the last fixture" PRs. `inputs/`, `examples/`, `tests/`, and
  `docs/` land together with the code change.
- Proposals that invoke "backwards compatibility" or "legacy support" as
  a justification for added complexity are rejected. Strip the shim; the
  PR gets smaller.

## Why

Every shim is permanent weight. The deprecation never lands, the alias
never gets removed, the feature flag never gets toggled off. They
accumulate until removing them is its own refactor. We have no users to
protect, so the cost is pure waste.

Move forwards. The current codebase is the only thing that has to keep
working.

## Reviewer checklist trigger

If a PR under review introduces any of the forbidden constructs above,
that is an automatic `REQUEST-CHANGES` blocker. Point at this file.
