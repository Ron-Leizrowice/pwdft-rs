---
id: ERR2
status: draft
priority: medium
complexity: medium
risk: low
depends_on: [ERRH]
blocks: []
---

# ERR2: Production-Code Panic/Unwrap Census + Enforcement

> **Scope reminder:** Proposal only, no code changes. The audit is a snapshot
> of `src/` against the standing standard "concrete custom errors via
> `thiserror` for every failure mode, no panics or unwraps outside tests."
> ERRH (completed 2026-04-17) fixed the initial 25 production sites. This
> proposal inventories what remains post-ERRH and post-CAST, and proposes a
> clippy-enforced posture so future regressions fire at CI time rather than
> in code review.

## Problem

The `PwdftError` enum (`src/error.rs`) has eight variants today:

```
Io, ConvergenceFailure, InvalidInput, MissingPseudopotential,
Parse, Eigensolver, Gpu
```

ERRH migrated 25 call sites away from panic-family macros; CAST added a
handful of release-mode `assert!` guards for cast-safety invariants
(documented with `reason = "…"`). No clippy lints currently forbid the
panic family, so:

1. There is no CI-level floor on unwrap/expect/panic counts — regressions
   are only caught when a reviewer notices.
2. The boundary between "legitimate invariant" and "should-have-been-Result"
   is fuzzy; some recent `expect` calls carry `BUG:` prefixes while an
   equivalent site in `scf/mixing/broyden.rs` (landed after ERRH in the
   BROY series) does not.
3. `PwdftError::InvalidInput(String)` has accreted twelve distinct failure
   kinds — stringly-typed error handling that defeats the value of the
   `#[error(...)]` annotation.

## Census (post-ERRH, post-CAST, pre-FGRD)

Counts below are **production-only** (excluding `#[cfg(test)]` and
`mod tests { … }` bodies). Numbers verified against every candidate
`grep` hit by ±10-line context read.

### Category A — unwrap / expect / panic

| File | `.unwrap()` | `.expect(…)` | `panic!` | Notes |
|---|---|---|---|---|
| `fft.rs` | 0 | 4 | 0 | 4× `BUG: Array3 should be contiguous` around the 2-buffer FFT plan (FFTB). Legitimate invariant. |
| `gpu/mod.rs` | 0 | 4 | 0 | 1× `BUG: g_squared buffer not allocated`; 3× `BUG: GPU readback channel closed / BUG: GPU buffer mapping failed`. |
| `ewald.rs` | 0 | 1 | 0 | `BUG: atom has no matching pseudopotential (should have been validated at startup)` — documented precondition via `ScfContext::new`. |
| `scf/initial_density.rs` | 0 | 1 | 0 | Same `ScfContext::new` precondition; same `BUG:` message pattern as `ewald.rs`. |
| `scf/mixing/broyden.rs` | 0 | 1 | 0 | `.expect("Broyden+Kerker mode requires g_squared")` — **missing the `BUG:` prefix** used everywhere else. Cosmetic inconsistency. |
| `scf/mixing/anderson.rs` | 0 | 1 | 0 | `BUG: Kerker mode requires g_squared to be provided by caller`. Matches ERRH pattern. |
| `scf/driver.rs` | 0 | 1 | 0 | `BUG: invalid progress bar template` — `indicatif::ProgressStyle::with_template` on a static literal. |
| `symmetry/density/mod.rs` | 0 | 2 | 0 | 2× `BUG: empty fixed-size array` on `[usize; 3].iter().max()`. Could be simplified via `array::IntoIter::max()` returning `usize` directly, but the panic is vacuous (not a user-reachable path). |
| **Production total** | **0** | **15** | **0** | — |

**Post-ERRH diff:** ERRH reported 25 sites fixed and 13 `expect()` surviving
with `BUG:` prefixes. Current count is 15 `expect()`, +2 since ERRH
(`scf/mixing/broyden.rs:66` from BROY and `scf/driver.rs:103` from the
driver split), 0 regressions to `unwrap()` or `panic!()`.

All 15 production `expect` calls are **legitimate invariants** — none would
benefit from `Result` propagation because:
- FFT `as_slice()` calls are guaranteed by row-major `Array3` layout
- GPU buffer channels are closed exactly once after `device.poll()`
- PP-lookup `expect`s are preceded by `ScfContext::new` validation
- `[usize; 3].iter().max()` is vacuously non-empty
- `indicatif::ProgressStyle::with_template` with a static string cannot fail

These stay. The only cleanup is standardizing the comment style
(§ Phase 0 below).

### Category B — `assert!` / `debug_assert!` (grey zone)

Production-only release-mode `assert!`/`assert_eq!`:

| File | Site | Classification | Rationale |
|---|---|---|---|
| `fft.rs:65, 90` | `assert_eq!(data.len(), nx*ny*nz)` | **Keep-as-assert** | Precondition on a `&mut [Complex64]` slice. Bug in caller, not user input. |
| `fft.rs:135` | `assert!(n_max >= 0)` in `fft_grid_size` | **Keep-as-assert** (CAST) | Prevents sign-loss hang; has `reason = "…"` allow on the downstream cast. |
| `numerics.rs:28` | `assert_eq!(n, rab.len())` in `simpson_integrate` | **Keep-as-assert** | Caller-contract precondition. Documented in `# Panics` doc section. |
| `symmetry/density/real_space.rs:51` | `assert_eq!(rho.len(), n_grid)` | **Keep-as-assert** | On a deprecated real-space symmetrizer. |
| `symmetry/density/g_space.rs:163, 164` | `assert_eq!(rho_r.len(), n_grid)`, FFT dims | **Keep-as-assert** | Caller contract; PCFX-landed. |
| `symmetry/operations.rs:90` | `assert!(d == 1 \|\| d == -1)` in `SpaceGroupOp::inverse` | **Keep-as-assert** | Detected-rotation determinant is ±1 by construction. Defense against a future non-rotation `SpaceGroupOp`. |
| `symmetry/detect.rs:33` | `assert!(ops.iter().any(\|op\| op.is_identity…))` | **Keep-as-assert** | Post-construction sanity: identity is always in the group. Bug in detect, not user input. |
| `potential/nonlocal.rs:132` | `assert!(proj.l >= 0)` in `NonlocalPotential::new` | **Keep-as-assert** (belt-and-suspenders post-UPFV) | UPFV has landed (commit `eb405f9`, PR #72) — `angular_momentum < 0` now rejected at parse time in `pseudopotential/upf/xml.rs:44-47`. This `assert!` is now defence-in-depth. Add a CAST-style `reason` comment cross-referencing the UPF validator. |
| `potential/nonlocal.rs:331, 332` | `debug_assert_eq!(h.nrows() == self.n_pw)` | **Keep-as-debug_assert** | Dev-only GEMM shape check. Correct classification. |
| `potential/nonlocal.rs:366` | `assert!(lmax >= 0)` in `real_sph_harmonics` | **Keep-as-assert** (CAST) | Release-mode guard with `reason` on downstream sign-loss cast. |
| `potential/nonlocal.rs:367` | `debug_assert_eq!(out.len(), (lmax+1)² as usize)` | **Keep-as-debug_assert** | Pairs with the release-mode l≥0 check. |
| `potential/nonlocal.rs:493` | `assert!(l >= 0)` in `spherical_bessel_j` | **Keep-as-assert** (CAST) | Same family as `real_sph_harmonics`. |
| `potential/xc.rs:238` | `debug_assert_eq!(n, rho_down_r.len())` | **Keep-as-debug_assert** | Spin-channel grid-size match. |
| `gpu/mod.rs:299, 300` | `assert_eq!(v_h.len(), n_grid)`, etc. | **Keep-as-assert** | `v_eff_assembly` caller contract. Documented in the `#[allow(cast_…)]` `reason`. |
| `scf/mixing/anderson.rs:175` | `assert!(m >= 1)` in `AndersonMixer::diis_step` | **Keep-as-assert** | Doc comment explicitly says "Panics if called before any history". Private API. |
| `scf/mixing/anderson.rs:279` | `assert!(period >= 1)` in `PeriodicPulayMixer::new` | **Keep-as-assert** (already guarded by `ScfParams::validate`) | `ScfParams::validate` at `src/scf/mod.rs:124-129` already rejects `period == 0` with `InvalidInput`. This is the belt-and-suspenders counterpart. Add a comment cross-referencing `ScfParams::validate`. |
| `scf/grid.rs:68` (FGRD, PR #73) | `assert!(dims[i] <= MAX_FFT_DIM)` in `FftGrid::new` | **Convert → Result** in Phase 2 | User-reachable via pathological `ecut` + tiny lattice; should be `PwdftError::InvalidParam` so users get a helpful message. Not urgent — no legitimate calculation reaches 1024³. |

**No urgent Result conversions.** UPFV has closed the only
"should-be-Result" candidate in this list (negative angular-momentum is
now a parse-time `PwdftError::Parse`).

**FGRD landed (PR #73).** `src/scf/grid.rs:68` now has
`assert!(dims[i] <= MAX_FFT_DIM)` with `MAX_FFT_DIM = 1024`, enforcing the
cast-safety invariant used throughout `symmetry/density/*` and
`scf/grid.rs` index-wrapping code. **This site is the one genuine
"convert-to-Result in Phase 2" candidate** in the whole census — a
pathological user input (very high ecut + tiny cell) could theoretically
trigger it, and `PwdftError::InvalidParam { name: "fft_grid", reason:
"dims exceed MAX_FFT_DIM=1024" }` is more useful than a release-mode
panic. Not urgent (1024 is huge; no legitimate calculation reaches it),
hence deferred to Phase 2 after the clippy posture is in place.

### Category C — `unreachable!()` / `todo!()`

| File | Site | Rationale |
|---|---|---|
| `scf/mixing/anderson.rs:87` | `MixingMode::Broyden \| MixingMode::PeriodicPulay => unreachable!(…)` inside `AndersonMixer::new` | Branch-pruning in a sum-type dispatch; `AndersonMixer::new` is constructed only from `Mixer::new` after a match that never forwards these variants. Comment is explicit. **Keep.** |

Zero `todo!()`, zero `unimplemented!()`, zero `FIXME` or `XXX` in live code
(one `FIXME(faer-upstream)` comment in `eigensolver/iterative.rs:17-20`
tracking an external issue, not a production panic).

### Category D — `PwdftError::InvalidInput(String)` catch-all

The `InvalidInput` variant currently absorbs **twelve distinct failure
kinds** (counting `src/` production hits; test-only uses excluded):

| # | Site | String | Proposed structured variant |
|---|---|---|---|
| 1 | `scf/mod.rs:98` | `"n_bands must be > 0"` | `InvalidParam { name: "n_bands", reason: "must be > 0" }` |
| 2 | `scf/mod.rs:101` | `"conv_threshold must be positive"` | `InvalidParam { name: "conv_threshold", reason: "must be positive" }` |
| 3 | `scf/mod.rs:104` | `"mixing_beta must be in (0, 1], got {x}"` | `InvalidParam` with value |
| 4 | `scf/mod.rs:109` | `"smearing_sigma must be non-negative"` | `InvalidParam` |
| 5 | `scf/mod.rs:112` | `"ecutrho_ratio must be >= 1, got {x}"` | `InvalidParam` with value |
| 6 | `scf/mod.rs:117` | `"nspin must be 1 or 2, got {x}"` | `InvalidParam` with value |
| 7 | `scf/mod.rs:127` | `"pulay_period must be >= 1 for PeriodicPulay mixing"` | `InvalidParam` |
| 8 | `scf/mod.rs:212` | `"at least one atom is required"` | `InvalidCrystal("empty atom list")` |
| 9 | `scf/mod.rs:215` | `"at least one k-point is required"` | `InvalidCrystal("empty k-point list")` |
| 10 | `scf/mod.rs:219` | `"lattice has zero or near-zero volume"` | `InvalidCrystal("degenerate lattice")` |
| 11 | `settings.rs:456` | `"unknown element: {sym}"` | `UnknownElement(String)` — this is definitely user-facing. |
| 12 | `main.rs:47` | `"band_path mode requires a band path definition"` | `MissingBandPath` or keep in `InvalidInput` (CLI-glue only). |

**Recommended consolidation** (not this proposal's implementation scope):

```rust
#[error("invalid parameter {name}: {reason}")]
InvalidParam { name: &'static str, reason: String },

#[error("invalid crystal input: {0}")]
InvalidCrystal(&'static str),

#[error("unknown element symbol: {0}")]
UnknownElement(String),
```

This collapses the 12 sites into 3 structured variants, shrinks the
remaining `InvalidInput(String)` to main-binary glue (CLI-mode errors),
and gives callers a machine-readable discrimination (useful for the
eventual language-binding / CLI exit-code layer).

`Parse(String)` has 15 sites (UPF XML parsing); no consolidation needed —
all parse-error messages already interpolate the failing tag/field, and
a single variant is idiomatic for parser errors.

## Proposed enforcement: clippy lints

Add to `Cargo.toml [lints.clippy]`:

```toml
# --- ERR2: panic-family enforcement in production code ---
# Deliberately `warn` not `deny` so a module-level `#[allow(... reason = "...")]`
# can be used at surviving sites following CAST's discipline. Clean CI is
# still enforced via `cargo clippy -q --all-targets -- -D warnings` in the
# pre-PR gate (see CLAUDE.md § Code Quality).
unwrap_used      = "warn"
expect_used      = "warn"
panic            = "warn"
unreachable      = "warn"

# `todo` is always a bug in production — no reason should ever survive review.
todo             = "deny"

# Rarely-triggered companions; flip these on opportunistically:
unimplemented    = "deny"
```

### Test-scope override

Rust 2024 does not yet have a first-class per-cfg lint-group toggle, so the
test-suppression approach is per-test-module:

```rust
#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: test bodies are allowed to panic — production is enforced by the same lints at `warn` level"
)]
mod tests { … }
```

…emitted once per `mod tests` block in `src/`. Count to touch: ~33 test
modules (confirmed via scan for `mod tests` in `src/`). Integration tests
in `tests/` and benches in `benches/` are unaffected in spirit — clippy
lint tables in `Cargo.toml` apply to the whole crate including tests, but
`cargo clippy --all-targets --features gpu` post-fix should still be
warning-free because every `tests/` file would need analogous
module-level `#[allow]`s. The concrete Phase 0 PR will apply the same
`#[allow]` banner at the top of each `tests/*.rs` and `benches/*.rs`
file, or (cleaner) use the file-level form `#![allow(...)]`.

## Phased migration plan

Four PRs, each mergeable independently and not blocking each other beyond
the stated order:

### Phase 0 — Turn the lights on

- Add the `[lints.clippy]` block above to `Cargo.toml`.
- Add file-level `#![allow(clippy::unwrap_used, clippy::expect_used,
  clippy::panic, reason = "ERR2 § Phase 0")]` to every
  `tests/*.rs` and `benches/*.rs`.
- Add module-level `#[allow(...)]` on every `mod tests` in `src/` (~33
  sites).
- For each of the 15 surviving production `expect(…)` calls, add a
  line-level `#[allow(clippy::expect_used, reason = "<specific
  invariant>")]` matching the existing `BUG:` comment. Standardize the
  comment format so `scf/mixing/broyden.rs:66` gains the `BUG:` prefix
  the rest of the codebase uses.
- For `scf/mixing/anderson.rs:87` (`unreachable!`), add
  `#[allow(clippy::unreachable, reason = "<rationale>")]`.
- Zero behavior change. CI turns green on the new lints; regressions
  now fire as warnings (and `-D warnings` in the clippy gate turns them
  into errors for PRs).

Acceptance: `cargo clippy -q --all-targets` and `cargo clippy -q
--all-targets --features gpu` both clean (zero new warnings). `cargo
test` unchanged.

### Phase 1 — Category D splitting (stringly-typed → structured)

- Add `InvalidParam`, `InvalidCrystal`, `UnknownElement` variants to
  `PwdftError`.
- Migrate the 12 `InvalidInput(String)` sites in the census table.
- Update match-arms in `src/scf/mod.rs::tests` (lines 248, 254, 267)
  and `src/settings.rs::tests` to match on the new variants.
- No clippy lint change; pure `PwdftError` enrichment.

Acceptance: `grep PwdftError::InvalidInput src/` returns ≤ 2 sites (main
CLI glue + a deliberate catch-all at the `Settings` level).

### Phase 2 — Assert → Result where genuinely user-reachable

The only concrete target is `src/scf/grid.rs:68` (FGRD, PR #73):

```rust
// BEFORE (release-mode panic):
assert!(dims[0] <= MAX_FFT_DIM && dims[1] <= MAX_FFT_DIM && dims[2] <= MAX_FFT_DIM, …);

// AFTER (user-facing error with actionable message):
if dims[0] > MAX_FFT_DIM || dims[1] > MAX_FFT_DIM || dims[2] > MAX_FFT_DIM {
    return Err(PwdftError::InvalidParam {
        name: "fft_grid",
        reason: format!(
            "dims {}x{}x{} exceed MAX_FFT_DIM={} (lower ecut, or set fft_grid: [nx, ny, nz] explicitly in YAML)",
            dims[0], dims[1], dims[2], MAX_FFT_DIM
        ),
    });
}
```

`FftGrid::new` becomes `fn new(…) -> Result<Self>`. Callers in `ScfContext::new`
already return `Result`, so the `?` propagates cleanly.

`potential/nonlocal.rs:132` already has UPFV as its upstream guard
(PR #72) — no downgrade needed. Just document the relationship in Phase 0.

Acceptance: Production `assert!` count ≤ Phase-0 count minus the number
of FGRD-introduced grid-size bounds (one or two).

### Phase 3 — Opportunistic tightening

- Review new additions every 3 months against the census baseline.
- Consider promoting `unwrap_used = "warn"` → `"deny"` after no
  regressions for 3 months.

Anti-scope reminder: Phase 2 does **not** remove legitimate invariant
checks (shape preconditions on `&mut [f64]`, lookups behind
`ScfContext::new` validation, etc.). Those are the kind of defensive code
the standard explicitly allows.

## Anti-scope

This proposal does **not**:

1. Remove or weaken `assert!` / `debug_assert!` checks that document
   caller invariants. The standard forbids unwraps on user input, not
   defensive programming against internal bugs.
2. Audit WGSL shaders (`src/gpu/shaders/*.wgsl`) — WGSL does not have
   Rust's panic family and is out of scope.
3. Touch test-directory code semantics (`tests/`, `benches/`). Tests are
   allowed to panic; Phase 0 only adds `#![allow]` banners so clippy
   remains silent there.
4. Modify `src/error.rs` during this audit. The variant additions in
   Phase 1 are a separate change tracked against this proposal.
5. Block on in-flight agent work (FGRD, UPFV). Flags are in § Category
   B for post-landing review.

## Flagged for follow-up (out of this proposal's scope)

- **FGRD (landed as PR #73):** `src/scf/grid.rs:68`
  `assert!(dims[i] <= MAX_FFT_DIM)` is user-reachable on pathological
  inputs. Convert to `PwdftError::InvalidParam` in Phase 2.
- **UPFV (landed as PR #72):** Parse-time rejection of
  `angular_momentum < 0` is live. `potential/nonlocal.rs:132` is now
  belt-and-suspenders and stays as `assert!`; in Phase 0 just add a
  `reason` comment cross-referencing the UPF validator.
- **`scf/mixing/broyden.rs:66`** (BROY spillover): add the `BUG:` prefix
  to the `expect` message. Trivially fixable in Phase 0.
- **`symmetry/density/mod.rs:69, 86`**: `[usize; 3].iter().max()` could
  be rewritten as `dims.into_iter().max().unwrap_or_default()` or just
  `dims[0].max(dims[1]).max(dims[2])`. Cosmetic, not a panic concern.
- **Standalone follow-up nit**: `basis.rs:171` has an `unwrap()` inside
  a test — confirmed test-scope via `mod tests` at line 120, so no
  production concern; noted here because a casual grep flags it.

## Acceptance criteria (Phase 0 only)

1. **Clippy clean** under both `cargo clippy -q --all-targets` and
   `cargo clippy -q --all-targets --features gpu`. Zero new warnings.
2. **Every surviving `expect`, `unreachable`, `panic!` in `src/`** carries
   a line-level `#[allow(clippy::X, reason = "…")]` whose reason text
   identifies the specific invariant — not a generic "safe" or "won't
   happen".
3. **No behavior change.** `cargo test` before and after shows identical
   test counts and pass/fail state.
4. **Regressions fire at CI.** A proof-of-concept new `.unwrap()` added
   anywhere outside a `mod tests` block is caught by clippy at `warn`
   level and (because of the pre-PR `-D warnings` gate) blocks the PR.
5. **`scf/mixing/broyden.rs:66`** picks up the `BUG:` prefix aligning it
   with the rest of the codebase.
