---
id: ERR2
status: active
priority: medium
complexity: small
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

```text
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

P0 (clippy lint floor) landed as PR #86 and ERR2-AX (TYPE-A narrowing
`expect` annotations) landed as PR #111. P1 covers the Category D
`InvalidInput(String)` variant split, fully scoped in the
`## P1 — InvalidInput variant split (2026-04-19 scoping pass)` section
below. The scoping pass supersedes the brief outline originally drafted
here.

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

## P1 — InvalidInput variant split (2026-04-19 scoping pass)

P0 landed as PR #86 (blanket `clippy::unwrap_used` / `expect_used` /
`panic` = warn across the crate, plus the surviving-expect annotation
pass). ERR2-AX landed as PR #111 (TYPE-A narrowing `expect` sites got
`#[expect(reason = "...")]`). TYPB #134 closed the PR #80 spillover. P1
is the remaining Category D work: split today's catch-all
`PwdftError::InvalidInput(String)` into structured variants so callers
can branch on error kind rather than string-match the payload.

### Ground truth re-scan

The original Category D table (drafted pre-P0) listed 12 production call
sites. A full `PwdftError::InvalidInput` scan on `origin/main` as of
2026-04-19 surfaces **15** production sites (three new sites have landed
since the original census: `gaussian_sigma` from CFGN Phase 1,
`resolve_ecutwfc` from the recommended-cutoff path, and the UPF
`angular_momentum < 0` validator from UPFV).

Current `PwdftError` enum (`src/error.rs`, unchanged since P0):

```rust
pub enum PwdftError {
    Io(#[from] std::io::Error),
    ConvergenceFailure { iterations: usize, delta: f64 },
    InvalidInput(String),
    MissingPseudopotential(String),
    Parse(String),
    Eigensolver { size: usize, detail: String },
    Gpu(String),
    NotImplemented { what: String },
}
```

Note the shape is `InvalidInput(String)` — a newtype tuple, not the
`{ what: String }` struct shape the task brief hypothesized. The P1
migration keeps tuple-style consistency where it fits ("unknown element
X" is fundamentally `String`) and switches to struct-style only for
variants that benefit from named fields (`{ name, reason }` parameter
errors).

### 1. Production call-site enumeration

Each row is a production call site (`#[cfg(test)]` bodies excluded).
The *cluster* column groups sites by the invariant kind; sites sharing
a cluster will migrate in the same PR and map to the same new variant.

| # | File:line | What it represents | Cluster | Proposed variant |
|---|---|---|---|---|
| 1 | `src/scf/mod.rs:158` | `ScfParams::n_bands == 0` | **PARAM** | `InvalidParam { name: "n_bands", reason }` |
| 2 | `src/scf/mod.rs:161` | `ScfParams::conv_threshold <= 0` | **PARAM** | `InvalidParam { name: "conv_threshold", reason }` |
| 3 | `src/scf/mod.rs:164` | `ScfParams::mixing_beta ∉ (0, 1]` | **PARAM** | `InvalidParam { name: "mixing_beta", reason }` — interpolates value |
| 4 | `src/scf/mod.rs:169` | `ScfParams::smearing_sigma < 0` | **PARAM** | `InvalidParam { name: "smearing_sigma", reason }` |
| 5 | `src/scf/mod.rs:172` | `ScfParams::ecutrho_ratio < 1` | **PARAM** | `InvalidParam { name: "ecutrho_ratio", reason }` — interpolates value |
| 6 | `src/scf/mod.rs:177` | `ScfParams::nspin ∉ {1, 2}` | **PARAM** | `InvalidParam { name: "nspin", reason }` — interpolates value |
| 7 | `src/scf/mod.rs:187` | `PeriodicPulay::period == 0` | **PARAM** | `InvalidParam { name: "pulay_period", reason }` |
| 8 | `src/scf/mod.rs:196` | `gaussian_sigma` non-finite / non-positive | **PARAM** | `InvalidParam { name: "gaussian_sigma", reason }` — interpolates value |
| 9 | `src/scf/mod.rs:329` | `crystal.atoms` is empty | **CRYSTAL** | `InvalidCrystal { reason: "at least one atom is required" }` |
| 10 | `src/scf/mod.rs:332` | `kpoints` is empty | **CRYSTAL** | `InvalidCrystal { reason: "at least one k-point is required" }` |
| 11 | `src/scf/mod.rs:336` | lattice volume < 1e-10 ų | **CRYSTAL** | `InvalidCrystal { reason: "lattice has zero or near-zero volume" }` |
| 12 | `src/settings.rs:539` | unknown element symbol in YAML | **ELEMENT** | `UnknownElement { symbol: String }` |
| 13 | `src/settings.rs:617` | `basis.ecutwfc` unset + no recommended cutoff tabulated for any species | **PARAM** | `InvalidParam { name: "basis.ecutwfc", reason }` — reason string describes the tabulated-cutoff fallback |
| 14 | `src/main.rs:47` | `BandPath` mode but no high-sym path configured | **CLI** | stays `InvalidInput(String)` (CLI glue only; one site) |
| 15 | `src/pseudopotential/upf/xml.rs:46` | `angular_momentum < 0` in PP_BETA | **UPF** | `InvalidPseudopotential { source: String, reason: String }` — see discussion below |

**Total**: 15 production sites. Proposed outcome:

- **PARAM cluster** (8 sites): `InvalidParam { name: &'static str, reason: String }`
- **CRYSTAL cluster** (3 sites): `InvalidCrystal { reason: &'static str }`
- **ELEMENT cluster** (1 site): `UnknownElement { symbol: String }`
- **UPF cluster** (1 site): `InvalidPseudopotential { source: String, reason: String }`
- **CLI catch-all** (1 site): keep as `InvalidInput(String)` — single CLI-glue call in `main.rs`, not worth a dedicated variant.
- **Post-P1 `InvalidInput` count: 1 / 15 = 6.7 %**, beating the < 20 % target.

### 2. Variant design

Final proposed enum (additions only; existing variants preserved):

```rust
#[derive(Debug, thiserror::Error)]
pub enum PwdftError {
    // --- existing variants unchanged ---
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("SCF did not converge after {iterations} iterations (delta = {delta:.2e})")]
    ConvergenceFailure { iterations: usize, delta: f64 },
    #[error("missing pseudopotential for element {0}")]
    MissingPseudopotential(String),
    #[error("parse error: {0}")]
    Parse(String),
    #[error("eigendecomposition failed for {size}x{size} matrix: {detail}")]
    Eigensolver { size: usize, detail: String },
    #[error("GPU error: {0}")]
    Gpu(String),
    #[error("{what} is not yet implemented")]
    NotImplemented { what: String },

    // --- P1 additions ---

    /// User-supplied parameter failed a range/shape check.
    /// `name` is a compile-time-constant parameter label (YAML key or
    /// `ScfParams` field name). `reason` interpolates the offending value
    /// when useful.
    #[error("invalid parameter {name}: {reason}")]
    InvalidParam { name: &'static str, reason: String },

    /// Crystal / k-point / lattice structural precondition failed.
    /// Used for whole-input shape errors (empty atom list, degenerate
    /// lattice) rather than per-parameter range checks.
    #[error("invalid crystal input: {reason}")]
    InvalidCrystal { reason: &'static str },

    /// YAML referenced an element symbol not in the `crate::atoms`
    /// periodic table.
    #[error("unknown element symbol: {symbol}")]
    UnknownElement { symbol: String },

    /// A pseudopotential file parsed structurally but failed a
    /// physical-validity check (e.g. negative angular momentum, missing
    /// required block). Distinct from `Parse(String)`, which is reserved
    /// for syntax-level UPF errors. `source` is the file path or tag
    /// context; `reason` is the specific invariant that was violated.
    #[error("invalid pseudopotential {source}: {reason}")]
    InvalidPseudopotential { source: String, reason: String },

    /// Catch-all for the one remaining CLI-glue call site in
    /// `src/main.rs`. Not to be used by new code — the next unknown
    /// error kind should become its own structured variant. Post-P1
    /// usage sites: 1 (audit pre-merge).
    #[error("invalid input: {0}")]
    InvalidInput(String),
}
```

Design notes:

- **Field naming.** `InvalidParam` uses `name: &'static str` (all param
  labels are compile-time constants), `reason: String` (values
  interpolated at runtime). `InvalidCrystal` uses `reason: &'static str`
  because all three current sites have static reason strings.
  `UnknownElement` and `InvalidPseudopotential` use owned `String`
  fields because the payload is user-supplied.
- **Why a new `InvalidPseudopotential` vs reusing `Parse`?** The UPF
  `angular_momentum < 0` check at `xml.rs:46` is structurally distinct
  from the tag/attribute parse errors in the same file — the XML
  parsed fine; the *value* was out of physical range. `Parse` already
  covers 15 syntax sites in `upf/xml.rs` and `upf/convert.rs` (all
  using `PwdftError::Parse(format!(...))`); keeping those separate
  means callers can distinguish "malformed file" from "file is
  syntactically valid but carries impossible physics."
- **`InvalidCrystal` scope.** Deliberately narrow — just the three
  pre-SCF shape checks in `scf::run_scf`. Future "crystal-ish"
  validation (e.g. overlapping atoms, inversion-symmetry check
  failure) can either extend the reason space or graduate to its own
  variant.
- **CLI catch-all discussion.** The one surviving `InvalidInput`
  site at `src/main.rs:47` is pure CLI glue — `band_path` mode
  requires a path definition and the user omitted it. The error text
  is already user-facing and a caller unwinding to exit(1) doesn't
  need to discriminate. Keeping `InvalidInput` as a tiny remainder
  plus a doc comment "use a dedicated variant for new error kinds" is
  cheaper than inventing a one-shot `MissingBandPath` variant.
  Optional future step: if a 16th `InvalidInput` site shows up in
  review, re-audit.
- **Visibility.** `InvalidInput` stays `pub`. Making it `pub(crate)`
  would break the existing downstream pattern where integration tests
  (and any future consumers of the `pwdft_rs::error` module) match on
  this variant. Documented as "catch-all for CLI glue; prefer a
  structured variant" in the doc comment is enough friction.

### 3. Migration ordering

Four small PRs, each ≤ 10 call sites, each mergeable independently. No
behavior change; only the error *type* surface widens.

#### P1.a — Add variants + formatting tests (0 call-site migrations)

**Scope:** `src/error.rs` only. Add `InvalidParam`, `InvalidCrystal`,
`UnknownElement`, `InvalidPseudopotential`. Existing `InvalidInput`
usages are unchanged; the 15 production sites still compile.

**Tests added** (in `src/error.rs`'s existing `#[cfg(test)] mod tests`
or new if none exists):

- `fn invalid_param_fmt()` — asserts
  `PwdftError::InvalidParam { name: "n_bands", reason: "must be > 0".into() }.to_string()`
  returns `"invalid parameter n_bands: must be > 0"`.
- `fn invalid_crystal_fmt()` — analogous.
- `fn unknown_element_fmt()` — asserts `{symbol}` interpolation.
- `fn invalid_pseudopotential_fmt()` — asserts `{source}` and
  `{reason}` both appear.

**Acceptance:** `cargo test` Tier-1 + `cargo clippy` clean. No src/
changes outside `error.rs` and its test module.

**Estimated size:** ~50 lines added; ~30 lines of tests.

#### P1.b — Migrate the CRYSTAL cluster (3 call sites)

**Scope:** Sites #9, #10, #11 in `src/scf/mod.rs:329,332,336`. Convert
each `PwdftError::InvalidInput("...".into())` to
`PwdftError::InvalidCrystal { reason: "..." }`.

**Test changes:**

- The integration test at `tests/` currently has no direct match on the
  empty-atoms error (grep confirmed: `PwdftError::InvalidInput` returns
  zero hits in `tests/`). So no test breakage from the migration.
- One *new* integration test `tests/err2_p1_crystal.rs` (or a
  unit test in `src/scf/mod.rs::tests`) asserts that calling
  `scf::run_scf` with an empty atom list returns
  `PwdftError::InvalidCrystal { reason }` containing "atom".

**Acceptance:** `cargo test` Tier-1 still passes (Tier 2 not needed —
no SCF physics touched). `grep 'PwdftError::InvalidInput' src/`
decreases by 3.

**Estimated size:** ~10 lines of production change + ~20 lines of test.

#### P1.c — Migrate the PARAM cluster (8 call sites)

**Scope:** Sites #1–#8 and #13 in `src/scf/mod.rs:158,161,164,169,172,
177,187,196` and `src/settings.rs:617`. Convert each
`PwdftError::InvalidInput(format!(...))` to
`PwdftError::InvalidParam { name: "<param>", reason: "..." }`.

Worked example for site #3 (`mixing_beta`):

```rust
// BEFORE:
if self.mixing_beta <= 0.0 || self.mixing_beta > 1.0 {
    return Err(PwdftError::InvalidInput(
        format!("mixing_beta must be in (0, 1], got {}", self.mixing_beta),
    ));
}

// AFTER:
if self.mixing_beta <= 0.0 || self.mixing_beta > 1.0 {
    return Err(PwdftError::InvalidParam {
        name: "mixing_beta",
        reason: format!("must be in (0, 1], got {}", self.mixing_beta),
    });
}
```

**Test changes required** — the existing tests that currently match
`PwdftError::InvalidInput(msg)` must be updated:

- `src/scf/mod.rs:369` — `validate_rejects_zero_pulay_period` matches
  `PwdftError::InvalidInput(msg)` for site #7 → update to match
  `PwdftError::InvalidParam { name: "pulay_period", reason }`.
- `src/scf/mod.rs:388` — same test, `matches!(…, Err(InvalidInput(_)))`
  → `Err(InvalidParam { name: "pulay_period", .. })`.
- `src/settings.rs:1280` — `scf_params_validate_rejects_non_positive_gaussian_sigma`
  matches `PwdftError::InvalidInput(msg)` for site #8 → update.

**Acceptance:** `cargo test` Tier-1 passes; the nine affected sites
still produce actionable error text (eyeball-verified in test
assertions).

**Estimated size:** ~30 lines of production change + ~30 lines of test
updates.

#### P1.d — Migrate ELEMENT + UPF clusters (2 call sites) + drop `InvalidInput` usage

**Scope:**

- Site #12 in `src/settings.rs:539` → `UnknownElement { symbol }`.
- Site #15 in `src/pseudopotential/upf/xml.rs:46` →
  `InvalidPseudopotential { source: tag.to_string(), reason: format!("angular_momentum must be non-negative (got {l})") }`.
- Site #14 in `src/main.rs:47` stays as `InvalidInput(String)` — the
  CLI catch-all.

**Test changes required:**

- `src/pseudopotential/upf/convert.rs:555,561,596,602` — four test
  match arms currently expect `PwdftError::InvalidInput(msg)` →
  update to `PwdftError::InvalidPseudopotential { source, reason }`
  and assert the `reason` still mentions `angular_momentum`.
- Add one test in `src/settings.rs::tests` for the unknown-element
  path (pattern-match on `PwdftError::UnknownElement { symbol }`).

**Acceptance:**

- `grep 'PwdftError::InvalidInput' src/` returns exactly 1 site
  (`src/main.rs:47`) — down from 15.
- CLAUDE.md § Code Quality (or equivalent) grows a one-line note:
  "`PwdftError::InvalidInput(String)` is reserved for CLI glue in
  `main.rs`. New error kinds get a dedicated variant — see the P1
  taxonomy in `proposals/ERR2-panic-free-production-audit.md`."

**Estimated size:** ~10 lines production change + ~20 lines test
updates + ~5 lines CLAUDE.md.

### 4. Test policy

Minimum coverage per new variant:

- **Formatting test** in `src/error.rs::tests`: one `#[test]` per new
  variant asserting `to_string()` produces the exact `#[error(...)]`
  template output. Landed in P1.a.
- **Integration test per migration PR**: one test per cluster that
  triggers the error path through a realistic caller entry point and
  asserts the new variant discriminant (not the string). These keep
  the migration from silently downgrading an error to
  `InvalidInput` or `Parse`.
- **No regression tests on the old `InvalidInput(String)` shape.**
  The migration deletes the test match-arm on `InvalidInput` for every
  site it moves. Leaving test expectations behind defeats the purpose
  of the split.

Total new test count across P1.a–d: approximately

- P1.a: 4 formatting tests (unit)
- P1.b: 1 integration or unit test
- P1.c: 3 test updates + 1 new test
- P1.d: 4 test updates + 1 new test
- **Total: ~13 test touchpoints, ~8 genuinely new tests.**

Test overhead is deliberately minimal — the P1 migration is structural
refactoring, not a new feature. The existing UPF / YAML / SCF
integration tests already exercise the error *paths*; P1 only changes
the *shape* of the returned discriminant.

### 5. Summary of counts

| Quantity | Pre-P1 | Post-P1 target |
|---|---|---|
| `PwdftError::InvalidInput` production sites | 15 | 1 |
| `PwdftError::InvalidInput` catch-all share | 100 % | 6.7 % |
| New variants added to `PwdftError` | 0 | 4 (`InvalidParam`, `InvalidCrystal`, `UnknownElement`, `InvalidPseudopotential`) |
| PRs | — | 4 (P1.a, P1.b, P1.c, P1.d) |
| Max call-site migrations per PR | — | 9 (P1.c) |
| New tests | — | ~8 new + ~5 updates |

### 6. Surprising findings from the re-scan

- **Site count went 12 → 15, not 12 → 12.** CFGN (gaussian_sigma),
  EAUT-style recommended-cutoff (resolve_ecutwfc), and UPFV
  (angular_momentum) all landed after the original census. Always
  re-scan at the start of any P1-style migration — the drift is one
  site per month on this codebase.
- **`Parse(String)` is stable at 15 sites.** A separate cluster, not
  conflated with `InvalidInput`; no consolidation needed. The task
  brief asked "are half the sites parser errors that should become a
  new `PwdftError::Parse` family?" — the answer is *no*. `Parse` is
  already its own variant and already has 15 sites inside the UPF
  parser. The UPFV site (#15 above) is the *only* `InvalidInput` hit
  in UPF code, and it's a value-range check, not a syntax parse — so
  it needs `InvalidPseudopotential`, not another `Parse`.
- **Existing test match-arms are the main friction.** Seven test
  locations (`scf/mod.rs:369,388`, `settings.rs:1280`,
  `pseudopotential/upf/convert.rs:555,561,596,602`) match on
  `PwdftError::InvalidInput(_)`. Each P1.b–d PR must update the
  specific test arms it touches; P1.a adds no test-arm debt.
- **`InvalidInput { what: String }` in the task brief is a typo** —
  the actual enum is the tuple shape `InvalidInput(String)`. The P1
  plan keeps `InvalidInput` tuple-shaped for backward compatibility
  with the surviving CLI site; new variants use struct fields where
  they help (`InvalidParam { name, reason }`) and tuple fields where
  they don't (`UnknownElement { symbol }` — could be tuple, but named
  for API consistency with sibling variants).

### Anti-scope of P1 (explicit)

- **Does not touch P2** (FFT grid `assert!` → `Result`). That's
  `scf/grid.rs:68` and lives in the separately-scoped Phase 2.
- **Does not tighten clippy lints** beyond the P0 floor. `unwrap_used`
  stays `warn`; no promotion to `deny`.
- **Does not modify `Parse(String)`, `MissingPseudopotential(String)`,
  `Gpu(String)`.** Those already have coherent single-purpose semantics.
- **Does not touch `tests/` integration tests that don't currently
  match on `InvalidInput`.** Net-new tests in P1.b–d are additive.
- **Does not rebase onto any in-flight branch** (GGAP-C, VGCH-1b,
  Al-QE-regen, TSPL, MIXA). The 15 sites touched are all in
  `src/scf/mod.rs`, `src/settings.rs`, `src/main.rs`,
  `src/pseudopotential/upf/xml.rs`, and `src/error.rs` — files not
  owned by any of those parallel proposals.

## 2026-04-21 refresh — post-P1.d census + stricter standard

### P1 completion status

P1.a–d all landed. `InvalidParam`, `InvalidCrystal`, `UnknownElement`, `InvalidPseudopotential` added (P1.a); all 14 production `InvalidInput` call sites migrated (P1.b–d). Sole surviving `InvalidInput` is `src/main.rs:47` (CLI glue, by design).

### User's stated standard (2026-04-21)

> "Only `.expect()` for extremely implausible failure modes, with Results for genuinely fallible codepaths."

Implausible = requires another thread to have panicked, or a caller to have violated a documented precondition that is enforced at construction time.
Fallible = anything a user input, a hardware condition, or a runtime state could legitimately trigger.

### P1.e — transplant.rs two-site mop-up (not in original P1 scope)

Two `InvalidInput` sites in `src/scf/transplant.rs` were missed in the original P1 scan — the file is a VGCH-2B diagnostic stub that didn't exist when the census was written.

| # | File:line | Current | Proposed |
|---|---|---|---|
| 16 | `scf/transplant.rs:121` | `InvalidInput("VGCH-2B transplant diagnostic only supports nspin=1")` | `InvalidParam { name: "nspin", reason: "transplant diagnostic only supports nspin=1" }` |
| 17 | `scf/transplant.rs:130` | `InvalidInput(format!("transplanted rho_g_fft has {} entries; expected n_grid={}", ...))` | `InvalidParam { name: "rho_g_fft", reason: format!(...) }` |

Scope: ~5 lines in one file. No Tier-2 needed (transplant.rs is a diagnostic stub, not SCF production code).

**P1.e acceptance:** `grep 'PwdftError::InvalidInput' src/` returns exactly 1 site (`src/main.rs:47`). `cargo test` Tier-1 green.

### Updated `.expect()` census (2026-04-21 scan, commit `e113379`)

Fresh production census: **18 `.expect()` sites, 0 `.unwrap()`.**

Three new sites vs the original 15-count census — all in `gpu/mod.rs`:

| New site | Pattern | Classification |
|---|---|---|
| `gpu/mod.rs:413` | `pool.scratch_f32.lock().expect("scratch_f32 poisoned")` | **Extremely implausible** — Mutex poison requires another thread to have panicked; unrecoverable state. Keep. |
| `gpu/mod.rs:636` | same pattern | **Extremely implausible.** Keep. |
| `gpu/mod.rs:832` | same pattern | **Extremely implausible.** Keep. |

All three are `scratch_f32.lock().expect(...)` on the GPU buffer pool's internal Mutex. Mutex poison is only triggered when the lock holder panicked — at that point the program state is already unrecoverable and converting to `Result` would just defer an unavoidable abort. All three stay as `expect`.

### Phase 2.5 — gpu/mod.rs:755 (new candidate under stricter standard)

`read_staging_buffer` in `gpu/mod.rs` (~line 755):

```rust
slice.map_async(wgpu::MapMode::Read, |_| {}).expect("BUG: GPU buffer mapping failed");
```

`wgpu`'s `BufferAsyncError` is a real hardware error (device lost, OOM, invalid buffer state). Under the stricter standard this is **genuinely fallible** — the GPU could legitimately fail to map a staging buffer due to device loss or OOM; it does not require a logic bug or mutex poison. The `BUG:` prefix here is a misnomer.

**Proposed Phase 2.5:** Convert `read_staging_buffer` from `-> Vec<f32>` to `-> Result<Vec<f32>, PwdftError>`, propagate via `PwdftError::Gpu(...)` at line 755. Callers in `scf/driver.rs` already return `Result`, so `?` propagates cleanly.

Scope: ~5 lines production + caller propagation.

### Next moves — P1.e + Phase 2 + Phase 2.5 as one small PR

After current agent wave (URES, SKPL, PZPW) clears, bundle the three remaining items:

1. **P1.e** — `scf/transplant.rs:121,130` → `InvalidParam` (~5 lines)
2. **Phase 2** — `scf/grid.rs:68` assert → `PwdftError::InvalidParam { name: "fft_grid", ... }` (~10 lines)
3. **Phase 2.5** — `gpu/mod.rs:755` `read_staging_buffer` → `Result<Vec<f32>, PwdftError>` (~5 lines + caller propagation)

Total: ~30–40 LOC, all non-overlapping files, single PR. Tier-2 required (gpu/ and scf/ both on the trigger list).
