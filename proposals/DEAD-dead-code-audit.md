---
id: DEAD
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# DEAD: Dead-Code Audit (post MODR / NCFX / PCFX / ITEV)

## Motivation

Ten PRs merged on 2026-04-18 (NCFX, PCFX, CCMX, PRPL, ITEV, MODR-A/C/D,
NLCC audit, DBGC). Aggressive physics refactoring + three rounds of module
splits leaves orphans: `pub` items whose last external caller moved into
the same module, helpers the new pipeline no longer exercises, and stale
proposal/doc references to deleted files.

This is an **investigation-only** audit — no source edits. Findings point
at roughly 10 small atomic PRs totalling ~150 LOC of removals + a handful
of `pub → pub(crate)` tightenings. None change behaviour; `cargo test` is
the only quality gate.

## Findings

### 1. Unused / over-scoped `pub` items

| Path:line | Symbol | Current vis | Suggested | External callers | Notes |
|---|---|---|---|---|---|
| `src/potential/hartree.rs:22` | `hartree_potential` | `pub fn` | **delete** | 0 prod, 1 self-test | GPU path uses `gpu.hartree_potential`; CPU path uses `scf::energy::hartree_on_fft_grid`. Module test exercises a now-orphan function. |
| `src/potential/hartree.rs:48` | `hartree_energy` | `pub fn` | **delete** | 0 | SCF uses `scf::energy::hartree_energy` (crate-internal). Superseded. |
| `src/potential/local.rs:93` | `v_local_matrix_element` | `pub fn` | **delete** | 0 | Comment: "For now, we use the direct structure-factor approach". Grid-based pipeline is now the only one. |
| `src/kpoints.rs:72` | `fcc_high_sym_points` | `pub fn` | **delete** | 0 | 25 LOC + 0 callers outside its own module. The YAML `band_structure` path builds `HighSymPoint`s from user input via `Settings::to_high_sym_path` (`src/settings.rs:500`). |
| `src/symmetry/density/mod.rs:61` | `compatible_grid_dims` | `pub fn` | **delete** (or `#[cfg(test)]`) | 0 | Researcher flagged this on 2026-04-17 (`logbooks/researcher.md:430`). Only rotation-compatible; post-PCFX the G-space symmetrizer handles any dims. |
| `src/symmetry/density/mod.rs:41` | `check_grid_compatibility` | `pub fn` | `pub(crate)` → `#[cfg(test)]` | 0 outside `symmetry::density::*` tests | Ditto. Called only by the deprecated real-space symmetrizer's tests. |
| `src/potential/xc.rs:36` | `XcPoint` | `pub struct` | `pub(crate)` | 0 | Returned from `lda_xc` (only exercised by internal tests + `lda_xc_grid`). |
| `src/potential/xc.rs:189` | `XcSpinPoint` | `pub struct` | `pub(crate)` | 0 | Same — only ever consumed inside `xc.rs`. |
| `src/scf/initial_density.rs:31,39` | `InitialDensityConfig` + impl | `pub struct` | `pub(crate)` | 0 | Constructed only inside `scf::mod.rs`. |
| `src/scf/density.rs:11` | `DensityGrid<'a>` | `pub struct` | `pub(crate)` | 0 | Same — constructed only in `scf::mod.rs`. |
| `src/bandstructure.rs:10,53` | `BandStructure` + impl | `pub struct` | keep | 1 (main.rs) | Legitimate library API. |

**Aggregate:** 5 deletable functions (~140 LOC with tests), 4 `pub → pub(crate)` tightenings.

### 2. Narrowest-scope candidates

Focused on the big offenders only, per scope:

- `src/symmetry/density/real_space.rs::symmetrize_density` is `pub` → re-exported by facade. All non-test callers were removed by PCFX. `#[deprecated]` is correct and retained for test pinning (see §5). No action beyond that.
- `src/eigensolver/iterative.rs:67,74` — `DEFAULT_TOL` / `DEFAULT_MAX_RESTARTS` are `pub const`; callers are only the module's own tests. Could be `pub(crate)`.

### 3. `#[allow(dead_code)]` / `#[allow(unused)]`

Exactly one hit crate-wide:

- `src/gpu/mod.rs:45` — `#[allow(dead_code)]` on `BufferPool { real_bufs, real_staging, … }`. The fields are **allocated but never read** (verified: `real_bufs` / `real_staging` appear only in `BufferPool::new` and in the struct declaration — no getter, no consumer). The comment says "reserved for LDA XC pooled path" but LDA XC on GPU uses ephemeral buffers (`src/gpu/mod.rs` line ~600 and the bench). **Smell, not legit.** Either wire the pool into `lda_xc_grid` on GPU (separate PR) or drop the fields + the allocation sites (simpler, -~40 LOC incl. alloc scaffold).

### 4. Orphan / stale module references

- None. All `pub mod` declarations in `src/lib.rs`, `src/scf/mod.rs`, `src/symmetry/mod.rs`, `src/potential/mod.rs` resolve to live files. The MODR-A/C/D splits are clean.
- **Doc-side orphan:** `proposals/HD5I-hdf5-io.md:141` still references `src/input.rs`, deleted by YAML migration (legacy #32). `INDEX.md:110` already flags this; HD5I body should be updated when the proposal unblocks, not now. Listed here for traceability only.

### 5. `#[deprecated]` items

- `src/symmetry/density/real_space.rs:39` — `symmetrize_density`. **Legitimate.** All live callers are under `#[cfg(test)]` (see `real_space.rs:113` `#[allow(deprecated)]`) and `g_space.rs:417,473` cross-check tests. Pins the short-circuit + symmorphic (τ=0) behaviour the G-space form must match. Keep.

### 6. Commented-out code

Grep for `^\s*//\s*(fn|pub fn|match|if let|for|let|use)\s` returns 5 hits, all descriptive prose (English sentences that happen to start with "for " or "let "). **No action.**

### 7. TODO / FIXME ages

| Location | Tag | Age | Verdict |
|---|---|---|---|
| `src/eigensolver/iterative.rs:17,20` | `FIXME(faer-upstream)` | 0 days (2026-04-18, ITEV PR #45) | Active; cross-linked to `tests/itev_iterative_eigensolver.rs:100` `#[ignore]`. Keep. |

No other `TODO(`/`FIXME(`/`XXX`/`HACK` in `src/` or `tests/`. Prior backlog fully drained.

### 8. Duplicated expressions flagged in MODR follow-up

`ctx.v_local_g0 * ctx.n_electrons` appears 6 times in `src/scf/mod.rs` (lines 503, 512, 572, 944, 962, 1028). Originally MODR flagged 4 sites at lines 415/418/808/820; line numbers shifted post-PCFX. Small helper (`ctx.local_g0_shift()`) would collapse to a single expression. **Low priority** — pure stylistic consolidation.

### 9. Config-gated dead code

`#[cfg(feature = "gpu")]` audit: every gated item has a non-gpu counterpart (CPU fallback path in `scf::mod` and `gpu::mod`'s `try_new() -> Option`). Clean.

### 10. Unused `pub use` re-exports

All five crate-level `pub use` statements have external callers:

- `src/atoms.rs:5` `pub use mendeleev::Element` — widely used
- `src/eigensolver/mod.rs:4` `pub use dense::EigenResult` — used by `scf::mod`
- `src/symmetry/mod.rs:8` `pub use operations::{SpaceGroupOp, SymmOp}` — used by `detect.rs`, `kpoints.rs`
- `src/symmetry/density/mod.rs:30,32` — live facade

### 11. Hand-rolled linear solver

`src/scf/mixing/linalg.rs::solve_linear_system` (66 LOC Gauss-elim w/ partial pivoting). **Not dead** — 2 callers (`anderson.rs`, `broyden.rs`). Flagged in file docstring as "future swap" for `faer::linalg::solvers::PartialPivLu`. For 2–8 × 2–8 systems the cost is negligible, so this is taste/hygiene, not dead code. **Out of DEAD scope.**

## Recommended cleanups (atomic, prioritised)

1. **Delete `v_local_matrix_element`** (`src/potential/local.rs:93-120` + its test). Zero callers; obsoleted by grid path.
   Impact: -~30 LOC, public surface -1. Risk: **low**.

2. **Delete `hartree_potential` + `hartree_energy` in `src/potential/hartree.rs`** + the module's two tests; delete `pub mod hartree;` in `src/potential/mod.rs`.
   Impact: -~100 LOC, public surface -2 fns, -1 module. Risk: **low**. (Both replaced by `scf::energy` + `gpu::hartree_potential`.)

3. **Delete `fcc_high_sym_points`** (`src/kpoints.rs:72-99`).
   Impact: -~28 LOC, public surface -1. Risk: **low**.

4. **Delete `compatible_grid_dims`** and move `check_grid_compatibility` to `#[cfg(test)] mod` inside `real_space.rs`.
   Impact: -~35 LOC, public surface -2. Risk: **low**.

5. **Tighten visibility to `pub(crate)`: `XcPoint`, `XcSpinPoint`, `InitialDensityConfig`, `DensityGrid<'a>`**.
   Impact: public surface -4, 0 LOC. Risk: **low**. (No behaviour change.)

6. **Remove `BufferPool::{real_bufs, real_staging}`** and their allocation in `BufferPool::new`, then drop the `#[allow(dead_code)]`.
   Impact: -~40 LOC, clears the last `allow(dead_code)` in the crate. Risk: **medium** (touches GPU scaffolding; confirm the buffers are truly unread before deletion).

7. **Introduce `ScfContext::local_g0_shift()` helper** replacing the 6-site `ctx.v_local_g0 * ctx.n_electrons` expression.
   Impact: 6 sites → 1 definition; -~5 LOC net. Risk: **low**.

8. **Tighten `DEFAULT_TOL` / `DEFAULT_MAX_RESTARTS` to `pub(crate)`** in `src/eigensolver/iterative.rs:67,74`.
   Impact: public surface -2. Risk: **low**.

9. **Leave `symmetrize_density` + `#[deprecated]` as-is.** Documented decision; pins G-space symmetrizer's symmorphic-τ behaviour.

10. **HD5I / INDEX doc fix** — out of DEAD scope; flag for Technical Writer when HD5I unblocks.

Items 1–5 can ship as one PR ("delete 5 dead items, tighten 4 visibilities"); 6–8 are independent follow-ups.

## What this is NOT

- **Not removing `#[deprecated] symmetrize_density`.** Test-pinning use is legitimate and documented in-file.
- **Not triaging old TODOs into proposals.** None are stale; the only `FIXME` is fresh (1 day).
- **Not renaming for taste.** No "unclear name" findings — only real dead-removals and visibility tightenings.
- **Not touching `src/scf/mixing/linalg.rs`.** Not dead; faer-LU swap is a separate (future) decision.
- **Not rewriting HD5I** or other doc orphans — flagged for the owning role.

## Open questions (for EM)

1. **`BufferPool::{real_bufs, real_staging}`** — was the GPU LDA XC pooled path intentionally deferred, or is this a planning artefact? If deferred-with-intent, keep the allocation + add a tracking proposal; if forgotten, delete.
2. **`fcc_high_sym_points`** — any intent to keep as a library convenience for downstream users building band paths programmatically? If yes, document + add a test; if no, delete.

## Flagged for follow-up

- `proposals/HD5I-hdf5-io.md:141` references `src/input.rs` (deleted). **Technical Writer** should retarget to YAML `Settings` when HD5I is picked up.
- `src/scf/mixing/linalg.rs::solve_linear_system` — 66-LOC hand-rolled Gauss elim with a "swap to faer LU" note. **Core Engineer / Performance Engineer** call: at 2–8×2–8 sizes the swap is probably wash, but worth measuring before/after once MXBA / BROY iterations land.
- `src/scf/mod.rs` — 6× duplicated `ctx.v_local_g0 * ctx.n_electrons` is a fine helper-extraction PR for **Core Engineer** (item #7 above).
