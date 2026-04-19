---
id: TYPB
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# TYPB: Integer type cleanup — revert premature narrowings, fix sign-carrying invariants, retire narrowing `expect` sites

## Problem

Earlier type work (PR #80, TYPE-A) traded `i32 → i16` for Miller indices and `i32 → i8` for crystallographic rotation entries on the hypothesis that smaller storage + tighter invariants was a net win. Looking back with a season of CAST / CLAU / ERR2-AX suppressions and `expect`-site annotations layered on top, several of those narrowings — and one orthogonal sign-as-runtime-assert pattern in `fft.rs` — have visible maintenance cost and no measurable speed or memory benefit at SCF-relevant problem sizes. The `expect` sites the narrowings introduced have been tracked separately in FLUP as `TYPE-AX`; this proposal folds that follow-up in, so there is one unified integer-type cleanup rather than two overlapping PRs.

Four distinct but related issues, in order of impact:

1. **`i16` Miller indices (PR #80 Phase A).** `src/basis.rs` stores G-vector Miller indices as `i16`. Every call site that reads them widens to `i32` before use (e.g. `src/scf/grid.rs:57`). Storage savings for a typical Si ecut=15 Ry cell are 2 000 × 3 × 2 bytes = 12 kB (i16) vs 24 kB (i32) — fits in L1 cache either way, and BLAS matrices in the same SCF step are megabytes. The memory win is rounding error; the cost is multiple `#[allow(clippy::cast_possible_truncation)]` plus widening casts at every use, plus one TYPE-AX-tracked `expect` site at `src/basis.rs:65`.

2. **`fft_grid_size(n_max: i32)` with runtime sign assertion.** `src/fft.rs:131` takes `i32` and `assert!(n_max >= 0, ...)` at line 135. This is exactly the code smell the type system is meant to prevent — the parameter's invariant ("non-negative") belongs in the type, not a runtime assert. Right type: `u32`. The callers (`src/scf/grid.rs`, `src/symmetry/density/mod.rs`) already have non-negative `usize` / `u32`-sourced values; they only pass `i32` because the signature says so.

3. **Narrowing `expect` sites from PR #80 (formerly FLUP → TYPE-AX).** PR #80 introduced five `try_from(...).expect(...)` sites guarded by structural bounds:

   - `src/basis.rs:65` — `i16::try_from(n).expect(...)` on a Miller index. **Auto-resolved by Part A** below (reverting to `i32` deletes the `try_from` entirely).
   - `src/symmetry/operations.rs:71` — `i8::try_from(v).expect("SymmOp::from_flat: rotation entry out of i8 range")`.
   - `src/symmetry/operations.rs:120` — `i8::try_from(v).expect("SymmOp::inverse: adjugate entry out of i8 range")`.
   - `src/symmetry/operations.rs:151` — `i8::try_from(v).expect("SymmOp::compose: product entry out of i8 range")`.
   - `src/symmetry/detect.rs:185` — `i8::try_from(v).expect("symmetry::detect: rotation entry exceeds i8 range")`.

   The four remaining `i8` rotation-entry sites have a crystallographically structural bound: rotation entries live in `{-2, -1, 0, 1, 2}` for all space groups, and adjugate / product entries stay within `{-6..6}` — well inside `i8`. This is a legitimate invariant; the sites are candidates for "add a `reason = "..."` citing the bound and move on", **not** candidates for `Result` propagation. Leaving them as bare `expect` hides a structural invariant from readers.

4. **`#[allow(clippy::cast_possible_truncation, reason = "asserted elsewhere")]` sprinkling.** Each such suppression is a place where the type system could have caught the issue but didn't because we narrowed too eagerly. Several reason strings cite the same (TYPE-A-era) bounds that Parts 1–3 above remove.

Related housekeeping: the separate `cast_lossless` proposal (CLSS) landed as PR #122 and is mechanical cast-style, not integer-width. It is out of scope here — included in this proposal's background only so a future reader doesn't re-consolidate them.

## Proposal

Four coordinated changes, none physics-affecting:

### Part A — Revert `i16` Miller → `i32`

- `src/basis.rs`: store Miller indices as `i32`. Delete the `i16::try_from(...).expect(...)` at `src/basis.rs:65` (one TYPE-AX site auto-closed).
- Delete widening casts and `#[allow(clippy::cast_possible_truncation)]` annotations at call sites that read Miller triples.
- Remove any TYPE-A `reason = "..."` markers that refer to the Miller width.
- Expected net: ~20 lines removed, 8–12 `#[allow]` suppressions removed, 1 `expect` site deleted.

### Part B — `fft_grid_size` takes `u32`

- Change `fft_grid_size(n_max: i32)` → `fft_grid_size(n_max: u32)`. Delete the runtime `assert!(n_max >= 0, ...)`.
- Callers: two sites in `src/scf/grid.rs`, three in `src/symmetry/density/mod.rs`. All pass values known to be non-negative; update cast expressions.
- Remove the corresponding `#[allow(clippy::cast_sign_loss)]` suppressions.

### Part C — Retire the 4 remaining PR #80 `expect` sites

For each of the four `i8` rotation-entry sites (`src/symmetry/operations.rs:71,120,151`, `src/symmetry/detect.rs:185`):

- Preferred form — add a `#[expect(clippy::missing_panics_doc, reason = "...")]` (or equivalent `// BUG:`/`reason = ...` comment matching the ERRH + FGRD pattern) citing the structural bound: rotation entries in `{-2..2}` for 230 space groups, adjugate / product entries in `{-6..6}`. The bound is grep-discoverable at a glance.
- Alternative (out of scope here, documented for the next reader) — if TYPB reviewers decide these are user-reachable at UPF-load-time, convert the enclosing function to return `PwdftError::InvalidInput`. Today's call sites all originate from post-validation space-group detection; the invariant is upheld internally.

### Part D — Sweep remaining narrow-int suppressions

Run `rg -n 'clippy::cast_possible_truncation|clippy::cast_possible_wrap|clippy::cast_sign_loss' src/` and triage each hit:

- **Keep** if the narrowing is load-bearing (e.g., a type-erased length that genuinely cannot be wider).
- **Revert** if it's a premature optimization whose invariant is a runtime assertion.
- **Re-annotate** if the `reason = "..."` string is stale after Parts A / B / C (e.g. still cites "Miller narrowing").

Out-of-scope: `benches/` and `tests/` — those have pre-existing bench-code `usize → i32 → usize` modular arithmetic that ERR2-AX flagged for a dedicated follow-up.

## Performance

Expected impact on SCF wall-time: **none measurable**. The `i16` Miller storage was never on a hot path; widening casts fold into loads on aarch64. The `i8` rotation entries stay `i8` (Part C is annotations-only) so the 4.7–5.4 % `symmetrize_density_g_n72` win from PR #80 is preserved.

Benchmark the Si 4×4×4 SCF and the `symmetrize_density_g_n{18,36,72}` bench pre- and post-change to confirm no regression; if any bench shows a real slowdown (> 2 %), that specific narrowing was actually load-bearing and stays.

## Risk

- **Low.** Mechanical refactor with annotations-only Part C. Existing tests catch any arithmetic regression.
- **Zero external API surface.** Miller-index getter currently returns a slice of `i16` — update to `i32`. One downstream consumer in `src/symmetry/kpoints.rs`; unchanged internal arithmetic.

## Non-goals

- Not proposing `num_traits::NumCast` or any generic-over-int-type abstraction. Stay concrete.
- Not revisiting `scale * n_max` style arithmetic in `src/scf/grid.rs`. That stays `i32` — it's a computation, not a storage cell.
- Not touching `i8` rotation storage itself. PR #80's 4.7 % `symmetrize_density_g_n72` win came from the rotation narrowing, not the Miller narrowing; that win stays.
- Not bundling `cast_lossless` or `missing_errors_doc` / `missing_panics_doc` work — those were CLSS (PR #122, landed).

## Acceptance

- `rg -n '\bi16\b' src/basis.rs src/scf/` returns zero hits (or documented exceptions).
- `fft_grid_size` signature is `fn fft_grid_size(n_max: u32) -> usize` with no runtime sign assertion.
- All five PR #80 `expect` sites either (a) carry a `reason = "..."` annotation citing the structural bound, or (b) are deleted by the Part A revert. No bare `expect` remains in the narrowing path.
- Clippy suppression count for `cast_possible_truncation | cast_possible_wrap | cast_sign_loss` in `src/` drops by ≥ 5.
- `cargo bench --bench scf_benchmarks` on Si/Fe shows no regression > 2 % on end-to-end SCF wall, and `symmetrize_density_g_n72_ops48` stays within 2 % of today's post-PR-#80 baseline (~17.70 ms).
- Full quality gate green (`cargo test`, both clippy invocations, `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps`).

## Provenance

- Consolidates TYPB (original narrow-int audit, 2026-04-19) + FLUP's TYPE-AX entry (PR #80 follow-up, 2026-04-19).
- Supersedes the TYPE-AX FLUP entry as of this consolidation. FLUP entry struck with a pointer here.
