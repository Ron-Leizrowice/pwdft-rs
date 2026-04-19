---
id: TYPB
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# TYPB: Narrow-int audit — revert premature narrowings; fix `fft_grid_size` sign

## Problem

The earlier TYPE-A work (PR #86 — integer narrowings) traded i32 → i16 for Miller indices and similar narrowings, on the hypothesis that smaller storage + tighter invariants was a net win. Looking back with a season of CAST / CLAU / ERR2-AX suppressions layered on top, several of those narrowings appear not to pay for themselves:

1. **i16 Miller indices.** `src/basis.rs` stores G-vector Miller indices as `i16`. Every call site that reads them widens to `i32` before use (e.g., `src/scf/grid.rs:57`). Storage savings for a typical Si ecut=15 Ry cell: 2000 × 3 × 2 bytes = 12 KB (i16) vs 24 KB (i32) — fits in L1 cache either way. BLAS matrices in the same SCF step are megabytes. The memory win is rounding error; the cost is multiple `#[allow(clippy::cast_possible_truncation)]` + widening casts at every use.

2. **`fft_grid_size(n_max: i32)` with runtime sign assertion.** `src/fft.rs:131` takes `i32` and `assert!(n_max >= 0, ...)` at line 135. This is exactly the code smell the type system is meant to prevent — the parameter's invariant ("non-negative") should be in the type, not the runtime. Right type: `u32` or `usize`. The callers (`src/scf/grid.rs`, `src/symmetry/density/mod.rs`) already have non-negative `usize` / `u32`-sourced values; they only pass `i32` because the signature says so.

3. **`#[allow(clippy::cast_possible_truncation, reason = "asserted elsewhere")]`** sprinkling. Each such suppression is a place where the type system could have caught the issue but didn't because we narrowed too eagerly.

## Proposal

Three changes, none physics-affecting:

### Part A — Revert i16 Miller → i32

- `src/basis.rs`: store Miller indices as `i32`. Delete widening casts at call sites.
- Remove `#[allow(clippy::cast_possible_truncation)]` annotations that guarded i32 → i16 narrowing.
- Remove TYPE-A `reason = "..."` markers that referred to the Miller width.
- Expected net: ~20 lines removed, 8-12 `#[allow]` suppressions removed.

### Part B — `fft_grid_size` takes `u32`

- Change `fft_grid_size(n_max: i32)` → `fft_grid_size(n_max: u32)`. Delete the runtime `assert!(n_max >= 0, ...)`.
- Callers: two sites in `src/scf/grid.rs`, three in `src/symmetry/density/mod.rs`. All pass values known to be non-negative; update cast expressions.
- Remove the corresponding `#[allow(clippy::cast_sign_loss)]` suppressions.

### Part C — Sweep remaining narrow-int suppressions

Run `rg -n 'clippy::cast_possible_truncation|clippy::cast_possible_wrap|clippy::cast_sign_loss' src/` and triage each hit:

- **Keep** if the narrowing is load-bearing (e.g., a type-erased length that genuinely cannot be wider).
- **Revert** if it's a premature optimization whose invariant is a runtime assertion.

Out-of-scope: `benches/` and `tests/` — those have pre-existing bench-code `usize → i32 → usize` modular arithmetic that ERR2-AX already flagged for a dedicated follow-up.

## Performance

Expected impact on SCF wall-time: **none measurable**. The i16 storage was never on a hot path; widening casts fold into loads on aarch64. Benchmark the Si 4×4×4 SCF pre- and post-change to confirm no regression; if the bench shows a real slowdown (> 2%) anywhere, that specific narrowing was actually load-bearing and stays.

## Risk

- **Low.** Mechanical refactor. Existing tests catch any arithmetic regression.
- **Zero to external API.** Miller-index getter currently returns a slice of `i16` — update to `i32`. One downstream consumer in `src/symmetry/kpoints.rs`; unchanged internal arithmetic.

## Non-goals

- Not proposing `num_traits::NumCast` or any generic-over-int-type abstraction. Stay concrete.
- Not revisiting `scale * n_max` style arithmetic in `src/scf/grid.rs`. That stays `i32` — it's a computation, not a storage cell.

## Acceptance

- `rg -n 'i16' src/basis.rs src/scf/` returns zero hits (or documented exceptions).
- `fft_grid_size` signature is `fn fft_grid_size(n_max: u32) -> usize` with no runtime sign assertion.
- Clippy suppression count for `cast_possible_truncation|cast_possible_wrap|cast_sign_loss` in `src/` drops by ≥ 5.
- `cargo bench --bench scf_benchmarks` on Si/Fe shows no regression > 2% on end-to-end SCF wall.
- Full gate green.
