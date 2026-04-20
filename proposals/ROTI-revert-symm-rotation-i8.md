---
id: ROTI
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# ROTI: Revert SpaceGroupOp rotation from i8 to i32 and fix the cache-line comment

## Problem

Two related issues introduced by TYPE Phase A (`operations.rs` commit, 2026-04-18):

### 1. Wrong cache-line count in the module docstring

`src/symmetry/operations.rs:17` claims:

> "packs all 48 ops of Fd-3m into ~2 cache lines"

The actual numbers:

| Field | Size per op | 48 ops |
|-------|-------------|--------|
| `rotation: [[i8;3];3]` | 9 B | 432 B |
| padding (to align `f64`) | 7 B | 336 B |
| `translation: [f64;3]` | 24 B | 1 152 B |
| **Total struct** | **40 B** | **1 920 B = 30 cache lines** |

Even the rotation fields alone (432 B) are ~7 cache lines. The "~2" figure is off by a factor of 15. It appears to conflate the rotation-only footprint with the full `Vec<SpaceGroupOp>`.

The comparison baseline (i32 rotation) is also understated: `[[i32;3];3]` = 36 B + 4 B padding + 24 B = **64 B per op → 48 cache lines**. So i8 saves ~18 cache lines (30 vs 48), not 25 (2 vs 27).

### 2. Complexity that outweighs the benefit

TYPE Phase A measured a **composite** win from three simultaneous changes: i8 rotation + i16 Miller + `index_map` deletion. TYPB subsequently isolated and reverted the i16 Miller change, showing **−0.1% noise** for that sub-change. The isolated contribution of i8 rotation is unmeasured.

The i8 rotation adds the following complexity to support a type whose values are always in `{−2..=2}`:

| Added construct | File | Lines |
|-----------------|------|-------|
| `rotation_i32()` widen helper | `operations.rs:227–235` | 9 |
| 3 × `to_i8` closures + `#[expect]` blocks | `operations.rs:139–144, 182–188, 77–79` | ~18 |
| 27 `i8::try_from()` calls across `from_flat`, `inverse`, `compose` | `operations.rs` | 27 |
| `narrow_rotation()` + `#[expect]` | `detect.rs:189–206` | 18 |
| `transpose_rotation()` (now does widen + transpose) | `g_space.rs:33–40` | 8 |
| 9 `i32::from()` calls in `transpose_rotation` | `g_space.rs:36–38` | 9 |
| 9 `i32::from()` widening calls in `trace()`, `rotation_i32()` | `operations.rs` | 9 |

Estimated **−80 to −90 lines** on revert.

Furthermore, the hot path in `symmetrize_density_g` already precomputes
`r_ts: Vec<[[i32;3];3]>` before the per-G-point loop — so i8 storage only
affects the single cold setup pass, not the actual inner loop, making the
cache-density argument weaker than it appears.

## Implementation

**Step 1 — Change field types in `operations.rs`.**

```rust
// Before
pub struct SymmOp {
    pub rotation: [[i8; 3]; 3],
}
pub struct SpaceGroupOp {
    pub rotation: [[i8; 3]; 3],
    pub translation: [f64; 3],
}

// After
pub struct SymmOp {
    pub rotation: [[i32; 3]; 3],
}
pub struct SpaceGroupOp {
    pub rotation: [[i32; 3]; 3],
    pub translation: [f64; 3],
}
```

**Step 2 — Fix the module docstring comment.**

Replace the "~2 cache lines" claim with accurate numbers: 30 cache lines
for 48 i32-rotation ops (64 B/struct) or a brief honest description of the
struct size.

**Step 3 — Simplify `SymmOp` methods in `operations.rs`.**

- `from_flat`: remove `to_i8` closure and `#[expect]`; assign i32 entries directly.
- `inverse`: remove `to_i8` + `#[expect]`; the adjugate result is already computed in i32, assign directly.
- `compose`: remove `to_i8` + `#[expect]`; matrix-multiply stays in i32, assign directly.
- `rotation_i32()`: delete entirely; call sites that used it now use `self.rotation` directly.
- `trace()`: remove 3 `i32::from(...)` wrappings; entries are already `i32`.
- `apply()`: remove 9 `f64::from(r[...])` wrappings; replace with `r[...] as f64`.

**Step 4 — Update `SpaceGroupOp::new`.**

Signature changes from `rotation: [[i8;3];3]` to `rotation: [[i32;3];3]`.
Update the one `detect.rs` call site and the test helpers in `operations.rs`.

**Step 5 — Delete `narrow_rotation()` from `detect.rs`.**

Three call sites: lines 45, 169, 215. Each currently does
`SymmOp::from_flat(...)` or `SpaceGroupOp { rotation: narrow_rotation(&r), ... }`.
With i32 storage, pass `r` directly. Remove the function and its `#[expect]`.

**Step 6 — Simplify `transpose_rotation()` in `g_space.rs`.**

The function widens i8 → i32 and transposes in one pass. With i32 storage,
it only transposes:

```rust
// Before: transpose + widen
fn transpose_rotation(r: &[[i8; 3]; 3]) -> [[i32; 3]; 3] {
    [
        [i32::from(r[0][0]), i32::from(r[1][0]), i32::from(r[2][0])],
        ...
    ]
}

// After: transpose only
fn transpose_rotation(r: &[[i32; 3]; 3]) -> [[i32; 3]; 3] {
    [
        [r[0][0], r[1][0], r[2][0]],
        [r[0][1], r[1][1], r[2][1]],
        [r[0][2], r[1][2], r[2][2]],
    ]
}
```

The docstring "performs both the transpose and the i8 → i32 widen" becomes
"transposes the rotation matrix".

**Step 7 — Update the `density/mod.rs` grid-compatibility check.**

`src/symmetry/density/mod.rs` has `i32::from(r[i][j])` calls for a
divisibility check. These become bare `r[i][j]` reads.

## Previous attempt (2026-04-20, aborted)

A core-engineer agent (ID `a7fe9ca547748e383`) implemented Steps 1–7
but aborted before `/pr-submit`. The work is preserved on branch
**`origin/ROTI/revert-symm-rotation-i8`** at commit `8b6fbbc`.

- **Diff**: +48 / −151 LOC across 4 files (`src/symmetry/operations.rs`,
  `src/symmetry/detect.rs`, `src/symmetry/density/g_space.rs`,
  `src/symmetry/density/mod.rs`). Matches the proposal's −80 to −90
  projection (actually slightly larger).
- **Not done**: `/quality-gate`, `/test --tier2` (symmetry is a Tier-2
  trigger), PR creation, session logbook.

**Recovery plan for the next agent.** Check out the preserved branch,
rebase on `origin/main`, run `/pr-draft` immediately as a safety
checkpoint, then finish `/quality-gate` + `/test --tier2` + `/pr-submit`.
Don't re-derive the diff — it's structurally complete.

```bash
git -C <MAIN> fetch origin
git worktree add .claude/worktrees/<new-agent> -b ROTI-resume/revert-symm-rotation-i8 origin/ROTI/revert-symm-rotation-i8
cd .claude/worktrees/<new-agent>
git rebase origin/main
/pr-draft "ROTI recovery pickup — Steps 1-7 inherited from 8b6fbbc"
# ...quality gate + tier-2...
/pr-submit
```
