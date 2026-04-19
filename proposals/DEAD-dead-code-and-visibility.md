---
id: DEAD
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# DEAD: Dead-code and visibility cleanup (DRSD + SMRT + DHPC)

Bundle of three small, mechanically independent reductions in source-tree
surface area. Bundled because each item is too small to warrant its own
proposal, and the natural review focus (visibility / unused public API /
misleading test) is shared.

## Problem

### DRSD — Real-space `symmetrize_density` is unreachable production code

`src/symmetry/density/real_space.rs` (264 LOC, including its test sub-mod)
is `#[deprecated]` (line 45) and has zero production callers. The PCFX
landing (G-space symmetrization, 2026-04-18) made it permanently
unreachable from any SCF path; it is retained only because three test
sites still cross-check against it:

| Caller                                                | Purpose |
|-------------------------------------------------------|---------|
| `src/symmetry/density/g_space.rs:452,511`             | Cross-check G-space ≡ real-space on symmorphic grids |
| `src/symmetry/density/real_space.rs:182-254`          | Self-tests pinning legacy behaviour |
| `src/symmetry/mod.rs:189`                             | "Identity-only is no-op" sanity test |

Every call site is `#[allow(deprecated)]`-decorated and sits inside a
`#[cfg(test)]` block — i.e. the code only exists to validate itself.
Maintenance cost: ~250 LOC, four `#[allow(deprecated)]` shims, and a
re-export in `src/symmetry/density/mod.rs:30` that newcomers' rust-analyzer
will eagerly suggest as a completion target.

The cross-check value can be preserved by inlining a 30-LOC symmorphic
reference implementation directly inside `g_space.rs`'s test module —
just enough to prove the G-space form behaves like a per-grid-point
average when τ=0 and the grid is compatible. That removes the temptation
for any future code to import `symmetrize_density` thinking it's a real
API.

### SMRT — `symmetry/mod.rs:178-199` test pins the deprecated path

`symmetrize_with_identity_only_is_noop` in `src/symmetry/mod.rs:179`
calls the deprecated real-space `symmetrize_density` to assert the
`n_ops <= 1` short-circuit. Redundant: `test_symmetrize_g_identity_only_is_noop`
in `src/symmetry/density/g_space.rs:470` asserts the same bit-identical
no-op property on the G-space form, which is what SCF actually
exercises. The legacy test is the only thing keeping the
`#[allow(deprecated)]` shim on `src/symmetry/mod.rs:141` alive — once
DRSD removes the deprecated function entirely, the test would no longer
compile anyway. Delete it.

### DHPC — `eigensolver::dense::diagonalize_hermitian` exposed broader than needed

`src/eigensolver/dense.rs:36` exposes `diagonalize_hermitian` as `pub`.
Production callers all go through wrappers:

- `diagonalize_lowest` (line 73) — used by SCF.
- `diagonalize_subspace` (line 220, WFRX warm-start) — used by SCF.

External `pub` callers are tests and benches only:

| Site                                  | What it does |
|---------------------------------------|--------------|
| `tests/free_electron_bands.rs:325`    | Validates all eigenvalues vs sorted diagonal — needs the full spectrum |
| `tests/free_electron_bands.rs:353`    | `H v = λ v` reconstruction — uses only the lowest 10 eigenvectors |
| `benches/scf_benchmarks.rs:80`        | Microbench of full diagonalization |

The integration-test sites are the only reason the function can't be
`pub(crate)`. Site 325 genuinely needs all eigenvalues; site 353 only
needs the first 10. Both can use `diagonalize_lowest` instead — line 75
of dense.rs (`let n = n_bands.min(full.eigenvalues.len())`) means
passing `basis.len()` returns the full spectrum, and passing 10 returns
exactly what site 353 iterates. No new public symbol needed; the
existing `diagonalize_lowest` wrapper is the public entry point.

## Implementation

### Step 1 — DRSD: delete the real-space module

1. Delete `src/symmetry/density/real_space.rs` entirely.
2. In `src/symmetry/density/mod.rs`:
   - Drop `mod real_space;` (line 25).
   - Drop the `#[allow(deprecated)]` + `pub use real_space::symmetrize_density;`
     re-export (lines 29-30).
   - Update the module docstring (lines 1-22) to no longer reference the
     real-space form. Replace the "Two forms are implemented" framing
     with single-form documentation focused on `symmetrize_density_g`.
3. In `src/symmetry/density/g_space.rs`:
   - Lines 437-467 (`test_symmetrize_g_matches_real_space_on_compatible_grid`):
     replace the `symmetrize_density(&mut rho_real, dims, &symmetry)`
     reference with a 30-LOC inline `nint`-based symmetrizer scoped to
     the test (only used as a fixed-point check; symmorphic-only).
   - Lines 491-525 (`test_symmetrize_g_projects_pre_symmetric_density`):
     same treatment — inline the symmorphic reference. The projector
     property `P · (P · ρ) = P · ρ` is what the test pins; the reference
     implementation just needs to produce *some* G-symmetric density to
     project.

### Step 2 — SMRT: delete the redundant test

In `src/symmetry/mod.rs`:
- Delete `symmetrize_with_identity_only_is_noop` (lines 178-199, the
  `#[test]` attribute and fn body).
- Delete the `#[allow(deprecated)]` blanket on line 141 (no longer
  needed after the test is gone).
- Update the comment on lines 137-140 to no longer reference the legacy
  module.
- Coverage is preserved by `test_symmetrize_g_identity_only_is_noop` in
  `src/symmetry/density/g_space.rs:470`.

### Step 3 — DHPC: narrow `diagonalize_hermitian` to `pub(crate)`

1. In `src/eigensolver/dense.rs:36`, change `pub fn diagonalize_hermitian`
   to `pub(crate) fn diagonalize_hermitian`.
2. In `tests/free_electron_bands.rs`:
   - Line 325: `dense::diagonalize_hermitian(&h)` → `dense::diagonalize_lowest(&h, basis.len())`.
   - Line 353: `dense::diagonalize_hermitian(&h)` → `dense::diagonalize_lowest(&h, 10)` (the test only iterates the first 10 eigenpairs anyway, line 357).
3. In `benches/scf_benchmarks.rs:80`: `dense::diagonalize_hermitian(&h)` → `dense::diagonalize_lowest(&h, h.nrows())`. The bench measures full-decomposition cost; passing `h.nrows()` keeps the measurement equivalent (the wrapper is one `subcols(0, n)` call away from the underlying `self_adjoint_eigen`).

### Step 4 — Verify the bundle hangs together

Run the full quality gate. Bundled landing means a single PR — easier to
revert if any one piece misbehaves.

## Verification

```bash
.claude/bin/machine-lock acquire "Core Engineer" "DEAD bundle validation"
cargo test                                              # 265 tests pass; symmetry tests intact
cargo test --features gpu                               # GPU tests still build (no symmetry coupling)
cargo clippy -q --all-targets                           # no new warnings
cargo clippy -q --all-targets --features gpu            # ditto
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps          # no broken intra-doc links
cargo bench --bench scf_benchmarks -- 'eigensolver/(faer|diagonalize)' --sample-size 10
                                                         # diag bench within noise of pre-DEAD
.claude/bin/machine-lock release
```

Acceptance:

- `src/symmetry/density/real_space.rs` no longer exists.
- `cargo test` passes; in particular `test_symmetrize_g_matches_real_space_on_compatible_grid`
  and `test_symmetrize_g_projects_pre_symmetric_density` still pin the
  G-space form against an inline symmorphic reference.
- `grep -r 'diagonalize_hermitian' src/` returns only the definition
  site (now `pub(crate)`) and internal callers; tests use
  `diagonalize_lowest`.
- `LOC delta` is approximately −250 (real_space.rs gone) +30 (inlined
  test reference) −20 (SMRT test gone) ≈ −240 net.

## Out of scope

- No changes to G-space symmetrization correctness. Bit-identical SCF
  results before and after.
- No changes to `check_grid_compatibility` or `compatible_grid_dims` —
  both stay in `src/symmetry/density/mod.rs` (they don't depend on the
  real-space form).
- No changes to `EigenResult`, `diagonalize_lowest`, or
  `diagonalize_subspace` semantics. Only the visibility of one entry
  point shifts and three external call sites swap to the wrapper.
