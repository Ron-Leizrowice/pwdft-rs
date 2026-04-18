---
id: FLUP
title: Follow-up backlog — seeded from 2026-04-18 code reviews and agent reports
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# FLUP — Follow-up backlog from 2026-04-18

Seven distinct follow-up items surfaced during today's 13-PR merge wave but
were deliberately not fixed inline — either they were out of their source
PR's scope, blocked on upstream, or better-owned by a different agent
role. Rather than leave them as comments in review replies (which
evaporate), this file records them with enough detail to spin each into
its own proposal and PR when the EM schedules it.

Each entry includes a suggested 4-letter ID, owner role, priority,
concrete file paths, and an acceptance criterion. When the EM activates
an entry, it gets promoted to a standalone `proposals/<ID>-<slug>.md`
and the FLUP entry is struck through (not deleted — history of what was
seeded when).

## Entries

### G0SH — DRY the V_local(G=0) shift expression

- **Role:** Core Engineer
- **Priority:** low, **Complexity:** small, **Risk:** low
- **Source:** MODR audit ("Flagged for follow-up" section of
  `proposals/completed/MODR-modular-refactor.md`); re-surfaced by MODR-B.

The expression `+ ctx.v_local_g0 * ctx.n_electrons` appears four times
around total-energy and Harris-Foulkes computations — originally at
`src/scf/mod.rs:415, 418, 808, 820` pre-MODR-B, now redistributed across
`src/scf/driver.rs` and `src/scf/driver_spin.rs`. Replace with a
`with_g0_shift(energy, ctx)` helper in `src/scf/energy.rs` and call it
from all four sites.

**Acceptance criterion:** every call site passes through the new helper;
all existing tests still pass bit-identical. No changes to `EnergyComponents`
semantics.

### GLUS — Replace hand-rolled Gauss-elim with `faer` LU

- **Role:** Performance Engineer
- **Priority:** low, **Complexity:** small, **Risk:** low
- **Source:** MODR audit; `src/scf/mixing/linalg.rs::solve_linear_system`.

`solve_linear_system` is a hand-rolled partial-pivoting Gauss-elimination
for the DIIS/Broyden normal-equations system (≤ `mixing_ndim` × `mixing_ndim`,
typically ≤ 8×8). `faer` is already linked project-wide. Swap to
`faer::linalg::lu::partial_pivoting::solve` or equivalent. The current
implementation is correct (17 unit tests in `mixing/`), but replacing
it removes a hand-maintained linalg primitive and small-matrix edge-case
risk surface.

**Acceptance criterion:** `solve_linear_system` becomes a thin wrapper
(or is deleted outright if the callers can inline the `faer` call);
all existing mixer tests pass; benchmarks show no regression (the matrices
are tiny, so wins are unlikely — the point is maintainability, not speed).

### VNLB — Block-wise `D·B^H` construction in V_NL assembly

- **Role:** Performance Engineer
- **Priority:** medium, **Complexity:** small, **Risk:** low
- **Source:** VNLM code review, PR #49 nit 3.

`src/potential/nonlocal.rs` currently builds `D·B^H` via scalar
`faer::Mat` indexing (lines ~252-281 on the landed commit). This is the
root cause of the 1.5× regression in `vnl_new` wall-time (43 → 78 ms at
n_pw=725) that pushes the VNLM GEMM lift's break-even point to ~2 SCF
iterations. A per-(atom, l) block-wise `faer::matmul` should recover
most of the one-time cost and move break-even below 1 iteration.

**Acceptance criterion:** `bench scf_benchmarks -- hamiltonian/vnl_new`
at n_pw=725 drops to ≤ 55 ms (break-even ≤ 1 SCF iter).
Correctness: all 265+ tests pass bit-identical.

### SYMP — Parallelize `symmetrize_density_g` with `par_iter_mut`

- **Role:** Performance Engineer
- **Priority:** low, **Complexity:** small, **Risk:** low
- **Source:** PCFX code review, PR #44 nit 3.

`src/symmetry/density/g_space.rs::symmetrize_density_g` (post-MODR-C
path) has a serial per-G-vector outer loop. For a 100³ grid with N_ops=48
(Fd-3m) this is ~5M MACs per SCF iteration — small in absolute terms
but a free win with `rayon::prelude::par_iter_mut`. Inner loop writes
into one G-point's slot and reads only rotated-source slots, so
parallelization is safe without atomics.

**Acceptance criterion:** `bench scf_benchmarks` shows ≥ 2× speedup
for the symmetrization step on large grids; all PCFX tests pass
bit-identical (they must — `par_iter` is a pure permutation of a
commutative sum).

### FDLT — Expose `ScfResult.final_delta`

- **Role:** Core Engineer
- **Priority:** low, **Complexity:** small, **Risk:** low
- **Source:** CCMX code review, PR #43 nit 1 + TACC finding (silent
  regression risk).

Add `pub final_delta: f64` (the last Δρ the SCF saw before convergence
or `max_iter`) to `ScfResult` in `src/scf/mod.rs`. Compute it in both
`src/scf/driver.rs` and `src/scf/driver_spin.rs` — both already track
the value in `last_delta` locally but discard it. Update
`tests/spin_polarization.rs::test_ccmx_fe_free_magnetization_converges`
to assert `result.final_delta < 1e-2` — the specific pathology the
CCMX limit cycle produced was Δρ pinned at 0.254, not iter-count, so
this is the pathology-specific regression guard the review flagged.

**Acceptance criterion:** `ScfResult` has the new field; CCMX Fe test
asserts on it; at least one more convergence test (Si Γ-only or the
BROY regression) also asserts a reasonable `final_delta` upper bound.

### VNMT — m-isolation defense-in-depth test for V_NL

- **Role:** Core Engineer
- **Priority:** low, **Complexity:** small, **Risk:** low
- **Source:** VNLM code review, Step 2 of the VNLM proposal's validation
  plan that never landed.

Add an integration or unit test that pins a single-(l, m)-channel V_NL
matrix element: build a test `NonlocalPotential` with `D_ij = 0` except
for one `(l=2, m=+2)` entry, then assemble H at two specific G-vectors
and compare to a hand-computed value using the known real-Y_lm formula.
Any future edit to `real_sph_harmonics` or the shift-and-flip matmul
builder that wrongly normalizes a single m-channel (e.g. √2 off on
m≠0) would be caught — the addition-theorem test only pins sums over m.

**Acceptance criterion:** new test passes post-landing; manually
breaking the `(l=2,m=+2)` normalization (e.g. multiplying by √2)
makes the test fail with an intelligible error message.

### ITVF — Flip ITEV default to `Iterative` post-faer-upstream-fix

- **Role:** Core Engineer (blocked on external faer fix)
- **Priority:** medium (when unblocked), **Complexity:** small, **Risk:** low
- **Source:** ITEV code review, PR #45 follow-ups.

Tracker proposal for the day faer upstream fixes the `iterate_lanczos`
Gram-Schmidt infinite-loop bug. When that lands:
1. Un-`#[ignore]` `tests/itev_iterative_eigensolver.rs::itev_iterative_matches_dense_si_total_energy` and remove the `FIXME(faer-upstream)` markers.
2. Add a proper `iterative_cold_n{89,259,725}` bench in `benches/scf_benchmarks.rs` (the current bench has only the dead-code `let _ = &iterative::DEFAULT_TOL;` import-retention hack).
3. Benchmark Dense vs Iterative SCF wall-time at production `n_pw` (ITEV proposal projected 2.5-4× speedup).
4. If speedup is confirmed, flip `EigensolverKind` default in `src/settings.rs` from `Dense` to `Iterative`.

**Acceptance criterion:** integration test un-ignored + green; bench
shows ≥ 2× speedup at n_pw ≥ 200; default flipped.

**Blocker:** faer upstream issue not yet filed. **File-it-now action**
(separate from this tracker): a Performance Engineer should
file the issue against faer 0.24 with the minimum reproducer pwdft-rs
already has (Si Γ-only at n_pw=89) and paste the issue URL into the
`FIXME(faer-upstream)` markers in `src/eigensolver/iterative.rs` and
`tests/itev_iterative_eigensolver.rs`. One-PR action, not worth a
proposal — just do it.

### RDOC — Clean up 15 pre-existing rustdoc warnings

- **Role:** Technical Writer
- **Priority:** low, **Complexity:** small, **Risk:** low
- **Source:** DOCS post-MODR sweep, PR #53 Technical Writer flag.

`cargo doc --no-deps` emits 15 warnings on main today, all pre-dating
the DOCS sweep. Three clusters:
- `src/pseudopotential/upf/mod.rs:7-9` — rustdoc links `[xml]` and
  `[convert]` point at private sub-modules; two warnings.
- `src/settings.rs:70` — lattice docstring `[[ax,ay,az], ...]` trips
  rustdoc's link resolver; three warnings.
- `src/symmetry/operations.rs:17` — `rotation[i][j]` in a `///`
  docstring.

Mix of escaping brackets (`\[foo\]` or backticks) and dropping
would-be-links to private items. No text rewrites, just mechanical fixes.

**Acceptance criterion:** `cargo doc --no-deps` emits zero warnings.

### UPFV — Parse-time validation for UPF `angular_momentum`

- **Role:** Researcher or Core Engineer
- **Priority:** low, **Complexity:** small, **Risk:** low
- **Source:** CAST code review, PR #56 nit 2.

`src/pseudopotential/upf/xml.rs::extract_beta_angular_momentum` returns
`am_str.trim().parse().ok()`, which happily accepts `-1`, `-2`, etc.
CAST added a runtime assertion in `NonlocalPotential::new` (real
belt-and-suspenders, not vacuous), but a parser-level
`PwdftError::InvalidInput("angular_momentum must be non-negative")`
would (a) fail fast at UPF load time instead of at projector-build
time and (b) give a clearer error message ("bad UPF file" vs "internal
assertion"). Add a regression test that synthesizes a malformed UPF
with `angular_momentum="-1"` and asserts the parse returns
`PwdftError::InvalidInput`.

**Acceptance criterion:** parse of a malformed UPF fails with a
domain-specific error string; `NonlocalPotential::new` runtime assert
remains as defense-in-depth.

### FGRD — Explicit FFT-grid upper-bound check

- **Role:** Core Engineer
- **Priority:** low, **Complexity:** trivial, **Risk:** low
- **Source:** CAST code review, PR #56 nit 3.

Several CAST `#[allow(clippy::cast_possible_truncation, reason = "...")]`
sites in `src/scf/grid.rs:80,106` and `src/symmetry/density/*.rs` cite
"FFT grid dims ≤ ~512 per axis in practice". In practice ecut + Miller
truncation + grid-size heuristics keep us well under 512, but the
bound is load-bearing on convention, not a checked invariant. Either
(a) add an explicit `debug_assert!(nx <= 1024, "FFT grid too large")`
at `FftGrid::new`, or (b) rewrite the reason strings to cite the `ecut`
bound that actually limits them — whichever makes the invariant
visible to a future grep.

**Acceptance criterion:** every CAST `#[allow]` that currently cites
"grid ≤ 512" either backs the bound with a `debug_assert!` or cites a
first-principles bound (ecut ≤ X Ry → grid ≤ Y).

## What this is NOT

- **Not an implementation plan.** Each entry needs to be promoted to
  its own proposal before coding starts (except the "file-it-now" ITEV
  upstream issue, which is a one-PR drive-by).
- **Not ordered by priority.** Ordering within this file is roughly
  "flagged first came first"; the EM picks a real order when activating.
- **Not a promise to land all seven.** Some may turn out to be not
  worth the round-trip cost — ITVF in particular is a tracker whose
  preconditions may never arrive (faer may deprecate the affected API
  path and we swap to LOBPCG, closing ITVF as obsolete).

## Proposal-file etiquette when activating

When the EM activates an entry:
1. Create `proposals/<ID>-<slug>.md` with proper frontmatter.
2. Copy the FLUP entry's body into the new proposal's **Motivation**
   and **Acceptance criterion** sections; expand into full Implementation
   plan + Verification sections.
3. Strike through the FLUP entry (prepend `~~` to each line) so the
   seeding history stays legible, but don't delete.
4. Add the new proposal's row to `proposals/INDEX.md`.
