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

## FLUP — Follow-up backlog from 2026-04-18

Follow-up items surfaced during 2026-04-18's 13-PR merge wave (and the
2026-04-19 follow-on reviews) but were deliberately not fixed inline —
either they were out of their source PR's scope, blocked on upstream, or
better-owned by a different agent role. Rather than leave them as comments
in review replies (which evaporate), this file records them with enough
detail to spin each into its own proposal and PR when the EM schedules it.

Each entry includes a suggested 4-letter ID, owner role, priority,
concrete file paths, and an acceptance criterion. When the EM activates
an entry, it gets promoted to a standalone `proposals/<ID>-<slug>.md`
and the FLUP entry is struck through (not deleted — history of what was
seeded when).

### Status summary (post-2026-04-19 sweep)

- **Landed:** G0SH, GLUS, SYMP, FDLT, VNMT, RDOC, UPFV, FGRD, MXB1, MXB3, VNLT, VNLB (struck by VNLT), DWGT, DFLT, G2ZT (PR #153).
- **Still live (drive-by):** none.
- **Still live (larger):** MXB2 (Fe CCMX retune, small-medium), ITVF (tracker, blocked on faer 0.25).
- **Still live (performance investigation):** EIGV, EIGW (bench-noise triage; may self-resolve on next bench pass).
- **Added 2026-04-19:** TYPE-AX (5 `try_from` expect sites flagged by ERR2 P0 report).
- **Added 2026-04-18 (TRV2-F3 wake):** FLP3 (NLCC ρ_core(G) parametric expansion — 60 unpinned PPs; defensive, ~1 day).

### Entries

#### ~~G0SH — DRY the V_local(G=0) shift expression~~

- **Role:** Core Engineer
- **Priority:** low, **Complexity:** small, **Risk:** low
- **Source:** MODR audit ("Flagged for follow-up" section of
  `proposals/completed/MODR-modular-refactor.md`); re-surfaced by MODR-B.
- **Status:** done — `scf::energy::with_g0_shift` helper landed; 4 arithmetic
  sites in `scf/driver.rs` and `scf/driver_spin.rs` now route through it.

The expression `+ ctx.v_local_g0 * ctx.n_electrons` appears four times
around total-energy and Harris-Foulkes computations — originally at
`src/scf/mod.rs:415, 418, 808, 820` pre-MODR-B, now redistributed across
`src/scf/driver.rs` and `src/scf/driver_spin.rs`. Replace with a
`with_g0_shift(energy, ctx)` helper in `src/scf/energy.rs` and call it
from all four sites.

**Acceptance criterion:** every call site passes through the new helper;
all existing tests still pass bit-identical. No changes to `EnergyComponents`
semantics.

#### ~~GLUS — Replace hand-rolled Gauss-elim with `faer` LU~~ (landed)

~~- **Role:** Performance Engineer~~
~~- **Priority:** low, **Complexity:** small, **Risk:** low~~
~~- **Source:** MODR audit; `src/scf/mixing/linalg.rs::solve_linear_system`.~~

~~`solve_linear_system` is a hand-rolled partial-pivoting Gauss-elimination~~
~~for the DIIS/Broyden normal-equations system (≤ `mixing_ndim` × `mixing_ndim`,~~
~~typically ≤ 8×8). `faer` is already linked project-wide. Swap to~~
~~`faer::linalg::lu::partial_pivoting::solve` or equivalent. The current~~
~~implementation is correct (17 unit tests in `mixing/`), but replacing~~
~~it removes a hand-maintained linalg primitive and small-matrix edge-case~~
~~risk surface.~~

~~**Acceptance criterion:** `solve_linear_system` becomes a thin wrapper~~
~~(or is deleted outright if the callers can inline the `faer` call);~~
~~all existing mixer tests pass; benchmarks show no regression (the matrices~~
~~are tiny, so wins are unlikely — the point is maintainability, not speed).~~

Landed as `faer::Mat::partial_piv_lu().solve_in_place(rhs)` thin wrapper
in `src/scf/mixing/linalg.rs`; the singular-pivot fallback (uniform
`1/(n+1)` coefficients on `|U[i,i]| < 1e-15`) is preserved by scanning
`lu.U()` diagonals after factoring. Kept the wrapper — callers in
`anderson.rs` and `broyden.rs` pass flat `Vec<f64>` + `n`, and inlining
would leak `faer::Mat` plumbing (and the per-call singular guard) into
two sites. All 39 `scf::mixing` tests pass bit-identical, full release
test suite green (225 lib + integration), both clippy invocations clean,
`RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` clean. See PR GLUS.

#### ~~VNLB — Block-wise `D·B^H` construction in V_NL assembly~~ (struck 2026-04-19 by VNLT)

- **Role:** Performance Engineer
- **Priority:** medium, **Complexity:** small, **Risk:** low
- **Source:** VNLM code review, PR #49 nit 3.
- **Status:** struck — the 1.5× regression that motivated VNLB does not
  exist. VNLT investigation bench-verified three clean runs of
  `hamiltonian/vnl_new_n725` at 44.37, 44.36, 44.23 ms on current main —
  matching the pre-VNLM baseline (42.33 ms) to within 5 %. The 78.5 ms
  in VNLM PR #49's table was a single-run criterion outlier. The bench
  harness (`benches/scf_benchmarks.rs`) is byte-identical across the
  VNLM merge and today's HEAD; only CAST `#[allow]` attributes and two
  cheap asserts have touched `src/potential/nonlocal.rs` since VNLM
  merged (commits `1291729`, `6840237`), and neither can explain a
  ~35 ms wall-time swing. See the "Note 2026-04-19" section of
  `proposals/VNLM-vnl-blocked-matmul.md` for the full audit.

~~`src/potential/nonlocal.rs` currently builds `D·B^H` via scalar
`faer::Mat` indexing (lines ~252-281 on the landed commit). This is the
root cause of the 1.5× regression in `vnl_new` wall-time (43 → 78 ms at
n_pw=725) that pushes the VNLM GEMM lift's break-even point to ~2 SCF
iterations. A per-(atom, l) block-wise `faer::matmul` should recover
most of the one-time cost and move break-even below 1 iteration.~~

~~**Acceptance criterion:** `bench scf_benchmarks -- hamiltonian/vnl_new`
at n_pw=725 drops to ≤ 55 ms (break-even ≤ 1 SCF iter).
Correctness: all 265+ tests pass bit-identical.~~

~~### SYMP — Parallelize `symmetrize_density_g` with `par_iter_mut`~~

~~- **Role:** Performance Engineer~~
~~- **Priority:** low, **Complexity:** small, **Risk:** low~~
~~- **Source:** PCFX code review, PR #44 nit 3.~~

~~`src/symmetry/density/g_space.rs::symmetrize_density_g` (post-MODR-C~~
~~path) has a serial per-G-vector outer loop. For a 100³ grid with N_ops=48~~
~~(Fd-3m) this is ~5M MACs per SCF iteration — small in absolute terms~~
~~but a free win with `rayon::prelude::par_iter_mut`. Inner loop writes~~
~~into one G-point's slot and reads only rotated-source slots, so~~
~~parallelization is safe without atomics.~~

~~**Acceptance criterion:** `bench scf_benchmarks` shows ≥ 2× speedup~~
~~for the symmetrization step on large grids; all PCFX tests pass~~
~~bit-identical (they must — `par_iter` is a pure permutation of a~~
~~commutative sum).~~

**Landed:** parallelized via `par_chunks_mut(ny·nz)` over destination
G-point xy-planes (not `par_iter_mut` — per-item rayon overhead eats
small-N_ops wins; xy-slab chunking amortizes scheduling across
`ny·nz` slots and matches cache locality). Inner per-slot reduction
stays serial → bit-identical result, all 10 PCFX tests pass. On a
72³ grid: ops=48 106.3 → 17.7 ms (6.0×), ops=8 25.8 → 10.7 ms (2.4×);
scaling holds at 36³ (ops=48: 5.2×). Crossover below 18³·8 (1.1×),
above which every tested config clears 2×. See PR SYMP.

#### ~~FDLT — Expose `ScfResult.final_delta`~~ (landed)

~~- **Role:** Core Engineer~~
~~- **Priority:** low, **Complexity:** small, **Risk:** low~~
~~- **Source:** CCMX code review, PR #43 nit 1 + TACC finding (silent~~
~~ regression risk).~~

~~Add `pub final_delta: f64` (the last Δρ the SCF saw before convergence~~
~~or `max_iter`) to `ScfResult` in `src/scf/mod.rs`. Compute it in both~~
~~`src/scf/driver.rs` and `src/scf/driver_spin.rs` — both already track~~
~~the value in `last_delta` locally but discard it. Update~~
~~`tests/spin_polarization.rs::test_ccmx_fe_free_magnetization_converges`~~
~~to assert `result.final_delta < 1e-2` — the specific pathology the~~
~~CCMX limit cycle produced was Δρ pinned at 0.254, not iter-count, so~~
~~this is the pathology-specific regression guard the review flagged.~~

~~**Acceptance criterion:** `ScfResult` has the new field; CCMX Fe test~~
~~asserts on it; at least one more convergence test (Si Γ-only or the~~
~~BROY regression) also asserts a reasonable `final_delta` upper bound.~~

**Landed:** `ScfResult.final_delta: f64` added to
`src/scf/mod.rs::ScfResult` and populated from `last_delta` in both
`scf::driver::run_scf_unpolarized` and `scf::driver_spin::run_scf_spin`.
Semantics: on success, the converged Δρ; on `max_iter` exhaustion the
driver returns `PwdftError::ConvergenceFailure { delta, .. }` so the
field only ever carries a converged value. Two regression guards
landed: CCMX Fe test asserts `final_delta < 1e-2` (pathology pin —
pre-CCMX Δρ was pinned at ≈0.254; post-CCMX observed 5.7e-4 at iter
14); BROY `test_broyden_vs_plain_scf_convergence` asserts both plain
and Broyden Si SCF runs produce `final_delta < 1e-4` (two orders below
the `conv_threshold=1e-6`, empirical ≤ 1e-7). See PR FDLT.

#### ~~VNMT — m-isolation defense-in-depth test for V_NL~~ (landed)

~~- **Role:** Core Engineer~~
~~- **Priority:** low, **Complexity:** small, **Risk:** low~~
~~- **Source:** VNLM code review, Step 2 of the VNLM proposal's validation~~
~~ plan that never landed.~~

~~Add an integration or unit test that pins a single-(l, m)-channel V_NL~~
~~matrix element: build a test `NonlocalPotential` with `D_ij = 0` except~~
~~for one `(l=2, m=+2)` entry, then assemble H at two specific G-vectors~~
~~and compare to a hand-computed value using the known real-Y_lm formula.~~
~~Any future edit to `real_sph_harmonics` or the shift-and-flip matmul~~
~~builder that wrongly normalizes a single m-channel (e.g. √2 off on~~
~~m≠0) would be caught — the addition-theorem test only pins sums over m.~~

~~**Acceptance criterion:** new test passes post-landing; manually~~
~~breaking the `(l=2,m=+2)` normalization (e.g. multiplying by √2)~~
~~makes the test fail with an intelligible error message.~~

Landed as `src/potential/nonlocal.rs::tests::test_single_channel_l2_m_isolation`
(multi-m pin with hand-computed per-m breakdown for three live channels —
m=0, +1, +2 — since Si ONCV PP's D_ij is diagonal in m by construction,
so a pure single-m-only test is architecturally unreachable without
synthesizing a non-physical projector). Manual-break sanity check
(multiply `Y_{2,+2}` by √2) triggered a 3.5e-7 residual vs 1e-10
tolerance, confirming detectability. See PR VNMT.

#### ITVF — Flip ITEV default to `Iterative` post-faer-upstream-fix

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

#### ~~RDOC — Clean up 15 pre-existing rustdoc warnings~~ (done)

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

#### ~~UPFV — Parse-time validation for UPF `angular_momentum`~~ (landed)

~~- **Role:** Researcher or Core Engineer~~
~~- **Priority:** low, **Complexity:** small, **Risk:** low~~
~~- **Source:** CAST code review, PR #56 nit 2.~~

~~`src/pseudopotential/upf/xml.rs::extract_beta_angular_momentum` returns~~
~~`am_str.trim().parse().ok()`, which happily accepts `-1`, `-2`, etc.~~
~~CAST added a runtime assertion in `NonlocalPotential::new` (real~~
~~belt-and-suspenders, not vacuous), but a parser-level~~
~~`PwdftError::InvalidInput("angular_momentum must be non-negative")`~~
~~would (a) fail fast at UPF load time instead of at projector-build~~
~~time and (b) give a clearer error message ("bad UPF file" vs "internal~~
~~assertion"). Add a regression test that synthesizes a malformed UPF~~
~~with `angular_momentum="-1"` and asserts the parse returns~~
~~`PwdftError::InvalidInput`.~~

~~**Acceptance criterion:** parse of a malformed UPF fails with a~~
~~domain-specific error string; `NonlocalPotential::new` runtime assert~~
~~remains as defense-in-depth.~~

Landed: `extract_beta_angular_momentum` now returns
`Result<i32, PwdftError>`; negative values map to
`PwdftError::InvalidInput`, missing / unparseable to `PwdftError::Parse`.
CAST's `NonlocalPotential::new` `assert!(proj.l >= 0)` stays as
defense-in-depth. Regression tests
`test_extract_beta_l_rejects_negative`,
`test_extract_beta_l_accepts_zero_and_positive`, and
`test_parse_rejects_negative_angular_momentum_in_upf` in
`src/pseudopotential/upf/convert.rs::tests` pin both the helper-level
rejection and the end-to-end UPF path (Si ONCV file with
`angular_momentum="0"` → `angular_momentum="-1"` substitution).

#### ~~FGRD — Explicit FFT-grid upper-bound check~~ (landed)

~~- **Role:** Core Engineer~~
~~- **Priority:** low, **Complexity:** trivial, **Risk:** low~~
~~- **Source:** CAST code review, PR #56 nit 3.~~

~~Several CAST `#[allow(clippy::cast_possible_truncation, reason = "...")]`~~
~~sites in `src/scf/grid.rs:80,106` and `src/symmetry/density/*.rs` cite~~
~~"FFT grid dims ≤ ~512 per axis in practice". In practice ecut + Miller~~
~~truncation + grid-size heuristics keep us well under 512, but the~~
~~bound is load-bearing on convention, not a checked invariant. Either~~
~~(a) add an explicit `debug_assert!(nx <= 1024, "FFT grid too large")`~~
~~at `FftGrid::new`, or (b) rewrite the reason strings to cite the `ecut`~~
~~bound that actually limits them — whichever makes the invariant~~
~~visible to a future grep.~~

~~**Acceptance criterion:** every CAST `#[allow]` that currently cites~~
~~"grid ≤ 512" either backs the bound with a `debug_assert!` or cites a~~
~~first-principles bound (ecut ≤ X Ry → grid ≤ Y).~~

Landed via option A: added `pub(crate) const MAX_FFT_DIM: usize = 1024`
and a runtime `assert!` at `src/scf/grid.rs::FftGrid::new`. The seven
`#[allow(cast_*)]` reason strings in `src/scf/grid.rs` and
`src/symmetry/density/{mod,real_space,g_space}.rs` now cite
"asserted <= MAX_FFT_DIM (1024) at `scf::grid::FftGrid::new`"
instead of "<= ~512 per axis in practice". See PR FGRD.

#### ~~G2ZT — Hoist bare `1e-12` `|G|=0` threshold into `consts`~~ (PR #153)

~~- **Role:** Code Reviewer~~
~~- **Priority:** trivial, **Complexity:** trivial, **Risk:** low~~
~~- **Source:** CFGN re-scope (PR #78) drive-by finding.~~

~~`src/pseudopotential/mod.rs:137` uses a bare `1e-12` literal for its~~
~~`|G|=0` check instead of `crate::consts::G2_ZERO_THRESHOLD` (or its~~
~~sqrt). One-line edit; makes the convention grep-discoverable.~~

~~**Acceptance criterion:** bare `1e-12` replaced with the existing~~
~~constant; `cargo test` bit-identical.~~

Landed as PR #153. The call site at `src/pseudopotential/mod.rs:159`
compares `g_norm` (i.e. `|G|`) directly against the threshold, not
`|G|²`, so reusing the existing `G2_ZERO_THRESHOLD` would have silently
rescaled the floor by 10⁶. Added a sibling `G_ZERO_THRESHOLD = 1e-12`
in `src/consts.rs` with a docstring flagging the distinction between
`|G|` and `|G|²` floors; the call site now reads
`if g_norm < G_ZERO_THRESHOLD`. Zero behavior change (bit-identical
Tier-1 pass at 330 tests, clippy unchanged at 18/24, rustdoc clean).

#### ~~DFLT — Document density-skip/normalization thresholds in `scf/density.rs`~~ (landed)

~~- **Role:** Code Reviewer~~
~~- **Priority:** trivial, **Complexity:** trivial, **Risk:** low~~
~~- **Source:** CFGN re-scope (PR #78) drive-by finding.~~

~~`src/scf/density.rs:61,92` has two `1e-15` literals for occupation-skip~~
~~and normalization-integrand safety. Semi-internal; a small~~
~~`const DENSITY_SKIP_THRESHOLD: f64 = 1e-15` (or two named constants at~~
~~module scope) would document intent for the next reader.~~

~~**Acceptance criterion:** both literals replaced by named constants~~
~~with a one-line docstring each explaining the physical/numerical~~
~~motivation.~~

Landed as two module-private constants in `src/scf/density.rs`:
`OCCUPATION_SKIP_THRESHOLD` (guards the forward FFT + |ψ(r)|² accumulation
from bands below eigensolver round-off) and `NORMALIZATION_INTEGRAL_FLOOR`
(guards the final ρ-rescale against a vanishing integral). Both sit at
`1e-15` — bit-identical SCF convergence preserved.

#### FLP3 — NLCC ρ_core(G) regression: parameterize over all 64 NLCC-active PPs

- **Role:** Code Reviewer (refactor) + Researcher (reference-value generation)
- **Priority:** medium (defensive — no known regression), **Complexity:** small, **Risk:** low
- **Source:** TRV2-F3 agent correction on PR #96 (landed 2026-04-18).

TRV2's Finding #3 originally estimated "~7 NLCC-active PPs" blind. The
actual count from `grep -l 'core_correction="T"' pseudopotentials/nc/lda/*.upf`
is **64**. PR #96 landed Cu + Mn regression coverage on top of the
existing Si + Fe pins, bringing the hand-coded total to 4 elements
(8 tests: {Si, Fe, Cu, Mn} × {G=0, first shell}). **60 NLCC-active
pseudopotentials remain unpinned.**

NCFX (PR #40) fixed a universal bug class in the NLCC path — the
unit conversion (e/Bohr³ → e/Å³) and the r²·4π radial weight — that
applies identically to all 64 NLCC pseudos. The regression risk for
the 60 unpinned elements is the same class the inline tests were
written to catch; blind coverage is precisely the situation NCFX was
created to close.

**Scope:**

1. Refactor the 8 hand-coded `test_{si,fe,cu,mn}_rho_core_of_g_{zero,first_shell}`
   tests at `src/pseudopotential/upf/convert.rs:305-554` into a single
   parameterized test driven by a static table of
   `[(element, cell_type, a_Å, ref_g0_e_per_ang3, ref_g_shell_1_e_per_ang3), ...]`
   rows. The Rust test iterates the table, loads the matching UPF,
   builds the appropriate cell (FCC for the Si/Cu family, BCC for the
   Fe/Mn family — same pattern PR #96 already uses), integrates
   `rho_core(G)` with the existing trapezoidal quadrature, and asserts
   each row at the Fe-matching 1e-4 e/Å³ tolerance.
2. Extend `scripts/validate/rho_core_g_reference.py` (currently emits
   rows for Si, Fe, Cu, Mn — 24 rows in `rho_core_g_reference.csv`) to
   iterate all 64 NLCC-active pseudos and emit reference rows for each.
   Cell convention: reuse the existing FCC/BCC auto-selection heuristic
   PR #96 established; pick a defensible lattice constant per element
   (experimental `a` from standard tables or any value that puts the
   first shell in a well-resolved G-range — the reference calc is
   self-consistent with whatever lattice the Rust test reads from the
   table).
3. The parameterized test should list all 64 element names in the
   table so a missing-row panic is loud if Python and Rust drift.

**Acceptance criterion:**

- Single parameterized `test_nlcc_rho_core_of_g_all_elements` (or
  similar) iterates all 64 NLCC-active PPs; each passes at ≤ 1e-4 e/Å³.
- `scripts/validate/rho_core_g_reference.csv` contains 64 × 2 = 128
  rows (G=0 + first shell per element).
- The 8 existing hand-coded tests are either deleted (subsumed by the
  parameterized path) or retained as the four most load-bearing rows
  (Si, Fe, Cu, Mn) with a comment pointing at the parameterized table
  for the rest.
- No change to `compute_core_density` semantics; all other NLCC tests
  (Si core-charge integral, NLCC SCF integration tests) remain
  bit-identical.

**Cost estimate:** ~0.5 CE-day for the Rust-side refactor and table
wiring; ~0.5 Researcher-day to run
`scripts/validate/rho_core_g_reference.py` end-to-end over all 64
elements once (output is checked-in CSV — one-time cost, then
deterministic).

**Why not a standalone proposal:** the work is bounded (< 1 CE-day),
mechanical, and closes a coverage gap from a just-landed PR rather than
introducing new physics. FLUP is the right home.

#### ~~MXB1 — Verify Eyert §3.3 vs §5 threshold constants~~ (struck 2026-04-18)

- **Role:** Researcher
- **Priority:** low, **Complexity:** small, **Risk:** low
- **Source:** MXBA code review, PR #57 flagged follow-up.

Resolved by citation-fix amendment in
`proposals/completed/MXBA-adaptive-mixing-beta.md` (Note, 2026-04-18).
The code implements the §3.3 residual-norm monitor (matching its own
docstrings); the §5 block + constants `(γ_up=1.2, γ_down=0.5, c_down=0.8,
c_up=1.0)` in the archived proposal were a draft residue, superseded
during implementation and never rewritten. Paper full text was not
accessible, so the exact equation number within §3.3 remains unverified
(noted in the amendment). The constants appear empirically tuned; a
sweep script is recommended as part of MXB2, not MXB1.

#### MXB2 — Re-diagnose MXBA Fe failure + retune

- **Role:** Core Engineer
- **Priority:** medium (blocks default-on), **Complexity:** small-medium, **Risk:** low
- **Source:** MXBA code review, PR #57 RCA correction.

The PR body says "flat residual damps β" — wrong. A flat ratio ≈ 1.0
is inside the hysteresis band `[0.5, 1.2]` and fires neither damp nor
restore. The `tests/mxba_adaptive_beta_fe.rs` trajectory shows β stays
at 0.3 for iters 1-9 (monitor silent) and only starts dropping iter 10+.
Real failure mode: the residual *oscillates past 1.2× often enough to
chain damps, while the 3-iter streak of <0.5× needed to restore is
unreachable once DIIS is starved*. The "DIIS warm-up window" the MXBA
PR flagged is therefore the wrong mitigation. Candidates to try instead:
(a) require a 2-iter *growth* streak before damping, (b) raise `β_min`
floor from `0.05·β_start` to e.g. `0.2·β_start`, (c) relax
`restore_threshold` from 0.5 to 0.8. Bench on Fe CCMX + C diamond at
30 Ry plain mixing.

**Acceptance criterion:** `tests/mxba_adaptive_beta_fe.rs` either
converges with adaptive=on or remains `#[ignore]`'d with an updated
reason explaining which tuning was tried and why it didn't work. If a
tuning makes Fe converge, flip `adaptive_beta` default to `true`.

#### ~~MXB3 — Direct `AdaptiveBeta::update` Fe-trajectory unit test~~ (landed 2026-04-18)

- **Role:** Code Reviewer (or Core Engineer)
- **Priority:** low, **Complexity:** trivial, **Risk:** low
- **Source:** MXBA code review, PR #57 companion-test recommendation.

~~Add a ~20-line direct unit test on `AdaptiveBeta::update` that feeds
the documented Fe-failure trajectory shape (initial flat plateau
followed by the 1.2×+ oscillation pattern) and asserts β floors at
`β_min` within the observed 80-iter envelope. Deterministic and
millisecond-cost — a cheap companion to the expensive `#[ignore]`'d
integration test in `tests/mxba_adaptive_beta_fe.rs`. If MXB2 changes
the failure mode, this test updates with it.~~

Landed as
`src/scf/mixing/mod.rs::adaptive_beta_tests::adaptive_beta_fe_failure_trajectory_floors_to_beta_min`.
Synthetic 80-iter trajectory: 10-iter flat plateau at 0.34 (ratio ≈ 1.0,
monitor silent) followed by 70 iters of `[1.3, 0.9]` oscillation on the
same base (ratios alternate 1.444 → damp, 0.692 → band; the 3-iter
streak of ratios < 0.5 that would restore β is architecturally
unreachable). β reaches β_min = 0.015 by iter ~27 and stays there
through iter 80. When MXB2 lands a fix, this test updates with it —
either asserts β recovers or gains an `#[ignore]` marker matching the
integration test.

#### ~~VNLT — Investigate non-reproducing VNLM `vnl_new` regression~~ (done 2026-04-19)

- **Role:** Performance Engineer
- **Priority:** medium (may retire VNLB), **Complexity:** small, **Risk:** low
- **Source:** PERF post-MXBA benchmark pass (PR #59), ANOM-3.
- **Resolution:** outcome (a) — the 78.5 ms `vnl_new_n725` in VNLM PR #49
  was a criterion outlier. Three clean runs on current main
  (Apple M3 Max, machine-locked, criterion `--measurement-time 6 s / 100 samples`,
  same config as PERF): 44.37 ms [44.28, 44.51] → 44.36 ms [44.30, 44.45]
  (p = 0.93 vs run 1) → 44.23 ms [44.19, 44.28]. All three land within
  ±0.3 % of each other and within 5 % of the pre-VNLM 42.33 ms baseline.
  Before-bench diff analysis showed `benches/scf_benchmarks.rs` is
  byte-identical across the VNLM merge and today's HEAD, and the only
  two commits to `src/potential/nonlocal.rs` since VNLM merged
  (`1291729` CAST, `6840237` RDOC) added zero-cost attributes and two
  cheap asserts — none can explain a ~35 ms wall-time swing. **VNLB
  struck** (see its entry). VNLM's archived proposal amended with a
  "Note 2026-04-19" paragraph.

VNLM (PR #49) reported a 1.5× `vnl_new_n725` regression (43.2 → 78.5 ms)
as the one-time cost of the GEMM lift, with "break-even at ~2 SCF iters"
amortization caveat. Today's PERF pass measured 43.4 ms — *matching the
pre-VNLM baseline*. Three possibilities:

1. The VNLM regression was real but something after VNLM inadvertently
   fixed it (MXBA / CAST / MODR rebases changed the code path or the
   SIMD vectorization pattern).
2. The VNLM bench was noisy and the 78.5 ms was an outlier.
3. The bench input or criterion config drifted between PR #49 and PR #59.

Clean re-bench with longer `--measurement-time`, compare to PR #49's
criterion baseline file if preserved in `target/criterion/`. If
confirmed that `vnl_new` is genuinely at 43 ms on current main, **VNLB
becomes moot** (no regression to recover) and VNLM becomes a pure
speedup with no amortization caveat — update VNLM's archived proposal
accordingly.

**Acceptance criterion:** either (a) confirm vnl_new ≈ 43 ms is stable
across 3 bench runs → strike VNLB from FLUP and amend VNLM notes, or
(b) reproduce the 78 ms regression → keep VNLB and add a note to this
entry explaining what toggled it.

#### EIGV — Investigate `faer_eigen_n259` +36% regression anomaly

- **Role:** Performance Engineer
- **Priority:** low (may be pure noise), **Complexity:** small, **Risk:** low
- **Source:** PERF post-MXBA benchmark pass (PR #59), ANOM-1.

PERF flagged `faer_eigen_n259` as +36% slower with an unusually wide
95% CI (±15%). No landing today claims the eigensolver at n=259.
Re-bench with longer `--measurement-time 15` (default is ~5s) and
`--sample-size 50` to tighten the CI. If the regression tightens,
bisect — most likely culprit is `Cargo.lock` drift in a faer-adjacent
dep that changed its internal SIMD tuning.

**Acceptance criterion:** either (a) regression disappears under longer
measurement → mark as noise, strike entry, or (b) regression confirmed
within ±5% CI → open a proposal to bisect the cause (likely a single
cargo update commit from this week).

#### EIGW — Investigate `faer_eigen_n725` unclaimed −16% win

- **Role:** Performance Engineer
- **Priority:** low (free win but needs understanding), **Complexity:** trivial, **Risk:** low
- **Source:** PERF post-MXBA benchmark pass (PR #59), ANOM-2.

PERF noted that `faer_eigen_n725` quietly improved by −16% with no
claimed source. Free wins are nice but unexplained ones are a smell —
may indicate a related regression is masked. Same re-bench protocol as
EIGV (longer measurement, tighter CI). If the win is real, figure out
which landing caused it (MODR? CAST? Some incidental inline change?) so
future bisections against this baseline have a reference. Combine
investigation with EIGV — same tool, same worktree, same afternoon.

**Acceptance criterion:** either (a) win disappears under longer
measurement → mark as noise, or (b) win confirmed → git bisect across
the week's landings to attribute and document.

#### ~~DWGT — Add `cargo doc` to the quality gate~~ (landed 2026-04-18, commit `5456c80`)

~~- **Role:** Technical Writer~~
~~- **Priority:** low, **Complexity:** trivial, **Risk:** low~~
~~- **Source:** RDOC agent suggestion, PR #61.~~

~~Today's RDOC cleanup cleared 17 pre-existing `cargo doc --no-deps`
warnings to 0. Without CI enforcement, regressions will accumulate
again. Add `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` to the "Code Quality"
section of `CLAUDE.md` alongside the existing two clippy invocations.
Also update each agent definition under `.claude/agents/*.md` whose
workflow mentions the quality gate.~~

~~**Acceptance criterion:** CLAUDE.md + agent defs updated; any new PR
touching docstrings that introduces a warning is blocked by the gate.~~

Landed: CLAUDE.md § Code Quality now lists `RUSTDOCFLAGS='-D warnings'
cargo doc --no-deps` alongside both clippy invocations (lines 47 and 52
of CLAUDE.md on today's main). Agent definitions updated:
`.claude/agents/core-engineer.md:53` chains the cargo-doc gate into the
quality-gate block, and `.claude/agents/technical-writer.md:64` documents
DWGT as mandatory and forbids `#[allow]` on rustdoc warnings. ERR2 P0
(PR #86) and subsequent merges have run clean against the gate. (DOCX
2026-04-18 corrected the historical `cargo doc --no-deps -- -D warnings`
recipe to the env-var form; current cargo rejects `-D warnings` passed
after `--`.)

#### ~~TYPE-AX — Decide on TYPE-A narrowing `expect` sites~~ (folded into TYPB, 2026-04-19)

Consolidated into `proposals/TYPB-narrow-int-audit.md` as Part C (the four
`i8` rotation-entry sites gain `reason = "..."` citing the crystallographic
bound) and Part A (the `src/basis.rs:65` Miller `i16` `expect` is auto-closed
by reverting the narrowing). The single integer-type cleanup PR now takes all
five sites in one sweep rather than splitting the decision across TYPB +
ERR2-P1. History of the original entry preserved below.

~~- **Role:** Core Engineer (ERR2 P1 owner)~~
~~- **Priority:** low, **Complexity:** trivial (decision + either `reason`
  comments or Result-returning refactor), **Risk:** low~~
~~- **Source:** ERR2 P0 post-landing signal (PR #86); five new
  `.expect(...)` sites introduced by TYPE-A (PR #80) when narrowing
  `i32 → i8` rotation entries and `i32 → i16` Miller indices.~~

~~The TYPE-A narrowing introduced five `expect` call sites guarded by
`try_from`:~~

- ~~`src/basis.rs:65` — `i16::try_from(n).expect(...)` on Miller indices.
  Comment (line 60) asserts "infallible under physically meaningful
  `ecut`".~~
- ~~`src/symmetry/operations.rs:71` — `i8::try_from(v).expect("SymmOp::from_flat: rotation entry out of i8 range")`.~~
- ~~`src/symmetry/operations.rs:120` — `i8::try_from(v).expect("SymmOp::inverse: adjugate entry out of i8 range")`.~~
- ~~`src/symmetry/operations.rs:151` — `i8::try_from(v).expect("SymmOp::compose: product entry out of i8 range")`.~~
- ~~`src/symmetry/detect.rs:185` — `i8::try_from(v).expect("symmetry::detect: rotation entry exceeds i8 range")`.~~

~~ERR2 P0 intentionally left all 15+ production `expect` sites as warnings
rather than fixing them; P1 is where each site gets decided. For TYPE-A's
five sites the bound is structural (crystallographic rotation entries are
in {-2..2}, Miller indices are bounded by `sqrt(ecut / HBAR2_OVER_2M)`),
so these are candidate "legitimate invariant, add `reason = "..."` and
move on" — **not** candidates for Result propagation. P1 should either
(a) add a `reason` comment citing the bound, matching the `BUG:` pattern
established by ERRH + FGRD, or (b) if treating them as user-reachable,
convert to `PwdftError::InvalidParam`. Option (a) is the sane default.~~

~~**Acceptance criterion:** each of the five sites either carries a
`reason = "..."` comment referencing the structural bound, or its
enclosing function returns `Result<_, PwdftError>`. No bare `expect`
remains in the narrowing path. Does not need its own proposal — fold
into ERR2 P1 when that starts.~~

### What this is NOT

- **Not an implementation plan.** Each entry needs to be promoted to
  its own proposal before coding starts (except the "file-it-now" ITEV
  upstream issue, which is a one-PR drive-by).
- **Not ordered by priority.** Ordering within this file is roughly
  "flagged first came first"; the EM picks a real order when activating.
- **Not a promise to land every entry.** Some may turn out to be not
  worth the round-trip cost — ITVF in particular is a tracker whose
  preconditions may never arrive (faer may deprecate the affected API
  path and we swap to LOBPCG, closing ITVF as obsolete). EIGV/EIGW may
  self-resolve under longer `--measurement-time` and just get struck as
  noise on the next bench pass.

### Proposal-file etiquette when activating

When the EM activates an entry:

1. Create `proposals/<ID>-<slug>.md` with proper frontmatter.
2. Copy the FLUP entry's body into the new proposal's **Motivation**
   and **Acceptance criterion** sections; expand into full Implementation
   plan + Verification sections.
3. Strike through the FLUP entry (prepend `~~` to each line) so the
   seeding history stays legible, but don't delete.
4. Add the new proposal's row to `proposals/INDEX.md`.
