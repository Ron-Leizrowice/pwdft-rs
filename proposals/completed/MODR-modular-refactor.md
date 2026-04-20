---
id: MODR
title: Modular refactor audit — split god-modules into orthogonal sub-modules
priority: medium
complexity: large
risk: medium
depends_on: []
blocks: []
status: proposed
---

## MODR — Modular refactor audit

### Motivation

Two drivers.

**1. Conflict surface.** With 5–6 concurrent agents, most active proposals
land in the same handful of files. Counting file-touches in the last six
weeks (`git log --since="6 weeks ago" --name-only -- src/ tests/`):

| Rank | File                          | Touches |
|------|-------------------------------|---------|
| 1    | `src/scf/mod.rs`              | 19      |
| 2    | `src/pseudopotential/mod.rs`  | 21 (mostly doc) |
| 3    | `src/scf/mixing.rs`           | 5       |
| 4    | `src/settings.rs`             | 6       |
| 5    | `tests/gpu_consistency.rs`    | 15      |
| 6    | `src/scf/initial_density.rs`  | 4       |
| 7    | `src/pseudopotential/upf.rs`  | 4       |
| 8    | `src/symmetry/density.rs`     | 2       |

Cross-referencing the active proposal bodies against these files:

- **`src/scf/mixing.rs`** — touched by PRPL (new variant + split `mix()`),
  MXBA (adaptive β hooks at L129 + L278), CCMX (if it grows a coupled-channel
  variant). Three proposals, one file, ~500 LOC.
- **`src/scf/mod.rs`** — touched by CCMX (`run_scf_spin` mixer wiring at
  ~L498 and ~L771), ITEV (eigensolver call sites at L299–307), WFRX, and
  every energy-accounting or logging change. Two of the in-flight branches
  already need a manual rebase in this file.
- **`src/symmetry/density.rs`** — PCFX will add a G-space path alongside the
  existing real-space one, doubling the file's scope.
- **`src/pseudopotential/upf.rs`** — next unit-conversion bug or UPF v1
  support lands here; NLCC and VERF both rewrote adjacent blocks.

**2. Code quality.** This is a research codebase — correctness over cleverness.
"High quality" for MODR means:

- single-responsibility modules,
- thin public API per module (prefer `pub(crate)`; module is at most 3–5
  re-exported symbols),
- explicit data flow (no globals, no hidden mutation),
- no speculative abstractions (traits and generics require ≥ 2 concrete
  present-day callers),
- tests co-located with the code they pin (unit tests next to the function,
  integration tests in `tests/`).

`src/scf/mod.rs` at 1073 lines combines driver, hot loop, energy-assembly
dispatch, and per-iteration logging — every change is a merge hazard. This
is the cheapest class of fix: move code around to reduce conflict surface
without touching physics.

### Current module map

LOC from `wc -l`; ★ marks god-modules. Production LOC (excluding
`#[cfg(test)] mod tests`) is shown in parens for the stars.

```text
src/
├── lib.rs                    19   crate root, re-exports
├── main.rs                  153   CLI + two calculation modes inlined
├── consts.rs                 42
├── error.rs                  27
├── atoms.rs                  71
├── crystal.rs               132
├── basis.rs                 179
├── kpoints.rs               226
├── hamiltonian.rs           121
├── fft.rs                   216
├── numerics.rs              153
├── ewald.rs                 271
├── bandstructure.rs         161
├── settings.rs              964 (484 prod)   ★ all knobs + YAML parse + tests
├── eigensolver/
│   ├── mod.rs                 1
│   └── dense.rs             153
├── gpu/
│   ├── mod.rs               629 (492 prod)   ★ device + 3 pipelines + pool + shaders
│   └── shaders/*.wgsl
├── pseudopotential/
│   ├── mod.rs               305   PP data + `v_local_of_g` + tests
│   └── upf.rs               300   ★ parse + unit convert + xml helpers + tests
├── potential/
│   ├── mod.rs                 4
│   ├── hartree.rs           100
│   ├── local.rs             177
│   ├── nonlocal.rs          451
│   └── xc.rs                525 (396 prod)   LDA + spin + parallel threshold
├── scf/
│   ├── mod.rs              1073 (967 prod)   ★★ driver + hot loop + spin variant + logging + tests
│   ├── mixing.rs            971 (503 prod)   ★ Anderson + Broyden + Kerker + solver + Mixer enum
│   ├── smearing.rs          431
│   ├── density.rs           100
│   ├── energy.rs            308
│   ├── initial_density.rs   305
│   ├── context.rs           146
│   ├── potentials.rs        166
│   └── grid.rs               90
└── symmetry/
    ├── mod.rs               247   SymmetryInfo + identity_only + is_trivial
    ├── operations.rs        439
    ├── detect.rs            347
    ├── kpoints.rs           325
    └── density.rs           287   ★ (small but doubling scope under PCFX)
```

### Proposed module map

Every new sub-module has a one-line purpose. "API:" lists the intended public
surface (3–5 symbols).

```text
src/
├── scf/
│   ├── mod.rs                 ~150   Re-exports + ScfParams + ScfResult + run_scf() dispatcher
│   │                          API: run_scf, ScfParams, ScfResult, EnergyComponents
│   ├── driver.rs              ~350   Non-spin SCF hot loop (moved from mod.rs)
│   │                          API: pub(crate) run_scf_unpolarized
│   ├── driver_spin.rs         ~400   Spin-polarized SCF hot loop (moved from run_scf_spin)
│   │                          API: pub(crate) run_scf_spin
│   ├── report.rs              ~120   Convergence logging + per-iteration progress + final summary
│   │                          API: pub(crate) IterationReport, log_iteration, log_converged
│   ├── mixing/
│   │   ├── mod.rs             ~80    MixingMode enum + Mixer dispatcher enum
│   │   │                      API: Mixer, MixingMode
│   │   ├── anderson.rs        ~200   AndersonMixer (moved as-is)
│   │   ├── broyden.rs         ~200   BroydenMixer (moved as-is)
│   │   ├── kerker.rs          ~60    precondition_residual + auto_q_tf_squared
│   │   │                      API: pub(crate) apply_kerker, pub(crate) auto_q_tf_sq
│   │   └── linalg.rs          ~80    solve_linear_system (module-local; tests in file)
│   │                          (no pub exports; pub(super) solve_linear_system)
│   ├── energy.rs              308    unchanged (already well-factored)
│   ├── density.rs             100    unchanged
│   ├── initial_density.rs     305    unchanged
│   ├── context.rs             146    unchanged
│   ├── potentials.rs          166    unchanged
│   ├── grid.rs                 90    unchanged
│   └── smearing.rs            431    unchanged
├── symmetry/
│   └── density/
│       ├── mod.rs             ~30    Public API facade + grid_compatibility helpers
│       │                      API: symmetrize_density, check_grid_compatibility, compatible_grid_dims
│       ├── real_space.rs      ~180   Current nint-based real-space symmetrizer
│       │                      API: pub(super) symmetrize_real_space
│       └── g_space.rs         ~180   (PCFX lands here) phase-factor G-space symmetrizer
│                              API: pub(super) symmetrize_g_space
├── pseudopotential/
│   ├── mod.rs                 305    unchanged (PseudopotentialData + v_local_of_g)
│   └── upf/
│       ├── mod.rs             ~40    Public `parse(&str)` entry point
│       │                      API: parse
│       ├── xml.rs             ~80    extract_attr, extract_data_block, extract_beta_angular_momentum
│       │                      (all pub(super); no leak outside upf/)
│       └── convert.rs         ~200   Unit conversion + rho_atom/NLCC blocks
│                              (pub(super) parse_body called by mod.rs)
├── settings.rs                964    unchanged (already sectioned; tests are half the LOC)
├── gpu/
│   └── mod.rs                 629    unchanged (deferred; see "What this is NOT")
└── main/                             (split main.rs)
    ├── main.rs                ~50    CLI parse + dispatch
    ├── modes/
    │   ├── scf.rs             ~80    MonkhorstPack → run_scf + pretty-print
    │   └── bands.rs           ~40    BandPath → compute_band_structure + TSV write
```

#### Consequence for in-flight proposals

- **PRPL** → `scf/mixing/anderson.rs` + one enum variant in `mixing/mod.rs`.
- **MXBA** → `scf/mixing/anderson.rs` + `scf/mixing/broyden.rs` (adaptive-β
  hook per mixer; different files from PRPL's edit).
- **CCMX** → `scf/driver_spin.rs` only (mixer wiring; the mixing module is
  untouched).
- **ITEV** → `scf/driver.rs` + `scf/driver_spin.rs` (two one-line call-site
  swaps). Only overlap with CCMX is the spin driver file, and the diffs
  don't touch the same lines.
- **PCFX** → brand-new `symmetry/density/g_space.rs`, zero overlap.

After Phases A+B+C, no two in-flight proposals share a file.

### Phased rollout

Ordered so earlier phases unblock the most concurrent work. Each phase is a
separate PR with its own implementation proposal (MODR-A, MODR-B, MODR-C,
MODR-D).

#### ✅ Phase A (PR #46) — Split `scf/mixing.rs` (unblocks PRPL + MXBA)

**Why first:** three active proposals target this file.

Move `AndersonMixer` → `mixing/anderson.rs`; `BroydenMixer` →
`mixing/broyden.rs`; `precondition_residual` + `auto_q_tf_squared` →
`mixing/kerker.rs`; `solve_linear_system` → `mixing/linalg.rs`.
`MixingMode` enum and the dispatcher `Mixer` stay in `mixing/mod.rs`.
Tests travel with their production code. **No behavior change.**
`cargo test` is the validator.

Size: ~2 h. Risk: low (pure move).

#### ✅ Phase B (PR #50) — Split `scf/mod.rs` (unblocks CCMX + ITEV + future energy work)

Moved `run_scf` → `scf/driver.rs::run_scf_unpolarized` (pub(crate));
`run_scf_spin` → `scf/driver_spin.rs::run_scf_spin` (pub(crate));
per-iteration progress and final-summary logging → `scf/report.rs`
(`IterationReport` + `log_iteration` + `log_convergence_summary` +
`log_entropy` + `log_components`, all pub(super)). `scf/mod.rs` now
carries only `ScfParams`, `ScfResult`, `run_scf` as a thin dispatcher,
and re-exports `EnergyComponents`. `EnergyComponents` itself moved
into `scf/energy.rs` along with the helper-pinning unit tests
(`test_real_to_g_space_*`, `test_assemble_v_eff_adds_correctly`,
`test_hartree_on_fft_grid_g0_zero`, `test_density_diff_*`) that
travel with those helpers (which already lived in `scf::energy`).

`mod driver`, `mod driver_spin`, `mod report` are all private — tighter
than Phase A's `pub mod mixing`, since no external caller needs them.
The shared helpers `diagonalize_dispatch`, `compute_occupations`,
`scf_progress_bar` live in `scf::driver` as `pub(super)` and are
imported sibling-to-sibling by `driver_spin`.

Size: ~3 h. Risk: low–medium (the hot loop is long but the split
follows a clean seam — spin vs non-spin vs logging).

#### ✅ Phase C (PR #48) — Split `symmetry/density.rs` (prepares PCFX)

Turn `symmetry/density.rs` into a folder. `real_space.rs` is the current
implementation verbatim. `g_space.rs` is a stub file (empty module) that
PCFX fills in. `mod.rs` exposes `symmetrize_density` and picks the
back-end (real-space today; PCFX will add a runtime/compile-time switch).
`check_grid_compatibility` and `compatible_grid_dims` stay in `mod.rs`
or move into `real_space.rs` — both become pub(super).

Size: ~1 h. Risk: low (pure move + tiny facade).

#### ✅ Phase D (PR #47) — Split `pseudopotential/upf.rs` (prepares future PP formats)

Turn `upf.rs` into a folder. `xml.rs` holds the three text helpers
(`extract_attr`, `extract_data_block`, `extract_beta_angular_momentum`).
`convert.rs` holds the unit-conversion body (Ry→eV, Bohr→Å, NLCC,
rho_atom). `mod.rs` is a 40-line `parse(content)` front door. All helpers
become `pub(super)`; nothing new leaks out.

Size: ~1 h. Risk: low. Future benefit: PSP8 support (deferred) drops in
as a sibling folder `pseudopotential/psp8/`.

**Not included in any phase:**

- `main.rs` split (153 LOC, low touch-rate, not worth the churn yet)
- `gpu/mod.rs` (rewriting behind CUCL proposal; splitting now would
  fight that rewrite)
- `settings.rs` (964 total but 484 prod; already well-sectioned with
  one struct per concern; splitting would add navigation cost without
  reducing conflict, since two agents rarely edit the same sub-struct)
- `potential/xc.rs`, `potential/nonlocal.rs`, `scf/smearing.rs` — not
  god-modules; cohesion is high and conflict rate is low.

### Risk + validation

The whole refactor is semantics-preserving. The validator is `cargo test`
plus both clippy invocations. Concerns per phase:

- **Phase A (mixing):** test coverage is strong (17 `#[test]` functions
  in `mixing.rs`, including the round-trip `test_broyden_vs_plain_scf_convergence`
  that MXBA is built against). Moves are safe.
- **Phase B (scf driver):** the driver itself is pinned by every integration
  test in `tests/` — `scf_convergence`, `parallel_consistency`,
  `spin_polarization`, `gpu_consistency`. That's the strongest possible
  pin: if the hot loop breaks, every one of these tests fails. The
  current `scf/mod.rs` unit-test block covers only three small helpers
  (`real_to_g_space`, `assemble_v_eff`, `density_diff`); those travel
  with the helpers' new homes if any.
- **Phase C (symmetry/density):** 6 unit tests in `symmetry/density.rs`
  plus `SYKP`-era integration coverage. Safe as a pure move.
- **Phase D (upf):** 10 tests in `pseudopotential/mod.rs` +
  `pseudopotential/upf.rs` parse Si/Fe/C PPs end-to-end. The three XML
  helpers are exercised by every parse path. Safe.

**Weak-coverage areas that MODR does NOT touch:** `main.rs` (integration
only, no unit tests) and `gpu/mod.rs` (dependent on GPU-enabled CI).

### What this is NOT

- **Not a trait-introduction refactor.** A `Mixer` trait object was
  considered and rejected: dispatch is once per SCF call, the current
  `enum Mixer` is zero-cost, and a trait object obscures data flow
  researchers need to read. Same argument rules out a generic
  parameter on the driver.
- **Not a settings split.** `settings.rs` is 964 LOC but 484 is tests,
  and it already has one struct per concern. Splitting forces two-file
  edits for a single knob.
- **Not a GPU-module split.** `gpu/mod.rs` has its own rewrite proposal
  (CUCL); splitting here would guarantee a conflict with CUCL.
- **Not a renaming pass.** Only location changes (same symbol name)
  so `git log --follow` stays clean. No type or function renames.
- **Not a public-API expansion.** All new sub-modules are `pub(crate)`
  or `pub(super)`. Nothing private becomes public.

### Open questions

1. **EnergyComponents home.** Currently in `scf/mod.rs`, produced only
   at the end of the drivers. Proposed: move into `scf/energy.rs`,
   re-export from `scf/mod.rs`.
2. **`scf/report.rs` vs inline logging.** Phase B extracts per-iteration
   `log!()` calls into a helper. Extraction wins because adding a new
   metric would touch one file instead of two.
3. **Phase C vs PCFX ordering.** Recommend landing MODR-C (1 h,
   zero-behavior-change) first so PCFX drops into a ready folder.

**Flagged for follow-up (not MODR scope):**

- `src/scf/driver.rs` (e_total/e_harris assembly) and
  `src/scf/driver_spin.rs` (same pair) — the G0-shift expression
  `+ ctx.v_local_g0 * ctx.n_electrons` is duplicated four times around
  `total_energy` / `harris_foulkes_energy`. A `with_g0_shift()` helper
  in `scf/energy.rs` would DRY it. **Core Engineer.** (Line numbers
  shifted after MODR-B; search for `v_local_g0 * ctx.n_electrons`.)
- `src/scf/mixing.rs:446–502` — `solve_linear_system` is a hand-rolled
  Gauss-elimination in a codebase that already links `faer`.
  **Performance Engineer** (swap to faer's LU when MXBA lands).
- `src/input.rs` stale reference — INDEX.md notes HD5I still points at
  the deleted path. **Technical Writer** / whoever refreshes HD5I.
