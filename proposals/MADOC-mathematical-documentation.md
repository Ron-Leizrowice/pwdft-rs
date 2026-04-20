---
id: MADOC
title: Systematic mathematical-content docstring push for physics-critical modules
priority: medium
complexity: large
risk: low
depends_on: []
blocks: []
status: active
author: Researcher
date: 2026-04-18
---

## MADOC — Systematic mathematical-content docstring push

### Problem

Three physics bugs caught in the last week (NCFX — NLCC unit convention;
PCFX — τ phase direction in real-space density symmetrization; CCMX —
(ρ_total, m) vs (ρ↑, ρ↓) mixer basis) each had one root cause in common:
**a docstring that paraphrased the signature instead of stating the
equation, convention, and units**. In each case a math-complete docstring
would have failed code review on inspection alone — the bug would never
have landed.

> *"We need robust mathematical documentation in all docstrings to better
> sanity check and confirm the mathematics in the underlying code."*
> — user, 2026-04-18

Current docstring quality is bimodal. Some modules are excellent (post-NLCC
`scf/energy.rs`; post-PCFX `symmetry/density/g_space.rs`; post-VNLM
`potential/nonlocal.rs` module header). Most physics-critical `pub` items
are one-line signature paraphrases ("Get V_local(G) for a given G-vector
index.") — true, useless.

This proposal sizes a systematic push to raise every physics-critical
`pub` / `pub(crate)` docstring to a **math-complete** standard: a
physicist reading only the docstring can reproduce the computation on
paper without opening the body.

### Relationship to DLNT / DWGT (lint gates)

The EM has drafted FLUP/DWGT (add `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` to
the quality gate) and is considering a follow-up DLNT (enable
`missing_docs` + `missing_errors_doc` + `missing_panics_doc` as `warn`).
**MADOC is the content half; DLNT is the lint half**. They are
complementary:

- DLNT without MADOC forces every new `pub` to have *some* docstring,
  which the author will satisfy with "Returns the Hartree potential.".
  The lint is green; the bug-catching power is zero.
- MADOC without DLNT is a one-shot content upgrade. Drift will
  re-accumulate as new `pub` items land without docs at all.

**Recommended ordering: MADOC first, then DLNT.** Rationale:

1. DLNT surfaces *missing* docstrings, not *shallow* ones. The bulk of
   today's problem is shallow-but-present; DLNT does not catch it.
2. Running DLNT first produces a list of hundreds of warnings, each of
   which the fixer has local context for "what should the docstring say
   here?". Without MADOC's target-shape contract written down, each
   fixer satisfies the lint at the cheapest level (one-liner) and we
   lock in exactly the paraphrase pattern we are trying to eliminate.
3. Post-MADOC, DLNT becomes a drift guard: once the critical-path
   docstrings are math-complete, DLNT ensures new additions include at
   least a docstring, and the MADOC contract (plus code review) ensures
   that docstring meets the math-complete bar.

DWGT (quality-gate flip for zero-warning `cargo doc`) is orthogonal to
both and can land at any point.

### Scope — module triage

Sampled one `pub` or `pub(crate)` item per module. Each item is
classified by the target-shape criteria in § "Target docstring shape"
below:

- **E**xcellent — equation + citation + variables + units + invariants
- **A**dequate — equation *or* citation, not both; units usually present
- **P**araphrase — restates the signature; no equation; no citation
- **M**issing — no docstring at all

#### Tier 1 — wrong equation produces wrong physics silently

| Module | Pub items | E | A | P | M | Notes |
|--------|----------:|--:|--:|--:|--:|-------|
| `potential/xc.rs` | 7 | 2 | 4 | 1 | 0 | `lda_xc_energy:99` E; `lda_xc_spin_grid:233` A (formula but no citation); `pz_correlation_spin:316` A; `XcPoint:36` P |
| `potential/nonlocal.rs` | 2 | 1 | 1 | 0 | 0 | Module-header post-VNLM E; `add_to_hamiltonian` A (no explicit `H_{GG'} = Σ_α` form in the body doc) |
| `potential/local.rs` | 4 | 1 | 1 | 2 | 0 | `LocalPotential::new:25` E; `v_of_g:69`, `as_slice:74` P; `LocalPotential:18` P |
| `scf/energy.rs` | 17 | 11 | 4 | 2 | 0 | Post-NLCC gold standard; outliers: `band_energy:45`, `density_diff:163`, `assemble_v_eff:193` missing unit annotations |
| `scf/initial_density.rs` | 2 | 0 | 2 | 0 | 0 | Module doc has SAD formula but `generate_initial_density:55` doesn't state the normalization equation (`∫ρ(r)d³r = N_el`) |
| `scf/smearing.rs` | 4 | 1 | 3 | 0 | 0 | `fermi_dirac_01:122` has formula, no Mermin citation; `find_fermi_energy:57` has bisection invariant but no convergence tol doc |
| `ewald.rs` | 1 | 1 | 0 | 0 | 0 | `ewald_energy:34` post-EWSF is excellent — all four terms with explicit formulas |
| `symmetry/density/g_space.rs` | 1 | 1 | 0 | 0 | 0 | `symmetrize_density_g:155` gold standard, post-PCFX |
| `symmetry/density/real_space.rs` | 1 | 0 | 1 | 0 | 0 | Deprecated path; should get a `#[deprecated]` note pointing at g_space, not a new math block |
| `pseudopotential/upf/convert.rs` | 0 pub | — | — | — | — | `parse_body` is `pub(super)`; inline comments are excellent (NCFX post-audit); no action |
| `pseudopotential/mod.rs` | 8 | 3 | 2 | 3 | 0 | `PseudopotentialData` fields post-NCFX E; `BetaProjector:57`, `load:66`, `find_for_atom:79` P — no docstring-level statement of what a projector *is* |

**Tier 1 total: ~47 `pub`/`pub(crate)` items. Distribution: ~21 E / ~18 A / ~8 P / 0 M.**

#### Tier 2 — feeds Tier 1 (numerics layer)

| Module | Pub items | E | A | P | M | Notes |
|--------|----------:|--:|--:|--:|--:|-------|
| `fft.rs` | 6 | 1 | 3 | 2 | 0 | Module header states sign + normalization convention; `forward:63`, `inverse:86` P |
| `eigensolver/dense.rs` | ~4 | 0 | 2 | 2 | 0 | Hermitian eigendecomp; generalized eigenvalue problem `Hψ = εSψ` not written out (S=I here, worth noting why) |
| `eigensolver/iterative.rs` | ~6 | 0 | 3 | 3 | 0 | Lanczos/Davidson algorithms; ITEV proposal has the math but the docstrings don't |
| `basis.rs` | 5 | 1 | 2 | 2 | 0 | `BasisSet::new:28` A — has ecut formula but no Parseval/normalization note |
| `kpoints.rs` | ~4 | 1 | 2 | 1 | 0 | `monkhorst_pack` has the fractional formula but not the IBZ weight sum convention |
| `symmetry/operations.rs` | 22 | 3 | 10 | 9 | 0 | `SymmOp:15` E; the cardinal-direction constructors `c2x`, `c3_111`, ... P — could cross-reference ITA Vol. A table |
| `symmetry/kpoints.rs` | 1 | 1 | 0 | 0 | 0 | `reduce_kpoints:30` post-SYKP gold standard |
| `symmetry/detect.rs` | 1 | 0 | 1 | 0 | 0 | |

**Tier 2 total: ~49 `pub`/`pub(crate)` items. Distribution: ~7 E / ~23 A / ~19 P / 0 M.**

#### Tier 3 — infrastructure (out of scope for MADOC)

`settings.rs`, `context.rs`, `scf/grid.rs`, `scf/report.rs`, `main.rs`,
`error.rs`, GPU scaffolding. Math content is minimal; a docstring push
here produces busy-work without bug-catching value. Defer.

#### Summary

**Top-tier (Tier 1) scope: ~47 items, ~8 paraphrase + ~18 adequate to upgrade
(roughly 26 items need real work; the 21 already-excellent ones get
spot-checks, not rewrites).** Tier 2 adds ~42 more items with paraphrase
or adequate status. Nothing is missing — every sampled item has at least
one sentence.

### Target docstring shape

A physicist reads the docstring and can reproduce the computation on
paper. Concretely, every **public** physics-function docstring must have:

1. **Equation in Unicode + ASCII math** (rustdoc doesn't render MathJax
   by default — treat `///` as plaintext plus Unicode `ρ Σ ∫ ∂ ε_xc`).
   Example: `` `E_xc = ∫ ε_xc(ρ↑, ρ↓) · (ρ↑ + ρ↓) d³r` `` rather than
   "compute the exchange-correlation energy".

2. **Citation** — paper + equation number *or* textbook + page.
   Prefer published-literature citations (e.g. *"Perdew & Zunger, *Phys.
   Rev. B* **23**, 5048 (1981), Eq. (C1)"*). **Do NOT add QE source-line
   cross-references in public rustdoc** (per DCLN — the crate's public
   documentation describes the physics directly rather than positioning
   itself as a QE derivative). QE cross-references belong in `//` dev
   comments or test-module docstrings only.

3. **Variable definitions** — every symbol in the equation gets a
   one-line gloss including units. *"`ρ` in e/Å³; `G` in 1/Å; `Ω` cell
   volume in Å³"*.

4. **Invariants / preconditions** — what the caller must guarantee.
   *"Caller must have FFT-normalized `rho_g` so that `rho_g[0].re`
   equals `n_electrons / Ω` on the shifted-MP convention"*.

5. **Returned units / shape** — mandatory for functions returning
   `Vec<f64>`, `Mat<Complex64>`, or other shape-opaque types.

#### Before / after example

**Before** (`src/potential/local.rs:67`, current):

```rust
/// Get V_local(G) for a given G-vector index.
#[must_use]
pub fn v_of_g(&self, ig: usize) -> Complex64 { ... }
```

**After** (target shape):

```rust
/// Local pseudopotential in reciprocal space at the `ig`-th basis
/// G-vector.
///
/// Returns
///   V_local(G) = Σ_α S_α(G) · v_local^(α)(|G|)
/// in eV, where
/// - `S_α(G) = exp(−iG·τ_α)` is the structure factor for atom α at
///   Cartesian position τ_α (Å);
/// - `v_local^(α)(|G|)` is the spherical Bessel transform
///   `(4π/Ω) ∫₀^∞ [r² v_local^(α)(r) − (−Z_α e²/|G|²) cancelled at G=0]
///   j₀(|G|r) dr` (Å³·eV) stored by
///   `PseudopotentialData::v_local_of_g`;
/// - the G=0 component is **zeroed at construction** and added back
///   as `V_local(G=0) · N_el` in the total-energy accounting (see
///   [`crate::scf::energy::with_g0_shift`] and proposal NCFX).
///
/// Reference: Kleinman & Bylander, *Phys. Rev. Lett.* **48**, 1425 (1982)
/// for the separable form; the spherical Bessel transform conventions
/// follow Martin, *Electronic Structure*, §11.4.
///
/// Panics if `ig >= n_basis` (plain slice-index OOB).
```

Roughly 4× longer. Every symbol has a unit; every convention has a
cross-reference; the "why is G=0 zero?" question is answered before it's
asked. This is the bar.

### Phased rollout

Full coverage of Tier 1 alone is ~47 items with ~26 needing rewrites.
That's too much for one PR. Six phases, each ≤ one afternoon of
Core-Engineer-or-Researcher time, each a separate PR:

1. **MADOC-A — `scf/energy.rs` + `scf/driver.rs` + `scf/driver_spin.rs`.** The total-energy
   path. Post-NLCC mostly excellent; this phase closes the 6 adequate
   items and the 2 paraphrase items. Highest familiarity (fresh from
   NCFX); quickest win. **First — warm-up + template.**
2. **MADOC-B — `potential/xc.rs` + `potential/nonlocal.rs` + `potential/local.rs`.** The
   physics kernel. ~13 items to touch. Spin LSDA interpolation formulas
   need explicit ζ(1±...) derivations. **Second — highest bug-catching
   ROI** (three of the last four physics bugs lived here).
3. **MADOC-C — `pseudopotential/upf/convert.rs` + `pseudopotential/mod.rs`.** The UPF
   unit-conversion sites. Post-NCFX the `PseudopotentialData` fields are
   already documented; `BetaProjector`, `load`, `find_for_atom` need
   upgrades. **Third** — same unit-convention territory as NCFX.
4. **MADOC-D — `symmetry/density/*.rs` + `symmetry/operations.rs` + `symmetry/kpoints.rs`.** PCFX
   territory. `g_space.rs` and `kpoints.rs` are already gold; this phase
   upgrades the `operations.rs` constructors (9 paraphrase items) and
   deprecates the real-space density path.
5. **MADOC-E — `ewald.rs`, `fft.rs`, `basis.rs`, `scf/smearing.rs`, `scf/initial_density.rs`.** Remaining
   Tier 1 + highest-math Tier 2. `ewald.rs` is already gold; `fft.rs`
   convention is the one every caller needs to know (sign of the forward
   FFT, normalization).
6. **MADOC-F — `eigensolver/*` + `scf/mixing/*`.** Numerics layer. Math
   is linear algebra, citations are straightforward (Saad, Eyert,
   Anderson original). Last because the bug-catching ROI is lowest —
   the algorithms are well-known and well-tested.

**Ordering justification**: A→B→C front-loads the modules where bugs
have actually been caught in the last week. D consolidates the PCFX/SYKP
work. E/F are maintenance-grade — fewer surprises, fewer citations to
hunt down. Each phase stands alone; no phase blocks another; EM can
parallelize across agents if capacity allows.

### Anti-scope

Explicitly **not** in MADOC:

- **No `/// # Examples` blocks.** Doctests are code; broken doctests break
  the gate. Defer to a dedicated DTST proposal if the need arises.
- **No prose-style rewrites.** Only math content. Do not fix comma
  splices, re-order sentences for "WHY vs WHAT", or touch non-physics
  prose. Keep the PR diffs reviewable.
- **No private functions.** Only `pub` and `pub(crate)`. Module-private
  helpers stay whatever they are; private implementation detail is
  allowed to be cryptic.
- **No GPU/WGSL shaders.** WGSL has its own documentation norms; a
  follow-up proposal (MADG?) can cover them.
- **No test docstrings.** Tests can stay cryptic — the test body is the
  documentation. A test that passes is self-validating.
- **No Tier 3 infrastructure.** `settings.rs`, `scf/context.rs`,
  `main.rs`, error types. Useful documentation, but not math.

### Acceptance criteria

Per phase:

- Every `pub` / `pub(crate)` item in the phase's module set has a
  docstring meeting the target shape (§ 2).
- At least one spot-check from § "Module triage" (the items quoted
  there) is verifiably upgraded — the PR description cites the
  before/after.
- `cargo doc --no-deps` emits zero warnings (compatible with DWGT when
  it lands).
- `cargo test` passes. (Mechanical rewrites of docstrings shouldn't
  break tests, but we verify.)
- `cargo clippy -q --all-targets` and `cargo clippy -q --all-targets
  --features gpu` both green.

Per whole proposal:

- Tier 1 paraphrase count drops from ~8 to 0.
- Tier 2 paraphrase count drops from ~19 to ≤ 5 (numerics-layer
  internals allowed to stay at adequate if the referenced algorithm
  paper is cited).
- Project-wide rustdoc warning count stays at 0 (DWGT compatibility).
- DLNT can be flipped to `warn` without producing a warning cascade
  on Tier 1 or Tier 2 modules.

### Risk

Low. Pure content changes; no behavior modified. The only way MADOC
*introduces* a bug is if a hand-written equation in a docstring contains
a typo (`ρ^(1/3)` written as `ρ^(1/2)`) that a future reader trusts over
the code. Mitigation: each phase PR is reviewed by Researcher against
the cited reference, same as a physics-code PR.

### Out of scope — flagged for follow-up

None. This proposal is self-contained.
