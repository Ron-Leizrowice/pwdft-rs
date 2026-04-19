---
id: MOAD
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# MOAD: Module-orientation docstrings (`//!` headers)

## Problem

Fourteen files in `src/` have no module-level (`//!`) docstring. The
omissions include the crate root and several large modules:

| File                         | LOC | What's in it (no docstring to say so) |
|------------------------------|-----|----------------------------------------|
| `src/lib.rs`                 |  ~60 | **Crate root** — empty cargo doc landing page |
| `src/atoms.rs`               | 200+ | Z=1–92 element table |
| `src/bandstructure.rs`       | 400+ | Band structure assembly + TSV writer |
| `src/basis.rs`               | 300+ | Plane-wave basis (G-vectors up to ecut) |
| `src/consts.rs`              | 100+ | Physical constants (HA_TO_EV, BOHR_TO_ANG, …) |
| `src/crystal.rs`             | 400+ | Lattice + atoms + Crystal aggregate |
| `src/error.rs`               | 100+ | `PwdftError`, `Result` |
| `src/hamiltonian.rs`         | 300+ | `build_hamiltonian` (kinetic + V_eff + V_NL) |
| `src/kpoints.rs`             | 600+ | Monkhorst-Pack + high-symmetry path |
| `src/scf/mod.rs`             | 400+ | `ScfParams`, `ScfResult`, `run_scf` dispatcher |
| `src/potential/mod.rs`       |   ~30 | XC + local + nonlocal facade |
| `src/eigensolver/mod.rs`     |   ~30 | Dense + iterative facade, `EigensolverKind` |
| `src/eigensolver/dense.rs`   | 500+ | `diagonalize_lowest`, `diagonalize_subspace`, WFRX subspace path |
| `src/pseudopotential/mod.rs` | 300+ | `PseudopotentialData`, `v_local_of_g`, `load` |

The contrast: `src/{ewald,fft,numerics,settings}.rs`,
`src/symmetry/mod.rs`, `src/symmetry/density/mod.rs`,
`src/scf/mixing/mod.rs`, and `src/gpu/mod.rs` already have well-written
`//!` headers covering "what's in here, what calls it, where to read
next." Those headers exist because someone needed them when the module
was reworked. The rest of the tree is a coverage gap — not a project
philosophy.

Concrete cost:

1. **`cargo doc --no-deps` lands on an empty page.** With `src/lib.rs`
   carrying no `//!`, `target/doc/pwdft_rs/index.html` shows the module
   tree but no orientation paragraph. Newcomers (and the user, when
   spot-checking) lose the chance to see "this is a plane-wave DFT
   solver; start at `scf::run_scf`" without opening source.
2. **Duplicate context across CLAUDE.md and code.** CLAUDE.md
   § Architecture lists the module groups; the modules themselves
   silently disagree (or rather, say nothing). When CLAUDE.md drifts
   (DOCX caught one such drift; UNTS another), there is no in-source
   anchor to disambiguate.
3. **Existing rustdoc is partial.** Function-level `///` docs are
   plentiful in these files — but a reader landing on
   `target/doc/pwdft_rs/scf/index.html` sees a list of items with no
   prose tying them together. This is exactly the orientation tier
   that MADOC (mathematical docs) and MAUD (correctness audit)
   explicitly do *not* cover.

This is not a "rewrite all the docs" proposal — MADOC owns mathematical
documentation, MAUD owns correctness review. MOAD is just the missing
1-paragraph headers on the 14 files above. Each header should answer:

- **What** is in this module (one sentence)?
- **Who calls it** from elsewhere in the crate (one sentence)?
- **Where to read next** if you need depth (link to the most-relevant
  sibling, or to a CLAUDE.md section)?

That's it. No deep prose, no math, no long examples. Aim for 4-8 lines
per header.

## Implementation

### Step 1 — Crate root (highest leverage)

Write `src/lib.rs`'s `//!` header first; everything else is downstream
of getting the cargo-doc landing page right. Suggested template:

```rust
//! Plane-wave Density Functional Theory (PWDFT) solver.
//!
//! This crate implements a self-consistent Kohn–Sham SCF loop on a
//! plane-wave basis with norm-conserving pseudopotentials. Targets LDA
//! today; PBE/hybrid functionals are tracked under the GGAP and HYBR
//! proposals. Optional Apple Metal GPU acceleration via the `gpu`
//! feature.
//!
//! Entry points:
//! - [`scf::run_scf`] — full self-consistent calculation.
//! - [`bandstructure::compute_band_structure`] — non-self-consistent
//!   eigenvalues along a k-path.
//!
//! Module groups (see `CLAUDE.md` § Architecture for the full map):
//! - **Crystal & basis:** [`crystal`], [`basis`], [`kpoints`], [`atoms`].
//! - **Pseudopotentials:** [`pseudopotential`].
//! - **Potentials & Hamiltonian:** [`potential`], [`hamiltonian`].
//! - **SCF loop:** [`scf`] (driver, mixing, energy, density, smearing).
//! - **Numerics:** [`fft`], [`eigensolver`], [`numerics`], [`ewald`].
//! - **Symmetry:** [`symmetry`].
//! - **GPU:** [`gpu`] (behind the `gpu` feature flag).
//!
//! Internal units are eV for energies, Å for lengths. Conversions live
//! in [`consts`]; pseudopotential parsing is the only Ry/Bohr boundary
//! (see [`pseudopotential::upf::convert`]).

#![cfg_attr(...)]  // existing crate attributes if any
```

### Step 2 — Module-mod.rs files (next-highest leverage)

`src/scf/mod.rs`, `src/potential/mod.rs`, `src/eigensolver/mod.rs`,
`src/pseudopotential/mod.rs` are the entry points users hit first when
clicking through the cargo-doc tree. Each gets a 4–8 line `//!`
covering: what's in the folder, what `pub` items the user typically
reaches for, what private siblings exist for context. Mirror the style
of the existing `src/symmetry/mod.rs` header.

### Step 3 — Single-file modules

The remaining 9 files (`atoms.rs`, `bandstructure.rs`, `basis.rs`,
`consts.rs`, `crystal.rs`, `error.rs`, `hamiltonian.rs`, `kpoints.rs`,
`eigensolver/dense.rs`) each get a 4–6 line `//!` header. Use the
function-level `///` docs that already exist as raw material — most of
the orientation prose is already written; it just lives at the wrong
granularity.

### Step 4 — Verify rustdoc cleanliness

Run `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps`. The DOCX gate
catches broken intra-doc links — every `[`fn`]` and `[`mod`]` reference
in the new headers must resolve. Spot-check the rendered HTML for the
crate root and `scf/`, `potential/`, `eigensolver/`, `pseudopotential/`
landing pages.

## Verification

```bash
.claude/bin/machine-lock acquire "Core Engineer" "MOAD validation"
cargo test                                              # unchanged behaviour
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps          # all intra-doc links resolve
cargo clippy -q --all-targets                           # no new warnings
.claude/bin/machine-lock release

# Manual spot-check:
open target/doc/pwdft_rs/index.html                     # crate root has prose
open target/doc/pwdft_rs/scf/index.html                 # scf has orientation
```

Acceptance:

- All 14 enumerated files start with `//!` (verifiable via the same
  loop the proposal-research used).
- The crate-root rustdoc page renders with the orientation paragraph,
  not just a module list.
- `cargo doc -- -D warnings` (via RUSTDOCFLAGS) is green.
- No source code changes outside `//!` blocks (i.e., zero behaviour
  change, zero risk of regression).

## Out of scope

- MADOC owns mathematical documentation (the dense math derivations,
  references to papers, equation rendering). MOAD writes "what" and
  "who calls it"; MADOC writes "and here's the derivation."
- MAUD owns correctness audit prose (e.g. "this implementation matches
  Eq. 12 of Marzari 1999"). MOAD does not duplicate that work.
- `src/main.rs` is excluded — binaries don't appear in `cargo doc`,
  and `main.rs` is short enough that a `//!` would be a header without
  a body.
- `src/symmetry/operations.rs`, `src/symmetry/detect.rs`,
  `src/symmetry/kpoints.rs`, and the `src/scf/{driver,driver_spin,
  context,density,energy,grid,initial_density,potentials,report,
  smearing}.rs` files: each may or may not have a header today; this
  proposal scopes only to the 14 enumerated above. A second sweep can
  cover them if MOAD lands cleanly.
