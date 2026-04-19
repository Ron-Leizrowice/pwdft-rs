---
id: DCLN
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: [MOAD]
---

# DCLN: Strip proposal IDs and QE references from public rustdoc

## Problem

Grep of `src/**/*.rs` `///` docstrings (2026-04-19):

- **57 proposal-ID tokens** (`CFGN`, `MPSH`, `ALOC`, `WFRX`, `MADOC`, `PCFX`, `NCFX`, `VGCH`, `CCMX`, `MXBA`, `MLFX`, `CAST`, `ERR2`, `TYPE`, `GGAP`, `HYBR`, `GOPT`, `CLNP`, `QLNT`, `RDOC`, `DOCX`, `UPFV`, `GLUS`, `FGRD`, `FDLT`, `VNLM`, `VNLT`, `VNMT`, `SYMP`, `G0SH`, `MODR`, `DEAD`, `XCNI`, `TRV2`, `VGCMP`) across 20 files. `settings.rs` alone has 11.
- **51 "QE" / "Quantum ESPRESSO" / "qe_validation" references** across 15 files.

Public rustdoc is the landing page that `cargo doc` builds, visible to anyone who depends on pwdft-rs or reads the API surface. These references leak engineering workflow and reference-implementation identity into what should be self-contained physics documentation. Two concrete harms:

1. **Audit-trail drift.** A docstring that says "PCFX moved density symmetrization to G-space" becomes meaningless after PCFX is six months in the git log. Future readers skim past it or are confused. Git history is the right home for "why did this change", not rustdoc.
2. **"We copy QE" framing.** Every mention of QE in user-facing docs reinforces the wrong reading: that this is a QE clone. It isn't — it's an independent implementation that happens to validate against QE. User directive 2026-04-19: *"QE is a reference for accuracy, we are not aiming to exactly copy it."*

## Proposal

Sweep `src/**/*.rs`, operating on `///` (rustdoc) lines only:

1. **Remove proposal-ID tokens.** Rewrite sentences that reference a proposal to describe the behavior or invariant directly. If the proposal ID is load-bearing (it explains a subtle choice), move the mention to a `//` developer comment above the docstring, not inside the docstring.

2. **Remove QE references from public-facing rustdoc.** Describe the physics (e.g., "Kleinman-Bylander separable non-local pseudopotential") instead of the reference-implementation comparison ("matches QE 7.5 non-local block"). When a field corresponds to a QE input variable, document the physics the field controls, not the QE name it matches.

## Scope — what stays

- **Test files** (`tests/`) and test fixtures (`qe_validation/*.toml`, `qe_validation/*.in|*.out`): keep QE references. Tests are engineering artifacts, not user-facing docs.
- **Developer `//` comments** (not `///`): keep proposal IDs where they explain a non-obvious invariant ("see ALOC F-5 for why we preallocate here").
- **`CLAUDE.md`, `proposals/**.md`, `.claude/**`**: unchanged.
- **Unit-conversion boundary in `src/pseudopotential/upf/convert.rs`**: mentions UPF fields by name (`PP_LOCAL`, `PP_DIJ`, …). These are file-format identifiers, not QE references — keep.
- **`README.md`**: separate sweep, out of scope here.

## Risk

- **Zero to runtime behavior.** Prose-only change.
- **Low to rustdoc gate.** `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` must still pass; rewriting prose may break intra-doc links ([`foo::Bar`] forms). Each rewritten docstring is re-checked.

## Non-goals

- Not proposing a new docstring *style guide*. Scope is deletion + minor rephrasing, not a policy document. A later proposal (MOAD + a style section in CLAUDE.md) can codify the principle.
- Not rewriting docstrings that are simply *short*. The goal is to remove leakage, not to expand coverage. MOAD handles missing headers.

## Acceptance

- `rg -p '///.*\b(CFGN|MPSH|ALOC|...etc-full-list...)\b' src/` returns zero hits.
- `rg -p '///.*\b(QE|Quantum ESPRESSO|qe_validation)\b' src/` returns zero hits.
- `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` green.
- Spot-check: `settings.rs`, `potential/xc.rs`, `scf/mod.rs`, `scf/energy.rs`, `scf/mixing/mod.rs`, `kpoints.rs`, `potential/nonlocal.rs` — all should read as self-contained physics documentation, understandable without knowing the project's engineering history.

## Blocks: MOAD

MOAD (PR #112) writes 14 new module-orientation headers. If DCLN lands first, MOAD writes them clean. If MOAD lands first, those 14 headers are likely to reference proposal IDs and need rewriting. Sequence: DCLN → MOAD.
