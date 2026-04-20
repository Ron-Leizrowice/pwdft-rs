# 2026-04-20 — VGCH follow-up proposal drafts (PZPW + CNLC + VNLM-CUD)

## Scope

Drafted three physics-hypothesis follow-up proposals against VGCH-MECH
taxonomy, promoting the INDEX-row stubs (PZPW / CNLC / VNLM-CUD) to full
proposal bodies ready for Core Engineer pickup. Proposal-only PR, no
src/ changes.

## Files created

- `proposals/PZPW-lda-pw92-correlation.md`
- `proposals/CNLC-c-diamond-nlcc-pin.md`
- `proposals/VNLM-CUD-cu-d-two-radial-projector-pin.md`

INDEX.md rows unchanged — the existing one-line summaries at lines
17–19 accurately describe the proposal bodies as drafted.

## Surprises found while grounding the drafts

### 1. PW92 unpolarized wrapper is `#[cfg(test)]`-gated

`pw92_correlation(rho)` at
`pwdft/pwdft-core/src/potential/xc.rs:893` is gated `#[cfg(test)]`:

```rust
#[inline]
#[cfg(test)]
fn pw92_correlation(rho: f64) -> (f64, f64) { ... }
```

This is material for PZPW's implementation shape — VGCH-2D's
proposal implied PW92 helpers were already "usable," and they are,
but the unpolarized eV/Å wrapper specifically is not compiled into
release builds. PZPW Phase 1 step 1 explicitly calls out un-gating
this helper (it's already tested; the helper just needs its
`#[cfg(test)]` attribute stripped).

**Why the gating exists:** PBE's correlation path goes through
`pw92_correlation_au` (the Hartree-AU primitive) directly, not
through the eV/Å wrapper, because PBE mixes PW92's (ε_c, v_c) with
its own gradient correction `h(ρ, t)` before unit conversion.
The eV/Å wrapper was introduced for PW92 LDA unit tests but never
had a live call site — hence `#[cfg(test)]`.

This is a genuinely clean implementation shape for PZPW: no
partial-wiring issue, just an attribute strip + new grid lifter.

### 2. Cu.upf and Si.upf have identical PP_DIJ shape

VGCH-2F § H-C5 evidence states this as the key structural fact,
but it is worth re-emphasizing from the code side: both UPFs
carry 6 β projectors in the same (l=0,0,1,1,2,2) layout with
strictly diagonal PP_DIJ. Si ONCV passes VNMT; the only new
structural surface Cu adds vs Si is the **non-trivial magnitude**
of both l=2 radial projectors (both D_44 and D_55 are material
on Cu, vs Si where the d-projectors are gap corrections of small
amplitude).

This strengthens the VNLM-CUD hypothesis: the Cu-vs-Si bug class
is not a novel code path that Si never exercises, but rather an
amplification of the same code path on a PP where both ν-channels
carry load-bearing amplitude. VNMT zeros D_55 to isolate one
ν; VNLM-CUD keeps both D_44 and D_55 alive.

### 3. Reference scripts migrated from `scripts/validate/` to `pwdft_validation` package

VGCH-2F's § Files section still references
`scripts/validate/rho_core_g_reference.py` (the session-2 artifact).
The actual current location is
`pwdft/pwdft-validation/pwdft_validation/reference/nlcc.py`; the CLI
entry point is `uv run pwdft-validate reference nlcc ...`. `scripts/`
no longer exists at repo root.

CNLC Phase 1 calls this out correctly — the C addition lands in the
Python package, regenerated via `uv run pwdft-validate`, not via the
stale script path.

### 4. VGCH-2D/2E/2F proposal files are still in `proposals/` (not `completed/`)

INDEX.md treats them as landed (PR #177/#178/#179 referenced in
the "most recent grooming pass" paragraph at line 7), but their
files haven't been moved to `proposals/completed/` yet. Not
surprising — EM moves those as part of `/proposal complete <id>`,
and GRM12 (yesterday's consolidation) may have focused the move
operation on other completed tickets. I linked the proposals at
their current (active) paths in each follow-up's Provenance
section, not at speculative `completed/` paths.

### 5. C.upf has a moderate Q_core (Z_core = 2, 1s²)

Before drafting CNLC I assumed C's NLCC effect would be near-zero
(light atom, should be 1s² in the core). Confirmed `core_correction="T"`
in `pseudopotentials/nc/lda/C.upf` by grep. The NLCC smearing of
1s² shows up as a Q_core ≈ 0.6 e (same order as Si's 0.74 e),
which explains why ΔE_xc = +0.43 eV on C at shared density is
physically reachable. Si and C both have small Q_core, similar
XC magnitude — making the NLCC-path attribution hypothesis
coherent.

## Acceptance of each proposal (inherited from VGCH-MECH)

- **PZPW:** Fe LDA transplant residual ≤ 2 eV under `LdaPw92` at
  nspin=1 confirms H-2D-A. Full scoreboard run recommended at
  Phase 2 to attribute Class A 1 eV floor.
- **CNLC:** C `ρ_core(G)` pinned at 1e-5 e/Å³; NLCC-off ablation
  either attributes ΔE_xc = +0.43 eV to NLCC path or moves
  suspect to V_NL l=1.
- **VNLM-CUD:** Cu two-radial-d per-m test passes with both
  mutation sanity checks firing; VQEF Cu pre/post numbers
  unchanged.

All three proposals are diagnostic-first. Each has a clean
success branch (attribution → follow-up fix) and a clean
refutation branch (escalate to next suspect).

## Grounding citations anchored in code

- PZPW line numbers verified by direct read at
  `potential/xc.rs:143,394,844,894,1111`.
- CNLC template tests at `pseudopotential/upf/convert.rs:303,329,
  362,387,430,454,489,513,565,620,674,728`.
- VNLM-CUD reference at `potential/nonlocal.rs:914` (VNMT) and
  `:111,195,364-365,412` (production V_NL assembly).
- C lattice constant `celldm(1) = 6.7409 Bohr` confirmed at
  `data/qe/c_diamond_scf.in`.
- QE PW label `'PW' = PW92` at
  `qe-7.5/XClib/qe_dft_list.f90:48,72–78`.
- QE `pw(rs, iflag=1)` interpolation at
  `qe-7.5/XClib/qe_funct_corr_lda_lsda.f90:378-389`.

## Flagged for follow-up

- **VGCH-2D/2E/2F archival:** the three proposal bodies are still
  in `proposals/` (active), not `proposals/completed/`, despite
  their PRs (#177/#178/#179) having merged. Worth asking EM to
  batch their `/proposal complete` moves in the next grooming
  pass to tidy `proposals/INDEX.md` Completed table and
  `proposals/completed/` directory.
- **Staleness in VGCH-2F § Files:** references
  `scripts/validate/rho_core_g_reference.py` — the file has moved
  to `pwdft/pwdft-validation/pwdft_validation/reference/nlcc.py`.
  Minor docs-drift; not blocking.
- **Si LDA tolerance pin post-PZPW promotion:** if PZPW-F flips
  the LDA default to PW92, Si LDA may shift by 30–60 meV (r_s
  sign-transition regime). The existing Si LDA VQEF tolerance
  (already GREEN at 0.023 meV post-SiEF-B1) may need reassessment.
  Not an issue for the current proposal (diagnostic only), but
  should be on the PZPW-F follow-up checklist.
- **Multi-functional LDA grid-lifter surface area:** after PZPW
  Phase 1, pwdft-rs has Pz + LdaPw92 lifters. If PZPW-F promotes
  LdaPw92 to default, retire the Pz grid lifters to reduce
  maintenance surface (keep `perdew_zunger_correlation` as a
  library helper for historical reference but delete the grid
  pointwise lift). Not in PZPW's scope.
- **CNLC `disable_nlcc` flag lifecycle:** the diagnostic-only
  runtime flag lands as `#[doc(hidden)]` in Phase 3. After CNLC
  closes (either attribution or refutation), the flag must be
  removed. Track this on the CNLC closure checklist.

## Time

~1.5 h wall to draft all three proposals (grounded reads +
authoring). Proposals are intentionally fuller than typical
Medium-priority drafts because each is the lead hypothesis for a
live VGCH mechanism class — Core Engineer should have enough to
land without re-researching citations.
