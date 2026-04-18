---
id: MPSH
title: Monkhorst-Pack shift alignment with QE convention
status: draft
priority: high
complexity: small-medium
risk: low
depends_on: []
blocks: [VQEF]
owner: core-engineer
---

# MPSH — Monkhorst-Pack shift alignment with QE convention

## TL;DR

`src/kpoints.rs::monkhorst_pack` hard-codes the **shifted** MP-1976
convention (`f_j = (2·i_j − N_j + 1)/(2·N_j)`). Our QE reference
calculations all use QE's **Γ-centered default** (`K_POINTS automatic /
Nx Ny Nz 0 0 0`), so the two codes sample **physically different
k-meshes** at identical grid density. This is the dominant residual on
every light-atom validation test in `tests/qe_validation.rs`:

| system | grid | pwdft-rs | QE | residual | test status |
|---|---|---|---|---|---|
| Si diamond | 4×4×4 | −231.865 eV | −231.610 eV | ~0.26 eV | `#[ignore]` |
| Al FCC     | 8×8×8 |  −64.197 eV |  −64.269 eV | ~73 meV  | `#[ignore]` |
| C  diamond | 4×4×4 |  stalls      |  converges  | SCF fails | `#[ignore]` |

SYKP (completed 2026-04-17) documented the mismatch and clarified
docstrings but deliberately deferred the code fix to a follow-up
proposal. This is that follow-up.

## Motivation

After NCFX closed the 13.43 eV Si gap to ~0.26 eV (see
`proposals/completed/NCFX-nlcc-core-density-fix.md`), the remaining
residual on every light-atom system (Z ≤ 14) in `qe_validation.rs`
tracks the MP-shift convention mismatch. Until pwdft-rs and QE sample
the **same** physical k-points, we cannot close this residual below the
50 meV tolerance that VQEF's validation matrix requires. Four
`#[ignore]` markers are attributable to this single issue:

- `test_si_diamond_vs_qe` — 0.26 eV residual, 23 meV Γ eigenvalue
  discrepancy
- `test_c_diamond_vs_qe` — SCF stalls at Δρ ≈ 4.1e-6 after 80 iters
  (the shifted grid misses the special Γ/X/L points that nucleate
  C diamond's bonding structure efficiently; QE's Γ-centered grid
  converges in 9 iters)
- `test_al_fcc_vs_qe` — 73 meV on 8×8×8 Al (just above the 50 meV
  tolerance)
- Every heavy-atom test in `qe_validation.rs` **also** carries this
  residual on top of the VGCH Z>14 V_local issue; closing MPSH is a
  prerequisite for cleanly diagnosing VGCH because the two residuals
  currently sum.

## Problem

### The convention asymmetry

Monkhorst & Pack (PRB **13**, 5188 (1976)) give two distinct formulas
for a uniform k-grid:

- **Eq. 4 (shifted):** `f_j = (2·i_j − N_j + 1)/(2·N_j)`,
  `i_j ∈ {1, …, N_j}`. For N=4: `{−3/8, −1/8, 1/8, 3/8}` — **no
  point at Γ**, symmetric about Γ. This was the original recipe,
  motivated by the observation that avoiding BZ-boundary points improves
  integration accuracy for many insulators and simple metals.
- **Γ-centered:** `f_j = (i_j − 1)/N_j`, `i_j ∈ {1, …, N_j}`. For N=4:
  `{0, 1/4, 1/2, 3/4}` — **Γ is included**, plus all high-symmetry
  points (X, L, K, W on the FCC BZ boundary) that coincide with the
  mesh.

For odd N both formulas reduce to the same Γ-centered grid modulo a
cyclic permutation. For even N they are genuinely distinct meshes:

- Shifted grid IBZ reduces differently from Γ-centered (e.g. Si FCC
  4×4×4: 10 vs 8 irreducible points).
- Symmetrically-equivalent on paper (both sample the reducible mesh
  uniformly), but at finite `ecut` the Γ-centered grid's high-symmetry
  points contribute eigenvalues that the shifted grid never samples.
  This shows up as a small but non-zero total-energy difference.
- For covalent systems with narrow bonding states around Γ (like
  diamond C), the mixer needs either many more iterations or a better
  initial sample to reach the same SCF tolerance on the shifted grid.
  This is why C diamond stalls.

### QE's convention

`qe-7.5/PW/src/kpoint_grid.f90:67-78` (QE 7.5 source):

```fortran
DO i=1,nk1
   DO j=1,nk2
      DO k=1,nk3
         n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
         xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
         xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
         xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
      ENDDO
   ENDDO
ENDDO
```

QE's formula is `xkg(α,n) = (i_α − 1)/N_α + k_α/(2·N_α)`, where
`k_α ∈ {0, 1}` is the user-supplied shift flag from the
`K_POINTS automatic / Nx Ny Nz k1 k2 k3` line. For `k1=k2=k3=0`
(the **default in every `qe_validation/*.in` file**) this is the
Γ-centered grid. For `k1=k2=k3=1` it matches pwdft-rs' shifted grid
modulo a cyclic reshuffle of the index ordering.

### What pwdft-rs currently does

`src/kpoints.rs:37-80` (pwdft-rs):

```rust
let f1 = f64::from(2 * i1 as i32 - n1 as i32 + 1) / (2.0 * f64::from(n1));
let f2 = f64::from(2 * i2 as i32 - n2 as i32 + 1) / (2.0 * f64::from(n2));
let f3 = f64::from(2 * i3 as i32 - n3 as i32 + 1) / (2.0 * f64::from(n3));
```

This is MP-1976 Eq. 4 directly — the shifted convention, with no
configurable shift. `src/symmetry/kpoints.rs:109-115` (`mp_fractional`)
uses the same formula so the IBZ reduction is consistent; but both
sides of the identity are equally non-QE.

### Why SYKP deferred the fix

SYKP (completed 2026-04-17) established that the convention mismatch
is real, the two IBZ reductions (10 vs 8 on Si 4×4×4) are both correct
for their respective grids, and the total-energy difference at finite
ecut is small (sub-meV expected from k-convergence theory). SYKP chose
to not ship a code change because at the time the Si gap was dominated
by the ~13.4 eV NLCC bug — a 23 meV shift correction was below the
signal threshold.

Post-NCFX (2026-04-18), the NLCC bug is fixed. The 23 meV Si residual
and 73 meV Al residual are now the **dominant** remaining signal on
light-atom systems, so MPSH is promoted from "documentation-only" to
a real code change.

## References

### Primary literature

- Monkhorst, H. J.; Pack, J. D. *Special points for Brillouin-zone
  integrations.* **Phys. Rev. B 13, 5188 (1976).** Eq. 4 is the
  shifted formula; §III discusses the equivalence of shifted and
  unshifted grids for integration of periodic functions.
- Pack, J. D.; Monkhorst, H. J. *"Special points for Brillouin-zone
  integrations" — A reply.* **Phys. Rev. B 16, 1748 (1977).**
  Clarifies that the Γ-centered grid is preferable when symmetry
  points lie on special lines in the BZ.

### QE source

- `qe-7.5/PW/src/kpoint_grid.f90:47-78` — the `kpoint_grid` subroutine;
  lines 67-78 contain the grid formula reproduced above.
- `qe-7.5/Modules/input_parameters.f90` — declares the `K_POINTS
  automatic` syntax including the `k1 k2 k3` shift integers.
- `qe-7.5/PW/src/setup.f90:673` — applies spin-degeneracy factor to
  weights (post-reduction); unrelated to shift but worth noting when
  comparing printed weight sums.

### pwdft-rs source

- `src/kpoints.rs:37-80` — `monkhorst_pack` (the target of this fix).
- `src/symmetry/kpoints.rs:109-115` — `mp_fractional` (the grid-index
  inverter used inside IBZ reduction; must stay in lockstep with
  `monkhorst_pack`).
- `src/settings.rs:108-128` — `KPointSettings::MonkhorstPack { grid }`
  (YAML-facing config; needs a `shift` field).
- `src/main.rs` — wire through YAML `shift` into the call to
  `monkhorst_pack`.

### Related proposals

- `proposals/completed/SYKP-symmetry-ibz-audit.md` — documented the
  mismatch, deferred the fix (§D2).
- `proposals/completed/NCFX-nlcc-core-density-fix.md` — closed the
  13.4 eV Si gap, leaving MPSH as the dominant residual.
- `proposals/VQEF-*` (in flight) — will cite this proposal as a
  dependency for its full-matrix acceptance.

## Implementation

Pure prose; no code diff in this proposal. The fix is mechanical and
should fit in 1-2 CE-days.

### Method (high-level)

1. **Add a `KGridShift` enum** with variants `Gamma`, `MP1976`, and
   `Custom([u32; 3])` (where the u32 shift components are 0 or 1, per
   QE's convention). `Custom` is forward-looking; the immediate need
   is just `Gamma` and `MP1976`.
2. **Change `monkhorst_pack` signature** from
   `monkhorst_pack(n1, n2, n3, lattice)` to
   `monkhorst_pack(n1, n2, n3, shift, lattice)`.
3. **Pick the default.** Two reasonable options:
   - **(a) Γ-centered as new default.** Matches QE out of the box;
     every existing validation test drops its `#[ignore]` after a
     single config change. Breaks the 10-IBZ-point behavior on
     4×4×4 Si that legacy tests pin, so pin updates are required in
     `src/symmetry/kpoints.rs::tests`, `tests/free_electron_bands.rs`
     (if any asserts on k-count), and any YAML `examples/`.
   - **(b) Keep shifted as default, opt-in Γ-centered via YAML.**
     Preserves all existing behavior; requires only the validation
     tests to opt into Γ-centered. Safer, but the default remains
     out-of-step with QE for new users.

   **Recommendation: (a).** The argument for matching QE out of the box
   is that (i) every paper-scale validation comparison currently costs
   the user a manual YAML flag, (ii) users migrating from QE will be
   surprised that the "MP 4×4×4" line gives different answers, and
   (iii) the shifted grid was historically motivated by a subtle
   integration-accuracy argument that modern codes have mostly dropped
   in favor of the simpler Γ-centered default (VASP, Abinit, and QE
   all default Γ-centered for the `automatic` specifier).

4. **Add `KPointSettings::MonkhorstPack { grid, shift }`** in
   `src/settings.rs`. YAML syntax:
   ```yaml
   kpoints:
     type: monkhorst_pack
     grid: [4, 4, 4]
     shift: gamma       # or mp1976, or [k1, k2, k3]
   ```
   Default the field via `#[serde(default = "...")]` to `gamma`.
5. **Update `mp_fractional`** in `src/symmetry/kpoints.rs` to accept
   the same `shift` parameter (or embed it in the grid state), so
   the IBZ reduction and the grid generator stay consistent.
6. **Drop `#[ignore]` on the 4 affected tests** (`test_si_diamond_vs_qe`,
   `test_c_diamond_vs_qe`, `test_al_fcc_vs_qe`, and the implicit
   "eigenvalues ≤ 10 meV" assertions across heavy-atom tests that
   currently inherit this residual). The 5 heavy-atom tests will
   **still** be `#[ignore]` pending VGCH.

### Risks

- **Pin churn.** Every integration test that asserts on k-point counts
  or on total energies computed from a specific MP grid will need its
  pins refreshed. Grep for `monkhorst_pack(` calls and audit; expect
  ~6-10 pin updates.
- **API break.** `monkhorst_pack` gains a parameter. Internal API only
  (no public API impact on users who drive via YAML), but all internal
  call sites need updating in lockstep with the signature change.
- **IBZ reduction.** `reduce_kpoints` in `src/symmetry/kpoints.rs` must
  correctly reduce the Γ-centered grid (which includes high-symmetry
  points). The existing reducer handles arbitrary fractional
  coordinates; the only pitfall is that Γ itself has a stabilizer
  equal to the full point group (weight = 1/N_tot after normalization),
  which is already what `reduce_kpoints` computes. Unit-test the
  4×4×4 Si case expecting 8 irreducible points (to match QE).

### What NOT to do

- Do **not** silently change the default under user YAML files that
  don't specify a `shift` field. If the decision is (a) Γ-centered
  default, that is a breaking change and must be announced; add a
  migration note to CHANGELOG / README. Alternative: default (b) plus
  an aggressive docstring.
- Do **not** skip updating `mp_fractional`. Its inverse-map
  (`frac_to_grid_index`) is called by the IBZ reducer; if only
  `monkhorst_pack` is updated the reducer will silently mis-index
  Γ-centered grids.

## Acceptance

Ship MPSH when **all** of the following hold:

1. **Si 4×4×4 Γ eigenvalues** match QE `qe_validation/si_scf.out`
   to ≤ 10 meV on every printed band at Γ (currently differ by ~1 eV
   on the topmost valence bands).
2. **Si 4×4×4 E_total** matches QE to ≤ 50 meV (currently 0.26 eV).
3. **Al 8×8×8 E_total** matches QE to ≤ 50 meV (currently 73 meV).
4. **C 4×4×4 SCF converges** within the default 80-iter budget
   (currently stalls at Δρ ≈ 4.1e-6).
5. **Fe BCC 8×8×8 E_F** matches QE to ≤ 10 meV. (The E_total gap is
   VGCH territory, but E_F should close once the k-sample matches.)
6. **IBZ reduction pins updated.** `src/symmetry/kpoints.rs::tests::
   test_si_4x4x4_reduces_to_8` passes with `== 8` (not `∈ [8, 10]`)
   for the Γ-centered branch.
7. `tests/qe_validation.rs` has `#[ignore]` dropped on
   `test_si_diamond_vs_qe`, `test_c_diamond_vs_qe`, and
   `test_al_fcc_vs_qe` (the 5 heavy-atom tests remain ignored pending
   VGCH — that is the correct residual state after MPSH alone).
8. Default behavior and any migration instructions documented in
   `CLAUDE.md` and `README.md`.

Pass-criterion tolerance rationale: 10 meV Γ-eigenvalue is the rough
floor set by the ONCV LDA PP's internal interpolation precision; 50 meV
total-energy is the VQEF matrix target for light-atom systems.

## Cost

~2 CE-days, broken down roughly:

- 0.5 d: signature change (`monkhorst_pack` + `mp_fractional`) plus
  all internal call-site updates.
- 0.5 d: YAML settings field + serde default + parser wiring.
- 0.5 d: pin updates across integration tests + IBZ reducer test
  for the new expected 8-point reduction.
- 0.5 d: validation run (the 4 light-atom `#[ignore]` drops) +
  regression check on full test suite + docstring / CHANGELOG updates.

No new physics; no new dependencies; pure plumbing plus a careful
choice of default. Risk is low.

## Non-goals

- **Shifted grid removal.** The shifted convention remains useful for
  Γ-phonon calculations and some response calculations; the enum
  keeps it as a first-class option.
- **Automatic shift selection.** Some codes pick shift based on
  symmetry (e.g. shifted for bcc, unshifted for fcc). Out of scope
  for MPSH — users can pick explicitly.
- **Unified shift + IBZ reducer path.** `reduce_kpoints` already
  handles arbitrary fractional coordinates; this proposal does not
  require any change to the reducer's core algorithm, only to the
  grid generator it consumes.

## Related

- SYKP (completed) — documented the mismatch; deferred the fix.
- VGCH (draft, sibling proposal) — the heavy-atom V_local residual
  that MPSH does not address. VGCH + MPSH together unblock VQEF's
  full 8-system validation matrix.
- VQEF (in flight) — cites MPSH as a dependency.
