---
id: VGCH
title: Heavy-atom V_local(G) residual — post-VGCMP continuation
status: active
priority: high
complexity: medium-large
risk: medium
depends_on: []
blocks: [VQEF]
owner: researcher
---

# VGCH — Heavy-atom V_local(G) residual (post-VGCMP continuation)

## Status (2026-04-19, post-Phase-1a)

Phase 1a diagnostic landed as **PR #139**. Key empirical findings on
Cu and Fe heavy-atom per-component audit:

- **One-electron sum: +37 eV too high** (on Cu).
- **Hartree: −24.6 eV too low** (on Cu).
- These partially cancel to the net +17 eV Cu residual.
- Signature: **SCF converges to a different density**, not a
  form-factor bug.

**Ruled out by Phase 1a:**

- V_local(G=0) Z-scaling — Python reference in
  `scripts/validate/vgch_vloc_heavy.py` matches Rust `v_local_of_g(0, Ω)`
  to all printed digits on every heavy-atom PP.
- Ewald Z² scaling — `test_fe_bcc_ewald_vs_qe` stays green at <0.01 eV;
  Cu Ewald Δ = 3·10⁻⁵ eV.

**Phase 1b scope (remaining work, ~1 CE-week):**

Investigate three hypotheses in order of prior probability:

1. **Non-local β_q projector form factors** on heavy species. If the
   pwdft-rs `radial_fourier_beta` handles semicore states or large-l
   projectors differently than QE's `init_us_1.f90`, the non-local
   projection would systematically shift the one-electron sum. Cross-
   check β_l(q) on Fe / Cu at production ecut.
2. **Initial density (SAD) for heavy atoms with semicore states.**
   If SAD mis-represents the semicore density, the SCF may converge
   to a density basin that QE avoids (QE initializes from atomic
   orbitals via `starting_wfc`). Test: initialize pwdft-rs from a
   density that matches QE's first-iteration density exactly and see
   if subsequent iters converge.
3. **Mixer basin / multi-basin SCF.** If defect (2) isn't the cause,
   the mixer might be stabilizing a different local minimum of the
   energy functional. Harder to diagnose cleanly; requires comparing
   occupation numbers and eigenvectors at each iter against QE's.

Phase 1a's `tests/vgch_per_component_heavy.rs` is the reusable
harness. Each hypothesis adds one arm to that test.

## TL;DR

VGCMP Phases 1–4 proved that the pseudopotential → Hamiltonian
assembly pipeline is bit-correct on Si: V_local(G) to 3×10⁻⁸ eV,
β_l(q) to 10⁻¹² Bohr^(3/2), D_ij to machine precision, assembled
H[G,G] to 4×10⁻⁹ eV. NCFX (post-VGCMP) closed the remaining 13.4 eV
Si gap via two NLCC Bessel-transform bugs. After both landed, light
atoms (Z ≤ 14) match QE to ~0.03 eV (residual = MPSH).

**But every Z > 14 system still carries a 7–34 eV total-energy gap
vs QE** that scales with species, cannot be explained by MPSH, and is
not closed by the Si-validated NLCC fix:

| system | species | pwdft-rs | QE | residual | NLCC? | semicore |
|---|---|---|---|---|---|---|
| Fe BCC     | Fe (Z=26, 16e⁻)        | −3050.80 eV | −3060.16 eV | ~9.5 eV  | yes | 3s/3p/3d |
| GaAs       | Ga (Z=31, 13e⁻) + As (Z=33, 15e⁻) | −4155.95 eV | −4189.59 eV | ~33.6 eV | yes+yes | 3d(Ga) / 3d(As) |
| Cu FCC     | Cu (Z=29, 19e⁻)        | −4837.47 eV | −4853.64 eV | ~16.2 eV | yes | 3s/3p/3d |
| NaCl       | Na (Z=11, 9e⁻)  + Cl (Z=17, 7e⁻)  | −1621.94 eV | −1629.69 eV | ~7.7 eV  | no+yes  | 2s/2p(Na) |
| MgO        | Mg (Z=12, 10e⁻) + O  (Z=8,  6e⁻)  | −1993.15 eV | −2003.24 eV | ~10.1 eV | no+yes  | 2s/2p(Mg) |

Observation: systems match the heavy-atom pattern even when some
species are nominally "light" (Na Z=11, Mg Z=12), because those PPs
include semicore states (Na 2s/2p, Mg 2s/2p) that raise
`z_valence` well above the valence-only count. **The common
signature is `z_valence ≥ 7` and/or the presence of semicore
states that extend the PP's radial support into the core region.**

## Motivation

5 of 8 `tests/qe_validation.rs` systems remain `#[ignore]`'d solely
on VGCMP-residual grounds. VQEF's full validation matrix cannot
close without VGCH landing. The Si-validated pipeline has been proven
bit-correct on one system — but that system's PP has
`z_valence=4.00`, l_max=2, no semicore. Extrapolating bit-correctness
from Si to Fe/Cu/Ga PPs is not physics; it's hope.

## Problem

### What VGCMP Phases 1–4 proved

On **Si only**, independent Python references confirmed agreement to
floating-point round-off:

- **Phase 1** — V_local(G), 20 G-shells, max |Δ| = 2.78×10⁻⁹ Ry.
- **Phase 2** — β_l(q), 6 projectors × 20 q, max |Δ| = 3.03×10⁻¹²
  Bohr^(3/2).
- **Phase 3** — D_ij, 6×6 diagonal, max |Δ| = 0.0 (bit-exact).
- **Phase 4** — assembled H[G,G] at Γ, first 5 shells, max |Δ| =
  2.9×10⁻¹⁰ Ry.

Then **VGC5** (per-component energy accounting) identified E_xc
(NLCC) as the 13.4 eV Si gap, and **NCFX** fixed two Bessel-transform
bugs: (i) PP_NLCC unit conversion `/BOHR_TO_ANG` should be
`/BOHR_TO_ANG³`, (ii) missing `r²` radial weight and `4π` prefactor in
`compute_core_density`. Post-NCFX Si gap = 0.26 eV (MPSH residual).

### What VGCMP Phases 1–4 did **not** prove

- Phase 1 checked Si.upf only. It did not check that
  `v_local_of_g` is bit-correct for PPs with **larger `z_valence`**
  (e.g. Cu 19e⁻, Fe 16e⁻), where the erf-subtraction term
  `Z·e²·erf(r)/r` is a ~5× larger number and any sign-error, missing
  constant, or radial-mesh edge artifact scales linearly with Z.
- Phase 2 checked 6 projectors at l=0,1,2. Heavy atoms in our UPFs
  go up to l=2 (l_max=2), but Fe has 6 projectors (2 per l), same
  as Si — so projector *count* scaling is already tested. What was
  **not** tested: projectors with **no node at large r** (transition
  metal d-orbitals decay slowly) or projectors whose radial support
  extends near the log-mesh outer boundary.
- Phase 3 checked Si's D_ij is strictly diagonal. For heavier atoms,
  the ONCVPSP diagonalization may or may not produce a strictly
  diagonal D_ij — the UPF format allows block-diagonal D_ij per l,
  and there is no unit test that pwdft-rs handles off-diagonal
  within-l blocks correctly. Neither Cu nor Fe have been spot-checked.
- Phase 4 checked H[G,G] at Γ. Off-diagonal H[G,G'] (which exercises
  the structure factor `S(G−G')` for non-zero argument) was never
  tested (Phase 4b was explicitly deferred).

### NCFX closed Si, but does it close heavy atoms?

NCFX changed two code paths:

1. `src/pseudopotential/upf.rs` — PP_NLCC unit conversion. Applies to
   **every** PP with `core_correction="T"`. So any NLCC-heavy element
   benefits.
2. `src/scf/potentials.rs::compute_core_density` — missing `r²` + `4π`
   in the Bessel transform. Applies to **every** PP with
   `core_correction="T"`.

Post-NCFX Fe 4×4×4 went from ΔE_xc = −48.85 eV to +0.69 eV — the
NLCC fix dominates. But Fe 8×8×8 still carries a 9.5 eV total-energy
residual. Since E_xc is now ~ correct and V_local(G) is bit-correct on
**Si**, the Fe 9.5 eV must live elsewhere.

**Critical observation:** NaCl (7.7 eV gap) has Na with NLCC=F
(no NLCC!) and Cl with NLCC=T. Mg (10.1 eV via MgO) has NLCC=F.
If the residual were purely NLCC-related, NaCl and MgO would have
smaller residuals than Fe/Cu/GaAs. They don't — they all sit in the
7–34 eV range. **Whatever causes the heavy-atom residual cannot be
purely an NLCC effect.**

### Candidate root causes

The observed residual energy scales loosely with `z_valence × n_atoms`,
ranging from ~0.5 eV/valence-electron to ~1 eV/valence-electron.
Candidates (not mutually exclusive):

**(a) V_local(G=0) compensating shift for large Z.**
`src/pseudopotential/mod.rs:137-148` computes the G=0 branch as
`(4π/Ω) ∫ r²·[V_loc(r) + Z·e²/r] dr` — the bracket is a
short-ranged `v_short`. The integrand scales with Z. Any sub-leading
error (e.g. radial mesh truncation at large r where `V_loc → −Z·e²/r`
only approximately cancels the Coulomb tail) would scale with Z. On Si
(Z=4) this error is ≤ meV; on Cu (Z=19) it could be ~1 eV.

**(b) NLCC radial extent and mesh truncation for transition metals.**
Fe/Cu 3d core densities extend further than Si 3s/3p. If the UPF
radial mesh truncates the core density before it decays to zero, the
missing tail contributes a systematic shift. This is tested on Si
(`test_si_core_charge_integrates_to_partial_core` checks ∫ρ_core = 0.74 e)
but **not** on Fe/Cu. QE's `rhoc_mod.f90:107-115` does the same
Simpson integration — if both codes truncate identically, the residual
is convergence-meaningful, not a bug.

**(c) Semicore state handling.**
Na (Z=11) has 1s²2s¹ → `z_valence=9` means the 2s²2p⁶ semicore
is explicitly included in the valence. Ga (Z=31) has 3d semicore;
`z_valence=13` (3d¹⁰4s²4p¹). These inner valence states have radial
nodes near the core, which pushes the required ecutrho higher than
the pure-valence cases. Our 4×4×4 Fe run uses `ecut=15 Ry` — a value
calibrated for Si valence states, potentially insufficient for Fe
3d semicore. QE adapts its internal FFT grid via `ecutrho`; pwdft-rs
uses `ecutrho_ratio=4` which may not resolve the semicore radial
oscillations. **Check: run Fe at ecutwfc = 40 Ry and see if the
residual shrinks.**

**(d) Non-local projector radial node truncation.**
ONCVPSP generates projectors with carefully-tuned cut-off radii.
For heavy atoms, the l=2 d-projectors can extend well past the
typical Si projector range. VGCMP Phase 2's 6-projector Si check
was at q_max = 7 Bohr⁻¹; for heavier atoms with finer radial detail,
the Bessel transform may need higher q_max or denser q-grid. Unlikely
to explain 10+ eV (VGCMP showed β_l(q) agrees to 10⁻¹² Bohr^(3/2) on
Si; the analogous check on Fe would confirm this is a non-issue).

**(e) Ewald for large ionic charges.**
Ewald energy scales as Z². Si Ewald is −228.5 eV; Fe Ewald is not
yet pinned in VGC5 on Fe 8×8×8. VGCMP Phase 5 (VGC5) pinned
E_ewald(Si) to QE within 0.011 eV, but Fe was only tested at 4×4×4
nspin=1 (Δ_ewald not reported separately). At large Z, the erfc
convergence parameter `η` and the cutoff choice become more sensitive.

**(f) Assembled H[G,G'] for G ≠ G' (Phase 4b).**
Phase 4 tested diagonal only. The structure factor `S(G−G')` at
non-zero argument exercises `exp(−i(G−G')·τ_atom)` with real phase;
a sign error here would be invisible on the Si Phase 4 test (which
pinned `S(0) = N_atom`) and invisible on the PCFX symmetrization test
(which is a self-consistent identity). For a 2-atom basis (Si, C,
GaAs, NaCl, MgO) the structure factor at G≠G' is `1 + exp(−iτ·(G−G'))`;
for 1-atom (Fe, Cu, Al) it is `exp(−iτ·(G−G')) = 1` (τ=0) — i.e. the
structure factor is trivial for 1-atom cells. But Fe **still** has a
9.5 eV gap, which rules out (f) as the Fe cause.

### Prior-probability ranking (before any new data)

1. **(c) Semicore/ecut convergence.** Highest prior. Cheap to test
   (rerun at ecut=40 Ry). If residual shrinks, this is the
   dominant effect and the fix is documentation ("heavy atoms
   need ecut ≥ 40 Ry"), not code. Explains why residual scales
   with z_valence.
2. **(a) V_local(G=0) scaling with Z.** Moderate prior. Testable by
   VGCMP-style Python cross-check on Fe/Cu `v_local_of_g`.
3. **(d) β_l(q) heavy-atom extrapolation.** Moderate prior. Testable
   by extending `beta_q_reference.py` to Fe/Cu.
4. **(e) Ewald for large Z.** Moderate prior. Testable by adding Fe,
   Cu, Ga/As to the VGC5 per-component pins.
5. **(b) NLCC mesh truncation.** Low prior given NaCl/Mg residuals
   (NaCl Na + MgO Mg have `core_correction="F"`).
6. **(f) Structure factor off-diagonal.** Low prior given 1-atom Fe's
   9.5 eV gap (τ=0 makes S(G−G') trivial).

## References

### Primary literature

- Louie, S. G.; Froyen, S.; Cohen, M. L. *Nonlinear ionic
  pseudopotentials in spin-density-functional calculations.*
  **Phys. Rev. B 26, 1738 (1982).** Original NLCC paper; relevant to
  (b).
- Hamann, D. R. *Optimized norm-conserving Vanderbilt pseudopotentials.*
  **Phys. Rev. B 88, 085117 (2013).** ONCVPSP construction — relevant
  to (c), (d).
- Kresse, G.; Furthmüller, J. *Efficient iterative schemes for ab
  initio total-energy calculations using a plane-wave basis set.*
  **Phys. Rev. B 54, 11169 (1996).** §III.D discusses semicore state
  handling and ecut requirements; relevant to (c).
- Martin, R. M. *Electronic Structure: Basic Theory and Practical
  Methods.* Cambridge 2004. §11–13 cover pseudopotential construction
  and the long-range/short-range splitting underlying (a).

### QE source

- `qe-7.5/upflib/vloc_mod.f90:100-175` — `init_tab_vloc`; lines
  130-165 are the erf-subtracted V_local(G) tabulator; line 148
  final multiply by `fpi/omega` and line 138 erf subtraction
  (relevant to (a)).
- `qe-7.5/upflib/vloc_mod.f90:179-225` — `interp_vloc`; 4-point
  Lagrange interpolation. Our `v_local_of_g` is direct-evaluated (no
  tabulation), so the interpolation path is not exercised — this is
  a known difference but shouldn't be a heavy-atom bug (Si Phase 1
  showed direct-eval matches QE's tabulated values to 10⁻⁹ Ry).
- `qe-7.5/upflib/rhoc_mod.f90:101-120` — `init_tab_rhc`; NCFX's
  reference for the `r²·4π` NLCC Bessel transform.
- `qe-7.5/PW/src/init_us_1.f90` — initializes D_ij for each PP
  (relevant to (g) — ONCVPSP diagonalization check for heavy atoms).
- `qe-7.5/upflib/init_us_1_base.f90` — baseline PP initialization.
- `qe-7.5/PW/src/set_rhoc.f90:29-125` — applies NLCC to FFT grid.

### pwdft-rs source (for each candidate)

- `src/pseudopotential/mod.rs:131-184` — `v_local_of_g` (candidate (a)).
- `src/pseudopotential/upf.rs` — PP_NLCC parse (candidate (b) via
  truncation).
- `src/scf/potentials.rs::compute_core_density` — NLCC Bessel
  transform (candidate (b)).
- `src/scf/context.rs:93-94` — V_local(G=0) zeroing and stash
  (candidate (a)).
- `src/ewald.rs` — Ewald sum (candidate (e)).
- `src/potential/nonlocal.rs::bessel_transform_projector` — β_l(q)
  (candidate (d)).

### Related proposals

- `proposals/VGCMP-vloc-g-cross-check.md` (Phases 1–4 done; kept as
  active for the Phase 4b off-diagonal follow-up deferred to here).
- `proposals/completed/VLQR-vloc-qe-reference-data.md` — closed the
  Cube-based V_local test; superseded by VGCMP Phase 1.
- `proposals/completed/VGC5-per-component-energy-accounting.md` —
  per-component infrastructure VGCH will reuse.
- `proposals/completed/NCFX-nlcc-core-density-fix.md` — closed the
  Si 13.4 eV gap; demonstrates the per-component diagnostic pattern
  VGCH will follow.

## Implementation

Phased investigation. Do **not** commit to a specific code fix up
front — VGCMP Phases 1–4 showed that the form factors are correct,
so the fix (if any) is most likely in setup/convergence land, not
raw numerics. Each phase is a day or two of work.

### Phase 0 — Relabel `#[ignore]` strings (mechanical)

**Motivation.** VGCMP Phases 1–4 closed the Si V_local(G) assembly
pipeline to bit-precision, and NCFX closed the remaining Si NLCC
gap. Despite that, five `#[ignore]` strings in `tests/qe_validation.rs`
still attribute the heavy-atom residual to "V_local" — an attribution
that VGCMP has actively *disproven* on Si and that VGCH's candidate
root-cause table (above) treats as only one of six hypotheses (and
not the most likely: (c) semicore/ecut convergence ranks highest).
The current strings misinform any reader (human or CI log scraper)
about where the bug lives. Relabel them to reflect VGCH ownership and
an explicitly **TBD** root cause; any future narrower attribution
should come from Phase 2's per-component diagnostic, not from
inertia on the old VGCMP-era prose.

**Exact edits in `tests/qe_validation.rs`** (documentation-only;
all inside `#[ignore = "..."]` attribute strings; no code path, no
test behavior, no pin value changes):

- **Fe BCC (`test_fe_bcc_fm_vs_qe`, line ~368):**
  - Before: `"CCMX fixes convergence (E = -3050.80 eV); ~9.5 eV gap vs QE -3060.16 eV blocked on VGCMP (heavy-atom V_loc)"`
  - After:  `"~9.5 eV gap vs QE (E_pwdft = -3050.80, E_qe = -3060.16 eV); root cause TBD, tracked in VGCH"`

- **GaAs (`test_gaas_zincblende_vs_qe`, line ~423):**
  - Before: `"VGCMP: heavy-atom V_loc residual ≈33.6 eV on GaAs (Z=31+33); pwdft-rs E = -4155.954 eV, QE = -4189.586 eV"`
  - After:  `"VGCH: heavy-atom residual ≈33.6 eV (root cause TBD) on GaAs (Z=31+33); pwdft-rs E = -4155.954 eV, QE = -4189.586 eV"`

- **Cu FCC (`test_cu_fcc_vs_qe`, line ~471):**
  - Before: `"VGCMP: heavy-atom V_loc residual ≈16.2 eV on Cu (Z=29, 3s/3p/3d semicore); pwdft-rs E = -4837.466 eV, QE = -4853.641 eV"`
  - After:  `"VGCH: heavy-atom residual (root cause TBD) ≈16.2 eV on Cu (Z=29, 3s/3p/3d semicore); pwdft-rs E = -4837.466 eV, QE = -4853.641 eV"`

- **NaCl (`test_nacl_rocksalt_vs_qe`, line ~514):**
  - Before: `"VGCMP: heavy-atom V_loc residual ≈7.7 eV on NaCl (Cl Z=17); pwdft-rs E = -1621.944 eV, QE = -1629.686 eV"`
  - After:  `"VGCH: heavy-atom residual (root cause TBD) ≈7.7 eV on NaCl (Cl Z=17); pwdft-rs E = -1621.944 eV, QE = -1629.686 eV"`

- **MgO (`test_mgo_rocksalt_vs_qe`, line ~562):**
  - Before: `"VGCMP: heavy-atom V_loc residual ≈10.1 eV on MgO (Mg semicore PP); pwdft-rs E = -1993.146 eV, QE = -2003.241 eV"`
  - After:  `"VGCH: heavy-atom residual (root cause TBD) ≈10.1 eV on MgO (Mg semicore PP); pwdft-rs E = -1993.146 eV, QE = -2003.241 eV"`

The module-level `//!` docstring at the top of `tests/qe_validation.rs`
mentions "VGCMP" several times in its "Why some tests are `#[ignore]`d"
section; a CE implementing Phase 0 should *also* sweep that docstring
so the narrative matches (ignore → VGCH, V_local attribution → root
cause TBD). This is a drive-by — no need to mint a separate proposal.

**Nature of the change.** Documentation-only, inside test-attribute
string literals and the module-level doc comment. No executable code
semantics change, no test behavior change, no pin-value change, no
quality-gate risk. `cargo test` output is unchanged (the same five
tests remain ignored, only their ignore *reason* string differs).
`cargo doc` is unaffected.

**Cost estimate.** ≤ 30 minutes of CE time (five string edits + one
docstring sweep + run `cargo test -- --list --ignored` to eyeball
the new ignore reasons).

**Acceptance.**

1. All five `#[ignore = "..."]` strings updated per the before/after
   list above.
2. Module-level `//!` docstring no longer attributes the heavy-atom
   residual to "V_local" without qualification; VGCMP references in
   the "Why some tests are `#[ignore]`d" section updated to VGCH
   where appropriate (VGCMP references to the Phase 1–4 Si work
   that *is* done stay as-is).
3. `grep -n "VGCMP:" tests/qe_validation.rs` returns nothing (the
   colon-form is only used in the now-replaced ignore strings;
   other VGCMP mentions use no colon).
4. `cargo test` passes with identical pass/fail/ignore counts as
   before the edit.

**Sequencing.** Phase 0 is the first planned CE task under VGCH but
will ship in a separate future PR alongside (or immediately before)
the Phase 1 ecut sweep. It is intentionally carved out as a
stand-alone mechanical sub-task so the semantic investigation phases
(1–5) aren't blocked on the relabel, and so the relabel itself can
be reviewed quickly without physics context.

### Phase 1 — Semicore / ecut convergence sweep (0.5 CE-day)

**Cheapest test with highest prior.** Pick Fe BCC 8×8×8 (the
smallest heavy-atom system). Run at ecutwfc = 15, 25, 40, 60 Ry
and plot |E_total − E_QE| vs ecut. Repeat with the **same** ecut
sweep in QE (re-run `qe_validation/fe_scf.in` at matching ecuts).

- **If both codes converge to the same limit at high ecut:** the
  9.5 eV gap at ecut=15 is **k-point-adjacent ecut underconvergence
  specific to pwdft-rs** (e.g. we're sampling the basis or the
  FFT grid differently than QE at low ecut). Fix: document min
  ecut for heavy atoms; possibly increase default `ecutrho_ratio`
  for PPs with `z_valence ≥ 10`.
- **If pwdft-rs plateau differs from QE plateau:** the residual is
  not ecut convergence — proceed to Phase 2.
- **If pwdft-rs doesn't plateau:** we have a deeper numerical
  problem (FFT grid, basis cutoff, or eigensolver at large
  n_pw) — escalate to a separate proposal.

**Deliverable:** `scripts/validate/heavy_atom_ecut_sweep.py` with
measured values; short writeup in researcher logbook.

### Phase 1a — Per-component diagnostic (landed in this PR)

Ran VGC5-style per-component accounting on Fe BCC (8×8×8 nspin=1,
ecut=15 Ry) and Cu FCC (4×4×4 nspin=1, ecut=25 Ry) — see new
`tests/vgch_per_component_heavy.rs` (Tier-2 `#[ignore]`'d; run via
`-- --ignored` to print the diagnostic). Observed (all in eV):

**Fe BCC (8×8×8 nspin=1):**

| term | pwdft | QE (nspin=2, collapsed) | Δ (ours−QE) |
|---|---|---|---|
| one-electron (sum) | -680.95 | -691.94 | **+10.99** |
| E_hartree | +361.77 | +362.47 | -0.70 |
| E_xc | -392.33 | -393.26 | +0.93 |
| E_ewald | -2337.17 | -2337.17 | +0.006 (clean) |
| **E_total** | **-3048.66** | **-3060.16** | **+11.50** |

**Cu FCC (4×4×4 nspin=1):**

| term | pwdft | QE (nspin=1, 8×8×8) | Δ (ours−QE) |
|---|---|---|---|
| one-electron (sum) | -1996.85 | -2033.88 | **+37.03** |
| E_hartree | +1015.94 | +1040.51 | **-24.57** |
| E_xc | -554.40 | -559.13 | +4.73 |
| E_ewald | -3301.02 | -3301.02 | +3e-5 (clean) |
| **E_total** | **-4836.33** | **-4853.64** | **+17.31** |

**Key observation:** the Cu residual is **not concentrated in a single
term**. The one-electron sum is +37 eV too high while Hartree is -24.6
eV too low; they partially cancel to +17.3 eV net. This is the
signature of a **self-consistent state on a different density** — if
the SCF converged on a more-spread-out density, E_H drops and E_kinetic
rises together. A single-term form-factor bug would concentrate the
residual; this distribution points at an SCF-convergence / mixer / or
density-initialization issue on heavy-atom cells.

Independent Python cross-check of V_local(G=0) via the bare-Coulomb
integrand `(4π/Ω) ∫ r²[V_loc(r)+Z·e²/r] dr` on all 11 heavy-atom PPs
(`scripts/validate/vgch_vloc_heavy.py`, CSV
`scripts/validate/vgch_vloc_heavy.csv`) **matches pwdft-rs's
`v_local_of_g(0, Ω)` to all printed digits** on every element. On Fe
with Z=16: Python = +5.1736 eV, pwdft-rs v_local_g0 (from the info
log) = +5.1736 eV. Candidate (a) — V_local(G=0) scaling with Z — is
ruled out.

Γ eigenvalue comparison (Fe 8×8×8, from `test_fe_bcc_fm_vs_qe` failure
dump) shows every eigenvalue shifted by ~5.1–5.3 eV vs QE, consistent
with the V_local(G=0) convention difference (pwdft zeros the G=0
component of H at Hamiltonian-assembly time; QE keeps it in `vltot`).
`E_local_g0_shift = V_loc(G=0)·N_el` in `EnergyComponents` compensates
exactly for this in the total energy sum, so the convention does not
produce a total-energy residual by itself. The +11.5 eV residual is
elsewhere.

### Phase 1b — β_l(q) form factors (CLEARED)

Landed in a separate PR. `scripts/validate/vgch_beta_l_heavy.py` + 
`tests/vgch_beta_l_heavy.rs` (Tier-2) cross-check the KB projector
Bessel transform against an independent QE-convention Simpson reference
over 11 elements × all projectors × 10 q-values = 590 rows. Max |Δ| =
3.17e-12 Bohr^{3/2} — 4 orders below the 1e-8 tolerance. **H1 cleared:
β_l(q) is bit-perfect on every VGCH heavy-atom PP including Fe/Cu
semicore 3d projectors.**

### Phase 1c — SAD initial-density (CLEARED)

Landed in this PR. `scripts/validate/vgch_sad_heavy.py` +
`tests/vgch_sad_heavy.rs` (Tier-2) cross-check pwdft-rs'
`generate_initial_density` against a QE-convention Python reference
that replicates `qe-7.5/PW/src/atomic_rho.f90` (Simpson Bessel
transform over the log mesh, structure-factor sum per species,
unnormalized IFFT, G=0 renormalization via
`qe-7.5/PW/src/potinit.f90:218-223`). Seven systems: C, Al, Fe, Cu,
GaAs, NaCl, MgO.

**Raw-sample point-wise ρ(r) diff at 200 grid points/system,
post-clamp + post-renorm:**

| system | max \|Δρ\| (e/Å³) | mean \|Δρ\| (e/Å³) |
|---|---|---|
| C diamond | 5.4e-11 | 1.0e-11 |
| Al FCC    | 5.5e-12 | 2.4e-12 |
| Fe BCC    | 4.0e-10 | 4.8e-11 |
| Cu FCC    | 4.0e-10 | 5.8e-11 |
| GaAs      | 1.2e-5  | 4.4e-7  |
| NaCl      | 8.6e-11 | 9.0e-12 |
| MgO       | 4.8e-10 | 3.8e-11 |

All seven at ≤ 1e-5 e/Å³; six at ≤ 1e-9 e/Å³. GaAs's 1.2e-5 outlier
localizes to pwdft-rs' negative-density clamp
(`src/scf/initial_density.rs:144-148`) zeroing ~2e-5 e of Gibbs-
ringing near the As core that QE keeps (see
`qe-7.5/PW/src/atomic_rho.f90:186-188` — QE explicitly comments
that clamping is "useless" because FFT round-trip makes negative
values reappear). The clamp effect is confined to ≤1 bin per atom
and accounts for a total energy shift of O(1e-5 eV) — too small to
explain GaAs's 33 eV residual.

**Critical C-diamond verdict:** pre-clamp max |Δρ| = 8.1e-7 e/Å³,
post-clamp max |Δρ| = 4.1e-11 e/Å³. C is **bit-perfect** against the
QE-convention reference. The 1.45 eV C_total residual does NOT live
in SAD.

**Clamp statistics table** (computed by
`build_sad_density_for_diagnostic_verbose`):

| system | ∫ρ pre-clamp | neg mass clamped (e) | renorm factor |
|---|---|---|---|
| C diamond | 7.999996 | 0     | 1.000000 |
| Al FCC    | 2.999999 | 0     | 1.000000 |
| Fe BCC    | 15.999998 | 0    | 1.000000 |
| Cu FCC    | 18.999999 | 0    | 1.000000 |
| GaAs      | 27.999996 | 2.1e-5 | 0.999999 |
| NaCl      | 15.999994 | 0    | 1.000000 |
| MgO       | 15.999994 | 0    | 1.000000 |

`renorm factor = 1` on six of seven confirms the clamp+renorm is a
no-op on those cells. GaAs alone sees a 1e-6 scale adjustment.

**H2 verdict: CLEARED.** The SAD pipeline (PP_RHOATOM unit
conversion, Simpson Bessel transform, per-species structure-factor
sum, IFFT) is bit-correct. No fix needed in
`src/scf/initial_density.rs` or `src/pseudopotential/upf/convert.rs`.
The remaining clamp-vs-QE pipeline delta is a fractional-eV effect at
most and cannot explain the multi-eV VGCH residuals.

### Phase 1d — Residual target (next session, scoped as VGCH-2)

Both H1 (β_l(q)) and H2 (SAD initial density) are cleared. The 7-34
eV residuals on heavy-atom cells must therefore live in either:

- **H3 — Total-energy assembly.** Specifically the V_loc(G=0)
  compensation `e_local_g0_shift = v_local_g0 · n_electrons` in
  `src/scf/energy.rs` and the `with_g0_shift` closure. Phase 1a's
  observation that Si compensates cleanly but C doesn't suggests a
  narrow bug in this branch — maybe related to how semicore PPs or
  multi-species cells aggregate the G=0 shift.
- **H4 — SCF mixer basin.** The mixer could be stabilizing a
  different local minimum of the energy functional. Diagnosing this
  cleanly needs a "transplant" experiment: seed pwdft-rs from QE's
  converged density (via QE's `save/charge-density.dat` or a
  postprocessing dump), run SCF, and see whether pwdft-rs stays
  there or drifts to its own fixed point.

The total-energy-assembly target is cheap to investigate and has
Phase 1a's fingerprint (opposite-sign one-electron vs. Hartree
partial cancellation). **VGCH-2** is spawned to take this up;
VGCH-1c recommends starting there before climbing the mixer tree.

### Phase 2 — Per-component energy diagnostic on heavy atoms (1 CE-day)

Extend the VGC5 infrastructure (`scripts/validate/vgc5_per_component.py`
and `tests/vgc5_per_component_si.rs`) to Fe, Cu, and NaCl. For
each system, tabulate:

| term | pwdft-rs | QE | Δ |
|---|---|---|---|
| E_kinetic | ? | ? | ? |
| E_local(G≠0) | ? | ? | ? |
| E_local(G=0)·N_el | ? | ? | ? |
| E_nonlocal | ? | ? | ? |
| E_hartree | ? | ? | ? |
| E_xc | ? | ? | ? |
| E_ewald | ? | ? | ? |
| E_total | ? | ? | ? |

This is the same diagnostic that localized Si's 13.4 eV gap to
E_xc (NLCC). For Fe the prediction (if ecut is not the cause) is
that the gap will concentrate in **one** term. Candidates:

- **E_local(G=0)·N_el** (large because N_el scales linearly with Z).
- **E_ewald** (scales as Z²).
- **E_nonlocal** (semicore d-projectors contribute more).

Whichever term holds the gap determines Phase 3's target.

**Deliverable:** `tests/vgch_per_component_fe.rs` and
`tests/vgch_per_component_cu.rs` (new integration tests with pins);
CSV references committed under `scripts/validate/`.

### Phase 3 — Targeted VGCMP-style cross-check on the suspect term (1 CE-day)

Based on Phase 2 output, pick the dominant-Δ term and repeat the
VGCMP methodology at the suspect numerical step:

- If **E_local(G=0)·N_el**: extend `scripts/validate/vloc_g_reference.py`
  to Fe, Cu, Ga, As, Na, Cl, Mg. Compare 20 shells per element.
  Cite heavy-atom equivalents of the VGCMP Phase 1 table.
- If **E_nonlocal**: extend `scripts/validate/beta_q_reference.py` to
  Fe, Cu. Cite heavy-atom equivalents of the VGCMP Phase 2 table.
- If **E_ewald**: add `scripts/validate/ewald_reference.py` using
  `qe-7.5/PW/src/ewald.f90` as the reference (or a direct
  screened-pair-sum Python).

**Deliverable:** new `tests/vgch_{vloc,beta_q,ewald}_heavy_cross_check.rs`
with element-by-element Python-vs-Rust pins, at the same 10⁻⁴ Ry
precision as VGCMP Phase 1. Either a clean bill of health (pushing
the hunt elsewhere) or a smoking gun (proceed to Phase 4).

### Phase 4 — Fix (1–5 CE-days, depends on what Phase 3 finds)

Scope reserved for the actual fix once Phase 2/3 isolates the root
cause. Could be as small as "bump default `ecutrho_ratio` from 4 to
6 for PPs with z_valence ≥ 10" (half a day) or as large as "rewrite
the V_local(G=0) branch to handle the long-range Coulomb correction
explicitly per-atom" (5 days).

### Phase 5 — Validate the full matrix (0.5 CE-day)

Once Phase 4 lands, rerun all 5 heavy-atom `qe_validation.rs` tests.
Drop `#[ignore]` on each whose residual closes below 50 meV/atom.
Update test pins. If any system remains above tolerance, file
a narrow follow-up proposal with the measured per-component
diagnostic and the specific hypothesis that proposal tests.

## Acceptance

Close VGCH when **all** of the following hold:

1. Fe BCC 8×8×8 total-energy residual vs QE ≤ 50 meV/atom
   (currently 9.5 eV).
2. GaAs 4×4×4 residual ≤ 50 meV/atom (currently 33.6 eV ≈ 16.8
   eV/atom; target close to QE).
3. Cu FCC 8×8×8 residual ≤ 50 meV/atom (currently 16.2 eV).
4. NaCl 4×4×4 residual ≤ 50 meV/atom (currently 7.7 eV).
5. MgO 4×4×4 residual ≤ 50 meV/atom (currently 10.1 eV).
6. `tests/qe_validation.rs` has `#[ignore]` dropped on
   `test_fe_bcc_fm_vs_qe`, `test_gaas_zincblende_vs_qe`,
   `test_cu_fcc_vs_qe`, `test_nacl_rocksalt_vs_qe`,
   `test_mgo_rocksalt_vs_qe`.
7. Per-component diagnostic pins (new `vgch_per_component_*.rs`)
   pass on Fe and Cu at minimum.

MPSH is a **co-requisite** for closing these — MPSH alone cannot
close 9.5 eV, and VGCH cannot isolate a signal smaller than the
~0.1 eV MPSH residual on light atoms. The right sequencing is:

1. Land MPSH first (matches grids between codes; light atoms
   close to ~10 meV).
2. Then run VGCH Phase 1 (ecut sweep) on a grid-aligned Fe reference.
3. Then proceed through VGCH Phases 2–5.

This proposal does NOT assume a pre-existing MPSH; each phase's
diagnostic number is robust to the ~50 meV MPSH noise floor because
the heavy-atom residuals are 100–600× larger.

## Cost

1–2 CE-weeks total. Breakdown assumes Phase 1 (ecut sweep) resolves
it: 1 CE-week. Phase 2+ (per-component diagnostic + root-cause fix):
up to 2 CE-weeks if the bug is non-trivial.

- Phase 0 (relabel obsolete VGCMP `#[ignore]` strings): ≤ 0.1 d.
- Phase 1 (ecut sweep): 0.5 d.
- Phase 2 (per-component diagnostic on Fe, Cu, NaCl): 1 d.
- Phase 3 (VGCMP-style heavy-atom cross-check on the suspect
  term): 1 d.
- Phase 4 (fix): 1–5 d, highly dependent on what Phase 3 finds.
- Phase 5 (validate full matrix, update pins, drop `#[ignore]`): 0.5 d.

The proposal is medium-large **complexity** (because the outcome is
unknown up front) but the risk is **medium**: VGCMP Phases 1–4 have
already ruled out the most common bugs, so the remaining suspects
are narrow and diagnosable.

## Non-goals

- **USPP / PAW support.** All systems here are norm-conserving;
  USPP would be a separate proposal.
- **Spin-orbit coupling.** Fe nspin=2 collapses to M=0 under this
  PP/cutoff per the current test comment; getting real magnetism
  would require either a better Fe PP or SOC support. Separate work.
- **Full-Z validation.** The 5 systems above are the existing
  `#[ignore]`'d VGCMP-blocked tests. Once those close, VGCH is
  done; any wider Z coverage (Ag, Au, Pt, W, etc.) is VQEF's
  territory.
- **Replacing VGCMP.** VGCMP's Phase 4b off-diagonal H[G,G']
  cross-check is explicitly folded into VGCH Phase 3 if Phase 2
  points at the local-potential term; otherwise Phase 4b remains
  optional and deferred to a future proposal.

## Related

- VGCMP Phases 1–4 (done) — established bit-correctness of the Si
  assembly pipeline.
- NCFX (done) — fixed two NLCC Bessel-transform bugs; closed Si.
  The same fix applies to all NLCC-enabled PPs but does not
  explain the full heavy-atom residual (NaCl/MgO have NLCC=F
  species and still show the residual).
- VGC5 (done) — per-component energy accounting infrastructure;
  VGCH Phase 2 extends it to heavy atoms.
- MPSH (draft, sibling proposal) — k-grid alignment; light-atom
  residual. MPSH + VGCH together unblock VQEF's full matrix.
- VQEF (in flight) — cites VGCH as a dependency.
