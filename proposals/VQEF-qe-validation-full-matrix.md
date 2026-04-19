---
id: VQEF
title: Full LDA+PBE QE validation matrix (8 systems × 2 functionals, gating)
priority: high
complexity: medium
risk: low
depends_on: [VGCMP, GGAP, QELK]
blocks: []
status: draft
author: Researcher
date: 2026-04-18
---

# VQEF — Full LDA+PBE QE validation matrix

## Problem

The QEVL suite (landed 2026-04-17) generated QE 7.5 reference data for
8 systems (Si, C, Al, Fe, GaAs, Cu, NaCl, MgO) at a single functional
(LDA/PZ) and wired up `tests/qe_validation.rs` with one `#[test]` per
system. Of those 8 tests, **7 are currently `#[ignore]`d** with specific
blockers (SYKP/MPSH grid-convention mismatch, or VGCMP heavy-atom
V_local(G) residual). The test file behaves as an aspirational stub
rather than a gate: `cargo test --test qe_validation` reports green,
but it is green only because almost every assertion is skipped.

The user directive (2026-04-18):

> We should be able to match QE values closely for Al, C, Cu, Fe (FM),
> GaAs, MgO, NaCl, and Si — with both LDA and PBE in the case that each
> can converge for that system.

This proposal is the **roadmap** to close that gap. It enumerates the
blockers across the two functional tracks, estimates cost, and proposes
a phase order so that each `#[ignore]` deletion has a clearly-owned
upstream proposal. **No code is written here.** VQEF is scheduling and
scoping for the other tracks that already own the physics work.

Scope caveat:

- VQEF does **not** re-scope GGAP (Phases A–F own PBE implementation)
  or HYBR (hybrids). VQEF lives downstream of those.
- VQEF does **not** propose new physics. It enumerates known blockers.
- VQEF does **not** fix SCF convergence of C diamond, Fe FM, or any
  other system — the fixes live in the blocker proposals; VQEF tracks
  the consequences.

## 1. Current-state inventory

### Table: 8 systems × 2 functionals

Status legend:

- **GREEN** — test passes without `#[ignore]`, assertion tolerance met.
- **YELLOW** — `#[ignore]`d with a documented blocker (tracked upstream).
- **RED** — no test exists, or no QE reference data exists.

Values below were captured 2026-04-19 after GGAP Phase F-light (all 8
PBE cells wired), and reflect the latest LDA/PBE residuals in
`tests/qe_validation.rs`.

| System | Z / semicore     | LDA state                                       | LDA blocker                       | PBE state                                                   | PBE blocker                                    |
|--------|------------------|-------------------------------------------------|-----------------------------------|-------------------------------------------------------------|------------------------------------------------|
| Si     | 14 / no          | YELLOW — residual ≈33 meV                       | SYKP/MPSH (post-MPSH basin)       | **GREEN** — residual 12 meV (20 meV tol)                    | —                                              |
| C      | 6 / no           | YELLOW — residual 1.45 eV                       | VGCH-2 / VGCH light-atom          | YELLOW — residual 0.32 eV (4.5× PBE improvement)            | VGCH-2 class (partially functional-sensitive)  |
| Al     | 13 / no          | YELLOW — residual 75 meV                        | VGCH light-atom (VQEF-AL)         | YELLOW — residual 108 meV (slightly WORSE than LDA)         | VGCH light-atom (functional-insensitive)       |
| Fe     | 26 / 3s3p (NLCC) | YELLOW — residual ≈11.5 eV, NM collapse         | VGCH-2 (heavy-atom)               | YELLOW — residual 1.97 eV, FM M≈2.16μB (6× LDA improvement) | VGCH-2 class (Fe FM cell, partially fn-sens.)  |
| GaAs   | 31+33 / 3d, 3d   | YELLOW — residual 33.6 eV                       | VGCH-2                            | YELLOW — residual 17.30 eV (2× PBE improvement)             | VGCH-2 class (partially functional-sensitive)  |
| Cu     | 29 / 3s3p3d      | YELLOW — residual 16.2 eV                       | VGCH-2                            | YELLOW — residual 10.06 eV (1.6× PBE improvement)           | VGCH-2 class (partially functional-sensitive)  |
| NaCl   | 11+17 / no       | YELLOW — residual 7.7 eV                        | VGCH-2                            | YELLOW — residual 4.86 eV (1.6× PBE improvement)            | VGCH-2 class (partially functional-sensitive)  |
| MgO    | 12+8 / 2s2p (Mg) | YELLOW — residual 10.1 eV                       | VGCH-2                            | YELLOW — residual 1.56 eV (6.5× PBE improvement, largest)   | VGCH-2 class (strongly functional-sensitive)   |

**Count:** 1 GREEN / 15 YELLOW / 0 RED out of 16 target cells.

Before GGAP Phase F-light the PBE column was 8×RED + 0×YELLOW + 0×GREEN
(PBE was entirely unwired apart from Si PBE; Fe PBE was added by Phase
D as YELLOW). This phase wires the 6 remaining PBE cells and produces
the first GREEN cell in the matrix (Si PBE, 12 meV). Pre-Phase-F-light
scoreboard was `1 GREEN (Si LDA energy arm) / 8 YELLOW / 8 RED`; after
this phase `1 GREEN (Si PBE) / 15 YELLOW / 0 RED`. The Si LDA energy
arm was only GREEN because of a VQEF-QC split — its Fermi arm remains
YELLOW under V_loc(G=0) convention (VGCH Phase 1b). Net movement: one
GREEN rotated from LDA → PBE, 8 RED cells promoted to YELLOW.

### Physics findings from GGAP Phase F-light

- **Al is functional-insensitive.** Al PBE (108 meV) is *worse* than
  Al LDA (75 meV). This rules out the Al VGCH light-atom class being
  an XC-functional artifact; root cause must be in the density basin,
  projector, or symmetry machinery.
- **Every heavy/semicore system shows PBE improvement over LDA** of
  1.6×–6.5×, except Al. MgO has the largest PBE improvement (6.5×, Mg
  2s/2p semicore); Cu / NaCl have ~1.6× and the III-V GaAs ~2×. This
  tells VGCH-2 Part A that the heavy-atom partial-cancellation
  residual is *partially* functional-sensitive — a significant
  component of it moves with the XC gradient term. Worth factoring
  into the VGCH-2 diagnosis workflow.
- **C diamond is partially functional-sensitive** (4.5× PBE improvement)
  — different from Al; the same light-atom VGCH-2 class contains both
  functional-dependent (C) and functional-insensitive (Al) residuals,
  so VGCH-2 cannot be one mechanism.

Non-gating defensive tests that **do** run green today in the same file:

- `test_fe_bcc_xc_nlcc_regression_guard` — NLCC regression sentinel on
  Fe, tolerance 1 eV.
- `test_fe_bcc_ewald_vs_qe` — Fe Ewald ion-ion vs QE, tolerance 0.01 eV.

These are not part of the 8 × 2 matrix; they protect specific code paths
against regression and should stay in place.

## 2. LDA blockers (close-the-known-residuals track)

There are exactly two LDA blockers across the eight systems:

### 2a. SYKP / MPSH — Monkhorst-Pack grid convention (Si, C, Al)

**Scope inheritance.** Si (Z=14), C (Z=6), and Al (Z=13) are all light-
to-intermediate-Z; none has VGCMP heavy-atom V_local issues. Their
residuals are fully attributed to the MP-shifted vs Γ-centered mismatch
documented in SYKP (completed; closed as docstring-only with the coded
behavior unchanged).

**Proposal status.** `SYKP` landed as a documentation-only clarification
(`proposals/completed/SYKP-symmetry-ibz-audit.md`, §D2, lines 175–186).
It explicitly recommended filing a follow-up `MPSH` proposal that adds
a shift parameter to `KPointSettings::MonkhorstPack` and exposes
`monkhorst_pack_shifted(n1, n2, n3, shift, lattice)` in
`src/kpoints.rs`. **MPSH is not filed.** This is a **direct dependency
for VQEF** and should be the first new-proposal work spawned from this
roadmap.

**Proposed MPSH shape** (for the EM to sign off; researcher/core-eng
owns the detailed proposal):

- Add `shift: [u32; 3]` (each ∈ {0, 1}) to `KPointSettings::MonkhorstPack`
  in `src/settings.rs`.
- Add `monkhorst_pack_with_shift` in `src/kpoints.rs` alongside the
  existing `monkhorst_pack` (which stays backward-compatible as the
  `shift = [1, 1, 1]` call site).
- YAML schema: `{ nk: [4, 4, 4], shift: [0, 0, 0] }` for Γ-centered
  (QE `4 4 4 0 0 0`), default stays shifted.
- All three Si/C/Al `qe_validation.rs` tests set `shift: [0, 0, 0]` in
  their `QeComparisonConfig` to match QE exactly.

**Expected residual after MPSH.** Based on the documented numbers:

- Si: 0.26 eV → expected < 5 meV (MP sampling is the only remaining
  residual; per-component VGC5 audit on Si matches QE to 10⁻¹¹ eV once
  the grid aligns).
- Al: 73 meV → expected < 10 meV (simple metal, nearly-free-electron).
- C: non-convergence → expected to converge cleanly in ≤15 iters at
  same ecut=30 Ry, 4×4×4 Γ-centered (matches QE's 9 iters). If it
  doesn't, spawn a mixer investigation (wider gap = harder density
  convergence; Kerker q_TF tuning is the suspect).

**Cost estimate.** MPSH: ~1 CE-day (parameter plumbing + 2 unit tests).
VQEF test-arm update: trivially ~1 hour after MPSH lands (add `shift:
[0, 0, 0]` to 3 configs, drop 3 `#[ignore]` markers, tighten tolerances
from 50 meV to 10 meV).

### 2b. VGCMP Phase 5 — heavy-atom V_local(G) residual (Fe, GaAs, Cu, NaCl, MgO)

**Scope inheritance.** Residuals range from ~7.7 eV (NaCl, Cl Z=17) to
~33.6 eV (GaAs, Ga Z=31 + As Z=33). Fe's 9.5 eV gap carries over; Cu's
16.2 eV is an intermediate case with semicore 3s3p3d. The MgO residual
is surprising (Mg is Z=12) but attributed to the Mg PseudoDojo LDA PP
including 2s/2p semicore — same form-factor issue.

**Proposal status.** `proposals/VGCMP-vloc-g-cross-check.md` is
**active, owned by Researcher**. Phases 1–4 (landed) proved bit-level
agreement on Si's PP→H assembly pipeline. Phase 5 (energy-component
audit on heavy atoms) is the remaining work. VGCMP §"2026-04-17 Phase
1 Result" (lines 125–160) documents Si V_local(G) agreement to
2.78×10⁻⁹ Ry — i.e., whatever causes the heavy-atom residual is **not**
in the Si V_local(G) path. It is something that turns on for Z>14 or
for semicore PPs.

**What Phase 5 needs to produce** (for VQEF to drop 5 `#[ignore]`s):

- Per-component audit on Fe (as VGC5 did for Si) so we know which of
  E_kinetic / E_Hartree / E_xc / E_local / E_nonlocal / E_ewald holds
  the 9.5 eV. The Fe Ewald pin (`test_fe_bcc_ewald_vs_qe`, green) rules
  out ion-ion. NLCC regression guard (`test_fe_bcc_xc_nlcc_regression_guard`,
  green) rules out the NLCC E_xc double-counting path.
- Extend the audit to Cu (semicore 3s3p3d) — if the residual lives
  there too, the bug is shared.
- Ship the Rust-side fix (exact location TBD by Phase 5 diagnosis).

**Cost estimate.** VGCMP Phase 5 diagnosis: ~2–3 researcher-days (it is
a lifted VGC5 with the same per-component workflow). Ship-the-fix:
unknown until diagnosed, but bounded above at ~3–5 CE-days based on
NCFX's shape (which turned out to be ~150 LoC once located). VQEF test-
arm updates: ~1 CE-hour for 5 systems after VGCMP Phase 5 lands
(drop 5 `#[ignore]`s, tighten tolerances from 0.1 eV to ≤10 meV).

### 2c. No other known LDA blockers

Double-check against the file: `grep '#\[ignore' tests/qe_validation.rs`
returns exactly 7 matches (all accounted for above, plus the Fe nspin=2
test which appears under both SYKP/MPSH and VGCMP Phase 5 attribution
because it is the nspin=2 version of Fe; same root-cause). `cargo test
--test qe_validation` currently reports 8 passing tests (the two
defensive guards and the `compile` probes in the helpers), not 8
out of 16 production assertions.

## 3. Fe ferromagnetism

**Current state.** Under PseudoDojo NC/LDA at ecut = 15 Ry, Fe BCC a=2.87 Å
**collapses to non-magnetic** in both QE and pwdft-rs. The
`reference_data.toml` entry `fe_bcc_fm` pins `total_magnetization_mub =
0.00`. The user directive says "Fe (FM)" — so an NM-only test does not
satisfy the spec, even if the nspin=2 machinery is exercised.

This is a **pseudopotential/cutoff limitation**, not a pwdft-rs bug.
Experimentally Fe BCC is ferromagnetic with M ≈ 2.22 μB. LDA
systematically under-binds the FM state vs NM; PBE corrects this. The
choice of PP makes a large difference too.

### Three possible paths

**Path A — Higher ecut with same PP.** Rerun both codes at ecut ∈
{30, 40, 60} Ry; see whether the NM collapse lifts. Based on general
experience, PseudoDojo NC/LDA Fe is known to be marginal for FM at
any reasonable cutoff. Cost: ~1 researcher-day to sweep (QE runs take
minutes; pwdft-rs runs take O(hours) at 8×8×8 + ecut 60 Ry, 15 bands).
**Likelihood of fixing it: LOW** — published benchmarks typically need
PBE or a different PP family.

**Path B — Different pseudopotential.** Two candidates already in-tree:

- `pseudopotentials/nc/lda/Fe_dalcorso.upf` (dal Corso Ultra-Soft
  library, LDA). **Blocker:** we don't support USPP yet
  (`PseudopotentialData` parses NC only). USPP is a large separate
  track. REJECT.
- `pseudopotentials/nc/pbe/Fe.upf` (PseudoDojo NC/PBE) — this is the
  PBE PP needed for track §4 anyway. **This is the right Fe PP for
  FM validation**, but it requires GGAP Phase D to consume.

**Path C — Accept NM for LDA, validate FM under PBE only.** Keep the
existing Fe LDA test `test_fe_bcc_fm_vs_qe` as-is (nspin=2 machinery
validation, NM collapse documented). Add a new `test_fe_bcc_pbe_fm_vs_qe`
once GGAP Phase D lands, using PseudoDojo NC/PBE Fe at ecut=50 Ry,
8×8×8, `starting_magnetization = 0.7`. Target: M = 2.22 μB ± 0.05 μB
(GGAP Phase D's own acceptance criterion).

**Recommendation: Path C.** The Fe LDA test validates machinery; the
Fe PBE test validates FM physics. Both tests co-exist under the "Fe (FM)"
row of the matrix. Cost: absorbed into GGAP Phase D + VQEF integration;
no extra work beyond §4 for the PBE leg.

### Fe LDA test's eventual acceptance criterion

Even under Path C, the Fe LDA test should eventually drop `#[ignore]`
once VGCMP Phase 5 closes the 9.5 eV V_local gap. At that point:

- Assertion target: `|E_pwdft − E_QE| < 20 meV` (same tolerance as Si
  post-MPSH).
- Magnetization assertion: `M = 0.00 μB ± 0.01 μB` (documented NM
  collapse; **not** an FM target).
- Docstring explicitly states: "This test validates nspin=2 machinery
  on a system where the LDA PP drives NM collapse. The FM physics test
  is `test_fe_bcc_pbe_fm_vs_qe`."

## 4. PBE expansion track

### Dependency on GGAP

GGAP Phase A (in flight on PR `GGAP/phase-a-dispatcher`) is dispatcher
refactor only — it routes `XcFunctional::Pbe` to a `NotImplemented`
error and keeps LDA bit-exact. Phase A alone does **not** enable any
PBE validation; it is purely scaffolding to keep HYBR compatibility.

PBE validation requires:

- Phase B (PBE exchange non-spin) + Phase C (PBE correlation non-spin
  + PW92) — enables Si, C, Al, GaAs, NaCl, MgO, Cu PBE validation.
- Phase D (spin-polarized PBE) — enables Fe PBE FM validation.
- Phase F (QE validation suite expansion) — explicitly in GGAP's plan,
  lists Si + Al + Fe as reference arms. **VQEF extends this to all 8
  systems.**

### Pseudopotential library

`pseudopotentials/nc/pbe/` **exists and contains all 8 target elements**
(verified at proposal time: `{Al, As, C, Cl, Cu, Fe, Ga, Mg, Na, O, Si}.upf`
all present; 72 PPs total). No external downloads required. This is
the same statement GGAP makes in its "Pseudopotentials" section.

A parallel `qe_validation/pseudo-pbe/` symlink should be added pointing
at `pseudopotentials/nc/pbe/` so QE inputs can reference
`pseudo_dir = './pseudo-pbe'` without duplicating content.
The LDA `qe_validation/pseudo/` symlink (pointing at the LDA library)
stays.

### Per-system PBE reference generation

For each of the 8 systems, once GGAP Phase B+C (or Phase D for Fe FM)
lands, the VQEF task is to:

1. Create `qe_validation/pbe/<system>_pbe.in` cloning the LDA input but
   with `input_dft = 'PBE'`, the matched PBE pseudo filename, and
   possibly a higher ecut (see per-system table below).
2. Run QE via the `qe-runner` skill (must hold machine lock — see §5)
   and capture the raw `pw.x` output at `qe_validation/pbe/<system>_pbe.out`.
3. Extract total energy (Ry), Fermi energy (eV), magnetization (if
   nspin=2), Γ-eigenvalues.
4. Append a section to `qe_validation/reference_data.toml`:
   ```
   [si_diamond_pbe]
   input_file        = "pbe/si_pbe.in"
   ecutwfc_ry        = 30.0
   ...
   ```
5. Add `test_<system>_pbe_vs_qe` to `tests/qe_validation.rs` with a
   matching assertion harness. If the same crystal builder is reused,
   only the pseudo path, `ecut_ry`, `XcFunctional` setting, and the
   toml key differ.

### Per-system PBE table

Suggested ecuts follow PseudoDojo PBE recommendations and PBE PP
convergence norms (generally +5 Ry over LDA for the same element).
Cutoffs for Fe and Cu in particular need to be raised to capture PBE's
gradient-sensitive behavior near transition-metal d-bands.

| System  | LDA ecut (Ry) | PBE ecut (Ry) | k-grid | nspin | starting_mag | PBE convergence expected? | Dependency (besides GGAP Phase B+C)  |
|---------|---------------|---------------|--------|-------|--------------|---------------------------|--------------------------------------|
| Si      | 15            | 30            | 4×4×4  | 1     | —            | YES (insulator)           | MPSH (§2a) for test tolerance        |
| C       | 30            | 30            | 4×4×4  | 1     | —            | Likely (mixer sensitive)  | MPSH (§2a) + mixer tune              |
| Al      | 15            | 20            | 8×8×8  | 1     | —            | YES (simple metal)        | MPSH (§2a)                           |
| Fe      | 15            | 50            | 8×8×8  | 2     | 0.7          | FM expected at 50 Ry      | GGAP Phase D + VGCMP Phase 5         |
| GaAs    | 20            | 30            | 4×4×4  | 1     | —            | YES                       | VGCMP Phase 5                        |
| Cu      | 25            | 40            | 8×8×8  | 1     | —            | YES (metallic)            | VGCMP Phase 5                        |
| NaCl    | 25            | 35            | 4×4×4  | 1     | —            | YES (ionic)               | VGCMP Phase 5                        |
| MgO     | 30            | 40            | 4×4×4  | 1     | —            | YES (ionic wide-gap)      | VGCMP Phase 5                        |

### Cost estimate

- QE reference generation: **~0.5 researcher-day** per system (includes
  input file, pw.x run under machine lock, extraction, toml commit,
  pw.x output commit). 8 systems → ~4 researcher-days total. Can be
  parallelized across sessions; QE runs are minutes each.
- Rust test arm per system: **~0.5 CE-hour** (paste-and-edit from LDA
  arm, swap toml key + pseudo path + functional knob). 8 systems →
  ~4 CE-hours total.
- Per-component (VGC5-style) cross-check for 1–2 flagship systems (Si,
  Fe): **~1 researcher-day** total; catches any PBE sign/factor bug
  before the full suite is trusted.
- Defensive regression guards (PBE-specific, e.g. gradient-consistency
  on Si): **~1 researcher-day**.

**Total PBE expansion cost: ~6 researcher-days + 4 CE-hours.** This is
dwarfed by GGAP's ~9–13 CE-days for the physics.

## 5. Machine-lock policy (dependency, not owned here)

VQEF's QE runs require the shared machine lock per `CLAUDE.md` (§
"Machine Coordination" and pending proposal **QELK** — QE machine lock
policy, in flight on `QELK/qe-machine-lock-policy`). The `qe-runner`
skill has been updated (or will be, per QELK) to assume the lock is
held for any `pw.x` invocation.

VQEF's runbook (prose only, codified when QELK lands) should state
explicitly that:

1. Every QE reference generation session acquires the lock under role
   "Researcher".
2. Lock acquisition description should read
   `VQEF: QE reference <system>_<lda|pbe>` for traceability.
3. The lock is released between systems (not held across the whole 8
   × 2 grid) so that concurrent CE work isn't blocked for hours.

**No VQEF code change is required for QELK.** This is a coordination
note.

## 5a. Band-sum identity gate (BSUM, 2026-04-19)

VQEF's scalar-level gate (E_total + E_F) is a **global** agreement check:
E_total sums over kinetic + local + non-local + Hartree + XC + Ewald,
and its residual can mask opposite-sign component drift (the VGCH
"different converged density" signature). The band-sum identity gate
(proposal BSUM, landed 2026-04-19) adds a second scalar-level gate per
cell that is **orthogonal** to E_total on the cancellation axis.

**Physical quantity.** QE reports in every `pw.x` output a line labeled
`one-electron contribution = eband + deband Ry`
(`PW/src/electrons.f90:1719`). `eband = Σ w_k · f_{ik} · ε_{ik}` is the
raw band sum; `deband = -Σ <ψ|V_H + V_xc|ψ>` is the double-counting
correction that turns `<ψ|H|ψ>` (which carries V_H + V_xc) into
`<ψ|T + V_ion|ψ>` (the pure one-electron piece). Their sum
`<ψ|T + V_ion|ψ>` is independent of the V_loc(G=0) convention (the
rigid-shift tracked under VGCH Phase 1b cancels once `deband`'s shift
counterpart is added), and it is a sharper density-drift indicator than
E_total because it isolates the Hamiltonian-level one-electron piece
from the Hartree / XC / Ewald terms where cancellation happens.

**What BSUM catches that E_total misses.** If `|ΔE_one-electron|` is
small while `|ΔE_total|` is small, the two codes agree on both the
density **and** the Hartree/XC/Ewald bookkeeping. If `|ΔE_one-electron|`
is large while `|ΔE_total|` is small, E_total's smallness is an
**accidental cancellation** across Hartree/XC/Ewald masking a real
density-basin disagreement — exactly the VGCH signature (`Δ one-e =
+1.76 eV, Δ E_H = -0.59 eV, Δ E_xc = +0.29 eV` on C diamond). If
`|ΔE_one-electron|` is small while `|ΔE_total|` is large, the
disagreement lives in Hartree/XC/Ewald (functional or NLCC or Ewald
bug), not in the converged density.

**Implementation.**
- Helper `assert_band_sum_matches_qe(label, result, qe_one_electron_ry,
  tolerance_ev)` in `tests/qe_validation.rs` compares
  `e_kinetic + e_local + e_local_g0_shift + e_nonlocal` (pwdft-rs) vs
  QE's `one-electron contribution`.
- `QeComparisonConfig::one_electron_qe_ry: Option<f64>` (default `None`,
  additive — does not disturb concurrent-agent test bodies). When
  populated, `run_qe_comparison` prints the diagnostic residual
  unconditionally (so heavy-atom cells owned by other tracks also emit
  machine-parsable numbers under `cargo test -- --ignored --nocapture`).
- `qe_validation/reference_data.toml` gains a `one_electron_ry` key per
  cell (all 16 entries; values harvested from the existing `*.out`
  files, no QE re-runs required).
- BSUM assertion wired into 6 cells (Si LDA E, Si PBE, Al LDA, Al PBE,
  C LDA, C PBE). VGCH-2B / Si-EF-owned heavy-atom cells keep
  `one_electron_qe_ry: None` and no BSUM assertion until those tracks
  opt in — the diagnostic print is still emitted from
  `run_qe_comparison`.

**Tolerance policy.**
- GREEN cells (Si LDA/PBE, Al LDA/PBE): tolerance in the 40–120 meV
  band, chosen as ~2×|observed ΔE_1e| + margin so the gate catches
  real regressions (2× observed) while not being so tight it fires on
  routine compiler / numerics drift.
- VGCH-YELLOW cells (C LDA/PBE): tolerance = round-up of observed
  `|ΔE_1e|`. C LDA `|ΔE_1e| = 1.72 eV` → tol 2.0 eV. C PBE `|ΔE_1e| =
  0.37 eV` → tol 0.6 eV. Per the brief, don't tighten beyond what
  E_total already indicates.
- Heavy-atom cells (all Z>14 LDA + Z>14 PBE): no BSUM assertion in the
  BSUM-landing PR; the diagnostic is printed so VGCH-2B / Si-EF can
  harvest residuals as independent evidence for the density-drift
  hypothesis without BSUM owning those test bodies.

**Measured residuals (16-cell table at BSUM landing, post-TSEN, 2026-04-19):**

| Cell     | \|ΔE_total\| (eV) | \|ΔE_1e\| (eV) | Ratio | Notes                          |
|----------|-------------------|----------------|-------|--------------------------------|
| Si LDA   | 0.0448            | 0.0173         | 0.39  | GREEN; E_1e ≈ E_total          |
| Si PBE   | 0.0024            | 0.0127         | 5.29  | GREEN; sub-meV floor           |
| C LDA    | 1.4500            | 1.7217         | 1.19  | VGCH light-atom, pinned        |
| C PBE    | 0.3217            | 0.3746         | 1.16  | VGCH light-atom, pinned        |
| Al LDA   | 0.0259            | 0.0030         | 0.12  | GREEN; sub-meV basis noise     |
| Al PBE   | 0.0081            | 0.0008         | 0.10  | GREEN; sub-meV basis noise     |
| Fe LDA   | 11.14             | 10.99          | 0.99  | VGCH heavy-atom (no assert)    |
| Fe PBE   | 1.70              | 5.51           | 3.25  | VGCH heavy-atom (no assert)    |
| GaAs LDA | 35.17             | 53.10          | 1.51  | VGCH heavy-atom (no assert)    |
| GaAs PBE | 17.28             | 28.52          | 1.65  | VGCH heavy-atom (no assert)    |
| Cu LDA   | 16.64             | 36.87          | 2.22  | VGCH heavy-atom (no assert)    |
| Cu PBE   | 9.97              | 17.68          | 1.77  | VGCH heavy-atom (no assert)    |
| NaCl LDA | 7.99              | 14.97          | 1.87  | VGCH heavy-atom (no assert)    |
| NaCl PBE | 4.86              | 9.71           | 2.00  | VGCH heavy-atom (no assert)    |
| MgO LDA  | 10.71             | 18.10          | 1.69  | VGCH heavy-atom (no assert)    |
| MgO PBE  | 1.56              | 2.87           | 1.84  | VGCH heavy-atom (no assert)    |

**Independent evidence for VGCH-2B.** The heavy-atom ratios sit in
**1.5–3.3×**. For ratio > 1, E_total's residual is carried
disproportionately on the **non-one-electron** side (Hartree + XC +
Ewald) in a partial-cancellation pattern with E_1e moving the opposite
way. Fe PBE's 3.25× is the most extreme: E_1e residual is 3× larger
than E_total's residual, confirming that the converged density differs
between pwdft-rs and QE more than E_total alone reveals. This is
independent support for VGCH-2B's transplant hypothesis (density drift,
not a Hamiltonian assembly bug). Conversely, ratios close to 1 (Fe LDA
at 0.99, light-atom cells) indicate the two codes converge to similar
densities and differ mostly on the accounting side.

## 6. Acceptance criteria

VQEF is complete when **all** of the following hold on `main`:

### 6a. LDA track (8 systems, non-ignored)

- `cargo test --test qe_validation` with no `--ignored` flag runs
  exactly **8 LDA `test_<system>_vs_qe`** tests (plus the existing
  defensive guards) and all pass.
- Tolerances on the 8 production LDA tests:
  - Light systems (Si, C, Al): `|ΔE| ≤ 10 meV`, `|ΔE_F| ≤ 20 meV`.
  - Heavy/semicore (Fe, GaAs, Cu, NaCl, MgO): `|ΔE| ≤ 20 meV`,
    `|ΔE_F| ≤ 50 meV`. Tighter targets (≤10 meV) deferred to a
    follow-up once VGCMP Phase 5 residuals are characterized in full.
- Each test reports the per-component breakdown on `eprintln!` (VGC5
  style) so a CI diff against `reference_data.toml` is human-readable.

### 6b. PBE track (8 systems, non-ignored where each can converge)

- `cargo test --test qe_validation` runs exactly **8 PBE
  `test_<system>_pbe_vs_qe`** tests.
- Acceptance depends on convergence:
  - If a system converges under GGAP-Phase-B+C (Phase D for Fe), the
    PBE test is non-`#[ignore]` and uses the same tolerance band as
    LDA (10/20 meV).
  - If a system fails to converge under our PBE implementation after
    reasonable mixer/ecut exploration, the test is `#[ignore]`d with
    a **clear, linked follow-up proposal** (not just a string). The
    PR body of the proposal-that-ignores must name the blocker.
- PBE track includes an explicit Fe FM test
  `test_fe_bcc_pbe_fm_vs_qe` with magnetization target
  `M = 2.22 μB ± 0.05 μB` (sets the user's "Fe (FM)" directive clearly
  under the PBE leg).

### 6c. Reference data hygiene

- `qe_validation/reference_data.toml` contains 16 named sections
  (8 LDA + 8 PBE), each with provenance (`generated_date`,
  `qe_version = "7.5"`, `mpi_ranks`, host architecture).
- `qe_validation/pbe/` folder exists and contains 8 `*_pbe.in` +
  `*_pbe.out` pairs.
- `qe_validation/pseudo/` (LDA) and `qe_validation/pseudo-pbe/` (PBE)
  symlinks both present and pointing at the canonical
  `pseudopotentials/nc/{lda,pbe}/` directories.
- `qe_validation/README.md` updated to describe the 16-cell matrix
  and which systems use which PP family.

### 6d. What is explicitly NOT a VQEF acceptance criterion

- Bit-exact (sub-meV) agreement. The Γ-centered MP grid and QE's
  internal FFT-grid rounding introduce sub-meV residuals that VQEF
  is not chasing.
- Band-structure agreement beyond Γ eigenvalues. VQEF pins the SCF
  loop only; dispersion relations at arbitrary k are a separate test
  track.
- USPP or PAW validation. This suite is NC only.
- Non-LDA/non-PBE functionals (PBEsol, SCAN, hybrids). Those go in
  their own proposals (HYBR for hybrids, FLUP for PBEsol/revPBE).

## 7. Phase ordering and cost estimate

### Critical path

The critical path is governed by two hard dependencies:

1. **MPSH must land** before Si/C/Al LDA tests can be made to pass.
2. **VGCMP Phase 5 must land** before Fe/GaAs/Cu/NaCl/MgO LDA tests
   can be made to pass.
3. **GGAP Phase B+C must land** before any PBE test can run.
4. **GGAP Phase D must land** before Fe FM PBE can run.

Since MPSH is small (~1 CE-day) and VGCMP Phase 5 is the only
physics-diagnosis step currently blocking 5 of 8 LDA cells, the
critical path is:

```
   ┌── MPSH ────────────────────────────┐
   │   (1 CE-day; code change)           │
   │                                     ▼
   │                              VQEF §2a LDA drop 3 ignores (Si, C, Al)
   ▼                                     │
VGCMP Phase 5 ───────────────────────────┤
(2–3 res-days diag + 3–5 CE-days fix)    │
                                         ▼
                                  VQEF §2b LDA drop 5 ignores (Fe, GaAs, Cu, NaCl, MgO)
                                         │
                                         ▼
                                  *** LDA GATE GREEN ***  (8/8 LDA cells)
                                         │
                                         │  GGAP is independent and can
                                         │  run in parallel
                                         ▼
GGAP Phase A (in flight) ──── Phase B (1d) ── Phase C (2-3d) ── Phase D (2d)
                                                  │                  │
                                                  ▼                  ▼
                                     VQEF §4 PBE non-spin (7 cells)  VQEF §4 Fe FM PBE (1 cell)
                                                  │                  │
                                                  └────────┬─────────┘
                                                           ▼
                                                  *** PBE GATE GREEN *** (8/8 PBE cells)
```

### Recommended phase order

Phase ordering is written to **maximize the earliest moment the matrix
goes half-green** (all 8 LDA cells), because that is the highest-value
single milestone for the user directive.

| Phase   | What                                                   | Blocking proposal | Cost (days) | Cumulative matrix state |
|---------|--------------------------------------------------------|-------------------|-------------|-------------------------|
| **P0**  | File MPSH proposal                                     | — (VQEF spawns)   | 0.25 res    | 0/16                    |
| **P1**  | MPSH implementation                                    | MPSH              | 1 CE        | 0/16                    |
| **P2**  | VQEF §2a test arms (Si/C/Al LDA)                       | P1 done           | 0.5 CE      | 3/16                    |
| **P3**  | VGCMP Phase 5 diagnosis                                | VGCMP             | 2–3 res     | 3/16                    |
| **P4**  | VGCMP Phase 5 fix (TBD)                                | VGCMP             | 3–5 CE      | 3/16                    |
| **P5**  | VQEF §2b test arms (Fe, GaAs, Cu, NaCl, MgO LDA)       | P4 done           | 0.5 CE      | **8/16** (LDA green)    |
| **P6**  | GGAP Phase A (dispatcher) — already in flight          | GGAP              | —           | 8/16                    |
| **P7**  | GGAP Phase B+C (PBE exchange + correlation non-spin)   | GGAP              | 3–4 CE      | 8/16                    |
| **P8**  | VQEF §4 PBE reference generation (7 non-Fe systems)    | P7 done           | 3 res + 3 CE-hrs | 15/16              |
| **P9**  | GGAP Phase D (spin-polarized PBE)                      | GGAP              | 2 CE        | 15/16                   |
| **P10** | VQEF §4 Fe FM PBE reference + test                     | P9 done           | 1 res + 1 CE-hr  | **16/16** (gate green) |

**Total cost (VQEF-owned work only, excluding MPSH, VGCMP, GGAP):**

- QE reference generation: ~4 researcher-days
- Rust test arm updates (all 16 cells): ~1 CE-day
- Per-component VGC5-style cross-checks on flagship systems: ~1–2
  researcher-days
- Defensive PBE regression guards: ~1 researcher-day

**≈ 7 researcher-days + 1 CE-day, scheduled as VQEF P2, P5, P8, P10.**

**Total cost including MPSH + VGCMP Phase 5 + GGAP (upstream tracks
that need to land regardless):**

- MPSH: ~1 CE-day
- VGCMP Phase 5: ~5–8 CE/res-days combined
- GGAP Phases B+C+D: ~7–9 CE-days (Phase A already in flight)

**Grand total ≈ 20–25 working days of CE+researcher effort.** If two
agents work in parallel (one on MPSH/VGCMP Phase 5, one on GGAP Phases
B→D) with VQEF test-arm work absorbed at the tail of each track, this
compresses to ~3 calendar weeks.

### Why this order

- **MPSH before VGCMP Phase 5:** MPSH is 1 day and unblocks 3 of 8
  LDA cells. Even if VGCMP Phase 5 takes a month, the matrix visibly
  moves from 0/16 to 3/16 after day 2.
- **VGCMP Phase 5 before GGAP PBE validation:** the heavy-atom V_local
  residual will also bite PBE (same V_local form factor code path,
  just with a different XC functional). Closing it under LDA first
  means we can cleanly attribute any subsequent PBE residual to GGAP
  bugs, not lingering V_local issues. Running GGAP against 5 systems
  that are already known-broken at 9+ eV would waste validation cycles.
  Note that GGAP Phases A-D themselves are independent of VGCMP and
  can proceed in parallel; it is only the VQEF PBE test-arms at P8/P10
  that want VGCMP Phase 5 landed first.
- **GGAP Phase B+C (7 non-Fe PBE cells) before Phase D (Fe PBE FM):**
  Phase D is spin-polarized PBE, the highest-risk PBE code path.
  Debugging it is easier with Phase B+C already validated on
  insulators/simple metals.
- **Fe FM PBE last:** it is the single cell with the most compounded
  dependencies (GGAP Phase D + VGCMP Phase 5 + starting_magnetization
  tuning). Saving it for last means every other input knob is pinned
  when we debug it.

### Risk registry

| Risk                                                                          | Impact                                                                   | Mitigation                                                                                                                                                                |
|-------------------------------------------------------------------------------|--------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| VGCMP Phase 5 diagnosis reveals multiple compounding bugs, not one            | Phase 5 cost blows from 5 days to 10–15                                  | Start with VGC5-style audit on Fe; if component-1 diff is O(eV) we likely have one bug, if it's distributed across components it's multiple.                              |
| MPSH's Γ-centered grid exposes a different mixer sensitivity on C diamond     | C diamond PBE still doesn't converge after MPSH                          | Mixer tune separately; already flagged in §2a. C LDA is the canary.                                                                                                       |
| PseudoDojo PBE PPs disagree with PBE theory at high ecut                      | GGAP validates against QE but residuals remain vs literature             | Accept QE as ground truth (CLAUDE.md directive); literature comparison is out of scope.                                                                                   |
| Fe FM under PBE still collapses to NM at ecut=50                              | VQEF §3 Path C fails; user directive unmet                               | Fall back to PseudoDojo PBE Fe_sv (semicore-valence) PP if available; failing that, document NM collapse under PBE too and file a separate investigation proposal.         |
| QE 7.5 Γ-centered + nk odd introduces edge-case IBZ reductions we haven't hit | Reference data generation fails for one system                           | Pin reference sub-cases individually; don't regenerate the 8-system table in one batch.                                                                                   |

## 8. Deliverables summary

VQEF delivers, across P2/P5/P8/P10:

- **16 entries in `qe_validation/reference_data.toml`** (8 LDA + 8 PBE).
- **16 test functions in `tests/qe_validation.rs`**, 0 `#[ignore]`
  unless a system genuinely cannot converge under PBE (each such
  `#[ignore]` cites its follow-up proposal).
- **8 new `qe_validation/pbe/*.in` + `.out` file pairs** (LDA files
  already exist).
- **`qe_validation/pseudo-pbe/` symlink.**
- **Updated `qe_validation/README.md`** describing the 16-cell matrix.
- **Defensive VGC5-style per-component cross-check tests** on Si LDA,
  Si PBE, Fe LDA, Fe PBE (4 component audits) to catch future PP/XC
  regressions.

## References

- `proposals/completed/QEVL-qe-validation-test-suite.md` — original
  QEVL proposal that landed the Tier-1+2 LDA reference data.
- `proposals/completed/QEDX-qe-energy-discrepancy.md` — historical
  Si 13.4 eV investigation (closed by NCFX).
- `proposals/completed/SYKP-symmetry-ibz-audit.md` §D2 — MPSH
  recommendation (this proposal's Phase P1 dependency).
- `proposals/VGCMP-vloc-g-cross-check.md` — heavy-atom V_local(G)
  audit (this proposal's Phase P3/P4 dependency).
- `proposals/GGAP-gga-pbe-functional.md` — PBE functional
  implementation (this proposal's Phase P6/P7/P9 dependency).
- `proposals/HYBR-hybrid-functional-support.md` — hybrid functionals
  (out of scope for VQEF, referenced for downstream context).
- `tests/qe_validation.rs` — live test harness at `origin/main`
  0520c8c (2026-04-18).
- `qe_validation/reference_data.toml` — current LDA reference data.
- User directive, 2026-04-18: "We should be able to match QE values
  closely for Al, C, Cu, Fe (FM), GaAs, MgO, NaCl, and Si — with both
  LDA and PBE in the case that each can converge for that system."
