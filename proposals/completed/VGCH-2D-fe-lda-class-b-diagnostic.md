---
id: VGCH-2D
title: VGCH-2D — Fe LDA Class B diagnostic (PZ vs PW92 vs spin-driver bias)
status: completed
priority: high
complexity: small
risk: low
depends_on: [VGCH-MECH, VGCH-2]
blocks: [VQEF]
owner: researcher
---

## VGCH-2D — Fe LDA Class B diagnostic

### Context

VGCH-MECH (PR #168) fractured the heavy-atom residual into three
mechanism classes. **Class B is a one-cell class**: Fe LDA, singled
out by BSUM (PR #165) because its `|ΔE_1e| / |ΔE_total|` ratio is 0.99×
— band-sum and E_total move in lockstep — while every other Class A
cell (including Fe PBE) sits at 1.5–3.3×, the classic partial-
cancellation signature of "different converged density." Fe LDA's
11.5 eV E_total residual therefore does not look like a density-basin
gap; it looks like a Hamiltonian-level d-manifold shift.

Two suspects under VGCH-MECH § Class B:

1. **PZ-81 vs PW92 LDA correlation.** pwdft-rs `XcEvaluator::Pz` calls
   `perdew_zunger_correlation` and `pz_correlation_spin` — PZ-81
   (Perdew & Zunger, Phys. Rev. B 23, 5048, Eq. C1 + §III
   spin-interpolation). QE's LDA reference deck on every validation
   system reports `Exchange-correlation= SLA  PW   NOGX NOGC`
   (verified: `grep "Exchange-correlation=" qe_validation/*.out` —
   all 8 LDA cells, Si/Al/C/Cu/Fe/GaAs/NaCl/MgO). `PW` is QE's
   `qe-7.5/XClib/qe_dft_list.f90:48` label for PW92 (Perdew & Wang,
   Phys. Rev. B 45, 13244 (1992)). **QE LDA = Slater + PW92**;
   pwdft-rs LDA = Slater + PZ-81. This is a global mismatch on the
   LDA path, not Fe-specific — but its spin-polarized impact on
   transition metals is larger than its nonmagnetic impact on Si.
2. **Spin-driver CCMX code path bias.** `src/scf/driver_spin.rs`
   switches to the coupled-channel (ρ_total, m) basis at every mixer
   step (CCMX, PR landed 2026-04-18). Fe LDA runs nspin=2; Al / C /
   Si LDA all run nspin=1 and go through `driver.rs`. If CCMX
   introduces a systematic bias in the spin-polarized Hartree or XC
   assembly that the nspin=1 path doesn't hit, Fe LDA would pick it
   up and Al / C / Si LDA would not.

This proposal designs a single diagnostic that separates the two. It
is diagnostic-only, ~1 CE-day, no code changes required.

### Key additional facts from validation data

- **QE Fe LDA collapses to M = 0** (verified:
  `qe_validation/fe_bcc_fm_scf.out`, last iteration
  `total magnetization = 0.00 Bohr mag/cell` with
  `starting_magnetization(1) = 0.5`). PZ and PW92 both underbind Fe
  ferromagnetism; both correctly predict a NM ground state at LDA.
  The **QE reference we compare against is spin-polarized-collapsed-
  to-NM at a converged nspin=2 density**, not a ferromagnetic LDA
  state.
- **QE Fe PBE stays ferromagnetic** at M = 2.34 μB
  (`fe_bcc_fm_scf_pbe.out`). The LDA→PBE transition on Fe therefore
  changes *both* the spin state (NM → FM) *and* the correlation
  functional (PW92 → PBE), which is why Fe PBE sits in Class A
  (different-density class) while Fe LDA sits in Class B (Hamiltonian-
  side class).
- **Everywhere else, the QE LDA validation deck uses SLA+PW92.** So
  if PZ vs PW92 is a large driver, we should already be seeing
  residual on Si / Al / C LDA too. Current numbers: Si LDA Δone-e =
  +17 meV, ΔE_xc = −59 meV; Al LDA +75 meV total; C LDA +1.45 eV
  total (Class C). The size of the LDA-only PZ↔PW92 disagreement
  on *insulators and simple metals* is therefore ≲ 1 eV. On Fe the
  disagreement is +11.5 eV — either PZ↔PW92 is anomalously large on
  transition metals, or it's riding on top of a spin-driver bug.

### Hypotheses

#### H-2D-A: PZ-vs-PW92 is the dominant driver

**Prediction if true.** Running pwdft-rs Fe with a PW92 LDA
correlation (replacing PZ) closes Fe LDA to a residual comparable to
Class A's Fe-PBE bucket (1.97 eV, VGCH-class). The residual would
move to ΔE_xc at shared density:

- Signature: at iter 1 transplant with ρ = ρ_QE, `|ΔE_xc| ≫
  |Δone-e|, |ΔE_H|, |ΔE_NL|`. The one-electron term should match QE
  to the same +0.26 eV/band Cu gauge tolerance (VGCH-2B) because the
  Hartree and non-local pieces are functional-independent; only the
  XC double-counting and the V_xc-on-diagonal piece shift.
- Magnitude prediction: at Fe's NM density, switching PZ→PW92 on
  unpolarized correlation shifts ε_c by ~5 mHa / electron
  (Ortiz-Ballone PRB 50, 1391, Table I at r_s ≈ 2–3). For a valence
  of 8 e/atom and a cell volume of 10.6 Å³ with ρ_val concentrated
  in the 3d manifold, the integrated ΔE_c is O(1 eV/atom). Fe's 11.5
  eV residual is 8–10× larger — **so PZ-vs-PW92 alone cannot account
  for Fe's residual** unless the spin-polarized PZ extrapolation to
  fully-polarized Fe adds substantially more. But QE reports M = 0 on
  Fe LDA, so the spin interpolation is barely exercised. Therefore
  **H-2D-A's quantitative prediction is a residual drop of O(1 eV),
  not 11 eV**. If we see an 11 eV drop, H-2D-A explains it; if we
  see a ≤ 2 eV drop, H-2D-A is partial and H-2D-B carries the rest.

**Falsifier:** if the Fe-LDA ΔE_xc-at-shared-ρ residual is < 500
meV and the Δone-e + ΔE_H + ΔE_NL residual is > 2 eV, PZ-vs-PW92 is
not the dominant driver.

#### H-2D-B: Spin-driver (CCMX / driver_spin.rs) introduces a Hamiltonian-side bias on nspin=2

**Prediction if true.** Running pwdft-rs Fe with nspin=1 (forced
non-magnetic) produces a residual much smaller than nspin=2. The nspin=1
path goes through `driver.rs` and uses the un-polarized XC; the nspin=2
path goes through `driver_spin.rs` + CCMX + `pz_correlation_spin` /
`assemble_v_eff_spin`. QE's M-collapsed Fe LDA is effectively an
nspin=1 state on an nspin=2 code path, which we can reproduce either
as pwdft-rs nspin=1 (forced NM) or pwdft-rs nspin=2 with starting
magnetization 0 (CCMX still active, but at ζ=0 throughout).

- Signature: residual appears only in the nspin=2 runs, not nspin=1.
- Specifically, Δone-e / Δρ at shared density in nspin=1 should be
  at or below Cu/GaAs's Class A magnitudes (Δone-e ~ O(10 eV)
  at ρ_QE transplant, carried over from Class A), while nspin=2
  adds an extra O(11 eV) of residual. The delta between the two
  paths isolates the spin-driver bias.

**Falsifier:** if nspin=1 Fe LDA and nspin=2 Fe LDA give residuals
within ≤ 500 meV of each other, the spin driver is NOT the dominant
driver.

#### H-2D-C (null): Neither — escalate

**Prediction if true.** All four (PZ, PW92) × (nspin=1, nspin=2)
runs give the same Fe LDA residual (±500 meV). Fe LDA's 11 eV
residual is then neither a correlation-functional choice nor a
spin-driver bug; it's a deeper Hamiltonian-construction issue shared
with the Class A cells but with a different fingerprint.

In that case, escalate to:

- Extending the Cu G=G' `D_ij · Σ β·β` cross-check (Class A H-C5) to
  Fe's d-projectors at shared density;
- Auditing the Fe PP's `PP_RHOATOM` SAD seed for consistency with
  its UPF-reported total valence;
- Checking whether Fe's NLCC core charge (`PP_NLCC`) has a
  different radial grid convention between pwdft-rs and QE. Fe's
  NLCC is already regression-pinned by
  `test_fe_bcc_xc_nlcc_regression_guard` on the PZ path; the new
  pin would have to re-run that guard on PW92.

### Experimental design

Reuse the VGCH-2B transplant infrastructure on Fe with a small
scope extension to cover all four PZ/PW92 × nspin combinations.
The key change from VGCH-2B: `src/scf/transplant.rs` currently hard-
errors on `nspin != 1`. The Fe diagnostic needs either
(a) an nspin=2 extension to `transplant.rs`
(`run_scf_iter1_from_rho_g_fft_spin`, mirroring `driver_spin.rs`'s
iter-0 step and taking `(rho_g_total, rho_g_mag)`), or
(b) running Fe with nspin=1 (forced NM) through the existing
transplant and comparing against QE's M-collapsed Fe LDA.

Path (b) is strictly smaller and is what I propose for session 1.
If session 1 yields the H-2D-A 11 eV drop, we're done. If it yields
a ≤ 2 eV drop, session 2 lands path (a) and runs the nspin=2 leg.

#### Matrix

Run all four configurations and parse per-term energies:

| # | functional path | nspin | driver | CCMX active |
|---|-----------------|-------|--------|-------------|
| 1 | PZ (status quo) | 1     | `driver.rs`        | no   |
| 2 | PZ (status quo) | 2     | `driver_spin.rs`   | yes  |
| 3 | PW92 (new)      | 1     | `driver.rs`        | no   |
| 4 | PW92 (new)      | 2     | `driver_spin.rs`   | yes  |

Run each against the same QE reference (M-collapsed Fe LDA
`fe_bcc_fm_scf.out`, which is already SLA+PW92, nspin=2 collapsed to
ζ=0). For the transplant at iter 1, the input density ρ_in is the
QE converged ρ (M=0, so ρ_up = ρ_down = ρ_total/2 — identical in
nspin=1 and nspin=2 at the transplant step).

**Infrastructure delta needed.**

- Add an `XcEvaluator::PzOrPw92` knob (or a sibling
  `XcFunctional::LdaPw92` value) so the diagnostic can switch LDA
  correlation without touching the PBE path. **Not a production
  feature** — gate the variant behind `#[doc(hidden)]` or
  `#[cfg(feature = "diagnostic")]` and make clear in the docstring
  that the supported LDA path stays PZ until a separate feature
  proposal (post-VGCH-2D) promotes PW92 to the default.
- Add an nspin=2 variant of `run_scf_iter1_from_rho_g_fft` that
  parses QE's spin-resolved charge-density.dat. The Python parser
  (`scripts/validate/vgch2_parse_qe_density.py`) already reads
  nspin=2; it writes `ρ_g[ngm, 2]`. The Rust side needs to consume
  that and initialize `ctx.rho_up / rho_down` (or `ρ_total / m` in
  CCMX units) separately.

Both pieces are small (<200 LOC) and self-contained to the
diagnostic surface. Gate both behind `#[doc(hidden)]` so production
users never see them.

#### Measured quantities (per configuration)

For each of the four runs, record at iter 1 transplanted from ρ_QE:

- `E_HF` (Harris-Foulkes, paired with ρ_in = ρ_QE)
- `E_total` (same as ScfResult.total_energy, paired with ρ_out from
  diagonalizing H[ρ_QE])
- Per-term decomposition: `E_kinetic`, `E_local`, `E_NL`, `E_H`,
  `E_xc`, `E_vxc`, `E_ewald`, `E_smearing`. All at ρ_out (consistent
  with VGCH-2B's table).
- `Δone-e = E_kinetic + E_local + E_NL − QE's one-electron`
- `ΔE_H`, `ΔE_xc`, `ΔE_ewald`: same as VGCH-2B.
- Fermi level, entropy `TS`, per-k Γ-point eigenvalue table.
- `Δρ(in, out) RMS` on the FFT grid.

#### Expected outcome interpretation table

| Measurement                       | H-2D-A dominant  | H-2D-B dominant      | Null    |
|-----------------------------------|------------------|----------------------|---------|
| PZ nspin=1 residual               | ≈ 11 eV         | ≈ 11 eV              | ≈ 11 eV |
| PZ nspin=2 residual               | ≈ 11 eV         | ≈ 11 eV              | ≈ 11 eV |
| PW92 nspin=1 residual             | ≤ 2 eV          | ≈ 11 eV              | ≈ 11 eV |
| PW92 nspin=2 residual             | ≤ 2 eV          | ≤ 2 eV if driver bias is only in ρ-reconstruction; ≈ 11 eV if in V_H / V_xc | ≈ 11 eV |
| Dominant Δ-per-term (PZ nspin=2)  | ΔE_xc           | Δone-e + ΔE_H (CCMX mixes via V_H) | (mixed) |
| `Δ(nspin=1 − nspin=2)` at same XC | O(0.1 eV) k-grid | > 5 eV                | O(0.1 eV) |
| Fermi M = 0 in pwdft-rs?          | yes (PZ and PW92 both predict NM Fe) | yes     | yes     |

If the **PW92 nspin=1 row drops below 2 eV**, H-2D-A is confirmed as
the dominant single driver → escalate fix to VGCH-2D-F (flip LDA to
SLA+PW92 as the default, port spin-polarized PW92 to the LDA path —
pwdft-rs already has `pw92_correlation_spin_au` in
`src/potential/xc.rs:1111` because PBE needed it). Expected fix
scope: ~100 LOC rewire in `XcEvaluator::Pz` + sibling
`XcEvaluator::PwLda` that reuses PBE's PW92 helpers.

If the **nspin=1 − nspin=2 gap at fixed XC is > 5 eV**, H-2D-B is
confirmed → escalate fix to a spin-driver audit:

- CCMX basis-change rounding error on ρ_up ↔ (ρ_total + m)/2 and
  back;
- `assemble_v_eff_spin` vs `assemble_v_eff` dispatch — do they route
  through the same `v_local_fft` (they should; if not, bug);
- Hartree energy assembly in `driver_spin.rs` vs `driver.rs`.

If the **null row holds** (all four configurations give ≈ 11 eV), Fe
LDA is not a correlation-functional or spin-driver issue. Escalate
to deeper Hamiltonian-side audit — see H-2D-C.

### Falsifiable predictions (summary for the logbook)

1. PW92 nspin=1 Fe iter-1 transplanted E_HF residual **< 2 eV**
   ⇒ H-2D-A confirmed.
2. PW92 nspin=1 Fe iter-1 transplanted E_HF residual **in [2, 10] eV**
   ⇒ H-2D-A partial, combined with either H-2D-B or H-2D-C.
3. PW92 nspin=1 Fe iter-1 transplanted E_HF residual **> 10 eV**
   AND nspin=1 vs nspin=2 delta **> 5 eV** (at either XC choice)
   ⇒ H-2D-B confirmed.
4. All four configurations residual **≈ 11 eV** within 1 eV
   ⇒ null; escalate to H-2D-C investigation.

### Deliverables

- Extend `src/scf/transplant.rs` with
  `run_scf_iter1_from_rho_g_fft_spin` (nspin=2 variant) — gated
  `#[doc(hidden)]`. **Implementation out of scope for this
  proposal**; but the Researcher should pin the required interface
  (ρ_total_g, m_g, initial magnetization seed, return
  `TransplantIter1SpinResult` with per-spin eigenvalues) so the
  Core Engineer can land it in a small follow-up PR.
- Add `XcFunctional::LdaPw92` as a **diagnostic-only** variant
  (gated `#[doc(hidden)]` / hidden from the settings parser). Reuses
  `pw92_correlation` and `pw92_correlation_spin_au` already in
  `src/potential/xc.rs`. No PBE path change. Implementation also
  out of scope for this proposal — spec only.
- `scripts/validate/vgch2d_fe_lda_transplant.py` — driver script that
  runs the four configurations, parses the Rust test output, and
  emits `vgch2d_fe_lda.csv` with per-term deltas for each of the
  four runs.
- `tests/vgch_transplant_fe.rs` — Tier-2 (`#[ignore]`) test with four
  sub-cases. Prints per-term deltas; does not assert (diagnostic-
  only, like `tests/vgch_transplant_cu.rs`).

### Acceptance

Close VGCH-2D when:

1. All four configurations have been run and per-term deltas
   recorded in `vgch2d_fe_lda.csv`.
2. One of the three outcome branches (H-2D-A, H-2D-B, H-2D-C) is
   confirmed by the falsifiers above.
3. A follow-up proposal (VGCH-2D-F fix, OR VGCH-2E deeper audit)
   has been drafted with an explicit scope, based on the branch
   confirmed.

### Cost

- **Diagnostic run + analysis:** ~1 CE-day.
- **Follow-up fix (conditional on H-2D-A):** ~1 CE-day — rewire
  LDA default to SLA+PW92, reusing PBE's `pw92_correlation` and
  `pw92_correlation_spin_au`. Rerun all 8 LDA QE-validation cells;
  expect small shifts on every cell (Si, Al, C, Cu, GaAs, NaCl,
  MgO) — some may close further, one or two may open slightly.
  **Important:** this fix also affects `test_fe_bcc_xc_nlcc_
  regression_guard` — the PZ-based pin would need to be replaced
  with a PW92-based pin or made functional-parametric.
- **Follow-up fix (conditional on H-2D-B):** 2–3 CE-days — audit
  CCMX basis-change + `assemble_v_eff_spin` + `driver_spin.rs`
  Hartree assembly against `driver.rs` side-by-side. No scope
  estimate beyond that without the evidence.
- **Follow-up (conditional on H-2D-C null):** open as a new proposal;
  scope 3–5 CE-days.

### Non-goals

- Does NOT fix Fe LDA in this proposal. Diagnostic-only.
- Does NOT modify the production `XcEvaluator::Pz` variant or the
  default XC on any LDA run. Everything behind `#[doc(hidden)]`.
- Does NOT run nspin=2 transplant infrastructure — that's session 2
  if session 1 doesn't close the question.
- Does NOT touch Fe PBE — it's Class A and handled by VGCH-2 Part C.
- Does NOT touch the other Class B candidates (none exist today;
  Fe LDA is the sole Class B cell per BSUM).

### Risks

- **PW92 spin-polarized helpers already exist but aren't tested at
  the Fe LSDA density regime.** `pw92_correlation_spin_au` lives in
  `src/potential/xc.rs:1111`, was added as part of GGAP Phase D (PBE
  correlation spin). It's pinned at ζ=0 (matches unpolarized PW92)
  and at ζ=1 (matches polarized branch) but has no pin at
  intermediate ζ on a transition-metal r_s ≈ 1.5–2.5 point. If the
  PW92 LSDA at Fe densities has a bug in pwdft-rs (unlikely — it's
  a line-for-line port of QE 7.5's `pw_spin` per the docstring),
  this diagnostic would mis-attribute the result to H-2D-B. Mitigation:
  add a canonical-point regression test on
  `pw92_correlation_spin_au(rs=1.8, zeta=0.3)` comparing against a
  hand-evaluated QE `pw_spin` reference, before running the Fe
  diagnostic. ~30 minutes of work; catches the silent-port-bug case.
- **QE Fe LDA may not be a "true" Kohn-Sham minimum.** Starting M =
  0.5 collapsed to M = 0 means the SCF sampled the spin sector
  briefly and fell out. If pwdft-rs starts from the same seed but
  keeps M ≠ 0 at convergence (different mixer dynamics), the two
  are converging to different spin states and the 11 eV delta is
  mostly spin-state-difference. Mitigation: log pwdft-rs's
  total_magnetization per iteration for the nspin=2 Fe LDA baseline;
  if it's not at 0 by iter 30, we're comparing apples and oranges
  and the diagnostic must hold magnetization to 0 (fixed-M run).
- **The diagnostic-only `XcFunctional::LdaPw92` variant leaks into
  production if not properly gated.** Mitigation: cfg-gate +
  `#[doc(hidden)]` + Code Reviewer sign-off before the follow-up
  proposal lands.

### Provenance

- VGCH-MECH Class B hypothesis (proposals/VGCH-MECH-mechanism-
  taxonomy.md:152–170).
- VGCH-2 Part B transplant infrastructure
  (`src/scf/transplant.rs`, `scripts/validate/vgch2_parse_qe_density.py`).
- BSUM Fe LDA 0.99× ratio signature (PR #165).
- QE LDA dispatch ground truth: `qe-7.5/XClib/qe_dft_list.f90:48`
  (`'PZ'`, `'PW'` labels) + `qe-7.5/XClib/qe_dft_list.f90:72–78`
  (`dft_full(1) = PZ, dft_full(2) = PW`). Every pwdft-rs QE
  validation deck in `qe_validation/*_scf.out` uses SLA+PW (PW92),
  verified by `grep "Exchange-correlation=" qe_validation/*.out`.
- Fe LDA magnetic collapse: `qe_validation/fe_bcc_fm_scf.out` final
  iteration reports `total magnetization = 0.00 Bohr mag/cell`.
- pwdft-rs LDA dispatch: `src/potential/xc.rs:143`
  (`perdew_zunger_correlation(rho)` in `lda_xc`) +
  `src/potential/xc.rs:394` (`pz_correlation_spin(rho_up, rho_down)`
  in `lda_xc_spin`). `XcEvaluator::Pz` calls both verbatim.
- PW92 code path already present: `src/potential/xc.rs:894`
  (`pw92_correlation`, test-gated), `src/potential/xc.rs:1111`
  (`pw92_correlation_spin_au`, used internally by PBE). Both
  bit-matched to QE's `pw` / `pw_spin` per their docstrings and
  GGAP Phase D pin tests.
