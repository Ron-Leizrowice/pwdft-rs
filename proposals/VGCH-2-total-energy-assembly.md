---
id: VGCH-2
title: VGCH follow-up — total-energy assembly (V_loc(G=0) compensation + Harris-Foulkes pairings)
status: active
priority: high
complexity: medium
risk: medium
depends_on: [VGCH]
blocks: [VQEF]
owner: researcher
---

# VGCH-2 — Total-energy assembly for the heavy-atom residual

## Context

Both hypotheses under VGCH Phase 1 are now cleared:

- **H1 — β_l(q) form factors.** `tests/vgch_beta_l_heavy.rs` +
  `scripts/validate/vgch_beta_l_heavy.py` pin KB projectors for
  11 elements × 10 q-values = 590 rows to max |Δ| = 3.17×10⁻¹²
  Bohr^{3/2}, 4 orders below the 1e-8 tolerance. Fe/Cu semicore d
  projectors included. **Not the bug.**
- **H2 — SAD initial density.** `tests/vgch_sad_heavy.rs` +
  `scripts/validate/vgch_sad_heavy.py` pin pwdft-rs'
  `generate_initial_density` against a QE-convention Python
  reference for all 7 VGCH systems; raw-sample point-wise |Δρ| ≤
  1.2×10⁻⁵ e/Å³ (GaAs outlier from the negative-density clamp,
  O(1e-5 eV) total-energy effect). **C diamond is bit-perfect
  pre-clamp and post-clamp.** Not the bug either.

The 7–34 eV residuals must therefore live in **one of two places**:

1. **Total-energy assembly** — the `total_energy` + `with_g0_shift`
   + `harris_foulkes_energy` closure in `src/scf/energy.rs`, or the
   interaction between the G=0 compensation (`v_local_g0 · n_electrons`)
   and the NLCC double-counting subtraction, or a subtle mis-pairing
   of input-vs-output densities in the double-counting terms for
   multi-species cells.
2. **SCF mixer basin** — the mixer could be stabilizing a different
   local minimum of the energy functional.

Phase 1a's observation strongly points at (1):

- Cu one-electron is +37 eV too high, Hartree is −24.6 eV too low,
  they partially cancel to +17.3 eV net. The *signature* of a
  different converged density — but now we know the *initial* density
  is bit-perfect (Phase 1c).
- Si's one-electron vs. Hartree cancellation closes within
  `MPSH + SYKP` noise (≤ 50 meV). So the assembly works for Si.
  **Something in the assembly fails on Z > 14 cells.**
- On C diamond, every Γ eigenvalue is uniformly shifted by −3.16 eV
  vs. QE, exactly matching `2 · V_loc(G=0) = 3.09 eV` — see
  researcher logbook 2026-04-19 H1 entry. Band *gaps* agree to <50
  meV. The shift is absolute-reference-only and should be cancelled
  by `e_local_g0_shift = v_local_g0 · n_electrons`. Si cancels
  cleanly; C doesn't. That's a narrow fingerprint on the V_loc(G=0)
  compensation branch.

## Scope

**Part A — Isolate the mispairing (1 CE-day). COMPLETE 2026-04-19.**

Add a per-term trace to `src/scf/energy.rs` so we can cross-check
each assembled term against QE's per-term output (`E_one_electron`,
`E_hartree`, `E_xc`, `E_ewald`) using pwdft-rs' input-vs-output
density state at the moment each term is computed. Focus on:

- `e_local_g0_shift = v_local_g0 · n_electrons`. For Si this is +8
  electrons × small number ≈ small. For Cu this is +19 electrons ×
  ~5 eV. Does the sign match QE's convention?
- `e_xc` double-counting: the `(e_xc − e_vxc)` piece. NLCC makes
  this subtle — the core-density contribution to e_xc must *not*
  be subtracted from e_vxc. Is the sign right on all NLCC cells?
- Harris-Foulkes pairing: `E_HF = E_band − E_H[ρ_in] + (E_xc[ρ_in]
  − E_vxc[ρ_in]) + E_ewald`. The driver uses `ρ_in` for the double
  counting; does it actually pass `ρ_in` and not `ρ_out`?

### Part A findings (2026-04-19)

Artifacts:

- `scripts/validate/vgch2_per_term_trace.py` — parses QE
  `one_electron / hartree / xc / ewald / smearing_mts / total / fermi`
  from `qe_validation/*.out` for 8 systems → `vgch2_per_term_trace.csv`.
- `scripts/validate/vgch2_join_trace.py` — joins QE CSV with pwdft-rs
  CSV emitted by the extended `tests/vgch_per_component_heavy.rs`,
  computes `delta_meV = pwdft − QE`, ranks by |delta|.
- `tests/vgch_per_component_heavy.rs` extended from 2 → 8 cases
  (Si, C, Al, Fe, Cu, GaAs, NaCl, MgO) at the QE-reference SCF config
  where tractable. All 8 tests pass the PCFX self-check
  (`|Σ(components) − E_total|` < 33 meV on Fe, < 1 µeV elsewhere).

**Per-term delta table (eV, ours − QE):**

| system    | Δone-e   | ΔE_H     | ΔE_xc    | ΔE_ewald | ΔE_total |
|-----------|----------|----------|----------|----------|----------|
| Si        | +0.017   | −0.015   | −0.059   | +0.011   | −0.033   |
| Al        | −0.003   | −0.000   | −0.030   | +0.007   | +0.075   |
| C diamond | +1.722   | −0.586   | +0.327   | −0.013   | +1.450   |
| Fe BCC    | +10.986  | −0.702   | +0.931   | +0.006   | +11.502  |
| Cu FCC    | +37.033  | −24.568  | +4.728   | +0.000   | +17.308  |
| NaCl      | +14.970  | −8.934   | +1.852   | +0.100   | +7.988   |
| MgO       | +18.099  | −9.895   | +2.353   | +0.150   | +10.706  |
| GaAs      | +53.095  | −23.545  | +5.692   | −0.070   | +35.242  |

**Primary suspect: NOT an assembly mispairing.** The per-term
fingerprint is identical across every heavy-atom cell and is the
unambiguous signature of converging to a **different self-consistent
density**, not of a term-assembly bug:

1. **The direct-sum identity holds.** `|Σ components − E_total| <
   1 µeV` on every cell (the one exception is Fe at 23 meV,
   attributable to the nspin=1 mismatch with QE's nspin=2 ref —
   ρ-symmetrization of a magnetization-zero density is exact but
   the post-symmetrize occupation weight has numerical noise at
   1e-5 of the band sum). If `total_energy` / `with_g0_shift` /
   `harris_foulkes_energy` were double-counting a term, this sum
   would diverge by the same O(eV) as the QE residual — it doesn't.

2. **`e_local_g0_shift = N_el · Σ_sp V_loc(G=0)(sp)` is bit-correct
   against the direct formula.** Cross-checked against
   `scripts/validate/vgch_vloc_heavy.csv` to < 0.001 eV on all 8
   systems. The G=0 compensation path is not the bug.

3. **The Δone-e / −ΔE_H ratio tracks linear-response from a density
   perturbation**, not the factor-of-2 double-counting signature of
   an assembly bug:

   | system    | Δone-e    | −ΔE_H     | ratio  |
   |-----------|-----------|-----------|--------|
   | Si        | +0.017    | +0.015    | 1.19   |
   | C diamond | +1.722    | +0.586    | 2.94   |
   | Fe BCC    | +10.986   | +0.702    | 15.64  |
   | Cu FCC    | +37.033   | +24.568   | 1.51   |
   | NaCl      | +14.970   | +8.934    | 1.68   |
   | MgO       | +18.099   | +9.895    | 1.83   |
   | GaAs      | +53.095   | +23.545   | 2.26   |

   For a δρ that changes at fixed V_ext, linear response gives
   `Δone-e ≈ 2·∫V_H[ρ]δρ = 2·ΔE_H`, so ratio ≈ 2. Si/Cu/NaCl/MgO/GaAs
   sit in [1.2, 2.9]; C is 2.94; Fe at 15.6 is the outlier. A
   double-counting bug in `total_energy` would give ratio = 1 or
   ratio = ∞ (sign-dependent), not a smooth band. **This is a
   *different ρ* signature, not an assembly bug.**

4. **Δone-e + ΔE_H + Δxc + Δewald ≈ ΔE_total to within `−(−TS)`.**
   The residual after summing is exactly the QE `smearing contrib.
   (-TS)` term for Al/Fe/Cu/GaAs (100-260 meV), which means pwdft-rs'
   reported `total_energy` does NOT include the `−TS` smearing
   contribution (QE's `!    total energy` is F = E − TS). This is a
   separate ≤ 260 meV effect that should be flagged as its own
   follow-up but is NOT the cause of the 1.5-35 eV residuals.

5. **Si structural match: YES, with a twist.** On Si the |delta_meV|
   ranking is (ΔE_xc, Δone-e, ΔE_H) at (−59, +17, −15) meV. On every
   heavy system it is (Δone-e, ΔE_H, ΔE_xc) at O(1-53) eV. **The
   sign pattern is the same — Δone-e > 0, ΔE_H < 0, ΔE_xc > 0 — but
   the magnitude diverges by 3-4 orders.** Si's residual is O(k-grid
   noise + NLCC round-off + MP-shift convention); heavy systems
   have a *density-level* disagreement on top of those light
   effects. The assembly is bit-correct on both.

**Part B scope revision.** Given findings 1–4 above, the transplant
experiment moves from "only if Part A clears" to **primary next
step**. Seed pwdft-rs' SCF from QE's converged density (via a new
parser for QE's XML charge-density.xml or `.save/charge-density.dat`
binary). At iteration 1 with `ρ_in = ρ_QE`:

- If one-e / Hartree / xc match to < 50 meV/atom, the driver correctly
  reproduces QE's decomposition at QE's fixed point. The bug is then
  in the SCF *dynamics* — H3b (mixer-basin), H5 (symmetrization, but
  pinned clean by MPSH), H6 (initial density beyond SAD — already
  bit-perfect per H2), or H7 (eigensolver drift) — with H3b as the
  leading candidate because every failing cell uses a Kerker-family
  mixer and the `E_H` sign pattern is consistent with charge
  sloshing.
- If they don't match at iter 1, the bug is a subtle `v_xc` or `v_H`
  assembly term that is only visible when the density structure is
  "heavy" (semicore PP, multi-species, or compact valence overlap).

**Part B code pointers (post-revision):**

- `src/scf/driver.rs:593-595` — `total_energy(e_band, e_H, e_xc_corr,
  e_ewald)` followed by `with_g0_shift`. The pairing looks correct
  (same `rho_out` everywhere); add a sanity print of all 4 components
  at iter 1 for the seeded-density test.
- `src/scf/energy.rs:133-151` — `xc_energy_corrected`. Verify
  `rho_xc = rho_val + rho_core` and that `e_vxc` subtrahend
  integrates against `rho_val` not `rho_xc`. Already checked at
  docstring level — add an NLCC-specific unit test under Part B.
- `src/scf/context.rs:124` — `v_local_g0 = v_local_fft[0].re`. This
  is a sum over species (`sum_over_species(v_loc_of_g=0)`). On
  single-species cells this is fine; on GaAs/NaCl/MgO the sum is
  verified bit-correct above. Not the bug.

**What VGCH-2 Part A rules out:**

- Assembly mispairing in `total_energy` / `harris_foulkes_energy` /
  `with_g0_shift` (ratios, sum-identity, explicit PCFX check all
  pass).
- V_loc(G=0) compensation sign or magnitude (matches closed-form
  formula `N_el · Σ_sp V_loc(G=0)`).
- NLCC double-counting direction (e_xc − e_vxc has the right sign
  on all 5 NLCC-active cells).

**What VGCH-2 Part A does NOT rule out:**

- ρ-level disagreement at self-consistency (the primary finding).
- The missing `−TS` in pwdft-rs' reported `total_energy` (up to
  260 meV on Fe, up to 100 meV on Al/Cu — separate follow-up; see
  `src/scf/driver.rs:593-595` where `total_energy` is assembled
  without a smearing contribution).
- Double-counting in the LSDA driver for Fe specifically
  (Δone-e/−ΔE_H = 15.6 is a distinct outlier; rerun with nspin=2
  as a cross-check under Part B).

**Part B — Transplant experiment (1 CE-day, only if Part A clears).**

Only run this if Part A pins the assembly to bit-correct. Seed
pwdft-rs' SCF from QE's converged density (via a new parser for
QE's `save/charge-density.dat` binary or by regenerating QE with
`disk_io = 'high'`). Run 1 SCF iteration. Check:

- At iter 1, the one-electron and Hartree components should match
  QE to the per-component tolerance (≤ 50 meV/atom).
- If SCF drifts away in subsequent iterations, the mixer is the bug
  (H4); if SCF stays at QE's fixed point with residuals below the
  per-component tolerance, the assembly was the bug and Part A
  didn't catch it.

**Part C — Fix (1–5 CE-days).**

Scope reserved for the actual fix once Part A or B localizes the
bug. Expected fix classes, in decreasing prior:

- A sign or dimensional error in the V_loc(G=0) compensation for
  multi-species or semicore PPs.
- A mis-pairing of input vs. output density in `harris_foulkes_energy`
  or one of its feeders (the Phase 1a fingerprint is consistent with
  `ρ_in` vs. `ρ_out` confusion).
- A missing renormalization in the NLCC double-counting path when
  ρ_core doesn't integrate exactly to `z_core` due to log-mesh
  truncation.

**Part D — Close the matrix (0.5 CE-day).**

Once Part C lands, rerun all 5 heavy-atom `qe_validation.rs` tests.
Drop `#[ignore]` on those that close below 50 meV/atom. Update test
pins. Remaining cells filed as narrow follow-ups.

## Deliverables

- `scripts/validate/vgch_energy_assembly_heavy.py` — QE per-term
  parser (reads QE stdout `E_one_electron`, `E_hartree`, `E_xc`,
  `E_ewald`, plus the `one electron contribution` sub-terms from
  `verbosity='high'`).
- `tests/vgch_energy_assembly_heavy.rs` — Tier-2 integration test
  asserting per-term agreement for C diamond, Fe BCC, Cu FCC, GaAs,
  NaCl, MgO. Failure mode: print per-term Δ vs QE in eV.
- Fix on `src/scf/energy.rs` (most likely) or `src/scf/context.rs`.
- 3–5 `#[ignore]` removals in `tests/qe_validation.rs`.

## Acceptance

Close VGCH-2 when:

1. Fe BCC 8×8×8 ≤ 50 meV/atom residual (currently 9.5 eV).
2. C diamond residual ≤ 50 meV/atom (currently 1.45 eV).
3. At least 3 of 5 heavy-atom `qe_validation.rs` tests have
   `#[ignore]` dropped.
4. `tests/vgch_energy_assembly_heavy.rs` pins per-term agreement at
   ≤ 50 meV/atom for all 7 VGCH-scope systems.
5. `tests/vgch_sad_heavy.rs` and `tests/vgch_beta_l_heavy.rs` remain
   green (no regression on the Phase 1b/1c bit-perfect baselines).

## Cost

1–2 CE-weeks, depending on whether Part A's trace pins the bug or
Part B's transplant experiment is needed.

- Part A (per-term trace): 1 CE-day.
- Part B (transplant experiment): 1 CE-day.
- Part C (fix + tests): 1–5 CE-days.
- Part D (close matrix): 0.5 CE-day.

## Non-goals

- Fixing the SAD clamp (Phase 1c showed the clamp effect is O(1e-5
  eV) on the worst cell, GaAs). Can be a separate polish ticket.
- USPP / PAW support.
- Spin-orbit coupling.
- Whole-new mixer — the mixer is cleared for Si/Al/C (MPSH +
  SYKP closes those). Only revisit if Part B's transplant stays at
  the QE fixed point with no drift, ruling out assembly and pointing
  at mixer.

## Related

- **VGCH** (parent) — Phase 1a diagnostic landed PR #139, Phase 1b
  H1 (β_l(q)) landed PR #148, Phase 1c H2 (SAD) in this PR.
- **VQEF** — blocked on VGCH-2 for the 5 heavy-atom PBE-ready cells.
- **MPSH** — landed PR #110, closed Si/Al/Cu/Fe shift-convention
  gaps. MPSH is independent of VGCH-2; together they unblock the
  full VQEF matrix.
