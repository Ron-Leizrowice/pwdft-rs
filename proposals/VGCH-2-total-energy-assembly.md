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

**Part A — Isolate the mispairing (1 CE-day).**

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
