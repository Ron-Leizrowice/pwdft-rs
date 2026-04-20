---
id: VGC5
status: completed
priority: critical
complexity: medium
risk: low
depends_on: [VGCMP]
blocks: []
owner: researcher
---

# VGC5: VGCMP Phase 5 — Per-Component Energy Accounting (Si vs QE)

## Origin

VGCMP Phases 1+2+3+4 (PR #29 merged 2026-04-17) established that the entire pseudopotential → Hamiltonian assembly pipeline is bit-correct against independent Python references:

- V_local(G): max |Δ| = 2.78e-9 Ry (Phase 1)
- β_l(q): bit-equality (Phase 2)
- D_ij: bit-equality (Phase 3)
- Diagonal H[G,G]: max |Δ| = 2.9e-10 Ry (Phase 4)

**The 13.4 eV Si gap (pwdft-rs −218.18 eV vs QE −231.61 eV) is OUTSIDE the matrix assembly.** Possible loci:

- Total-energy term composition (kinetic, local, non-local, Hartree, XC, Ewald summation in `total_energy()`)
- The V_local(G=0) compensating background shift (likely missing — see Prime Suspect below)
- Ewald summation for the 2-atom Si diamond primitive vs the 1-atom Fe BCC cell
- SAD initial density pathology for covalent Si

Fe BCC matches QE to 0.02 eV. Si misses by 13.4 eV. The Fe-vs-Si asymmetry is the strongest signal: **whatever's wrong is geometry-dependent**, scaling with atoms/cell or basis covalency.

## Plan

Open branch `VGC5/per-component-energy-accounting`. Steps:

1. **Run Si SCF in pwdft-rs at ecut=30 Ry** (matches `qe_validation/si_scf.in`) with verbose energy logging:

   ```text
   E_band     = ...
   E_kinetic  = ...
   E_local    = ...
   E_local(G=0) compensating shift = ...
   E_nonlocal = ...
   E_Hartree  = ...
   E_xc       = ...
   E_ewald    = ...
   E_total    = sum
   ```

2. **Extract the same components from QE's `si_scf.out`** (QE prints all of these in its standard output).

3. **Tabulate side-by-side.** The 13.4 eV must localize to one specific term (or split across two).

4. **For the prime suspect (V_local(G=0)):** add a unit test asserting `total_energy()` includes the `V_local(G=0) · N_el` background shift. If it doesn't, write a fix (Researcher proposes; Core Engineer implements).

5. **Cross-check against Fe BCC:** run the same per-component report for Fe and confirm the discrepancy is Si-specific (not a systematic offset).

## Prime suspect

`src/scf/context.rs:93-94` zeroes `v_local_fft[0]` and stashes `v_local_g0` separately. The standard pseudopotential treatment requires adding back a constant background `V_local(G=0) · N_el / Ω` to the total energy to cancel the singular G=0 Coulomb piece. If `total_energy()` (in `src/scf/energy.rs` or wherever the per-term sum lives) skips this, the missing term scales with electron count and unit cell volume — explaining why Si (8 valence electrons / cell) misses by ~13 eV while Fe (1 atom / cell, 8 or 16 valence) matches to 0.02 eV (the Fe geometry happens to make the missing piece small or cancel).

To check: read `src/scf/context.rs:93-94`, then `src/scf/mod.rs::total_energy` (or wherever), and grep for `v_local_g0` to see whether it's ever consumed. If not consumed in the energy expression, that's the bug.

## Other suspects (rank-ordered)

1. **V_local(G=0) compensating shift missing in `total_energy()`** — most likely.
2. **Ewald for diamond primitive** — 2-atom Si vs 1-atom Fe; pair-sum double-counting plausible. Cross-check the lattice sum convention against QE's `electrons.f90` Ewald call.
3. **SAD initial density** — for covalent Si the SAD overlap may put the SCF in a basin that QE's QM-style start avoids. Likely produces convergence drift, not a 13 eV offset, so lower priority.
4. **Kinetic operator unit convention** — already cleared by Phase 4 to 7e−11 rel; would have to be in the energy aggregation, not the matrix element.

## Verification

The 13.4 eV gap should localize to one term. Once located:

- Write a unit test that pins the per-component value against an analytical expression or QE reference.
- Implement the fix in a separate proposal (e.g., `VGFX-vloc-g0-shift` or `EWFX-ewald-diamond-fix`).
- Re-run the Si SCF; confirm |E_total - E_QE| < 0.1 eV.

## References

- `proposals/VGCMP-vloc-g-cross-check.md` (Phase 5 recommendation section)
- `tests/vgcmp_assembled_h_cross_check.rs` (Phase 4 test, baseline for Phase 5)
- `scripts/validate/vgcmp_phase4_assembled_h.py` (template for Phase 5 Python reference)
- `qe_validation/si_scf.in` and `si_scf.out` (QE reference run)
- `.claude/logbooks/researcher.md` 2026-04-17 Phase 4 entry — full context
