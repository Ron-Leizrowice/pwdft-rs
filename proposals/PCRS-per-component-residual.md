---
id: PCRS
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# PCRS: Per-Component Energy Residual Investigation

## Origin

VGC5 Phase 5 (PR #34, 2026-04-17) instrumented per-component SCF energies and observed that for Si diamond, the identity

```
E_total = E_kinetic + E_local + E_nonlocal + E_Hartree + E_xc + E_ewald
```

closes to only **1.20 eV** at `conv_threshold = 1e-8`. Fe BCC shows 0.045 eV residual under the same conditions.

**Order-of-magnitude analysis:** with converged density RMS ~1e-8 e/Å³, V_xc scale ~10 eV, and ~8 valence electrons, the expected residual from incomplete self-consistency is ~1e-6 eV — six orders of magnitude below what we observe.

The discrepancy is structural, not a self-consistency artifact. VGC5 Code Reviewer's spot-analysis suggests the residual is geometry-specific (Si much worse than Fe) and possibly tied to `src/scf/density.rs:90-97` (rescale-to-N_electrons step) or the density symmetrization path. Not confirmed.

## Proposed investigation

1. **Tighten `conv_threshold` progressively** (1e-9, 1e-10, 1e-11). Does the residual shrink proportionally, or does it plateau? If it plateaus, the bug is not SCF noise — it's a bookkeeping inconsistency.

2. **Audit the decomposition identity.** For a fully-converged density, `Σ(components)` should equal `E_total` to machine precision. Read `src/scf/mod.rs::total_energy` and the VGC5 `EnergyComponents` computation side-by-side:
   - Are the components using the same density as `total_energy()`?
   - Does the rescale-to-N_electrons step (`src/scf/density.rs:90-97`) happen between component computation and total-energy computation?
   - Is density symmetrization applied once or twice in either path?

3. **Compare QE's behavior.** QE prints `total energy = -231.610 eV` plus per-component values. Do they close to the same precision we expect, or does QE exhibit a similar residual? (The VGC5 `scripts/validate/vgc5_per_component.py` already parses QE output — extending it to compute the QE-side identity check is ~10 lines.)

4. **Repeat for Fe BCC.** Is Fe's 0.045 eV residual the same structural issue in miniature, or a different class of discrepancy?

## Why it matters

Once NCFX (the NLCC core-density fix, new proposal from VGC5) lands and closes the 13.4 eV E_xc gap, the next source of truth for "does our total energy match QE?" will be the per-component identity. A 1.2 eV self-inconsistency between our own component computation and our own total-energy computation means we can't trust either number in isolation — we'd always need the sum. That's a regression-detection weakness; tightening it makes the test suite stronger.

## Scope

Investigation-only. Output is either:
- **No bug:** a writeup explaining why the residual is expected (e.g., known SCF normalization step), with the bound quantified analytically.
- **Bug:** a targeted fix proposal (e.g., `PCFX-per-component-density-consistency`).

Do NOT modify energy computations without a clear plan — this is measurement + analysis.

## Verification

Deliverable: a logbook entry or follow-up proposal documenting:
1. Residual vs `conv_threshold` curve for Si.
2. Whether QE shows the same.
3. Identified source of the 1.2 eV (or proof it's expected).

## References

- `proposals/VGC5-per-component-energy-accounting.md` — parent.
- `tests/vgc5_per_component_si.rs:128-135` — current observation comment.
- `src/scf/mod.rs` — `EnergyComponents` population.
- `src/scf/density.rs:90-97` — rescale-to-N_electrons step (possible culprit).
- `scripts/validate/vgc5_per_component.py` — QE output parser.
