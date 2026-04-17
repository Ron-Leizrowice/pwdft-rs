---
id: VGCMP
status: active
priority: critical
complexity: medium
risk: low
depends_on: [VERF]
blocks: [QEDX, QEVL]
owner: researcher
---

# VGCMP: V_local(G) and KB Projector Cross-Check vs QE 7.5

> **Context:** After SIMP (Simpson's rule) and VERF (erf subtraction), Si diamond total energy is still 13.43 eV above QE's reference (-218.18 vs -231.61 eV). Fe BCC matches to 0.02 eV. VERF alone is numerically equivalent to the bare-Coulomb form on our log mesh, so the Si discrepancy must come from either the reciprocal-space V_local(G) values, the Kleinman–Bylander non-local projectors, or their G-by-G assembly. This proposal nails down where.

## Problem

We do not currently know, to numerical precision, whether pwdft-rs and QE 7.5 agree on the individual form factors that enter the Hamiltonian:

1. **V_local(G)** — the Fourier transform of the local pseudopotential at each |G| shell.
2. **β_l(q)** — the non-local KB projector form factors at each |k+G|.
3. **D_ij** — the coupling matrix used inside the KB matrix element.

Any one of these can be off by a sign, a factor of 4π, a Rydberg↔Hartree, a Bohr↔Å, or a log-mesh interpolation step, and produce a multi-eV error in the total energy without breaking internal consistency tests. The 13.43 eV Si offset is the size of an *angular-momentum channel*, which points at the non-local machinery specifically.

## Goal

Produce a side-by-side numerical table for Si (Z=14, nc LDA UPF) of:

- `V_local(G)` for the first 20 distinct |G| shells (0, √3, √8, √11, √16, √19, …, in units of 2π/a)
- `β_l(q)` for each l-channel at 10 representative q values covering 0 → q_max
- `D_ij` — the full matrix as extracted from the UPF versus as actually used in the Hamiltonian

with tolerance **< 1e-4 eV absolute** (matches QE's UPF interpolation precision).

If all three quantities agree, the 13.43 eV Si error is elsewhere (Ewald, kinetic, symmetry) — that narrows the hunt by orders of magnitude. If any disagrees, we have the smoking gun.

## Reference QE machinery

- **V_local(G):** `qe-7.5/upflib/vloc_mod.f90:136-148` (erf-subtracted integrand, Simpson's rule on log mesh) → `tab_vloc(iq,nt)` at `nqx` equally-spaced q values → `interp_vloc` with 4-point Lagrange interpolation. Then `vloc_of_g` re-adds `−4π·Z·e²·exp(−G²/4)/(Ω·G²)` for the analytic Coulomb tail.
- **β_l(q):** `qe-7.5/upflib/init_tab_beta.f90` — same strategy, Simpson's rule in radial space, tabulated on equally-spaced q mesh.
- **D_ij:** UPF parser normalizes `dion` to Ry·e units; PP_DIJ block already sits in this space in UPF v2.

### How to extract QE reference values

QE does not expose `tab_vloc` directly on disk. Three options, in decreasing order of cost:

1. **Patch QE to dump `tab_vloc` and `tab_beta`** after `init_tab_vloc`. One-line `WRITE` statement, recompile, re-run Si SCF. Most rigorous.
2. **Call QE's routines from a small Fortran driver** that reads Si.upf, calls `init_tab_vloc` / `init_tab_beta`, and writes an ASCII table. Low-risk, one afternoon of Fortran.
3. **Re-derive V_local(G) from `upf%vloc(r)`** using an independent Python implementation with Simpson's rule. Redundant with our Rust code but provides a second witness.

Recommended: start with **option 3** (Python) in one morning, because it's the fastest signal. If Python matches QE's final energy (reconstructed from Si.upf + manual Ewald + manual KB assembly), then we know the Si.upf data is fine and the error is in our Rust assembly pipeline. If Python also disagrees with QE's final energy, it indicates a UPF parsing issue and option 1/2 is needed to isolate it.

## Methodology

### Phase 1: V_local(G) — Si FCC, 20 shells (1 day)

1. Write a Python reference implementation (`scripts/validate/vloc_g_reference.py`):
   - Parse Si.upf manually (extract `PP_LOCAL`, `PP_R`, `PP_RAB`, `Z_valence`)
   - Apply QE's erf-subtracted integrand (r in Bohr, erf(r) with Bohr argument, `exp(−G²·tpiba2/4)`)
   - Use `scipy.integrate.simpson(f, x=r)` on the log mesh
   - Output `V_local(G)` in Ry for the first 20 |G| shells, as |G|² = 3, 8, 11, 16, 19, 20, 27, 32, 33, 35, 36, 40, 41, 43, 44, 48, 49, 51, 52, 56 (in units of (2π/a)²)
2. Compare against pwdft-rs `PseudopotentialData::v_local_of_g(g_norm, omega)` at the same |G|, after converting:
   - |G| from (2π/a units, Å⁻¹) — our native units
   - Result from eV → Ry for direct comparison
3. **Pass criterion:** all 20 shells agree to **< 1e-4 Ry** (equivalent to 1e-3 eV per shell, which would be < 0.1 eV total for the ~100 V_local(G) terms in a tight Si SCF).
4. **If disagreement appears at specific shells (e.g. small |G|):** suspect the G=0 branch or the analytic Coulomb correction sign/factor. If it grows with |G|: suspect Simpson vs trapezoidal on the log mesh, or a missing upper-bound cutoff.

### Phase 2: KB β_l(q) form factors (1 day)

1. Extend `scripts/validate/` with `beta_q_reference.py`: parse `PP_NONLOCAL/PP_BETA.i` for each projector, compute `F_l(q) = 4π ∫₀^∞ χ(r) j_l(qr) r dr` using Simpson's rule, in QE's convention (χ stores r·β(r)).
2. Evaluate on a q-grid covering [0, q_max] where q_max = √(2·ecut) in Ry — for ecut=25 Ry this is q_max ≈ 7 Bohr⁻¹ ≈ 13 Å⁻¹.
3. Compare with `NonlocalPotential::F_l(q)` (or equivalent internal) at matched q values.
4. **Pass criterion:** < 1e-4 Å^(3/2) absolute per projector per q.
5. **If disagreement:** the trapezoidal-Simpson KBTF test 07 was relaxed (ratio 0.12 instead of 0.10) for the HGH l=1 projector. That was a symptom — confirm with this cross-check whether HGH is actually OK, or whether the issue is general.

### Phase 3: D_ij sanity check (half day)

For Si's nc LDA UPF, D_ij is diagonal (l=0 and l=1 channels don't mix). Just print the `dij` array from `PseudopotentialData` and compare with the `<PP_DIJ>` XML block, accounting for the RY_TO_EV conversion. This is a spot-check — if Phases 1 and 2 agree, D_ij is almost certainly fine.

### Phase 4: Assembled Hamiltonian element at a single k-point (half day)

At the Γ point, pick a pair of G-vectors and manually compute:
  - Kinetic: `|k+G|² · ħ²/(2m)`  (independent)
  - Local: `V_local(G−G') · S(G−G')`  (uses Phase 1)
  - Non-local: KB sum over projectors  (uses Phase 2 + Phase 3)

Evaluate `H_{GG'}` and compare with our internal matrix element. A single assembled Hamiltonian entry matching to 1e-4 eV is strong evidence the per-term agreement propagates correctly.

## Deliverables

- `scripts/validate/vloc_g_reference.py` — Python reference (runs in `uv` env)
- `scripts/validate/beta_q_reference.py` — Python reference for β_l(q)
- `tests/qe_numerical_cross_check.rs` — integration test calling the Python scripts and comparing against pwdft-rs (or alternatively, golden-file CSV tables generated once and checked in)
- A short writeup in `.claude/logbooks/researcher.md` with the table of discrepancies (if any) or a clean bill of health

## Success criteria

1. **If all three quantities match QE to < 1e-4 Ry:** declare the form-factor machinery correct. The 13.43 eV Si error then lies in Ewald, symmetry, mixing, or some other non-form-factor site — open a follow-up proposal to investigate those.
2. **If V_local(G) disagrees:** file a targeted fix proposal. Likely culprits (in order): wrong convention for erf Gaussian width (1 Bohr vs 1 Å), missing QE-style tabulation + Lagrange interpolation (vs our direct evaluation), sign error on the analytic Coulomb correction at small |G|.
3. **If β_l(q) disagrees:** file a fix proposal for the Bessel-transform code in `NonlocalPotential` or the UPF projector unit conversion (Bohr^(-1/2) → Å^(-1/2)).
4. **If D_ij is wrong:** fix the UPF parser unit conversion.

## Rationale for priority

This is the **bottleneck** proposal on the validation track: QEDX, QEVL, and any future Si-related physics work are all gated on understanding the 13.43 eV Si gap. The current state is uncomfortable — Fe passes, Si fails by an amount much larger than any known source of numerical error. Without VGCMP, we cannot trust any Si (or Si-like, including C) result from pwdft-rs.

## Estimated effort

3–4 days for a researcher who's comfortable with Python, UPF XML, and QE source. Phase 1 alone (V_local(G) cross-check) can be done in a day and will almost certainly isolate the issue if it lives there.

## Files touched

- New: `scripts/validate/vloc_g_reference.py`
- New: `scripts/validate/beta_q_reference.py`
- New: `tests/qe_numerical_cross_check.rs` (or `scripts/validate/compare.py` + golden CSV)
- Possibly: bugfix in `src/pseudopotential/mod.rs` or `src/potential/nonlocal.rs` depending on Phase 1/2 results

## Related

- VERF (completed) — established that the erf vs bare-Coulomb distinction is cosmetic on our mesh
- SIMP (completed) — established Simpson's rule; closed Fe gap but not Si gap
- QEDX (tracking) — will be archived once VGCMP lands or hands off the remaining Si error to a successor proposal
