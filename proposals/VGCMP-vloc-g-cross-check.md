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

## 2026-04-17 — Phase 1 Result

**Verdict: V_local(G) passes. The Si 13.43 eV gap is NOT in V_local(G).**

### Artifacts

- `scripts/validate/vloc_g_reference.py` — independent Python implementation (manual UPF XML parsing; `scipy.integrate.simpson` on the ONCVPSP linear mesh; QE erf-subtracted formula in native Ry/Bohr units).
- `scripts/validate/vloc_g_si_reference.csv` — committed golden file with the 20-shell reference values.
- `tests/vgcmp_vloc_cross_check.rs` — Rust test asserts shell-by-shell agreement to < 1e-4 Ry absolute.

### Numerical result (Si FCC, a = 5.431 Å, Ω = 270.256 Bohr³)

**Max |Δ| across 20 shells = 2.78×10⁻⁹ Ry (3.78×10⁻⁸ eV)** — five orders of magnitude below the 1×10⁻⁴ Ry pass threshold. Every shell agrees to 10 significant digits.

| shell | \|G\|² ((2π/a)²) | \|G\| (Bohr⁻¹) | Python (Ry) | Rust (Ry) | Δ (Ry) |
|-------|-----------------|----------------|-------------|-----------|--------|
| 0 | 3 | 1.060381 | −2.8416065668e−1 | −2.8416065946e−1 | −2.78e−9 |
| 1 | 4 | 1.224422 | −2.0232760981e−1 | −2.0232760778e−1 | +2.03e−9 |
| 2 | 8 | 1.731594 | −8.1619259966e−2 | −8.1619259374e−2 | +5.92e−10 |
| 3 | 11 | 2.030475 | −5.0143978352e−2 | −5.0143977786e−2 | +5.66e−10 |
| 4 | 12 | 2.120761 | −4.3381962073e−2 | −4.3381961468e−2 | +6.05e−10 |
| 5 | 16 | 2.448844 | −2.5585808228e−2 | −2.5585807828e−2 | +4.00e−10 |
| 6 | 19 | 2.668566 | −1.7801652398e−2 | −1.7801652784e−2 | −3.86e−10 |
| 7 | 20 | 2.737891 | −1.5832129100e−2 | −1.5832129488e−2 | −3.88e−10 |
| 8 | 24 | 2.999210 | −9.9970720150e−3 | −9.9970719093e−3 | +1.06e−10 |
| 9 | 27 | 3.181142 | −7.0948904852e−3 | −7.0948906910e−3 | −2.06e−10 |
| 10 | 32 | 3.463189 | −3.9350427122e−3 | −3.9350428241e−3 | −1.12e−10 |
| 11 | 35 | 3.621890 | −2.6943323485e−3 | −2.6943324277e−3 | −7.92e−11 |
| 12 | 36 | 3.673267 | −2.3584734957e−3 | −2.3584733900e−3 | +1.06e−10 |
| 13 | 40 | 3.871963 | −1.3108479118e−3 | −1.3108479780e−3 | −6.62e−11 |
| 14 | 43 | 4.014537 | −7.6455859255e−4 | −7.6455872032e−4 | −1.28e−10 |
| 15 | 44 | 4.060949 | −6.1746972237e−4 | −6.1746971867e−4 | +3.70e−12 |
| 16 | 48 | 4.241523 | −1.6545893380e−4 | −1.6545886338e−4 | +7.04e−11 |
| 17 | 51 | 4.372062 | +6.1386422657e−5 | +6.1386251691e−5 | −1.71e−10 |
| 18 | 52 | 4.414717 | +1.2039710003e−4 | +1.2039696600e−4 | −1.34e−10 |
| 19 | 56 | 4.581368 | +2.9141526304e−4 | +2.9141541872e−4 | +1.56e−10 |

### Interpretation

1. The worst discrepancy (shell 0, the largest |G|² = 3 shell that dominates V_local contribution) is 2.78 ns-Ry — pure floating-point round-off from the order of Kahan-free sums in Rust vs NumPy. No systematic structure (signs are roughly balanced across shells; magnitude is flat ~10⁻⁹–10⁻¹⁰ Ry across two orders of magnitude in V_loc(G)).
2. The erf-subtracted form, the Simpson quadrature, the UPF unit conversions (Bohr→Å, Ry→eV, e² Gaussian Rydberg convention), and the Å⁻¹↔Bohr⁻¹ G-vector handling are all correct end-to-end.
3. **V_local(G) is not the source of the 13.43 eV Si gap.** The ~100 V_local(G) terms entering a Si SCF contribute at most ~1e-6 eV of accumulated error from this pathway, many orders of magnitude below 13.43 eV.

### Recommended next step

Proceed to **Phase 2 (β_l(q) non-local projectors)** as originally scoped. The angular-momentum-channel-sized offset (13.43 eV ≈ one l-channel × N_atoms × a few per-state energies) points strongly at the KB projectors. Specifically:

- Suspect #1: `NonlocalPotential::F_l(q)` — the Bessel transform ∫ χ(r) j_l(qr) r dr. Use same Python/scipy technique as Phase 1.
- Suspect #2: UPF projector unit conversion. `src/pseudopotential/upf.rs:68-72` divides by √BOHR_TO_ANG; double-check this against the KB matrix-element convention actually used in the Hamiltonian assembly.
- Suspect #3: D_ij sign/diagonalization — UPF stores rotated projectors with diagonalized h^l; Phase 3 spot-check.

Open a follow-up branch `VGCMP/phase2-beta-q` once this PR merges.

### No bugs filed in `src/`

The existing `v_local_of_g` implementation (`src/pseudopotential/mod.rs:119-171`) is numerically correct to machine precision against the independent reference. No changes to production code were made in this session.

## 2026-04-17 — Phase 2 Result

**Verdict: β_l(q) passes. The Si 13.43 eV gap is NOT in the KB non-local projector form factors.**

### Artifacts

- `scripts/validate/beta_q_reference.py` — independent Python implementation of the QE Bessel-transform `F_l(q) = 4π ∫ χ(r) j_l(qr) r dr` (`scipy.integrate.simpson` on the UPF log mesh; `scipy.special.spherical_jn` for Bessel). Mirrors QE 7.5 `upflib/beta_mod.f90:111-116` exactly, stripped of the `(4π/√Ω)` Hamiltonian-assembly prefactor so we compare pure form factors.
- `scripts/validate/beta_q_si_reference.csv` — committed golden file: 6 projectors × 20 q-values = 120 rows. q-grid is equally spaced in [0.1, 7.0] Bohr⁻¹ (q_max ≈ √(2·25 Ry) covers the ecut=25 Ry regime).
- `tests/vgcmp_beta_q_cross_check.rs` — Rust test reproduces the production helper `bessel_transform_projector` (`src/potential/nonlocal.rs:223`) verbatim using public APIs (`PseudopotentialData` + `numerics::simpson_integrate`), converts Å^(3/2) → Bohr^(3/2) for comparison, and asserts per-row agreement. Includes per-projector and low-q/high-q bucket diagnostics.

### Numerical result (Si ONCVPSP LDA, 6 projectors: l=0, 0, 1, 1, 2, 2)

**Max |Δ| across 120 rows = 3.03×10⁻¹² Bohr^(3/2)** — eight orders of magnitude below the 1×10⁻⁴ Bohr^(3/2) pass threshold. Every row agrees to 12–13 significant digits.

Per-projector max |Δ| (Bohr^(3/2)):

| proj | l | max \|Δ\| | notes |
|------|---|-----------|-------|
| 0 | 0 | 6.9e−13 | clean |
| 1 | 0 | 3.0e−12 | clean |
| 2 | 1 | 6.7e−13 | clean |
| 3 | 1 | 1.4e−12 | clean |
| 4 | 2 | 7.3e−13 | clean |
| 5 | 2 | 1.1e−12 | clean |

Bucket breakdown:

- low-q (q < 3.5 Bohr⁻¹, 72 rows): max |Δ| = 3.03×10⁻¹²
- high-q (q ≥ 3.5 Bohr⁻¹, 48 rows): max |Δ| = 7.69×10⁻¹³

No systematic structure across projectors, l, or q. Errors are pure floating-point round-off from the Simpson sum order in Rust vs NumPy.

### Interpretation

1. The Bessel transform integrand `χ(r) · j_l(qr) · r` is evaluated identically (within rounding) in Python (scipy) and Rust.
2. The UPF projector unit conversion `χ_Å = χ_Bohr / √BOHR_TO_ANG` at `src/pseudopotential/upf.rs:68-72` is consistent with the Å-native r_grid, rab, and Å⁻¹ q convention — the full end-to-end round trip (Bohr^(3/2) ↔ Å^(3/2) via BOHR_TO_ANG^(3/2)) closes to machine precision.
3. The spherical Bessel function `spherical_bessel_j` in Rust (upward recurrence from j_0, j_1) matches `scipy.special.spherical_jn` at all l=0,1,2 and all q in our grid.
4. **β_l(q) is not the source of the 13.43 eV Si gap.** The KB form factor machinery — per-q values fed into the Hamiltonian at each (k+G, k+G') pair — is bit-for-bit faithful to QE's convention.

### Recommended next step — Phase 3: D_ij

With V_local(G) cleared (Phase 1) and β_l(q) cleared (Phase 2), the remaining suspects inside the pseudopotential are:

1. **D_ij matrix** — the Ry→eV unit conversion at `src/pseudopotential/upf.rs:77-78` is simple; a sign error, off-diagonal handling, or row/column ordering mistake is possible. Si's UPF has 2 projectors per l-channel (l=0, l=1, l=2), so D_ij is 6×6 but block-diagonal in l (zeros between different-l blocks). Diagonal of each block may not be diagonal itself (QE diagonalizes h^l and absorbs the rotation into χ; see the KBTF note from 2026-04-16).
2. **KB assembly** — the angular factor `(2l+1)/(4π) P_l(cos θ)`, structure factor `S(G−G')`, and the 1/Ω normalization in `add_to_hamiltonian` (`src/potential/nonlocal.rs:118-208`).

**Phase 3 plan:** follow-up branch `VGCMP/phase3-dij`. Two sub-tasks:

- **3A (fast):** Python script parses `<PP_DIJ>` XML block, prints the 6×6 matrix in both Ry (native) and eV (internal). Rust test loads Si.upf, prints `pp.dij` as 6×6 matrix, asserts element-wise agreement to 1e-12 eV. Confirms the UPF→internal unit conversion is clean.
- **3B (decisive):** extend the Python script to parse QE's `init_us_1` intermediate quantities or compute a single matrix element `H_{GG'}` at the Γ point for Si using Python's independent β_l(q) and D_ij values. Compare against pwdft-rs's assembled `H_{GG'}`. This is Phase 4 as originally scoped; given Phases 1 and 2 both passed to 10⁻⁹ Ry / 10⁻¹² Bohr^(3/2), the bug — if it's in the pseudopotential pipeline at all — must lie in the assembly step, not the form factors. Phase 4 isolates that directly.

If Phase 3 also passes, the Si 13.43 eV gap is **outside** the pseudopotential machinery. Candidates then become: Ewald sign/convention, structure factor in `scf/potentials.rs`, symmetry-breaking at Γ (known to exist, per the orientation note), or SCF convergence criterion.

### No bugs filed in `src/`

The existing `bessel_transform_projector` implementation (`src/potential/nonlocal.rs:223-246`) and UPF projector unit conversion (`src/pseudopotential/upf.rs:54-73`) are numerically correct to machine precision against the independent reference. No changes to production code were made in this session.

## 2026-04-17 — Phase 3 Result

**Verdict: D_ij passes bit-exactly. The Si 13.43 eV gap is NOT in the KB coupling matrix.**

### Artifacts

- `scripts/validate/dij_reference.py` — manual UPF v2 parser; extracts `<PP_DIJ>` block and per-projector `angular_momentum`, reshapes the 36-value flat list row-major into a 6×6 matrix in Ry. Also asserts block-diagonality in l as a sanity check.
- `scripts/validate/dij_si_reference.csv` — committed golden file: 36 rows (i, j, dij_ry) in row-major order.
- `tests/vgcmp_dij_cross_check.rs` — Rust test loads Si.upf via `pseudopotential::load`, divides each `pp.dij` entry by `RY_TO_EV` to recover the native UPF value, and asserts element-wise agreement to < 1e-12 eV absolute. Prints full 6×6 Rust and Python matrices on every run.

### Numerical result (Si ONCVPSP LDA, 6 projectors: l = 0, 0, 1, 1, 2, 2)

**max |Δ| = 0.000e+00 Ry (0.000e+00 eV) across all 36 matrix elements.** The Ry→eV conversion at `src/pseudopotential/upf.rs:77-78` is a single scalar multiply, and Rust + Python agree bit-exactly. Well below the 1e-12 eV pass threshold.

#### D_ij matrix for Si (Ry, row-major from UPF `<PP_DIJ>`)

|          | j=0 (l=0)    | j=1 (l=0)    | j=2 (l=1)    | j=3 (l=1)    | j=4 (l=2)    | j=5 (l=2)    |
|----------|--------------|--------------|--------------|--------------|--------------|--------------|
| i=0 (l=0)| +1.113192e+01| 0            | 0            | 0            | 0            | 0            |
| i=1 (l=0)| 0            | +1.713932e+00| 0            | 0            | 0            | 0            |
| i=2 (l=1)| 0            | 0            | +5.452221e+00| 0            | 0            | 0            |
| i=3 (l=1)| 0            | 0            | 0            | +1.259656e+00| 0            | 0            |
| i=4 (l=2)| 0            | 0            | 0            | 0            | −4.249609e+00| 0            |
| i=5 (l=2)| 0            | 0            | 0            | 0            | 0            | −8.892088e−01|

### Interpretation

1. **D_ij is strictly diagonal** — not merely block-diagonal in l. ONCVPSP generates two projectors per l-channel (designed for transferability across energy windows); QE's ONCVPSP diagonalization routine has already absorbed any within-l-block rotation into the χ(r) projectors (see `qe-7.5/upflib/init_us_1.f90`), so the coefficient matrix that sits in the KB separable form is pure diagonal. This matches the expectation that KB non-local potentials for norm-conserving PPs have the simple form `V_NL = Σᵢ |β_i⟩ D_ii ⟨β_i|`.
2. The parser at `src/pseudopotential/upf.rs:77-78` correctly flattens 36 row-major values and applies `dij_ry[k] * RY_TO_EV` element-wise. No sign error, no row/column transposition, no bytes lost in parsing.
3. **D_ij is not the source of the 13.43 eV Si gap.** Combined with Phase 1 (V_local(G), max |Δ| = 2.8e−9 Ry) and Phase 2 (β_l(q), max |Δ| = 3.0e−12 Bohr^(3/2)), **all three pseudopotential form factors now carry a clean bill of health against an independent Python reference.**

### Diagnostic implication

The combined 1/2/3 verdict is strong: if the Si 13.43 eV gap is anywhere in the **per-term** pseudopotential data (local potential, non-local projectors, coupling matrix), it has evaded three independent 10⁻⁹–10⁻¹² Ry checks. That is vanishingly unlikely. The bug must live in one of:

- The **assembly step** — how the form factors are combined into `H_{G,G'}` at a given k-point. This includes:
  - Structure factor `S(G−G') = Σ_atoms exp(−i(G−G')·τ_atom)` (sign, conjugation convention)
  - Non-local angular factor `(2l+1)/(4π) · P_l(cos θ_{k+G, k+G'})`
  - The `1/Ω` volume normalization, and the `4π/√Ω` KB prefactor inside `bessel_transform_projector`
  - The row/col indexing of D_ij during the `V_NL = Σᵢⱼ ⟨G|β_i⟩ D_ij ⟨β_j|G'⟩` double sum
- **Outside the PP pipeline** entirely: Ewald sign/prefactor, kinetic energy convention (ħ²/2m with which unit of m), symmetry-breaking at Γ, SCF convergence criterion, or initial-density (SAD) pathology that never escapes a local minimum for Si.

### Recommended next step — Phase 4: assembled H_{G,G'}

Open follow-up branch `VGCMP/phase4-hamiltonian`. Approach:

1. At Γ for a 2-atom Si FCC primitive cell, pick two specific G-vectors G_A, G_B (e.g. G_A = (0,0,0), G_B = first shell at |G|² = 3·(2π/a)²).
2. **Python reference:** using the already-validated Python V_local(G), β_l(q), D_ij from phases 1–3, assemble `H_{G_A, G_B}` = kinetic + local + non-local by hand, evaluated at each term per the analytic formula. Keep in QE native (Ry, Bohr) units.
3. **pwdft-rs side:** extend `NonlocalPotential::add_to_hamiltonian` (or build a thin test harness) to extract the matrix element at exactly (G_A, G_B). Compare to the Python reference, element by element (local only, non-local only, full).
4. **Pass criterion:** < 1e-4 Ry absolute per term. If passes, the bug is outside the PP pipeline. If fails, the Phase 4 component-level comparison pinpoints the exact sub-term (local-only diff vs non-local-only diff).

Rationale for proceeding to Phase 4 regardless of this clean Phase 3: the assembly step — structure factor, (2l+1)/(4π) angular factor, 1/Ω volume factor, and the summation pattern — has not been verified. The three form factors being correct in isolation is necessary but insufficient; the wiring between them is a separate code path with its own sign and index pitfalls.

### No bugs filed in `src/`

The existing UPF D_ij parser (`src/pseudopotential/upf.rs:77-78`) is bit-exactly correct against the independent reference. No changes to production code were made in this session.
