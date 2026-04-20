---
id: VGCH-MECH
title: VGCH mechanism taxonomy — separating density-basin / Hamiltonian-side / light-atom classes
status: active
priority: high
complexity: medium
risk: medium
depends_on: [VGCH-2, BSUM]
blocks: [VQEF]
author: Engineering Manager (post 2026-04-19 investigation wave)
---

# VGCH-MECH — Mechanism taxonomy for the heavy-atom residual

> **Scope**. After today's 56-PR investigation wave (VGCH Phase 1a/b/c, VGCH-2 Part A, TSEN, SiEF-B1, BSUM), the "VGCH heavy-atom residual" is no longer a single 7-system bucket. It fractures into three distinct mechanism classes with distinct physics and distinct fix-site candidates. This proposal documents the taxonomy, assigns each class to a targeted investigation, and sets acceptance criteria per class. VGCH-2 Part B is the active investigation for Class A and is already in flight; Classes B and C are new.

## Background: what today's wave established

### What was ruled out

1. **β_l(q) projector form factors** (VGCH Phase 1a, #139; Phase 1b, #148) — bit-perfect to 3 × 10⁻¹² Bohr^{3/2} across 11 elements and 10 q-values. The Kleinman-Bylander separable projectors carry no detectable error vs QE's `init_us_1.f90`.
2. **Superposition-of-atomic-densities (SAD) initial density** (VGCH Phase 1c, #156) — max |Δρ(r)| ≤ 5.5 × 10⁻¹⁰ e/Å³ across seven systems; bit-perfect on C diamond. The only outlier (GaAs at 1.2 × 10⁻⁵) sits in a negative-density clamp band and accounts for ≤ 10⁻⁵ eV.
3. **V_local(G=0) absolute-reference gauge** (Si E_F investigation #164, SiEF-B1 #166) — the 1.35 eV rigid offset on Si E_F and the uniform Γ-eigenvalue shifts on C (−3.16 eV), Fe (−5.17 eV), Cu (−7.75 eV), NaCl (−1.67 eV), MgO (−3.35 eV) were a pure gauge-convention mismatch. Fixed by keeping V_loc(G=0) on the Hamiltonian diagonal (matching QE) and retiring the `with_g0_shift` compensation path. Si E_F now 6.4 meV; C E_F 68.2 meV; Al E_F 5.3 meV — all GREEN.
4. **Smearing entropy reporting** (TSEN #162) — the missing `−TS` term in pwdft-rs's reported `total_energy` accounted for 70 meV (GaAs) to 260 meV (Fe) of the apparent residual on every metal. Fixed; closed Al LDA (75 → 26 meV GREEN) and Al PBE (108 → 8 meV GREEN). The `−TS` values themselves match QE to 1-5 meV.

### What was confirmed

1. **Per-term energy assembly is bit-correct** (VGCH-2 Part A, #160). Direct-sum identity holds to sub-μeV; `N_el · Σ V_loc(G=0)` matches QE exactly. The residual is not in the arithmetic.
2. **Different converged densities** on 6 of 7 heavy-atom systems. The per-term fingerprint — `|Δone-e| ∈ [1.7, 53] eV`, `|ΔE_H| ∈ [0.6, 25] eV` with opposite sign, ratio `|Δone-e|/|−ΔE_H| ∈ [1.2, 2.9]` — is a classic linear-response δρ signature. QE and pwdft-rs each converge to a variational minimum within their own code's numerical basin, but the two basins are not the same density.
3. **Band-sum vs E_total ratio structure** (BSUM #165). `|ΔE_1e| / |ΔE_total|` splits the remaining cells cleanly:
   - **Light atoms (Si, Al, C):** ratio 0.1–5.3×. E_total and band-sum track within ~20%.
   - **Heavy atoms except Fe LDA:** ratio **1.5–3.3×** uniformly. Classic partial-cancellation between (T + V_ion) and (E_H + E_xc) of opposite sign — the density-basin signature.
   - **Fe LDA is the exception:** ratio 0.99× (vs Fe PBE at 3.25×). Band-sum and E_total move together → not the density-basin class.

## The three classes

### Class A — Heavy-atom energy-functional-at-same-density gap

**Members:** Cu LDA+PBE, GaAs LDA+PBE, NaCl LDA+PBE, MgO LDA+PBE, Fe PBE. (8 of 12 remaining YELLOW E_total cells.)

**Signature:**
- `|ΔE_1e| / |ΔE_total|` ∈ [1.5, 3.3] (partial cancellation between (T + V_ion) and (E_H + E_xc)).
- `|Δone-e| / |−ΔE_H|` ∈ [1.2, 2.9] across systems.
- PBE reduces the residual 1.6–6.5× over LDA (MgO: 10.7 → 1.6 eV, the largest improvement). Semicore systems show the largest improvement.

**Refutation of the "density-basin" hypothesis (VGCH-2 Part B, #167).** The original hypothesis — that the mixer converges to a different self-consistent density than QE's — was tested by seeding pwdft-rs with QE's converged ρ_QE and evaluating one SCF iteration. Result on Cu FCC:

- `E_HF` (Harris-Foulkes, gauge-invariant, evaluated at ρ_in = ρ_QE) = **−4837.30 eV** vs QE's −4853.64 eV → **+16.34 eV gap at the SAME density**.
- Per-term at shared density: Δone-e = +62.87 eV, ΔE_H = −54.74 eV, ΔE_xc = +8.87 eV. The partial-cancellation pattern persists even with density held fixed.
- Γ eigenvalue offset = +0.26 ± 0.03 eV/band (post-SiEF-B1; was −7.47 eV pre-gauge-fix).
- **E_F mismatch at shared density = +2.08 eV.** With eigenvalues agreeing to 0.26 eV per band, 1.82 eV of the Fermi-level disagreement is pure DOS/occupation integration, not eigenvalue shift.

Because `E_HF` uses ρ_in for all DC terms and is invariant to mixer trajectory, the 16.34 eV gap **cannot be a mixer-basin effect**. The codes disagree on the energy-functional value at a shared density.

**Revised physics hypothesis — Fermi-finder / smearing on dense DOS.** Leading candidate, from VGCH-2B's diagnosis:

1. **`smearing::find_fermi_energy` convergence on dense d-manifold DOS.** Cu 3d bands are tightly clustered; QE's `ef.f90` bisection may have different convergence tolerances or subdivision granularity than pwdft-rs's. A 1.8 eV E_F disagreement at shared density means the two codes integrate occupations over slightly different ε-ranges, producing different f_ik per band, producing different densities after one SCF step, and producing different per-term energies.
2. **`n_bands` margin above E_F.** Cu at `n_bands = 14` has only 4.5 unoccupied bands above E_F. Fermi tails on dense DOS need more margin; QE likely computes more bands by default and the pwdft-rs cell settings may be under-specifying.
3. **Smearing function vs QE's.** pwdft-rs defaults to Fermi-Dirac; QE's reference runs use Methfessel-Paxton or cold smearing depending on the input deck. The two functions integrate differently at the band-edge discontinuity. If the QE reference runs use MP and pwdft-rs uses FD, that's a functional mismatch masquerading as "different density."
4. **NLCC ρ_core(G) cross-check** — Fe ρ_core(G) is pinned by `test_fe_bcc_xc_nlcc_regression_guard`, but Cu / GaAs / MgO have no equivalent. If Cu's NLCC unit-conversion regressed silently, it would shift E_xc at shared density by exactly the observed scale.
5. **D_ij · Σ_lm β·β non-local projector contraction** — individual β_l(q) projectors are bit-perfect (VGCH-1b) but the double-sum assembly in `NonlocalPotential::add_to_hamiltonian` could differ from QE's `vnlocal.f90::vloc_psi` in a way only a full single-channel cross-check would catch.

**Active investigation (Part C — new scope):**
1. Reproduce pwdft-rs's E_F on QE's converged Cu eigenvalues in Python (SciPy or numpy), matching QE's smearing function exactly. If the Python reference gives QE's E_F, pwdft-rs's `find_fermi_energy` has a bug; if it gives pwdft-rs's E_F, the smearing function mismatch is the issue.
2. Run Cu at `n_bands = 24` or 32 (margin = 14-22 unoccupied bands) and see if the shared-density E_HF gap shrinks.
3. Triangulate the same transplant on C (no semicore, insulator) and Fe (spin-polarized, same 3d semicore structure). If C's transplant is bit-correct and Fe's shows the same 16 eV gap, the problem is specifically in semicore/dense-DOS systems.
4. Pin Cu / GaAs / MgO ρ_core(G) with the same regression guard Fe has.

### Part C findings (2026-04-20) — H-C1 and H-C3 CLEARED

Artifacts:

- `scripts/validate/vgch2c_fermi_reference.py` — QE-convention Fermi
  finder in Python. Parses QE's final `bands (ev)` + `wk = ...` blocks
  from a converged pw.x stdout, applies the three QE smearing functions
  (Fermi-Dirac, Gaussian, Methfessel-Paxton order-1) and runs two
  bisection variants — pwdft-style (bracket-width ≤ 1e-14 eV) and
  QE-style (|N(E_F) − N_target| < 1e-10 electrons, matching
  `efermig.f90:47,289-307`).
- `scripts/validate/vgch2c_fermi_reference.csv` — 5 systems × 3
  smearings × 2 bisections = 30 rows.

**H-C1 verdict: CLEARED.** With QE's converged eigenvalues + weights
fed into the Python Fermi-Dirac finder, E_F matches QE's reported
value to **6 μeV on Cu**, **15 μeV on Fe FM (nspin=2)**, **10 μeV on
C diamond**, **66 μeV on MgO**, **0.036–0.18 meV on NaCl**. The
bisection algorithm itself is not the bug. Both
`bracket-width < 1e-14 eV` and `|ΔN| < 1e-10 electrons` convergence
criteria give identical results to the μeV. pwdft-rs's
`smearing::find_fermi_energy` is structurally identical to QE's
`efermig.f90` (sign-convention check: pwdft-rs has
`x = (ε − E_F)/σ; f = 1/(1+exp(x))`; QE has
`wgauss((E_F − ε)/σ, −99) = 1/(1+exp(−(E_F−ε)/σ)) = 1/(1+exp((ε−E_F)/σ))`.
Identical.).

**H-C3 verdict: CLEARED on Cu.** Every QE reference input deck under
`qe_validation/*.in` uses the smearing that matches pwdft-rs's runtime
smearing (both F-D for Cu; see `qe_validation/cu_fcc_scf.in:18`).
On Cu the cross-smearing-function gap is only **−0.08 eV (Gaussian)**
/ **−0.10 eV (MP)** below the F-D answer — far smaller than the
observed 1.82 eV DOS-origin gap in the Cu transplant. Even a
hypothetical smearing mismatch cannot explain the Cu residual.

On wide-gap insulators (C, NaCl, MgO) the cross-smearing gap is
1.1–2.2 eV — large by construction because the band edge has no DOS
in the gap, so E_F is weakly constrained. This is expected and NOT a
bug; for those systems pwdft-rs and QE use the same smearing (F-D)
and match to <1 meV.

**Key convention landmine documented (for future Rust↔Python
validation of occupations):** QE's `sumkg.f90` integrates
`Σ_ik wk(ik) · Σ_b wgauss(...)` where `wk` is pre-multiplied by
`degspin` in `setup.f90:673` — so `Σ wk = 2` for nspin=1, and
`sumkg` does NOT apply an additional spin factor. pwdft-rs's
`find_fermi_energy` takes `kpoint_weights` summing to 1 plus an
explicit `spin_factor = 2.0/nspin`, and integrates
`Σ_ik w_k · spin_factor · f_ik`. The products are identical, but a
naive Python reference that uses QE's `wk` (sum = 2) PLUS an
explicit `spin_factor = 2` double-counts by a factor of 2. Any
future cross-check must verify Σ wk at the interface — done in
`vgch2c_fermi_reference.py::total_n_electrons`.

**Remaining suspects (for Part D or a follow-up):**

- **H-C2 (n_bands margin)** — not closed in this session. Would
  require running 4 Cu SCFs at `n_bands ∈ {14, 18, 24, 32}`, a
  ~100-minute machine-locked experiment. On the Cu transplant
  diagnostic, Γ has 14 bands reaching 45.2 eV (≈ 26 eV above E_F)
  and other k-points reach 30+ eV — the margin should be enough for
  F-D tails at σ = 0.272 eV. Priority LOWER than H-C4/H-C5.
- **H-C4 (Cu/GaAs/MgO ρ_core(G) unpinned)** — not closed in this
  session. Fe's `test_fe_bcc_xc_nlcc_regression_guard` catches
  unit-conversion errors on a Z=26 NLCC profile; Cu (Z=29), Ga
  (Z=31), As (Z=33), Mg (Z=12) have no equivalent. A silent
  unit-conversion bug would affect E_xc at shared density by
  exactly the observed order (Cu transplant ΔE_xc = +8.87 eV).
- **H-C5 (V_NL `D_ij · Σ_lm β·β` contraction on Cu d-projectors)** —
  not closed in this session. Individual β_l(q) are bit-perfect
  per VGCH-1b, but the double-sum assembly is untested on Cu-sized
  cells.

**Revised leading suspect ordering post-Part C:**
1. H-C4 (ρ_core Fourier transform on heavy semicore elements) —
   highest prior; fingerprint matches ΔE_xc on transplant.
2. H-C5 (V_NL assembly) — middle prior; would explain Δone-e at
   transplant.
3. H-C2 (n_bands margin) — lower prior; Cu already has 26 eV of
   headroom above E_F at Γ.

**Acceptance criterion (revised):** Part C closes at least 4 of the 8 Class A cells to ≤ 500 meV at production ecut via either (a) Fermi-finder fix, (b) `n_bands` margin increase as a VQEF-config change, or (c) smearing function alignment. Remaining cells become Class A sub-classes for Part D.

### Class B — Fe LDA Hamiltonian-side (new, identified by BSUM)

**Members:** Fe LDA (1 cell). Fe PBE behaves like Class A (ratio 3.25×); Fe LDA is the outlier (ratio 0.99×).

**Signature:**
- `|ΔE_1e| / |ΔE_total|` ≈ 1.0 — band-sum and total-energy residuals move in lockstep. Not the partial-cancellation pattern of Class A.
- The 11 eV Fe LDA residual is consistent with a Hamiltonian-level difference in the d-manifold eigenvalues, not a density-basin drift.

**Physics hypothesis:**
1. **LDA (PZ) exchange-correlation on the Fe 3d manifold** is known to be stiffer than PBE for transition metals. If pwdft-rs's Perdew-Zunger implementation uses the original 1981 parametrization while QE uses the Perdew-Wang-92-adapted form even for the LDA path on metals, there would be a per-band d-manifold shift that E_total picks up directly.
2. **Non-local KB projector for the 3d semicore manifold** — projector-l ≥ 2 cases were tested bit-perfect in VGCH-1b's single-channel tests (l=2, m isolation), but Fe LDA specifically may route d-manifold states through a different code path than Fe PBE.
3. **Convention mismatch in the non-spin-polarized LDA-only path** that is bypassed when nspin=2 (Fe PBE stays ferromagnetic per GGAP Phase D #158). Fe LDA is run nspin=1 collapsed to non-magnetic — a specific code path.

**Proposed investigation (Class B diagnostic, Researcher, ~1 CE-day):**
1. **LDA vs PW92 parametrization check.** Read `src/potential/xc.rs::perdew_zunger_correlation` and compare against QE's `qe-7.5/XClib/qe_funct_corr_lda.f90::pw` + `qe-7.5/XClib/qe_funct_corr_lda.f90::pz`. QE uses the **PW92** correlation for its LDA reference calculations by default (not PZ 1981) — if pwdft-rs is using PZ-1981, we're comparing two different functionals.
2. **Compare Fe LDA at nspin=1 vs forced nspin=2 with M=0** — if the residual vanishes in the nspin=2-collapsed-to-NM path, the bug is in the non-spin-only code path, which the spin driver bypasses.
3. **Audit `assemble_v_eff` for LDA vs PBE dispatch** — verify the Hamiltonian assembly is functionally identical between the two XC paths, with only the V_xc-per-point values differing.

**Acceptance criterion:** Fe LDA residual either (a) drops to ≤ 200 meV after one of the three hypotheses closes (GREEN or approach-GREEN), or (b) closes to Class A's expected 1.5-3.3× band-sum ratio (reclassifying Fe LDA into Class A and matching Fe PBE's trajectory). Either outcome clears the "Fe LDA is anomalous" flag.

### Class C — Light-atom C diamond (functional-sensitive, mixer-stall-related)

**Members:** C LDA (1.45 eV), C PBE (0.32 eV). (2 cells.)

**Signature:**
- C PBE is 4.5× closer to QE than C LDA → partially functional-sensitive like Class A heavy-atoms, but at light-atom magnitude.
- `|ΔE_1e| / |ΔE_total|` ≈ 1.2 (near-unity; does NOT show the 1.5-3.3× Class A cancellation pattern).
- **Anderson mixer stalls** on plain C diamond — pinned as a negative regression by MIXA #149. VQEF-QC #144 found Broyden+Kerker converges cleanly in 10-15 iters matching QE's iter count, but to an E_total off by 1.45 eV with the opposite-sign one-e / Hartree split characteristic of Class A.

**Physics hypothesis:**
1. C diamond is a wide-gap covalent insulator. Anderson's fixed-point structure is sensitive to the occupation-discontinuity at the band edge. When it stalls (MIXA pins Δρ ≈ 1.75 × 10⁻⁸ at 150 iters with `MixingMode::Plain`), it's converging to an intermediate density that isn't the variational minimum.
2. With Broyden+Kerker, it does converge — but Kerker's q_TF for an insulator needs to be tuned differently than for a metal (default q_TF auto-estimate assumes a Thomas-Fermi metallic screening that C doesn't have).
3. The residual's reduction under PBE parallels Class A's: gradient-correction-sensitive, suggesting the C density basin QE and pwdft-rs reach differs in the gradient-sensitive sublattice (C's strong sp³ directional bonding).

**Proposed investigation (Class C diagnostic, Researcher, ~0.5 CE-day):**
1. **Kerker q_TF sweep for C.** Try q_TF ∈ {0.2, 0.5, 1.0, 2.0} Å⁻¹ and measure the converged C LDA total energy vs QE. If one value closes to ≤ 200 meV, it's a q_TF-tuning issue; document the insulator-appropriate default and add to CFGN knobs.
2. **Density-transplant on C** (after VGCH-2B lands its Cu infrastructure — reuse the QE ρ parser). If C converges to QE's density basin when transplant-seeded, it's the same Class A mechanism at a smaller magnitude.
3. **Check MIXA Plain-Anderson pathology interaction** — if Anderson's stall at Δρ ≈ 10⁻⁸ is a partial local minimum, Broyden may be settling in the same basin just more cleanly.

**Acceptance criterion:** C LDA residual drops below 500 meV (from 1.45 eV) by a combination of mixer tuning + density-basin fix; C PBE drops below 100 meV (from 0.32 eV). If neither closes, C joins Class A's fate — waiting for VGCH-2B's Part C deeper work.

## Cross-class sequencing (revised post-VGCH-2B)

```
[Class A]  VGCH-2 Part B (LANDED) — H3 mixer-basin CLEARED
              |
              v  +16.34 eV gap at shared density on Cu → functional disagreement
              |
              v  VGCH-2 Part C (NEW)
           Fermi-finder / smearing / n_bands / NLCC ρ_core / V_NL contraction
              |
              v  fix proposal (VGCH-2D, class-split TBD)
           Target: close 4-8 of 8 Class A cells to ≤ 500 meV

[Class B]  Fe LDA diagnostic (NEW)
              |
              v
           Finding: PZ-vs-PW92 / spin-path / dispatch
              |
              v  fix proposal (VGCH-MECH-B1)
           Closes Fe LDA (1 cell); Fe PBE separately handled by Class A

[Class C]  C diamond diagnostic (NEW, after VGCH-2B's QE ρ parser lands — done #167)
              |
              v
           Finding: transplant iter-1 match? / Kerker q_TF? / same as Class A?
              |
              v  fix (Class A-dependent OR separate CFGN knob)
           Closes 2 cells (C LDA, C PBE)
```

## Why this matters

Today's wave shifted the VQEF scoreboard from `1 GREEN / 15 YELLOW / 0 RED` to **`4 GREEN / 12 YELLOW / 0 RED` on E_total + 3 GREEN on E_F** in < 12 hours. The 12 remaining YELLOW cells were previously treated as a monolithic "VGCH heavy-atom residual" requiring a single breakthrough. The mechanism taxonomy reveals that:

- **Class A (8 cells)** needs one breakthrough: a fix to the Kerker / symmetrization / occupation path that makes pwdft-rs converge to QE's density basin.
- **Class B (1 cell, Fe LDA)** likely has an entirely different root cause — a small diagnostic may close it with ~1 CE-day of work.
- **Class C (2 cells, C)** may be a mixer-tuning issue independent of the heavy-atom problem.

Of the remaining 12 YELLOWs, **3 cells can be attacked right now** (Class B + Class C diagnostic) while VGCH-2B works the heavy-atom class. Without the taxonomy, these 3 would have waited behind Class A's breakthrough.

## Acceptance criteria (overall)

1. **Class A:** ≥ 4 of 8 cells close to ≤ 100 meV via VGCH-2B Part B + C (successful transplant + fix).
2. **Class B:** Fe LDA closes to ≤ 200 meV via one of the three hypotheses, OR reclassifies into Class A (band-sum ratio → 1.5-3.3×).
3. **Class C:** C LDA + C PBE close to ≤ 500 / 100 meV via mixer tuning OR density-basin fix.
4. **Post-VGCH-MECH scoreboard**: at least **7 of 12 YELLOW cells flip GREEN**, leaving ≤ 5 YELLOW for a VGCH-MECH-B2 follow-up.

## Out of scope

- GPU PBE shader (GGAP Phase E) — independent track.
- ITEV default flip (Phase-5 step-4 bench) — performance track.
- BSUM per-k extension to catch band-specific residuals — follow-up once scalar gate is stable.
- Clippy baseline drift (18/24 → 21/28 during today's merges) — handled in a separate grooming PR.

## Risks

- **Class A Part C may not close any cells.** If the Fermi-finder / smearing / n_bands triage doesn't localize the 16 eV shared-density gap, the root cause may be in the Hamiltonian construction itself despite VGCH-1b's β_l(q) checks — e.g., the `D_ij · Σ_lm β·β` contraction (the full V_NL operator, not individual projectors) differing from QE's `vnlocal.f90::vloc_psi`. In that case we escalate to (a) single-channel V_NL matrix element cross-check in G-space, or (b) ψ_nk overlap diagnostics (parse QE's wavefunction blocks + compute `⟨ψ^pwdft | ψ^QE⟩` at shared k-point).
- **Class B may reveal a larger rewrite.** If the LDA path has diverged from QE's convention not just in PZ-vs-PW92 but in larger structural ways (e.g., LDA XC applied at a different point in `assemble_v_eff`), the fix could be multi-file.
- **Class C may be a distraction.** If C LDA's 1.45 eV residual is just Class A at smaller magnitude, the "mixer tuning" angle burns Researcher time without closing cells. Mitigation: the transplant on C is a 30-minute reuse of VGCH-2B's Cu infrastructure — do it first; if C's shared-density gap is < 100 meV, C is NOT Class A and mixer tuning becomes worthwhile; if > 1 eV, C is Class A and waits.

## Provenance

- Diagnostic data: VGCH Phase 1a (#139), Phase 1b (#148), Phase 1c (#156), VGCH-2 Part A (#160), BSUM (#165), SiEF-B1 (#166), TSEN (#162), **VGCH-2 Part B (#167 — refuted density-basin hypothesis, identified Fermi-finder/smearing as new suspect)**.
- Eliminated hypotheses: H1 (β_q), H2 (SAD), V_loc(G=0) gauge, smearing reporting, **H3 mixer-basin** — collectively closed 5 potential root causes.
- Scoreboard evolution today: `0 GREEN / 8 YELLOW / 8 RED` → `4 GREEN / 12 YELLOW / 0 RED` on E_total + 3 new GREEN E_F cells. **+16.34 eV shared-density gap on Cu** is the new quantitative target for Part C.
- Active investigation as of this proposal: VGCH-2 Part C (Fermi-finder / smearing / n_bands triage).
