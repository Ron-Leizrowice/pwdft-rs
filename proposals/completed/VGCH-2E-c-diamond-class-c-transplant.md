---
id: VGCH-2E
title: VGCH-2E — C diamond Class C transplant diagnostic (mechanism localization)
status: completed
priority: high
complexity: small
risk: low
depends_on: [VGCH-MECH, VGCH-2, VGCH-2B]
blocks: [VQEF]
author: Researcher (2026-04-20, VGCH-MECH Class C diagnostic)
---

## VGCH-2E — C diamond Class C transplant diagnostic

### One-line finding

**C diamond is Class A at light-atom magnitude.** The VGCH-2B transplant
harness seeded pwdft-rs with QE's converged C-diamond density ρ_QE and
ran one SCF iteration; the per-term residuals show the same opposite-sign
one-electron / Hartree partial cancellation fingerprint as Cu, Fe PBE,
GaAs, NaCl, MgO — just 10× smaller.

### What was run

- QE regen: 4×4×4 Γ-centered, ecut=30 Ry, Fermi-Dirac σ=0.01 Ry,
  `disk_io='medium'` (the committed `qe_validation/c_diamond_scf.in`
  uses `low` and does not write `charge-density.dat`). Input staged at
  `/tmp/vgch2e_c/c.in`, run through the machine lock with 8 MPI ranks,
  9 iters to `conv_thr = 1e-8` — matches the committed reference
  (ΔE_total = 0.03 meV, same one-e/Hartree/XC/Ewald per-term to 6
  digits).
- Density parse:
  `scripts/validate/vgch2_parse_qe_density.py --rho
  /tmp/vgch2e_c/tmp/c.save/charge-density.dat --out
  /tmp/vgch2e_c/c_rho_qe.bin --verbose`.
  Output: ngm=1687, Miller range ±8, ρ(G=0)·Ω_Bohr = 0.10447·76.576 =
  **8.000 electrons** (spot-on; C has 8 valence electrons).
- Test harness: `tests/vgch_transplant_c.rs` (new; Tier-2, mirrors
  `vgch_transplant_cu.rs` byte-for-byte except crystal/PP/FFT-grid).
- Diagnostic wall: 1 min compile + 0.25 s run.

### Per-component residuals at shared density ρ_QE (iter-1)

Units: eV. `Δ = pwdft_iter1 − QE_converged`.

| term         | pwdft (iter-1) | QE (converged) | Δ (ours − QE) |
|--------------|---------------:|---------------:|--------------:|
| one-electron |      +117.5531 |      +115.6875 |      **+1.866** |
| Hartree      |       +24.0465 |       +24.8764 |      **−0.830** |
| XC (bare)    |      −116.6182 |      −117.0471 |      **+0.429** |
| Ewald        |      −347.9360 |      −347.9233 |        −0.013 |
| **Total**    |      −322.1303 |      −324.4065 |      **+2.276** |

One-electron breakdown (eV):

- E_kin = +216.2223
- E_loc (G≠0) = −91.7539
- E_loc (G=0)·N_el = 0 (post-VGCH-SiEF-B1 gauge)
- E_NL = −6.9152

Γ eigenvalues (eV):

| band | pwdft iter-1 | QE conv | Δ (meV) |
|------|-------------:|--------:|--------:|
| 0    |      −8.1393 | −8.1456 |    +6.3 |
| 1    |     +14.0367 | +14.0232 |   +13.5 |
| 2    |     +14.0367 | +14.0232 |   +13.5 |
| 3    |     +14.0367 | +14.0232 |   +13.5 |
| 4    |     +19.3667 | +19.3568 |    +9.9 |
| 5    |     +19.3667 | +19.3568 |    +9.9 |
| 6    |     +19.3667 | +19.3568 |    +9.9 |
| 7    |     +27.2652 | +27.2505 |   +14.7 |

Fermi energy: pwdft 15.8987 vs QE 15.8873 → **ΔE_F = +11.4 meV**.
Δρ_rms (in-out) = 2.16e-2 e/Å³. `|E_HF − E_KS|` = 0.830 eV (the E_HF
pair-gauge mismatch scales with the Hartree residual, as expected).

### Classification: Class A mechanism at light-atom magnitude

Two observations pin this:

#### 1. The opposite-sign ΔE_1e / −ΔE_H signature is present

`|ΔE_1e| / |−ΔE_H|` = **1.866 / 0.830 = 2.25** — inside the Class A
range [1.2, 2.9] documented in VGCH-MECH (Cu 2.9×, MgO 2.2×, NaCl 1.8×,
GaAs 1.5×, etc.). This is the quantitative fingerprint of "functional
disagrees at the same density" that VGCH-2B identified on Cu.

#### 2. The gap is 2.28 eV at shared ρ, but only 1.45 eV at independent convergence

This is a new constraint. Reported in the `test_c_diamond_vs_qe`
docstring, the SCF-converged per-term split is:
`Δ one-e = +1.76, Δ E_H = −0.59, Δ E_xc = +0.29, ΔE_total = +1.45` eV.
At shared ρ_QE the total gap is **1.57× larger** (+2.28 eV). In other
words: the mixer does NOT drive pwdft-rs toward QE's density basin —
pwdft-rs's self-consistent density is *further* from ρ_QE than ρ_QE is
from pwdft-rs's fixed-point-evaluated functional value. This is the
opposite sign from what a pure "mixer basin" hypothesis (Class C
original) would predict: if pwdft-rs were hung up in a local minimum,
seeding ρ_QE should pull the per-component gap *down*, not up.

#### 3. E_F is essentially correct — Class C is NOT Fermi-finder

11.4 meV matches Fermi-finder noise. Compare Cu shared-density E_F gap
of 2.08 eV (1.82 eV of it pure DOS/occupation origin per VGCH Part C
session-1). C diamond's E_F at ρ_QE is bit-correct to smearing noise.
Removes one of three hypotheses I came in with.

#### 4. Anderson mixer stall on Plain is irrelevant

VQEF-QC (#144) found Plain Anderson stalls at Δρ≈1.75e-8 at 150 iters;
Broyden+Kerker converges in 12 iters to the same 1.45 eV residual as
Kerker / Broyden / PeriodicPulay. This transplant uses Broyden+Kerker
and reaches the shared-density point in one iteration with
Δρ = 2.16e-2 e/Å³ — the mixer-stall angle was a red herring (MIXA
pinned it as a conditioning regression, not a physics residual).

### Revised mechanism hypothesis

**Primary (most likely):** C diamond carries the same mechanism as
Class A heavy-atom cells, just at light-atom magnitude. The specific
candidates that now apply to BOTH Class A and Class C:

1. **ρ_core(G) Bessel-transform (NLCC path).** C LDA's PP has
   `core_correction="T"` and uses the same
   `scf::potentials::compute_core_density` path as Cu/Fe/GaAs/Mg. If
   that helper has a silent unit or r²-weight regression that survived
   NCFX (#40) for some heavy-semicore subset, it should hit C too — Z=6
   just makes the magnitude smaller. ΔE_xc = +0.429 eV is consistent
   with a small-magnitude NLCC leak (compare Cu's +8.87 eV, Fe LDA's
   +0.93 eV pre-NCFX).
2. **V_NL `D_ij · Σ_lm β·β` contraction.** C has s+p semilocal
   projectors only (l=0,0,1,1 per the QE output). If the double-sum
   assembly differs from QE's `vnlocal.f90::vloc_psi` on the l=1
   manifold in a way that's absorbed into Si's smaller E_NL magnitude,
   C would be the cleanest light-atom witness. C's E_NL = −6.92 eV
   (pwdft) is 36× smaller than Cu's E_NL magnitude but the 1.87 eV
   one-electron residual means any percentage-wise drift in E_NL would
   be readable.
3. **Occupation / band-edge integration on an insulator.** C is a
   wide-gap insulator; QE and pwdft-rs both use Fermi-Dirac at σ=0.01
   Ry on a 4×4×4 grid. At σ=0.01 Ry (0.136 eV) the band-edge is
   sharply resolved; if the per-k occupations at the valence top (the
   triply-degenerate 14.02 eV cluster) differ by O(10⁻⁴) e between
   codes, that's O(meV) on one-electron — not enough to account for
   1.87 eV. Downgraded to tertiary.

**Refuted or downgraded by this transplant:**

- **Mixer basin (original Class C H3):** contradicted by the
  2.28 > 1.45 eV gap inversion. Not a local-minimum pathology.
- **Kerker q_TF default on an insulator:** the transplant uses
  Broyden+Kerker default (q_TF auto) and reaches the correct ρ_QE
  seed; the one-iter response is a Hamiltonian-at-shared-ρ result, not
  a mixer-trajectory result. Kerker tuning cannot fix a 2.28 eV gap
  that's already locked in at iter-0.
- **Fermi-finder / smearing function:** 11 meV ΔE_F is well below
  even Cu's resolved finder noise (VGCH-2C ~μeV on converged
  eigenvalues + weights).

### Proposed follow-up experiment to confirm NLCC hypothesis

The one-line test: **rerun C transplant with the same harness but at
ecut=60 Ry AND n_bands=16 AND NLCC cross-pinned against QE.**

Concretely:

1. **Pin ρ_core(G) for C against a SciPy reference** — add C to the
   same `test_fe_bcc_xc_nlcc_regression_guard` battery Fe has. The
   UPF `PP_NLCC` for C reads `size="1234"` on the same log mesh as the
   β_l(q) projectors; apply
   `scripts/validate/rho_core_g_reference.py` at 10 q-values in
   [0, 7] Bohr⁻¹ and compare against pwdft-rs'
   `compute_core_density` output on the 20×20×20 grid.
   If |Δρ_core(G)| > 1e-8 e/Bohr³ at any G, that's the bug.
2. **If (1) is bit-perfect**, run the C transplant at the same
   geometry but disable NLCC: patch `PP_NLCC` in a copy of the C UPF
   file to constant-zero (no valence impact since C is a light atom).
   If the shared-density residual drops from 2.28 eV to under 500 meV,
   the NLCC path is guilty and we escalate to auditing the Bessel
   transform of ρ_core at the specific C log-mesh shape.
3. **If NLCC is cleared**, the only remaining Class C suspect is
   Sub-hypothesis 2 (V_NL l=1 double-sum). That would require a
   single-channel V_NL matrix-element cross-check on C at Γ
   (paralleling the VGCMP Phase 4 Si test but for the l=1 block of
   the C projector set).

Total: ~0.5 CE-day for steps 1+2. Step 3 is ~1.5 CE-days if needed.

### Acceptance criteria

1. Per-component residual table committed to this proposal (done — see
   above).
2. Class assignment documented: C diamond is Class A at light-atom
   magnitude — the MECH taxonomy's three-class split now reduces to
   **Class A (9 cells, +C LDA absorbed)** + **Class B (Fe LDA)**. C PBE
   (0.32 eV) likely follows the same reclassification — to be confirmed
   by a parallel PBE transplant (same harness, same infrastructure;
   ~15 min additional run).
3. If the follow-up (NLCC pin + ablation, step 1+2) closes C LDA to
   ≤ 500 meV, that same fix-site is the prime candidate for the
   remaining 8 Class A cells — resolving C becomes a scoping lever
   for Class A.
4. If step 1+2 does not close C, the `D_ij · β·β` V_NL contraction
   hypothesis inherits the case and VGCH-2 Part C's ordering needs
   re-shuffling (H-C5 moves above H-C4).

### Artifacts (this proposal)

- `tests/vgch_transplant_c.rs` — Tier-2 test, 280 LOC, mirrors the
  VGCH-2B Cu harness. Gated on QE density bundle at
  `/tmp/vgch2e_c/c_rho_qe.bin` or `qe_validation/c_rho_qe.bin`.
- `qe_validation/c_rho_qe.bin` — VGCH2BIN flat binary, 47329 bytes,
  ngm=1687, nspin=1, C LDA converged density.
- `/tmp/vgch2e_c/c.in` + `/tmp/vgch2e_c/c.out` — QE regen for
  reproducibility (not committed; recipe in
  `tests/vgch_transplant_c.rs` header).

### Out of scope

- C PBE transplant — a 15-minute rerun with `pp_c = load_pp` pointed
  at `pseudopotentials/nc/pbe/C.upf` and `XcFunctional::Pbe` on the
  params. Expected to show ΔE_total at shared ρ of ~500 meV
  (C LDA 2.28 eV × C PBE/LDA SCF-converged ratio 0.22 = 0.5 eV).
  Defer to a VGCH-2E-PBE follow-up if the ratio needs confirmation.
- Implementation of the NLCC pin or ablation — separate CE-proposal
  once this diagnostic's hypothesis is confirmed.
- Any changes to `src/scf/transplant.rs` — the existing harness was
  designed to be geometry-agnostic and handled C as-is.

### Provenance

- VGCH-MECH (#168) — Class C scope carved out, "might be Class A at
  smaller magnitude" flagged as the mitigation risk.
- VGCH-2B (#167) — Cu transplant infrastructure; this proposal reuses
  it unchanged.
- VGCH-2 Part C (#174) — Fermi-finder refutation methodology carried
  over; ΔE_F on C confirms the same finder is correct.
- QE regeneration run 2026-04-20, machine-locked, 8s wall.
- Transplant iter-1 wall: 0.25 s after compile on M3 Max.
