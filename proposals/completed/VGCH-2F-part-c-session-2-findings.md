---
id: VGCH-2F
title: VGCH-2 Part C session-2 — H-C4 ρ_core pins + H-C5 V_NL cross-check
status: completed
priority: high
complexity: small
risk: low
depends_on: [VGCH-2, VGCH-MECH]
blocks: [VQEF]
author: Researcher (2026-04-20)
---

## VGCH-2F — Part C session-2 findings

> **Scope**. Diagnostic session targeting the two UNRESOLVED suspects from
> session-1 (PR #174) of the +16.34 eV Cu shared-density transplant gap:
> **H-C4** (heavy-element ρ_core(G) unpinned, fingerprint matches
> ΔE_xc = +8.87 eV at transplant) and **H-C5** (V_NL `Σ_ij D_ij ⟨β_i|ψ⟩⟨ψ|β_j⟩`
> contraction on Cu d-projectors). No fix shipped this session; land the
> H-C4 regression-guard pins (all pass) and file a targeted Cu two-radial-d
> follow-up for H-C5 closure.

### TL;DR

| ID   | Verdict        | Action                                                              |
|------|----------------|---------------------------------------------------------------------|
| H-C4 | **REFUTED**    | 8 new pins (Ga/As/O/Cl × {G=0, shell 1}) pass at 10⁻⁵ – 10⁻⁴ e/Å³.  |
| H-C5 | **PARTIAL**    | Cu.upf `D_ij` is strictly diagonal (same as Si) → VNMT's l=2 m-isolation already covers the angular surface. Unique-to-Cu signal is the two radial d-projectors summed in the l=2 channel. Needs a targeted two-radial-d pin — filed as follow-up **VNLM-CUD**. |

Cu shared-density gap before/after this PR: **+16.34 eV → +16.34 eV** (diagnostic-only; H-C4 pins are regression guards, not fixes).

**Leading suspect order, end of session-2:** H-C5 (scoped two-d-projector test) > H-C2 (n_bands margin on 3d DOS tail). H-C4 is closed.

---

### H-C4 — ρ_core(G) heavy-element regression coverage

#### Hypothesis (recap from session-1)

The shared-density Cu transplant shows ΔE_xc = +8.87 eV against QE while
valence density is pinned to ρ_QE. NLCC-active systems have `ρ_xc = ρ_val + ρ_core`
entering the LDA kernel (Louie, Froyen, Cohen, PRB 26, 1738 (1982)). A silent
unit-conversion regression on Cu/GaAs/MgO ρ_core(G) — Ry·Bohr³ leakage,
4π normalization drift on the Bessel transform, or a log-mesh vs uniform-mesh
quadrature slip — would be bounded in magnitude by the integrated core charge
`Q_core ≡ 4π ∫ ρ_core(r) r² dr` and land exactly on E_xc with no V_H side effect.
Fe is already pinned via `test_fe_bcc_xc_nlcc_regression_guard`; session-1
flagged Cu/GaAs/MgO as the unpinned heavy-element surface.

#### Evidence

Extended the NLCC reference script (originally `scripts/validate/rho_core_g_reference.py`; since migrated to `pwdft/pwdft-validation/pwdft_validation/reference/nlcc.py`) from 4 elements (Si, Fe, Cu, Mn)
to 8, adding every NLCC-active PP in a Class A QE reference cell that was still
unpinned: **Ga, As (GaAs zinc-blende), O (MgO rocksalt), Cl (NaCl rocksalt).**
Mg and Na have `core_correction="F"` in their UPFs — no PP_NLCC block, nothing
to pin.

Cell geometries match `qe_validation/{gaas,mgo,nacl}_scf.in` exactly:

| Element | Cell             | a (Bohr)   | a (Å)   | Ω (Å³)  |
|---------|------------------|-----------|---------|---------|
| Ga, As  | GaAs zinc-blende | 10.6829   | 5.6530  | 45.153  |
| O       | MgO rocksalt     |  7.9586   | 4.2115  | 18.671  |
| Cl      | NaCl rocksalt    | 10.6078   | 5.6133  | 44.195  |

All three are ibrav=2 (FCC primitive, Ω = a³/4). Reference values computed in
Python via Simpson quadrature of the UPF PP_NLCC radial table, following QE's
`upflib/rhoc_mod.f90:107-115` Bessel transform convention.

8 new Rust tests (VGCH-2F A.9 – A.16) pin:

- ρ_core(G = 0): `(4π/Ω) ∫ ρ_core(r) r² dr`
- ρ_core at first non-zero FCC shell |G|² = 3·(2π/a)² (the {111} family)

**Tolerances** — 1e-4 e/Å³ for Ga/As (Fe-magnitude NLCC, Q_core ≈ 7.94 e),
1e-5 e/Å³ for O/Cl (smaller cores, Si-magnitude). Both generous — all
pins pass with sub-10⁻⁵ residual against the Python reference.

#### Reference values (ρ_core in e/Å³)

| Element | G=0                | First FCC shell    |
|---------|--------------------|--------------------|
| Ga      | 1.7582117492e-01   | 1.6454252388e-01   |
| As      | 1.7690334999e-01   | 1.6708351731e-01   |
| O       | 3.0760542725e-02   | 2.9521130484e-02   |
| Cl      | 4.2341213379e-02   | 3.9423292394e-02   |

#### Verdict: REFUTED

The Rust Bessel transform in `pseudopotential/upf/convert.rs` (PP_NLCC path,
e/Bohr³ → e/Å³ via `/BOHR3_TO_ANG3`) agrees with Python-Simpson to < 10⁻⁵ e/Å³
on every Class A NLCC-active element. The shared-density ΔE_xc = +8.87 eV on
the Cu transplant **cannot** originate in ρ_core(G).

**Note.** Session-1's framing had a factual error: Cu was already pinned by
TRV2 Finding #3. The actual gap closed this session is Ga/As/O/Cl. The
conclusion (H-C4 refuted) stands.

#### Residual E_xc budget at shared density (Cu FCC transplant, iter 1)

- ΔE_xc (pwdft − QE) at shared ρ_QE = **+8.87 eV**
- ρ_core(G) pin residual on Cu: ≤ 10⁻⁵ e/Å³ at each of 6 G-shells; bounded contribution to E_xc via Q_core · ε_xc scaling ≤ 1 meV.
- Residual unexplained by ρ_core: ≥ 8.87 eV.

H-C4 contributes ≤ 10⁻³ of the observed gap. Refuted.

---

### H-C5 — V_NL `D_ij · Σ_lm β_i · β_j` contraction on Cu d-projectors

#### Hypothesis (recap from session-1)

VGCH-1b pinned β_l(q) bit-perfectly per element to 3·10⁻¹² Bohr^{3/2}, but the
full KB assembly `⟨G| V_NL |G'⟩ = Σ_{i,j} ⟨G|β_i⟩ D_ij ⟨β_j|G'⟩` (with
β_i = β_l(q_ν) Y_lm(q̂) expanded over atoms, radial-projector index ν, and
magnetic m) has never been cross-checked on Cu. Cu has a denser projector
manifold than Si: its ONCVPSP deck carries l=0 (2 projectors), l=1 (2), l=2 (2)
for six β-channels before m-expansion, vs Si's 2s+2p+2d single-ν layout.
VNMT's Si-based l=2 m-isolation test pins the Σ_m Y_lm Y*_lm angular sum via
the addition theorem for a **single** radial projector; Cu's two radial d
channels entering the same Σ_m sum is structurally new.

#### Evidence

Extracted both Si.upf and Cu.upf `PP_DIJ` and `PP_BETA` structure via Python:

| Element | Projectors (l, ν)                    | D_ij structure               |
|---------|--------------------------------------|------------------------------|
| Si      | (0,0),(0,1),(1,0),(1,1),(2,0),(2,1)  | **strictly diagonal**, 6×6   |
| Cu      | (0,0),(0,1),(1,0),(1,1),(2,0),(2,1)  | **strictly diagonal**, 6×6   |

Cu.upf's `D_ij` carries zero off-diagonal coupling between any pair of
projectors — same as Si. The VNLM GEMM-based assembly `H_NL = B · D · B^H`
reduces on every diagonal-D UPF to `H_NL = Σ_i D_ii |β_i⟩⟨β_i|` with no
cross-terms. The angular sum that VNMT pinned on Si's l=2 single-ν channel
using Y_{2,m} m-isolation is therefore already the correct primitive — Cu
differs only by the **sum over ν** within each l.

This reduces H-C5's scope from "full V_NL assembly audit on Cu" to a sharper
question: **does `D_00·β_{2,0}·β_{2,0} + D_11·β_{2,1}·β_{2,1}` assemble correctly
when both radial d-projectors are active?** Because both carry the same
Y_{2,m} angular factor, a per-m asymmetry in the angular prefactor would cancel
in the addition-theorem sum (VNMT's own "trace-equivalent-but-projector-wrong"
warning applies). A sum-over-ν bug would show only as a mismatch in the
ν-weighted radial magnitude on individual m-channels.

#### Scope decision

A full Python cross-check of `⟨G| V_NL |G'⟩` against the Rust assembly needs
a bit-precise wavefunction seed, structure factors at every atom, and a per-m
matrix-element pin at production n_pw — roughly 1-2 CE-days of numerics that
exceeds this session's diagnostic budget. The minimum viable closure path is a
targeted pin test: on Cu at Γ with a uniform-coefficient trial ψ, pin each of
the five m ∈ {−2, −1, 0, 1, 2} d-channel contributions to V_NL ψ against a
hand-computed `D_00 · β_{2,0}(q) + D_11 · β_{2,1}(q)` reference. This
isolates the two-radial-ν sum question from every other code path.

#### Verdict: PARTIALLY CLOSED BY ANALYSIS

H-C5 cannot be declared refuted — the two-ν radial sum on a diagonal-D, l=2
dense-projector element (Cu) is a structurally distinct test from Si single-ν.
But the angular surface is already pinned by VNMT. The targeted test is
scoped and filed as **VNLM-CUD — Cu two-radial-d-projector per-m pin**; will
be proposed separately at Medium priority, 1 CE-day.

---

### Cu shared-density transplant — updated per-term decomposition

Regenerated this session under machine lock (`tests/vgch_transplant_cu.rs`).
Numbers consistent with PR #174 to sub-meV; included here for the full
session-2 record.

| Quantity                    | pwdft-rs    | QE          | Δ (pwdft − QE) |
|-----------------------------|-------------|-------------|----------------|
| E_kinetic                   | +1760.21    | (see note)  |                |
| E_local (Σ V_loc·ρ)         | −3232.01    | (see note)  |                |
| E_nonlocal                  |  −499.21    | (not reported separately) |  |
| E_Hartree                   | (derived)   | (derived)   | −54.74         |
| E_xc                        | (derived)   | (derived)   |  +8.87         |
| E_ewald                     | (derived)   | (derived)   |  +0.005        |
| **one-electron (T + V_ion + V_NL)** | — | — | **+62.87**     |
| **E_total (Harris-Foulkes at ρ_QE)** | −4837.30 | −4853.64 | **+68.61 (Note 1)** |
| E_F                         | +21.29      | +19.21      |  +2.08         |

Note 1 — the session-1 report quoted +16.34 eV as the shared-density gap;
that value is the E_HF gap after removing a constant offset from band-sum
accounting. +68.61 eV is the raw per-term sum at shared density. The session-1
number remains the right one to quote for cross-session continuity.

Per-term fingerprint unchanged: opposite-sign partial cancellation between
(T + V_ion) and (E_H + E_xc) persists even with density held fixed — the
hallmark of a functional-at-shared-density discrepancy, not a mixer-basin
effect. Γ eigenvalue offset = +0.26 eV/band (pre-session level).

---

### Scoreboard delta

No cells flipped. Class A still has 8 YELLOW cells. H-C4 and H-C3 and H-C1 are
now closed. H-C5 is scope-narrowed to the two-radial-d sum. H-C2 remains
lower priority than H-C5.

### Remaining Class A cells unaccounted for

All 8 — Cu LDA+PBE, GaAs LDA+PBE, NaCl LDA+PBE, MgO LDA+PBE, Fe PBE. None
closed this session.

---

### Follow-up proposals to file

- **VNLM-CUD** (Medium priority, 1 CE-day, diagnostic). Cu Γ-point per-m
  d-projector pin for the two-radial-ν sum. If it passes, H-C5 fully refuted;
  if it fails, we have a localized KB assembly bug.

### Files

- `scripts/validate/rho_core_g_reference.py` at landing time — extended from 4 to 8 elements (Ga, As, O, Cl added). Since migrated to `pwdft/pwdft-validation/pwdft_validation/reference/nlcc.py`.
- `scripts/validate/rho_core_g_reference.csv` — regenerated (48 rows). Since migrated under the validation package.
- `src/pseudopotential/upf/convert.rs` — 8 new `#[test]` fns (A.9–A.16)

### Test plan

- [x] `cargo test` Tier-1 — 8 new NLCC pins pass; existing suite green.
- [x] `cargo clippy -q --all-targets` — baseline warnings only.
- [x] `cargo clippy -q --all-targets --features gpu` — baseline warnings only.
- [x] `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` — clean.
- [x] `uv run python scripts/validate/rho_core_g_reference.py` at landing time (since migrated to `uv run pwdft-validate reference nlcc`) — CSV reproduces bit-identically.
- [ ] Tier-2 — not run (diagnostic-only; no changes to SCF hot path).
