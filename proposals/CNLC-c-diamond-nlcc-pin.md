---
id: CNLC
title: CNLC — C diamond NLCC ρ_core(G) pin + NLCC-off ablation diagnostic
status: active
priority: high
complexity: small
risk: low
depends_on: [VGCH-2E]
blocks: [VQEF]
owner: core-engineer
author: Researcher (2026-04-20)
---

## CNLC — C diamond NLCC pin + NLCC-off ablation

### Problem

VGCH-2E (PR #178, `proposals/VGCH-2E-c-diamond-class-c-transplant.md`)
reclassified C diamond from Class C (mixer-basin suspect) to **Class A
at light-atom magnitude**: the shared-density transplant shows the same
opposite-sign ΔE_1e / ΔE_H partial-cancellation fingerprint as Cu / Fe /
GaAs / NaCl / MgO, just 10× smaller. Per-term residuals at ρ = ρ_QE,
iter-1:

| term         | Δ = pwdft − QE |
|--------------|---------------:|
| one-electron |    +1.866 eV   |
| Hartree      |    −0.830 eV   |
| **XC (bare)**|    **+0.429 eV** |
| Ewald        |    −0.013 eV   |
| **Total**    |    **+2.276 eV** |

`|ΔE_1e| / |−ΔE_H| = 2.25` — inside the Class A range [1.2, 2.9].

The primary suspect per VGCH-2E § Revised mechanism hypothesis is the
**NLCC ρ_core(G) Bessel transform** on C's UPF, which declares
`core_correction="T"` (verified via
`grep core_correction pseudopotentials/nc/lda/C.upf` → `core_correction="T"`).
C is the only Class A light-atom cell where NLCC is active — Si LDA /
PBE, Al LDA / PBE, NaCl LDA / PBE (Na and Cl both have
`core_correction="F"`), and MgO LDA (Mg `core_correction="F"`) are
all NLCC-inactive on the light-atom side of their cell.

**Why this is plausible.** NCFX (PR #40, completed) closed the 13.4 eV
Si gap by fixing unit conversion + `r²·4π` weighting in the ρ_core(r→G)
Bessel transform. That fix was validated by regression pins for Si and
Fe (`data/csv/rho_core_g_reference.csv` + the `test_fe_bcc_xc_nlcc_
regression_guard` E_xc guard). TRV2 (PR #98) extended the ρ_core(G)
pin to Cu and Mn. VGCH-2F Part C session-2 (PR #179) extended it to
Ga, As, O, Cl — 8 elements now pinned at sub-1e-4 e/Å³.

**C is the last NLCC-active UPF in the VQEF matrix that is unpinned.**
A small systematic regression of the same bug class that NCFX closed —
e.g. a log-mesh-specific quadrature slip on C's specific `PP_NLCC`
mesh shape — would fit the +0.429 eV ΔE_xc fingerprint.

**Magnitude estimate.** C's integrated core charge Q_core is
small (Z_core = 2 for 1s²; Q_core from `PP_NLCC` ≈ 0.6 e after
the pseudized-core smearing). The NLCC entry to LDA XC is
`ε_xc[ρ_val + ρ_core]` (Louie, Froyen, Cohen, *Phys. Rev. B*
**26**, 1738 (1982), Eq. 3–4). For a small ρ_core(G) error δρ_c(G),
the ΔE_xc perturbation scales as ≈ (dV_xc/dρ) · δρ_c · Ω. At C's
valence density a ΔE_xc = 0.43 eV corresponds to δρ_c ≈ O(10⁻³) e/Å³
integrated magnitude — well within the error budget of a silent
quadrature-slip class of bug and below the current pinned-element
tolerance floor.

### Research

#### What's already landed

`pwdft/pwdft-core/src/pseudopotential/upf/convert.rs::tests`:

| fn (line) | element | shells pinned | tolerance |
|-----------|---------|---------------|-----------|
| `test_si_rho_core_of_g_zero` (303) | Si | G=0 | 1e-5 e/Å³ |
| `test_si_rho_core_of_g_first_shell` (329) | Si | {111} FCC | 1e-5 |
| `test_fe_rho_core_of_g_zero` (362) | Fe | G=0 | 1e-4 |
| `test_fe_rho_core_of_g_first_shell` (387) | Fe | {110} BCC | 1e-4 |
| `test_cu_rho_core_of_g_{zero,first_shell}` (430, 454) | Cu | G=0, {111} | 1e-4 |
| `test_mn_rho_core_of_g_{zero,first_shell}` (489, 513) | Mn | G=0, {110} | 1e-4 |
| `test_ga_*`, `test_as_*` (565–653) | Ga, As | G=0, {111} | 1e-4 |
| `test_o_*`, `test_cl_*` (674–760) | O, Cl | G=0, {111} | 1e-5 |

All pins pass. The Rust helper `rho_core_of_g_ang(pp, g_norm, omega)`
(convert.rs:279–292) is the production path used by
`scf::potentials::compute_core_density` (`scf/potentials.rs:77`).

Python reference generator:
`pwdft/pwdft-validation/pwdft_validation/reference/nlcc.py` — a
Simpson quadrature over the UPF radial grid, following QE's
`upflib/rhoc_mod.f90:107-115` convention. Data output lands in
`data/csv/rho_core_g_reference.csv`. Extending this to C is an
8-line addition to `_SYSTEMS` (nlcc.py:64–82).

#### C.upf structural facts

- Element header: `element="C "` (nlcc.py's dataclass maps
  `element = "C"`).
- `core_correction="T"` — confirmed.
- Cell: diamond (FCC primitive), `ibrav=2`, `celldm(1) = 6.7409 Bohr
  = 3.5672 Å` (from `data/qe/c_diamond_scf.in`), Ω = a³/4 = 11.34 Å³.
- UPF radial mesh: log mesh (standard ONCVPSP shape). `size`
  attribute from `PP_NLCC` matches `PP_MESH/PP_R` `size`.

The FCC first non-zero |G| shell is the {111} family at
|G| = 2π/a · √3 = 3.0516 Å⁻¹. This is larger than Si's first shell
(a = 5.431 Å, |G| = 2.004 Å⁻¹) because C's diamond lattice constant
is ~66% of Si's — the Bessel-transform is probed at a tighter
q-range, which is a *different* part of the integrand than the 8
already-pinned elements. This is material: a log-mesh quadrature
regression tends to bite more at high q (small real-space r), so
C's {111} pin is in the more-sensitive half of the integrand space.

### Implementation

Three phases, landing in order.

#### Phase 1 — add C ρ_core(G) pin (half day)

##### Python reference side

In `pwdft/pwdft-validation/pwdft_validation/reference/nlcc.py`, add
C to `_SYSTEMS` (line 82, after Cl):

```python
# CNLC — C diamond NLCC pin closes the last NLCC-active UPF
# in the VQEF matrix. Primary suspect for VGCH-2E Class A
# light-atom residual (ΔE_xc = +0.43 eV at shared ρ).
# Lattice from data/qe/c_diamond_scf.in: celldm(1) = 6.7409 Bohr.
_NlccSystem("c", "C", 6.7409 * BOHR_TO_ANG, "fcc"),
```

Regenerate the CSV via
`uv run pwdft-validate reference nlcc --pseudo-dir pseudopotentials/nc/lda`
(invocation pattern per existing CSV regeneration in VGCH-2F §
Files). Commit the regenerated `data/csv/rho_core_g_reference.csv`.

##### Rust pin side

In `pwdft/pwdft-core/src/pseudopotential/upf/convert.rs` (around
line 760, after the Cl block), add:

```rust
fn c_content() -> String {
    std::fs::read_to_string(
        std::path::PathBuf::from(env!("CARGO_WORKSPACE_DIR"))
            .join("pseudopotentials/nc/lda/C.upf"),
    ).unwrap()
}

/// CNLC A.1: pin C ρ_core(G=0) — last NLCC-active UPF in VQEF.
/// FCC diamond a = 3.5672 Å (celldm(1) = 6.7409 Bohr in
/// data/qe/c_diamond_scf.in), Ω = a³/4 = 11.34 Å³.
#[test]
fn test_c_rho_core_of_g_zero() {
    let pp = parse(&c_content()).unwrap();
    assert!(pp.has_nlcc());
    let a = 6.7409 * crate::consts::BOHR_TO_ANG;
    let omega = a * a * a / 4.0;
    let rho_g0 = rho_core_of_g_ang(&pp, 0.0, omega);
    let expected = <TBD-from-CSV-row-c-shell-0>;
    assert!(
        (rho_g0 - expected).abs() < 1.0e-5,
        "C ρ_core(G=0) = {rho_g0:.8e}, expected {expected:.8e} e/Å³"
    );
}

/// CNLC A.2: pin C ρ_core(G≠0) at |G|² = 3·(2π/a)² ({111}).
#[test]
fn test_c_rho_core_of_g_first_shell() {
    let pp = parse(&c_content()).unwrap();
    assert!(pp.has_nlcc());
    let a = 6.7409 * crate::consts::BOHR_TO_ANG;
    let omega = a * a * a / 4.0;
    let g_norm = 2.0 * std::f64::consts::PI / a * 3.0_f64.sqrt();
    let rho_g = rho_core_of_g_ang(&pp, g_norm, omega);
    let expected = <TBD-from-CSV-row-c-shell-1>;
    assert!(
        (rho_g - expected).abs() < 1.0e-5,
        "C ρ_core({{111}}) = {rho_g:.8e}, expected {expected:.8e} e/Å³"
    );
}
```

The `<TBD>` expected values come from the regenerated CSV. Tolerance
1e-5 matches the O/Cl pin (C's Q_core is on the same small order as
O's light-atom core).

**Acceptance for Phase 1:** both tests pass in Tier-1. If either
fails at 1e-5, inspect ΔE_xc from the diagnostic (Phase 3) — that's
the direct NLCC-path attribution.

#### Phase 2 — extend ρ_core(G) coverage to multiple FCC shells on C (optional, half day)

Same rationale as the VGCH-2F regression-guard philosophy: a single
{111} shell pin might miss a high-q quadrature slip. Add a 6-shell
pin test parametric over `(l² = 3, 4, 8, 11, 12, 16)` |G|² values
(the first six FCC shells). Tolerance 1e-5 across all shells.

This step is optional; land it if (Phase 1 passes cleanly) and
(Phase 3 shows the residual is NOT NLCC). It adds confidence that
C's NLCC path is definitively clear.

#### Phase 3 — NLCC-off ablation on C (half day)

A direct attribution experiment. Two implementation options; option
A is preferred (no UPF mutation, stays in-tree):

##### Option A — temporary `rho_core_g.fill(0.0)` branch

Add a **diagnostic-only** knob to `ScfParams` (e.g.
`scf.disable_nlcc: Option<bool>`, `#[doc(hidden)]`, defaulting to
false). When true, `scf::potentials::compute_core_density`
(`scf/potentials.rs:77`) returns a zero-filled ρ_core(r) on the
FFT grid. All other NLCC bookkeeping (double-counting in E_xc,
V_xc assembly) remains wired; zeroing ρ_core at the Bessel-transform
entry point is the cleanest place to isolate the effect.

Run the C LDA QE-validation test case twice:

1. `disable_nlcc = false` (production path) — reproduces PZPW-era
   VQEF C LDA residual (currently +1.45 eV E_total at independent
   convergence, +2.28 eV at shared-density transplant).
2. `disable_nlcc = true` — same C LDA deck, NLCC suppressed.
   Compare QE reference is the same (NLCC is baked into QE's
   output); we're measuring pwdft-rs's internal sensitivity.

**Expected outcomes.**

- If disabling NLCC collapses the shared-density ΔE_xc from +0.43 eV
  to < 50 meV: the NLCC-on residual is confirmed as the mechanism.
  File a CNLC-F follow-up to deep-dive the C-specific mesh
  sensitivity (candidate: log-mesh Bessel quadrature at the
  {111}–{222} shell range). Do NOT promote `disable_nlcc` to a
  production knob.
- If disabling NLCC moves ΔE_xc by < 100 meV: the NLCC path is
  cleared. C Class A residual lives elsewhere — re-escalate to
  the V_NL `D_ij · β·β` l=1 hypothesis (VGCH-2E §
  Sub-hypothesis 2), which for C's s+p semilocal PP is the next
  candidate. Same outcome shape as PZPW's refutation branch:
  CNLC closes with "NLCC ruled out," moves the suspect along.
- If disabling NLCC changes ΔE_xc by ~0.43 eV **and** total energy
  changes similarly but with opposite sign on ΔE_1e: both NLCC
  and basis-density coupling are material; the real mechanism is
  the NLCC cross-coupling in the E_xc double-counting term. This
  is a subtler outcome and would need VGCH-2-ification.

##### Option B (rejected) — patched C.upf with zero `PP_NLCC`

Alternative: copy `pseudopotentials/nc/lda/C.upf` to
`pseudopotentials/nc/lda/C_noNLCC.upf`, set `core_correction="F"`,
remove the `PP_NLCC` block. Run the C LDA deck against the patched
UPF. Rejected because it mutates the pseudopotential library surface
and introduces a non-self-consistent PP (QE's own deck won't validate
against it). Option A's runtime switch is isolated and reversible.

### Verification

#### Phase 1

- `uv run pwdft-validate reference nlcc ...` regenerates CSV with
  new C rows.
- Two new Rust tests pass at 1e-5 tolerance.
- `/test` Tier-1 green.

#### Phase 3

- `disable_nlcc` flag round-trips through `ScfParams` (serde
  hidden, default false).
- The C diagnostic run emits an explicit logbook entry with:
  - NLCC-on ΔE_total vs QE at shared ρ.
  - NLCC-off ΔE_total vs QE at shared ρ.
  - Attribution verdict (NLCC guilty / cleared / partial).

#### Tier-2 PR policy

Touches `src/pseudopotential/upf/convert.rs` (Phase 1) and
`src/scf/potentials.rs` (Phase 3 flag). **Tier-2 required** per
CLAUDE.md. Expect no new Tier-2 regressions — the flag defaults to
false; production path is unchanged.

### Non-goals

- **Does NOT fix the C residual in this PR.** Diagnostic-only.
- **Does NOT touch the ρ_core(G→r) Bessel transform code** even if
  the diagnostic pins it as the culprit. That fix is CNLC-F scope.
- **Does NOT extend to C PBE** in the same PR; a parallel C-PBE
  transplant is a 15-min rerun (per VGCH-2E § Out of scope) and can
  land as a separate one-PR follow-up if the outcome shape needs
  confirmation.
- **Does NOT touch the NLCC double-counting subtraction in E_xc**
  (`scf::energy::compute_xc_energy_nlcc_correction` or equivalent —
  see Louie-Froyen-Cohen 1982 Eq. 4). That's a separate audit
  under VGCH-2.
- **Does NOT remove or modify existing ρ_core(G) pins.** All 8
  current pins must stay green.

### Risks

1. **Tolerance choice for C's 1e-5 pin.** C's Q_core ≈ 0.6 e —
   smaller than O's (~2 e) where 1e-5 worked. If C's magnitude at
   G=0 is on the order of ~0.05 e/Å³, a 1e-5 absolute tolerance is
   ≈ 2·10⁻⁴ relative, consistent with the Simpson-vs-trapezoidal
   residual seen on O/Cl. If the first pin fails at 1e-5 on an
   otherwise-valid transform, bump to 5e-5 and record the
   magnitude in the docstring. Document as pinned sensitivity in
   the logbook.
2. **The `disable_nlcc` flag stays in production.** Mitigated by
   `#[doc(hidden)]`, default false, YAML parser must not emit it
   (match PZPW's pattern of not adding a serde rename, or gate
   behind `#[cfg(any(test, feature = "diagnostic"))]` if that's
   too invasive). Remove after CNLC closes.
3. **Option A's zero-ρ_core path silently corrupts E_xc
   double-counting.** The E_xc correction term
   (Louie-Froyen-Cohen Eq. 4) uses ρ_core in the subtraction. If we
   zero ρ_core at one place but not another, ΔE_xc measurement
   becomes contaminated. Mitigation: zero ρ_core **at the source**
   in `compute_core_density` and verify all downstream consumers
   read from the grid (grep for `core_charge`, `rho_core`,
   `rho_core_g` in `src/scf/`). There should be exactly one
   producer and N consumers; audit once, land once.

### Provenance

- VGCH-2E: `proposals/VGCH-2E-c-diamond-class-c-transplant.md`
  (PR #178), Class A light-atom reclassification + NLCC hypothesis.
- VGCH-2F: `proposals/VGCH-2F-part-c-session-2-findings.md`
  (PR #179), 4 new NLCC pins (Ga/As/O/Cl) as the template.
- NCFX: `proposals/completed/NCFX-*.md` — original NLCC unit fix
  that closed the 13.4 eV Si gap; the sensitivity-class bug that
  CNLC Phase 3 is looking for would be a residual of that same class.
- Rust pin template: `pwdft/pwdft-core/src/pseudopotential/upf/convert.rs`
  lines 279 (`rho_core_of_g_ang`), 303 (Si test), 674 (O test),
  728 (Cl test).
- Python reference template:
  `pwdft/pwdft-validation/pwdft_validation/reference/nlcc.py:64–82`
  (`_SYSTEMS`).
- C.upf ground truth: `pseudopotentials/nc/lda/C.upf`
  (`core_correction="T"`, ONCVPSP LDA).
- C diamond cell ground truth: `data/qe/c_diamond_scf.in`
  (`ibrav=2, celldm(1)=6.7409` Bohr, a = 3.5672 Å).
- Paper citations:
  - Louie, Froyen, Cohen, *Phys. Rev. B* **26**, 1738 (1982),
    Eq. 3–4 (NLCC double-counting identity).
  - QE NLCC convention:
    `qe-7.5/upflib/rhoc_mod.f90:107-115` (Bessel transform),
    `qe-7.5/Modules/io_base.f90` (charge-density.dat schema used
    by the transplant side).

### Cost

- Phase 1: half CE-day (Python row + CSV regen + two Rust pins).
- Phase 2: half CE-day (optional multi-shell coverage).
- Phase 3: half CE-day (ablation flag + C diagnostic run +
  attribution logbook).

Total: **~1 CE-day** for Phases 1 + 3 (minimum viable CNLC).
Phase 2 optional.
