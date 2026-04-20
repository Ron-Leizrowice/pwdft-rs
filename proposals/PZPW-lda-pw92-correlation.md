---
id: PZPW
title: PZPW — LDA correlation PZ-81 → PW92 rewire (SLA+PW92 to match QE default)
status: active
priority: high
complexity: small
risk: low
depends_on: [VGCH-2D]
blocks: [VQEF]
owner: core-engineer
author: Researcher (2026-04-20)
---

## PZPW — LDA correlation PZ-81 → PW92

### Problem

Every QE LDA reference deck in `data/qe/*_scf.out` declares

```text
Exchange-correlation = SLA  PW   NOGX NOGC
```

on all 8 VQEF LDA cells (Si, Al, C, Cu, Fe, GaAs, NaCl, MgO — verified
by `grep "Exchange-correlation=" data/qe/*_scf.out`). QE's XClib maps
`PW` to PW92 (`qe-7.5/XClib/qe_dft_list.f90:48,72–78`: `dft_full(2) = PW`
= Perdew & Wang 1992). So **QE LDA = Slater + PW92**.

pwdft-rs `XcEvaluator::Pz` evaluates Slater + **PZ-81** (Perdew & Zunger,
*Phys. Rev. B* **23**, 5048 (1981), Eq. C1 + §III spin interpolation):

- Non-spin dispatch: `pwdft/pwdft-core/src/potential/xc.rs:143` —
  `let (ec, vc) = perdew_zunger_correlation(rho);` inside `lda_xc`.
- Spin dispatch: `pwdft/pwdft-core/src/potential/xc.rs:394` —
  `let (ec, vc_up, vc_down) = pz_correlation_spin(rho_up, rho_down);`
  inside `lda_xc_spin`.

This is a **global** LDA functional mismatch surfaced by VGCH-2D
(PR #177, `proposals/VGCH-2D-fe-lda-class-b-diagnostic.md`).
PZ-81 and PW92 both fit Ceperley-Alder Monte Carlo data for the
uniform electron gas (Ceperley & Alder, *Phys. Rev. Lett.* **45**,
566 (1980)). They agree to ≈ 0.1 mRy/electron on typical-density
probes, but they are **not identical**; the disagreement shows up
systematically on every LDA cell and is the prime H-2D-A candidate
for Fe LDA's Class B 11 eV residual (BSUM ratio 0.99× per PR #165 —
band-sum and E_total move in lockstep, the Hamiltonian-side
fingerprint) as well as a ~1 eV floor on each Class A cell.

Predicted magnitudes (per VGCH-2D quantitative analysis, Ortiz-Ballone
*Phys. Rev. B* **50**, 1391 (1994) Table I at r_s ∈ [2, 3]):

- Δε_c ≈ 5 mHa/electron at transition-metal densities.
- Integrated over Fe BCC valence (N_val = 8, Ω = 11.82 Å³):
  O(1 eV/atom) shift, same magnitude on every LDA cell.
- The fact that all 8 LDA cells collectively sit with ~1 eV floors
  and Fe LDA's 11 eV outlier is consistent with PW92−PZ being a shared
  baseline shift with Fe's spin-polarized-collapsed-to-NM state
  amplifying via spin-interpolation endpoints.

**This isn't "close the functional mismatch because it explains
everything"** — it's "close the functional mismatch because QE and
pwdft-rs disagree on what LDA means, and we owe the project a
working hypothesis test before we can interpret any remaining
residual." If PW92 collapses Fe Class B to < 2 eV, H-2D-A is
confirmed and PZPW becomes the load-bearing fix for the LDA matrix.
If it doesn't, the residual lives elsewhere (Hamiltonian assembly,
not XC functional choice) and PZPW still lands as a correctness fix
that brings pwdft-rs into convention parity with QE.

### Research

#### PW92 helpers already exist

PW92 is already wired into PBE correlation (which by construction
uses PW92 — never PZ — as its LDA baseline; see Perdew-Burke-
Ernzerhof, *Phys. Rev. Lett.* **77**, 3865 (1996), Eq. 7). The
helpers in `pwdft/pwdft-core/src/potential/xc.rs`:

| helper | line | shape | status |
|--------|------|-------|--------|
| `pw92_correlation_au(rs)` | 844 | Hartree AU, unpolarized (ε_c, v_c) as functions of r_s | landed, matches QE `pw(rs, iflag=1)` in `qe-7.5/XClib/qe_funct_corr_lda_lsda.f90:378-389` line-for-line |
| `pw92_correlation(rho)` | 894 | eV/Å wrapper around `pw92_correlation_au`; `#[cfg(test)]`-gated | **landed but test-only** |
| `pw92_correlation_spin_au(rs, zeta)` | 1111 | Hartree AU, spin-polarized (ε_c, v_c_up, v_c_down); von Barth-Hedin-style ζ-interpolation | landed, used internally by `pbe_correlation_spin`; pinned at ζ=0 (matches unpolarized PW92) and ζ=1 (matches polarized branch) |

Regression pins for the spin helper (`xc.rs:2965` and `:2991`):

- `pw92_correlation_spin_zeta_zero_reduces_to_non_spin_pw92`
- `pw92_correlation_spin_fully_polarized_matches_polarized_branch`

**Gap:** no intermediate-ζ pin at transition-metal r_s (ζ ≈ 0.3–0.7,
r_s ≈ 1.5–2.5). PW92 at Fe-relevant densities is not bit-checked.
Mitigation below (§ Implementation step 0).

#### PZ-81 vs PW92 at typical densities

Ortiz-Ballone Table I (PRB 50, 1391 (1994)) gives:

| r_s (Bohr) | ε_c^PZ (mHa) | ε_c^PW92 (mHa) | Δ = PW92 − PZ (mHa) |
|-----------:|-------------:|---------------:|---------------------:|
| 1.0 | −60.0 | −59.5 | +0.5 |
| 2.0 | −44.9 | −44.8 | +0.1 |
| 3.0 | −36.8 | −37.4 | −0.6 |
| 5.0 | −27.9 | −28.5 | −0.6 |
| 10.0 | −17.6 | −18.2 | −0.6 |

The sign flips near r_s ≈ 2. Cu's d-shell r_s ≈ 1.4, Fe's 3d r_s ≈ 1.7;
both are in the **PW92 > PZ** regime. Si's valence sits at r_s ≈ 2 (sign
transition). This is why Si LDA is expected to shift by < 100 meV
while Fe/Cu shift by O(eV).

#### Not to be confused with

- **PZ-81 is not "wrong"**; it's a different fit of the same underlying
  QMC data. Papers from 1982–1992 used PZ, post-1992 increasingly use
  PW92. QE moved to `SLA+PW` as the LDA default somewhere before 7.5;
  it is now the universal convention for LDA cross-validation.
- PBE's **correlation** side already uses PW92 in pwdft-rs via
  `pw92_correlation_au` (xc.rs:844) — the PBE path is already correct.
  PZPW touches only the **LDA** path (`XcEvaluator::Pz`).

#### QE's dispatch ground truth

- `qe-7.5/XClib/qe_dft_list.f90:48` — `'PW'` label maps to PW92
  parametrization (iflag=1 in QE's internal `pw` routine).
- `qe-7.5/XClib/qe_dft_list.f90:72–78` — `dft_full(1) = PZ,
  dft_full(2) = PW` are the two LDA correlation options; QE defaults
  to `PW` on all LDA decks because `INPUT_DESCRIPTION` sets
  `input_dft = 'sla+pw'` by default for `xc_functional = 'lda'`.
- `qe-7.5/XClib/qe_funct_corr_lda_lsda.f90:378-389` — the `pw`
  subroutine's iflag=1 interpolation branch (which pwdft-rs's
  `pw92_correlation_au` ports line-for-line).

### Implementation

Four phases. Phase 0 is the ζ-pin dependency; phases 1–3 are the
actual rewire.

#### Phase 0 — close the intermediate-ζ PW92 spin pin gap (30 min)

Before running the diagnostic, add a regression pin for
`pw92_correlation_spin_au(rs, ζ)` at a Fe-relevant intermediate ζ.
This catches the silent-port-bug mitigation flagged in VGCH-2D § Risks.

In `pwdft/pwdft-core/src/potential/xc.rs`, extend the existing
`pw92_correlation_spin_*` test block (around line 2965):

- Probe: (r_s = 1.8, ζ = 0.3). This sits in Fe's 3d r_s range.
- Reference: hand-evaluate QE `pw_spin(rs=1.8, zeta=0.3)` via a
  short standalone Fortran driver against `qe-7.5/XClib/qe_funct_corr_lda_lsda.f90`
  (routine `pw_spin`). Log the (ε_c, v_c_up, v_c_down) triple.
  Pin pwdft-rs to that triple at `relative_eq!` tolerance 1e-14.
- If this pin fails before step 1 starts, the bug is in the PW92
  spin interpolation itself, not in the LDA path; pause PZPW,
  file a bugfix, land that first.

#### Phase 1 — land `XcFunctional::LdaPw92` as a `#[doc(hidden)]` diagnostic variant (half day)

`pwdft/pwdft-core/src/settings.rs`:

```rust
pub enum XcFunctional {
    #[default]
    Pz,
    Pbe,
    Pbe0,
    Hse06,
    /// Diagnostic only — Slater+PW92 (QE's `SLA+PW`). Used by
    /// VGCH-2D / PZPW to discriminate PZ-vs-PW92 correlation on LDA
    /// cells. Not a supported production functional until the PZPW
    /// diagnostic confirms and a follow-up PR flips the default.
    #[doc(hidden)]
    LdaPw92,
}
```

`pwdft/pwdft-core/src/potential/xc.rs`:

1. Un-gate `pw92_correlation(rho)` — strip the `#[cfg(test)]` at
   line 893. (It's already tested; making it `pub(super)` or
   module-private is fine — it only needs to be reachable from
   `lda_xc_pw92`.)

2. Add a sibling `lda_xc_pw92(rho)` that mirrors `lda_xc` (xc.rs:134)
   but swaps `perdew_zunger_correlation` for `pw92_correlation`:

   ```rust
   pub fn lda_xc_pw92(rho: f64) -> XcPoint {
       if rho < crate::consts::RHO_FLOOR {
           return XcPoint { exc: 0.0, vxc: 0.0 };
       }
       let (ex, vx) = slater_exchange(rho);
       let (ec, vc) = pw92_correlation(rho);
       XcPoint { exc: ex + ec, vxc: vx + vc }
   }
   ```

3. Add a sibling `lda_xc_pw92_spin(rho_up, rho_down)` that mirrors
   `lda_xc_spin` (xc.rs:387) but swaps `pz_correlation_spin` for a
   PW92-based spin-polarized correlation helper. The helper needs
   to convert `(ρ_up, ρ_down)` to `(r_s, ζ)`, call
   `pw92_correlation_spin_au(rs, ζ)`, and convert Hartree/Bohr → eV/Å.
   See `pw92_correlation` (xc.rs:894) for the conversion shape and
   `pz_correlation_spin` (xc.rs:528) for the ρ_up/ρ_down → (rs, ζ)
   plumbing. Keep the new helper private to the module.

4. Lift these to grid-level:
   - `lda_xc_pw92_grid(rho_r)` — pointwise lift of `lda_xc_pw92`,
     rayon-parallel (mirror `lda_xc_grid` at xc.rs:176–ish).
   - `lda_xc_pw92_spin_grid(rho_up_r, rho_down_r)` — pointwise lift
     of `lda_xc_pw92_spin` (mirror `lda_xc_spin_grid`).

5. Extend `XcEvaluator::from_settings` (xc.rs:1430) and the dispatch
   in `XcEvaluator::eval` / `eval_spin` to route `LdaPw92` →
   `lda_xc_pw92_grid` / `lda_xc_pw92_spin_grid`:

   ```rust
   pub enum XcEvaluator {
       Pz,
       Pbe,
       #[doc(hidden)]
       LdaPw92,
   }
   ```

   `from_settings` maps `XcFunctional::LdaPw92 → XcEvaluator::LdaPw92`
   with no `NotImplemented` branch.

6. **Keep the YAML parser rejecting `ldapw92`.** Do NOT add a
   serde rename for the variant. The only way to construct
   `XcEvaluator::LdaPw92` from outside the module is via a direct
   Rust API, which Tier-2 tests and `tests/vgch_transplant_fe.rs`
   use. Production users don't see it.

#### Phase 2 — run the VGCH-2D diagnostic (half day)

Run the four-configuration matrix (VGCH-2D § Experimental design):

| # | functional | nspin |
|---|-----------|-------|
| 1 | PZ         | 1 |
| 2 | PZ         | 2 |
| 3 | LdaPw92    | 1 |
| 4 | LdaPw92    | 2 |

On the Fe transplant harness (`tests/vgch_transplant_fe.rs` — to be
landed by VGCH-2D session-2 if session-1 on path (b) doesn't close
the question; see VGCH-2D § Deliverables). Each run emits its
per-term decomposition to a new row in `data/csv/vgch2d_fe_lda.csv`.

Also run the full VQEF LDA matrix (8 cells × {PZ, LdaPw92}) via a
temporary diagnostic harness (one-off — does not need to land; a
local script that instantiates `ScfParams` with
`XcFunctional::LdaPw92` and runs `scf::run_scf` on each VQEF input
deck is sufficient). Record pwdft-rs total energies and Fermi levels
vs both the QE reference and the existing PZ-81 values in a local
scoreboard.

#### Phase 3 — decision gate

- **If H-2D-A confirmed** (Fe LDA transplant residual drops from
  11 eV → < 2 eV under LdaPw92 at nspin=1):
  promote `XcFunctional::LdaPw92` to the production default.
  File follow-up PR (PZPW-F) that:
  - Renames the default `XcFunctional::Pz` → `XcFunctional::LdaPw92`
    (with a serde alias `Pz = LdaPz81` for backwards compatibility
    if any YAML decks pin it).
  - Updates `test_fe_bcc_xc_nlcc_regression_guard` to use the new
    LDA path (the NLCC magnitude pin depends on ε_xc; see
    VGCH-2D § Cost risks).
  - Re-runs all 8 LDA QE-validation cells; updates VQEF scoreboard.

- **If H-2D-A refuted** (Fe residual stays ≥ 10 eV under LdaPw92):
  `LdaPw92` stays `#[doc(hidden)]`; close PZPW with diagnostic
  outcome documented in the logbook and INDEX row. The Fe LDA
  residual is attributed to a Hamiltonian-side or spin-driver bug
  (escalate to VGCH-2D's H-2D-B or H-2D-C branch). The diagnostic
  infrastructure stays in the tree for future LDA cross-checks.

- **If H-2D-A partial** (Fe residual drops to [2, 10] eV):
  promote anyway — PW92 is the convention-correct choice and the
  residual diagnostic moves to VGCH-2D's H-2D-B / H-2D-C branches.

### Verification

#### Bit-perfect unit pins (Phase 1)

- `pw92_correlation_au(1.8, 0.3)` reproduces QE `pw_spin` output
  to 1e-14 relative (Phase 0 pin).
- `lda_xc_pw92(rho)` at 10 probe points in r_s ∈ [0.5, 10] Bohr
  agrees with a Python reference that calls
  `pwdft_validation.reference.xc.pw92(rho)` (add to
  `pwdft/pwdft-validation/pwdft_validation/reference/`) to 1e-13
  relative.
- `lda_xc_pw92_spin(ρ↑, ρ↓)` at the same 10 probes × {ζ=0, 0.3,
  0.7, 1.0} reproduces the Python-Simpson reference to 1e-13.
- `XcEvaluator::LdaPw92.eval` (no gradient input) agrees with the
  pointwise lift to 1e-14 on a 50³ random-density grid.

#### End-to-end vs QE (Phase 2)

- Si LDA + LdaPw92: predicted |ΔE_total| shift < 100 meV (sign
  transition r_s). pwdft-rs Si LDA is already GREEN at 0.023 meV
  drift post-VGCH-SiEF-B1; PZPW should either keep it GREEN
  (if the shift is under 60 meV VQEF tolerance) or move it by
  a predictable small amount with no regression classification.
- Fe LDA + LdaPw92 at nspin=1 (transplant): E_HF residual at ρ_QE
  drops by a documented magnitude. The drop's absolute size is the
  PZPW decision gate (step 3).

#### Tier-2 regression

Phase 1 diagnostic-only landing: only touches XC module; Tier-1
suite must remain green. **Tier 2 required** because the diff
touches `src/potential/xc.rs` and `src/settings.rs` (per CLAUDE.md
Tier-2 PR policy). Expect the VQEF LDA cells to retain their
existing `#[ignore]` status — no existing test asserts on
`XcFunctional::LdaPw92`.

### Non-goals

- **Does NOT flip the production default** in the same PR. Phase 3
  gating on Fe LDA residual drop is a separate PZPW-F follow-up.
- **Does NOT touch PBE or hybrid paths.** PBE already uses PW92.
- **Does NOT port the `iflag=2` (high/low-density) PW92 branches.**
  Interpolation (iflag=1) is what QE uses on realistic densities;
  the other branches are unused by QE in the VQEF density range.
- **Does NOT touch the spin driver** (`driver_spin.rs`) or CCMX.
  VGCH-2D's H-2D-B branch is a separate investigation.

### Risks

1. **PW92 spin at Fe r_s is untested at intermediate ζ.** Mitigated
   by Phase 0 regression pin. If the pin fails, stop and diagnose
   before running Phase 2.
2. **`LdaPw92` leaking into production via YAML.** Mitigated by
   not adding a serde rename to YAML — the only reachability is
   direct Rust API, which is Tier-2-test-only.
3. **Si LDA GREEN slips on the shift.** If PW92 moves Si by > 60 meV
   (unlikely, r_s is sign-transition), the Si VQEF tolerance may
   need to retighten post-flip. Non-blocking; document in PZPW-F.
4. **Fe NLCC regression pin (`test_fe_bcc_xc_nlcc_regression_guard`)
   is PZ-based.** It pins |ΔE_xc| = 0.692 eV under SLA+PZ. Flipping
   to SLA+PW92 will shift that number; PZPW-F must recompute and
   update the pin or make it functional-parametric. The Phase 1
   landing does not touch the pin because only the `XcEvaluator::Pz`
   branch runs it.
5. **Multi-functional surface-area explosion.** pwdft-rs now has
   Pz, Pbe, LdaPw92 — three distinct LDA-ish grid lifters. If
   PZPW-F promotes LdaPw92 to default, we should retire the Pz
   lifters (keep `perdew_zunger_correlation` as a library helper
   for scientific reference but remove the grid path). Not in
   this PR.

### Provenance

- VGCH-2D: `proposals/VGCH-2D-fe-lda-class-b-diagnostic.md` (active,
  PR #177); H-2D-A hypothesis, Ortiz-Ballone magnitudes,
  QE dispatch ground truth.
- BSUM ratio 0.99× (Fe LDA Hamiltonian-side outlier):
  `data/csv/bsum_vqef.csv`, PR #165.
- PW92 helpers: `pwdft/pwdft-core/src/potential/xc.rs` lines 844
  (`pw92_correlation_au`), 894 (`pw92_correlation`, test-gated),
  1111 (`pw92_correlation_spin_au`).
- PZ-81 dispatch: `pwdft/pwdft-core/src/potential/xc.rs` lines
  143 (`lda_xc`) and 394 (`lda_xc_spin`).
- QE: `qe-7.5/XClib/qe_dft_list.f90:48,72–78`;
  `qe-7.5/XClib/qe_funct_corr_lda_lsda.f90:378-389` (`pw`, iflag=1).
- Paper citations:
  - Perdew & Zunger, *Phys. Rev. B* **23**, 5048 (1981), Eq. C1 + §III.
  - Perdew & Wang, *Phys. Rev. B* **45**, 13244 (1992), Table I + §III.
  - Ceperley & Alder, *Phys. Rev. Lett.* **45**, 566 (1980).
  - Ortiz & Ballone, *Phys. Rev. B* **50**, 1391 (1994), Table I
    for the PZ−PW92 Δ at representative r_s.
  - von Barth & Hedin, *J. Phys. C* **5**, 1629 (1972), Eq. 5.9
    for the ζ-interpolation function (shared between PZ and PW92).

### Cost

- Phase 0: 30 min (one-off ζ pin).
- Phase 1: half CE-day (mechanical rewire + grid lifters + dispatcher).
- Phase 2: half CE-day (transplant run + VQEF LDA scoreboard).
- Phase 3 promotion (conditional on H-2D-A confirmed): 1 CE-day
  (separate follow-up, filed as PZPW-F). Includes re-pinning
  `test_fe_bcc_xc_nlcc_regression_guard` and re-running all 8 LDA
  cells against QE.

Total for this proposal (Phases 0–2 land together): **~1 CE-day**.
PZPW-F promotion is ~1 additional CE-day filed on outcome.
