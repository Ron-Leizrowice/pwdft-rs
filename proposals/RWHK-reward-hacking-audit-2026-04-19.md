---
id: RWHK
status: draft
priority: high
complexity: small
risk: low
depends_on: []
blocks: [VQEF, VGCH-2, VGCH-MECH]
---

# RWHK: Reward-Hacking Audit (2026-04-19)

> Read-only adversarial audit of the validation machinery asking: "Is the
> 4 GREEN / 12 YELLOW / 0 RED scoreboard produced by honest physics
> agreement with QE, or by mechanisms that hide failures?"
>
> No `src/`, `tests/`, or `qe_validation/` changes in this PR. EM triage
> every Critical / Major finding below before merging downstream work.

## Scoreboard

- **Critical:** 1
- **Major:** 4
- **Minor:** 6
- **False-positive-ruled-out:** ~12 (listed inline; biggest ones noted)

## Executive summary

The scoreboard is **mostly honest but contains one reward-hack that
directly poisons the GPU-vs-CPU consistency signal** (Critical C1:
`test_gpu_vs_cpu_scf_direct_comparison` runs GPU on **both** branches —
the "CPU" branch still has the `gpu` feature enabled at compile time and
the GPU accelerator still initialises). The four GREEN E_total cells
(Si-LDA, Si-PBE, Al-LDA, Al-PBE) are all traceable to real QE
`pw.x` output files in `qe_validation/*.out`, with tolerances set at
observed + 20–60% margin and no self-referential pins that would mask a
physics regression. The 12 YELLOW cells correctly document
measured residuals in their `#[ignore]` reason strings and do not claim
QE agreement. XC constants (PZ, PW92, PBE β/κ/γ) cross-check against
QE 7.5 Fortran literals verbatim (non-issue). The `VGCH-SiEF-B1` gauge
fix was algebraically identity-preserving (0.023 meV drift), not
compensating-errors. Main scoreboard risk is **structural, not
fraudulent**: one GPU test mislabeled, one Fe-NLCC guard still pinned on
a pre-PCFX MP-1976 grid with 71× looser tolerance than its actual
post-MPSH measurement would permit, a handful of regression pins that
correctly disclaim QE agreement but still look like "checks" to a
drive-by reader, and a BSUM tolerance extrapolated without ever being
measured.

## Category A — Tolerance inflation

### A1 (False-positive-ruled-out) — LDA E_total GREEN tolerances are honest

Si-LDA 60 meV (observed 44.8 meV), Al-LDA 90 meV (observed 25.9 meV),
both documented in the test docstring with the observed-plus-margin
arithmetic shown. Same pattern for Si-PBE 20 meV (observed 12.4 meV)
and Al-PBE 20 meV (observed 8.1 meV). These are tight enough to catch
meaningful regressions (any E_total drift by 2×).

### A2 (Minor) — C-PBE BSUM tolerance is **extrapolated, never measured**

- **Location:** `tests/qe_validation.rs:1399-1403`
- **Evidence:** The comment reads "extrapolating from C-LDA's
  1.76/1.45 E_1e/E_total ratio suggests `|ΔE_1e| ≈ 400 meV` here, so
  600 meV leaves a small margin." But no one *measured* the C-PBE
  one-electron residual and printed it. The 400 meV is an inference.
- **Counterfactual:** Run the test once and record the observed
  `|ΔE_1e|`, then use observed + margin (like every other BSUM arm).
  The value is already printed by `run_qe_comparison`'s unconditional
  `E_1e^pwdft` log — just pull the actual residual from the eprintln.
- **Recommendation:** Add regression test — run the ignored test under
  `cargo test -- --ignored`, harvest the actual one-electron residual
  from the logs, update the tolerance to observed + ~30% margin.

### A3 (Minor) — VGC5 per-component pins are explicitly self-referential (acceptable, but visually identical to a QE-match claim)

- **Location:** `tests/vgc5_per_component_si.rs:280-317`, `:386-420`
- **Evidence:** The `pin()` helper takes an `expected` that was captured
  from pwdft-rs itself post-PCFX/post-B1, not from QE. This is correctly
  disclaimed in the module docstring ("Regression pins (pwdft-rs, NOT
  QE-match)"), so it is not a reward hack. But a drive-by reviewer
  skimming 10 `pin("E_local", c.e_local, -52.9772)` lines would see the
  structural shape of a ground-truth check.
- **Counterfactual:** Rename `pin` → `pin_regression`, and add an
  explicit `pin_vs_qe` helper used only for the one line
  (`E_ewald`) that IS lattice-only and QE-comparable.
- **Recommendation:** Non-blocking cosmetic improvement; file as
  Technical Writer drive-by when the pin values next change.

### A4 (Minor) — `test_fe_bcc_xc_nlcc_regression_guard` tolerance has a 71× headroom gap vs the Γ-centered observation

- **Location:** `tests/qe_validation.rs:1043-1110`
- **Evidence:** The test docstring explicitly says "switching to
  Γ-centered pushes `|Δ_xc|` to ≈ 1.41 eV" but the pinned ceiling is
  1.0 eV, so the test **must** stay on MP-1976 to pass. This is a
  regression guard pinned at the MP-1976 baseline value (0.69 eV) to
  catch the pre-NCFX 49 eV pathology — the 1 eV ceiling is the
  0.69 eV residual plus ~45% margin, which is fine for catching a
  50× regression. However, the test is now **locked** to MP-1976
  forever because its baseline was captured there; a future audit of
  grid-shift-vs-k-mesh-sensitivity in NLCC is silently blocked by
  this.
- **Counterfactual:** Track the Γ-centered observation as a sibling
  test (`test_fe_bcc_xc_nlcc_regression_guard_gamma`) with a
  tolerance of 2.0 eV (its 1.41 eV observed + ~40% margin), so
  NCFX-regression coverage exists on both k-grid conventions.
- **Recommendation:** Add as a follow-up proposal once VGCH Phase
  1a/b resolves the heavy-atom residual story (the MP-vs-Γ gap
  documentation is cleaner once the 0.7 eV vs 1.4 eV difference has a
  physical explanation).

### A5 (Minor) — GaAs-LDA tolerance 0.1 eV vs measured 35.2 eV residual

- **Location:** `tests/qe_validation.rs:856`
- **Evidence:** `assert_energy_matches_qe("GaAs", &result, -307.928_895_02, 0.1);`
  is masked behind `#[ignore]`. If the test is ever un-ignored without
  updating the tolerance, it fails (correct) — but the 0.1 eV claim
  looks like a live assertion.
- **Counterfactual:** Either delete the assertion body until the
  residual closes, or update the `tolerance` argument to match what the
  `#[ignore]` reason string records (35.2 eV). Same pattern on Cu-LDA
  (:904), NaCl-LDA (:956), MgO-LDA (:1009), Fe-LDA (:799) —
  the 0.1 eV / 0.05 eV claims in `assert_energy_matches_qe` do not
  match the `#[ignore]` residual disclosures.
- **Recommendation:** Pin each assertion at the observed residual + a
  sensible margin (e.g. Cu-LDA 20 eV, GaAs-LDA 40 eV). Prevents the
  "drop `#[ignore]` and see what happens" anti-pattern.

## Category B — Reference-value hardcoding / self-reference

### B1 (False-positive-ruled-out) — Every `total_energy_ry` in `reference_data.toml` matches its sibling `.out` file

- **Evidence:** Manual cross-check of every `! total energy =` line in
  `qe_validation/*.out` vs every `total_energy_ry` key in
  `reference_data.toml`. All 16 values agree to the last digit QE
  prints. Same for the PBE set. Category B on totals is clean.

### B2 (Minor) — `one_electron_ry` harvest is plausibly-correct but not independently verified

- **Evidence:** `reference_data.toml` header says "`one_electron_ry`
  values added 2026-04-19 under BSUM (harvested from the existing
  *.out files; no QE re-runs required)." No bulk-verification script
  exists in `scripts/validate/` to spot-check that the harvest
  actually pulls QE's `one-electron contribution` line (vs `eband` or
  `deband` or `hartree` mislabelled). The BSUM test does cross-check
  Si/Al on a passing assertion, so Si/Al are de-facto verified; the
  14 remaining values are one grep-and-typo away from being wrong.
- **Counterfactual:** Add a Python script (or inline `grep | awk`)
  in `scripts/validate/` that reparses every `one-electron
  contribution` line and compares it against `one_electron_ry` in
  TOML. Run once in CI as a data-integrity smoke test.
- **Recommendation:** File as Researcher follow-up — belt-and-braces
  on a dataset that is now baked into 16 tests.

### B3 (False-positive-ruled-out) — XC parameter literals (PZ, PW92, PBE) are canonical

- **Evidence:** `/src/potential/xc.rs:240-267` (PZ) and `:799-804`
  (PW92) and `:621-624,890-891` (PBE β/κ/γ) were cross-checked
  line-by-line against
  `qe-7.5/XClib/qe_funct_corr_lda_lsda.f90:43-47,355-356,897-898`
  and `qe_funct_corr_gga.f90:214-217`. All values are the published
  PZ/PW92/PBE paper constants, identical in pwdft-rs and QE. Not
  self-referential.

### B4 (Minor) — `test_si_total_energy_bit_identity_post_siefb1` pre-B1 pin `-231.6544` is self-captured and carries no physics reference

- **Location:** `tests/qe_validation.rs:519`
- **Evidence:** The test pins `E_pwdft = −231.6544 eV` as the
  "pre-B1 pin" with 1 meV drift tolerance. This value is
  self-referential — captured from pwdft-rs against itself, not from
  QE. Correctly disclaimed as an "algebraic-identity invariant"
  regression, not a physics agreement claim. Same risk as A3 — it
  looks like a physics check to a drive-by reviewer.
- **Counterfactual:** Rename to `test_si_total_energy_algebraic_identity_pin`
  and update the failure message to emphasize "algebraic identity, not
  QE match".
- **Recommendation:** Non-blocking; cosmetic improvement.

## Category C — Ignored tests masking bugs

### C1 (Critical) — `test_gpu_vs_cpu_scf_direct_comparison` does NOT compare GPU vs CPU

- **Location:** `tests/gpu_consistency.rs:486-612`
- **Evidence:** The test body comment at `:527-535` explicitly admits:
  > "Even with gpu feature, single-threaded rayon doesn't affect GPU
  > init. But the GPU will still be used in this path. To truly force
  > CPU-only, we'd need a runtime flag. For now, we verify
  > convergence consistency."
  Both branches call the same `pwdft_rs::scf::run_scf()` entry with
  the same inputs; the only difference is the surrounding rayon thread
  count. Under `#[cfg(feature = "gpu")]` the GPU path is unconditionally
  active in *both* branches. The test's `energy_diff < 0.1` assertion
  therefore verifies *GPU vs GPU* (trivially passes) rather than the
  f32/f64 boundary the test name implies.
- **Impact:** The "GPU consistency with CPU" signal is
  **structurally absent from the test suite**. Any f32 precision
  regression introduced by a WGSL kernel change (Hartree, XC, V_eff)
  will pass this test. The surviving f32-vs-f64 signal lives only in
  `test_gpu_scf_basic_convergence`'s ±0.1 eV pin against a
  self-captured `SI_REFERENCE_TOTAL_EV = -213.0283`, which is also
  GPU-vs-self.
- **Counterfactual:** Either (a) add a `ScfParams::force_cpu: bool`
  flag that the GPU driver respects and skip the GPU init in the CPU
  branch, or (b) move the CPU-side leg to a `#[cfg(not(feature =
  "gpu"))]` test so it literally cannot touch the GPU path. Option (b)
  is simpler but requires two `cargo test` invocations (with and
  without `--features gpu`) for full coverage — which is already the
  case per CLAUDE.md § Code Quality.
- **Recommendation:** **Fix immediately**. This is the single audit
  finding that the EM should treat as a blocker. Every GPU-feature PR
  claiming "consistency tests pass" has been silently
  passing via a GPU-vs-GPU tautology. Coordinate with Performance
  Engineer (owner of GPU path) to implement option (a) or (b).

### C2 (Minor) — `test_mxba_fe_documents_adaptive_failure` ignore has no deadline

- **Location:** `tests/mxba_adaptive_beta_fe.rs:56`
- **Evidence:** Ignore reason says "documents known adaptive-β failure
  mode on Fe BCC CCMX" but cites no tracking proposal ID for when the
  mode is fixed. MXB2 is in the FLUP backlog per INDEX.md but is not
  named in the ignore string.
- **Recommendation:** Update ignore reason to reference
  `proposals/FLUP-followup-backlog-seeding.md::MXB2` so a grep for
  MXB2 surfaces this test.

### C3 (Minor) — VGCH heavy-atom `#[ignore]` reasons are stale as of the PR body they reference

- **Location:** `tests/qe_validation.rs:829,882,929,982`
- **Evidence:** Five `#[ignore]` strings say "VGCH: heavy-atom
  residual (root cause TBD)" but VGCH Phase 1c (SAD cleared, PR #156)
  and VGCH-2 Part B (mixer cleared, PR #167) landed on 2026-04-19.
  The reason strings still say "root cause TBD" — which is *technically
  still true* (Class A, B, C classification from VGCH-MECH #168) but
  is less precise than it should be.
- **Counterfactual:** Update each ignore reason to reference its
  VGCH-MECH class (A / B / C). CLNP did this for other VGCH strings on
  2026-04-19 (PR #94); five arms were missed.
- **Recommendation:** Single-line ignore-string update follow-up —
  file under VGCH-MECH implementation work.

### C4 (False-positive-ruled-out) — Every other `#[ignore]` reviewed has a valid Tier-2 (TSPL) or VGCH-tracked reason

- **Evidence:** Of the ~45 `#[ignore]`s in `tests/`:
  - 28 are `TSPL Tier-2:` (valid — heavy SCF gated).
  - 10 are `VGCH Phase N` / `VGCH-2 class` (valid — tracked to
    proposal).
  - 3 are `VGCH-2 Part B` transplant diagnostics with explicit
    regeneration instructions.
  - 1 is MXBA `documents known adaptive-β failure` (C2 above).
  - 2 are diagnostic sweeps in `eigensolver/iterative.rs` with
    `"diagnostic sweep, not a gate"` (valid).
  None are bare `#[ignore]` with no reason.

## Category D — Feature-flag bypass

### D1 (False-positive-ruled-out) — `src/scf/transplant.rs::run_scf_iter1_from_rho_g_fft`

- **Evidence:** Reviewed line-by-line against `src/scf/driver.rs` iter
  0. Uses *identical* helpers: `add_core_density`,
  `hartree_on_fft_grid`, `assemble_v_eff`, `fill_hamiltonian_with_v_eff`,
  `diagonalize_dispatch`, `symmetrize_density_g`, `find_fermi_energy`,
  `compute_occupations`, `hartree_energy`, `xc_energy_corrected`,
  `kinetic_expectation`, `nonlocal_expectation`, `local_pp_energy_grid`.
  Identical symmetrization (G-space PCFX) and identical NLCC handling.
  Module explicitly marks itself `#[doc(hidden)]` and says "the
  diagnostic measures the same code path production SCF executes." Not
  a reward-hack path.

### D2 (False-positive-ruled-out) — No silent LDA fallback in `XcEvaluator::from_settings`

- **Evidence:** `src/potential/xc.rs:1401-1408`: exhaustive match on
  `XcFunctional::{Pz, Pbe, Pbe0, Hse06}` with `Pbe0`/`Hse06` returning
  `PwdftError::NotImplemented` at construction time (before any SCF
  work). The `scf::run_scf` dispatcher propagates this via `?`.

## Category E — Approximation that surface-matches

### E1 (False-positive-ruled-out) — PBE constants verbatim from QE

- **Evidence:** `PBE_KAPPA = 0.804`, `PBE_MU = 0.21951492776…`,
  `PBE_GAMMA = 0.0310906908696548950`, `PBE_BETA = 0.06672455060314922`
  cross-checked against `qe-7.5/XClib/qe_funct_exch_gga.f90` and
  `qe_funct_corr_gga.f90:214-217`. These are the canonical PBE paper
  values (Perdew, Burke, Ernzerhof, *Phys. Rev. Lett.* **77**, 3865
  (1996)), identical in every correct PBE implementation, not
  bit-tuned to QE.

### E2 (Major) — `krylov_max_dim` heuristic `max(128, max(2·n_request, n_pw/2))` was chosen to pass tuning tests

- **Location:** `src/eigensolver/iterative.rs:297`, tests `:558-581`
- **Evidence:** The heuristic was tuned by ITEV2 (PR #140) to close
  defect 1 on Si ecut=100 (|ΔE|=4.52e-12 eV). The ceiling comment
  cites "Lehoucq & Sorensen §3.2" but the specific floor (128) and
  ratios (2·n_request, n_pw/2) are empirical. Validated on
  `itev2_tune_n_request_across_ecuts` which sweeps Si/Fe/Cu at three
  ecuts. **No validation on GaAs, NaCl, MgO, C diamond, Al, or any
  system with pathological Hamiltonian condition number**. The
  default flip is explicitly gated behind an unstated additional
  check (CLAUDE.md says "Iterative is opt-in, not default") —
  which is the safety net.
- **Counterfactual:** (a) Broaden the coverage sweep to the full VQEF
  matrix before any default flip; (b) document the `max_dim_ceiling`
  pressure-release so operators can override when the heuristic
  underprovisions.
- **Recommendation:** Not-a-reward-hack for LDA, but **file as ITEV
  follow-up**: the empirical tuning was on a narrow system set and
  the default flip (ITEV Phase-5 step-4) should require a broader
  coverage pass first.

### E3 (False-positive-ruled-out) — `compute_density_gradient` Nyquist zeroing

- **Evidence:** Brief grep shows Nyquist zeroing is documented per
  Boyd §3.5 and matches QE's `fft_gga.f90` approach. Not checked
  exhaustively in this audit; flagging for Researcher review as a
  cross-referenced physics item (see Flagged for follow-up).

## Category F — Compensating-errors pattern

### F1 (False-positive-ruled-out) — VGCH-SiEF-B1 gauge fix is algebraically identity-preserving

- **Evidence:** `tests/qe_validation.rs::test_si_total_energy_bit_identity_post_siefb1`
  measured 0.023 meV drift (< 1 meV tolerance). The old code
  `total_energy + V_loc(G=0)·N_el` compensation was removed at the same
  time `e_band` gained the matching `V_loc(G=0)·N_el` piece via the
  Hamiltonian diagonal. `src/scf/energy.rs:585-600` documents the
  algebraic identity:
  `E_local_new = E_local_old + V_loc(G=0)·N_el;
  e_local_g0_shift_new = 0 (was V_loc(G=0)·N_el);
  e_band_new = e_band_old + V_loc(G=0)·N_el;
  E_total unchanged.`
  This is a gauge change, not a compensating error.

### F2 (False-positive-ruled-out) — TSEN `-TS` omission was not masked by another term

- **Evidence:** TSEN PR #162 landed Al GREEN 25.9 meV (observed delta
  for TS alone was ~110 meV on Al, matching QE within 1 meV). If a
  compensating term had been present, closing `-TS` would have
  **worsened** Al by 110 meV, not improved it by 82 meV. Simple
  before/after arithmetic confirms no compensating bug.

### F3 (Major) — NLCC XC double-counting: `xc_energy_corrected` uses `rho_val` for the vxc subtraction but `rho_xc = rho_val + rho_core` for the XC energy

- **Location:** `src/scf/energy.rs:136-154`
- **Evidence:** The implementation does:

  ```text
  e_xc  = ∫ (ρ_val + ρ_core) · ε_xc(ρ_val + ρ_core) dV
  e_vxc = ∫ ρ_val · V_xc(ρ_val + ρ_core) dV
  return e_xc - e_vxc
  ```

  Cross-checked against QE `PW/src/v_of_rho.f90::v_xc` — QE computes
  `etxc` on `rho_core_total = rho + rho_core` and `vtxc = Σ v_xc · rho`
  (valence-only) the same way. See
  `qe-7.5/PW/src/v_of_rho.f90:474-511`. Convention matches.
- **Concern:** Not a compensating-errors bug per se, but **this is the
  exact convention the Fe NLCC regression guard (Category A4 above)
  locks in via its 1.0 eV ceiling**. A sign flip in the convention
  (trivial one-liner) would still pass the Fe NLCC guard at 1.0 eV
  tolerance because the residual baseline was captured under the same
  (correct) convention. Defense-in-depth: add a *directional* regression
  test that swaps sign and asserts the test fails.
- **Recommendation:** File as Researcher follow-up to add a
  "wrong-sign NLCC" negative regression test in `tests/qe_validation.rs`.

## Category G — Stub / NotImplemented silently succeeding

### G1 (False-positive-ruled-out) — Every `NotImplemented` path is test-exercised

- **Evidence:** `grep -n 'PwdftError::NotImplemented' src/`:
  - `xc_functional 'pbe0'` / `'hse06'`: exercised by
    `src/potential/xc.rs:2007-2015` asserting error returned from
    `XcEvaluator::from_settings`. No PBE0 SCF test exists that could
    falsely pass.
  - `"pbe.eval requires rho_grad_r"`: internal guard that the driver
    path cannot trigger (the driver always supplies `rho_grad_r` when
    `needs_gradient() == true`). If driver is ever refactored wrong,
    the error fires loudly.
- **No `todo!()` or `unimplemented!()` in production paths** (0 hits on
  `grep 'todo!\|unimplemented!' src/`).

## Category H — Test coverage that doesn't exercise the physics it claims

### H1 (Major) — `test_si_pbe_non_spin_vs_qe` nominally exercises PBE end-to-end; does it really?

- **Location:** `tests/qe_validation.rs:1172-1202`
- **Evidence:** The test sets `xc_functional: XcFunctional::Pbe` and
  closes at 12.4 meV vs QE PBE reference. But the SCF-side confirmation
  that the evaluated XC was PBE (not LDA-silently-fallback) is
  indirect — the PBE-vs-LDA difference on Si is typically 100+ meV, so
  if LDA had silently been called the test would fail. Still, a
  positive assertion that `XcEvaluator::Pbe::eval` was *actually*
  invoked during the SCF would be stronger coverage.
- **Counterfactual:** Emit a debug log on `XcEvaluator::Pbe::eval`
  entry with the functional name, and assert in the test body that
  the log was observed. Alternatively add a `std::sync::AtomicUsize`
  counter inside `xc::Pbe::eval` and assert it bumped.
- **Recommendation:** File as TRV2 follow-up; soft assertion gap,
  physics-wise clean.

### H2 (False-positive-ruled-out) — Fe BCC PBE (`test_fe_bcc_fm_pbe_vs_qe`) exercises genuine nspin=2

- **Evidence:** Test sets `nspin: 2, starting_magnetization Fe=0.5`,
  then asserts `result.magnetization > 2.0 μB`. `magnetization` is
  computed from `∫|ρ↑ − ρ↓| dV` per
  `src/scf/driver_spin.rs::compute_magnetization` — a scalar that
  can only be nonzero if the SCF actually split the two spin densities.
  Genuine LSDA path. GGAP Phase D is legitimate coverage.

### H3 (Minor) — `assert_band_sum_matches_qe` uses k=0 of `result.eigenvalues` only (implicitly via `one-electron contribution`, not per-band)

- **Location:** `tests/qe_validation.rs:354-374`
- **Evidence:** BSUM does NOT do a per-band-per-k eigenvalue
  comparison; it compares the scalar composite `e_kinetic + e_local +
  e_local_g0_shift + e_nonlocal` against QE's `eband + deband` scalar.
  This is correct (and documented) — the docstring at `:321-353`
  explains why it is the shift-compensated quantity. But the phrase
  "band sum" in the function name misleads; a future reader
  could expect it to check `result.eigenvalues` per-k.
- **Counterfactual:** Rename `assert_band_sum_matches_qe` →
  `assert_one_electron_sum_matches_qe` (the function's actual
  semantics), and reserve `band_sum` for a future per-band pin if
  one ever lands.
- **Recommendation:** Cosmetic — file as Technical Writer follow-up.

## Category I — Silent fallbacks

### I1 (False-positive-ruled-out) — No default arm in XC dispatch

- **Evidence:** `src/potential/xc.rs:1460-1532` (eval) and
  `:1580-1670` (eval_spin) are exhaustive `match` on `XcEvaluator::{Pz,
  Pbe}`; no `_ =>` catch-all. `rustc` enforces exhaustiveness.

## Category J — Machine-precision asymmetry

### J1 (implied by C1) — GPU-path f32 tolerance not validated against CPU f64

- As noted in C1: the entire `test_gpu_vs_cpu_scf_direct_comparison`
  path is a GPU-vs-GPU comparison under `#[cfg(feature = "gpu")]`.
  See C1 for the full finding.

### J2 (False-positive-ruled-out) — No bit-identity assertions accidentally running through GPU f32

- **Evidence:** `test_si_total_energy_bit_identity_post_siefb1` tests
  1 meV drift (F4), far above f32 noise (~1e-3 eV on
  20-iter SCF); even if it were accidentally on the GPU path it would
  pass safely. Same for other `drift < 0.001` pins.

## Flagged for follow-up

- **GPU consistency test name vs implementation mismatch (C1)** —
  Performance Engineer fixes the GPU/CPU split with a `force_cpu`
  flag or splits the test across `#[cfg]` gates.
- **NLCC wrong-sign negative regression test (F3)** — Researcher
  verifies the convention and adds an inverted test that asserts the
  wrong-sign convention fails to match QE.
- **ITEV `krylov_max_dim` broader coverage (E2)** — Performance
  Engineer / Researcher sweep the heuristic across every VQEF cell
  before any `default: EigensolverKind::Iterative` flip.
- **Gradient Nyquist zeroing cross-check (E3)** — Researcher confirms
  `compute_density_gradient`'s Nyquist zeroing matches
  `qe-7.5/FFTXlib/fft_gga.f90`.
- **`one_electron_ry` harvest smoke test (B2)** — Researcher writes
  a Python script that reparses every `.out` and compares against TOML.
- **VGCH ignore-string refresh (C3)** — Technical Writer updates the
  5 stale "root cause TBD" strings to reference VGCH-MECH classes.
- **Regression-pin rename `pin` → `pin_regression` (A3, B4)** —
  Technical Writer cosmetic refactor the next time VGC5 pin values
  change.
- **GaAs/Cu/NaCl/MgO/Fe LDA bogus tolerances in `#[ignore]`d
  assertions (A5)** — Technical Writer updates each to match what the
  `#[ignore]` reason string documents.
- **C-PBE BSUM tolerance measured (A2)** — Code Reviewer or
  Researcher runs `cargo test -- --ignored test_c_diamond_pbe` to
  capture the real one-electron residual.
