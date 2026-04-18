# Researcher Logbook

Entries: date, what was validated, discrepancies found (with numbers), references used. Physics findings only — not code quality or docs.

## 2026-04-18 — CFGN re-scope (PR #78) — proposal only

Walked `src/` against `origin/main` @ `0be9290`. Census of user-facing numeric knobs shrinks from original ~25 to 10–12 after honoring NCFX (4 `1e-30` XC floors → `RHO_FLOOR` done), MODR (every file path drifted), CAST (invariants reinforce existing Settings shape — 0 new knobs), MXBA (AdaptiveBeta tunables = don't expose; Eyert paper defaults with test-encoded interactions).

**Top 3 highest-value exposures:** `electrons.fermi_search.{bounds_factor,tol,max_iter}` at `smearing.rs:73-92`; `scf.iterative_eigensolver.{tol,max_restarts}` at `eigensolver/iterative.rs:67,74`; `initial_density.gaussian_sigma` at `initial_density.rs:28` (already has runtime `InitialDensityConfig.gaussian_sigma`, just unwired from Settings).

**Two material errors in original proposal:** (a) "q-norm threshold 1e-12" claim — actual is `1e-9` at `nonlocal.rs:377`, branch on `|q|` (real-space Y_lm) not `cos(θ)` (VNLM rewrote this); (b) phantom `src/potential/hartree.rs` and phantom GPU "real buffer pool = 3" — neither exists.

Phase 1 = smearing/eigensolver/initial_density; Phase 2 = Ewald/RHO_FLOOR/G2_ZERO_THRESHOLD. Independent PRs.

## 2026-04-18 — NLCC audit (PR #42) — tests+docs only

No production code change (NCFX #40 was the fix). Added Python/SciPy ρ_core(G) pin at Si/Fe reference points (< 1e-7 residual), regression test on Fe E_xc (post-NCFX 0.69 eV vs QE — pre-NCFX was 48.85 eV, trips by >49× if regressed), and LFC (Louie/Froyen/Cohen PRB 26 1738 (1982)) citations wherever NLCC is mentioned.

**Fe regression test needed looser conv_thr on GPU** — f32 precision floor stalls at Δρ≈9e-8. Platform-specific tolerance for GPU XC tests is now a documented pattern.

## 2026-04-17 — PCRS: Si per-component 1.204 eV residual root-caused (→ PCFX)

**Plateau is structural, not SCF noise:**

| conv_thr | sym ON | sym OFF |
|----------|--------|---------|
| 1e-6     | 1.204  | 1.07e-07 |
| 1e-8     | 1.204  | 1.03e-08 |

sym OFF scales O(Δρ) (expected); sym ON is a hard plateau. QE's identity closes bit-exact.

**Root cause:** `symmetrize_density` applies `{R|τ}` via `nint(n·τ_i)` rounding. Si Fd-3m has τ=(¼,¼,¼); FFT grid 18 is not divisible by 4 → 4.5 rounds to 5, inverse op rounds back to a different grid point → density smeared. `check_grid_compatibility` only checks rotations, not translations. Fe Im-3m is symmorphic (τ=0) → bug absent.

E_KS itself is also biased ~17 meV (ρ_sym ≠ ρ_ψ leaks into `-e_H + (e_xc − e_vxc)` double-counting). Not just a diagnostic issue. Fix: G-space symmetrization via phase factors (QE `symme.f90::sym_rho`). Filed as PCFX.

**Tangential (preserved from session):**
- `compatible_grid_dims` is dead code outside tests.
- FFT grid 20/24 (both div-by-4) → Si Anderson mixer loses conditioning and oscillates. Separate mixer concern.

## 2026-04-17 — VGC5: per-component energy accounting — XC is the Si culprit

Added `ScfResult.components` (`EnergyComponents`) populated at convergence. Found:

| term        | pwdft-rs    | QE          | Δ (ours−QE) |
|-------------|-------------|-------------|-------------|
| one-electron| +68.608     | +66.225     | +2.38       |
| E_hartree   | +13.593     | +15.104     | −1.51       |
| **E_xc**    | **−70.658** | **−84.396** | **+13.74**  |
| E_ewald     | −228.519    | −228.530    | +0.011      |
| **E_total** | **−218.181**| **−231.610**| **+13.43**  |

Fe nspin=1: ΔE_xc = −48.85 eV. Scales with NLCC magnitude → NLCC core-density FT is the bug.

**Two compounding bugs in NLCC:**
1. `upf.rs:95-105` divides PP_NLCC by `BOHR_TO_ANG` — PP_NLCC is bare ρ_core(r) in e/Bohr³, correct divisor is BOHR_TO_ANG³.
2. `potentials.rs:92-107` integrates `ρ_c · j₀(Gr)` without r² weight and without 4π. QE's `rhoc_mod.f90:107` uses `ρ·r²·j₀` then `fpi/Ω`.

Filed as NCFX.

**V_loc(G=0) shift rules itself out as prime suspect.** `src/scf/mod.rs:415,424` adds `ctx.v_local_g0 * ctx.n_electrons` to both E_KS and E_HF — formula was right, docstring was incomplete (as the 2026-04-16 orientation note said).

## 2026-04-17 — VGCMP Phases 1-4: entire PP→H pipeline clean

All four cross-checks (V_local(G), β_l(q), D_ij, H_diag[G,G] at k=Γ) agree with independent Python (scipy simpson + spherical_jn) to 1e-9 to 1e-14 — 4-8 orders below tolerance.

| Phase | Quantity                      | max \|Δ\|                      |
|-------|-------------------------------|--------------------------------|
| 1     | V_local(G)                    | 2.8e−9 Ry                      |
| 2     | β_l(q)                        | 3.0e−12 Bohr^(3/2)             |
| 3     | D_ij                          | 0.0 Ry (bit-exact)             |
| 4     | H_diag[G,G] (kinetic + V_NL)  | 2.9e−10 Ry                     |

**Key finding:** Si.upf D_ij is **strictly diagonal** (not merely block-diagonal in l) — QE absorbs the within-l-block rotation into χ(r). Diagonal values (Ry): +11.132, +1.714, +5.452, +1.260, −4.250, −0.889 for (l=0,0,1,1,2,2).

Phase 4 kinetic residual 7e-11 relative is CODATA-vs-SI drift between `HBAR2_OVER_2M` (3.8099821159) and `RY_TO_EV·BOHR_TO_ANG²` (3.8099821161). Not a bug, but worth knowing — pin one convention if Phase 4 ever tightens.

The Si 13.4 eV gap is outside the entire PP→H assembly pipeline. Phase 5 pivoted to energy-component audit → VGC5 → NCFX (above).

## 2026-04-17 — SYKP: Si 4×4×4 IBZ 10 vs QE 8 is convention mismatch

pwdft-rs' `monkhorst_pack` hard-codes shifted MP-1976 (`frac = (2i−N+1)/(2N)`); QE's `si_scf.in` uses `4 4 4 0 0 0` which is Γ-centered unshifted. Both correct for their convention — not a bug. k-sampling convergence error at 4×4×4 is <10 meV for either grid; cannot produce the 13,400 meV Si gap.

**Dead wrong comment to fix:** `src/symmetry/kpoints.rs:148-151` claims "10 due to incomplete boundary handling" — replace with convention note. Logged as SYKP D1 (done via XCLN).

**Follow-ups:** MPSH proposal (add MP-shift parameter for byte-identical QE comparison) if anyone wants it.

## 2026-04-17 — QEVL + SPXC seed investigations

**QEVL:** QE reference data generated for 8 Tier 1+2 systems; tests in `qe_validation.rs` with `#[ignore]` reasons. All energies in Ry pinned: Si -17.02299344, C -23.84343910, Al -4.72371790, Fe -224.91744934, GaAs -307.92889502, Cu -356.73602869, NaCl -119.77970303, MgO -147.23547768. Tier 3 convergence studies deferred.

**SPXC:** confirmed bug in `run_scf_spin` E_KS: `exc_r` from INPUT mixed with `rho_xc_total` from OUTPUT. Non-spin path recomputes XC from OUTPUT (correct). Impact: O(Δρ) at convergence — sub-meV on total E but spoils quadratic |E_HF−E_KS| → linear, triggering false HF-KS warnings. QE never recomputes XC from output; uses `etxc`/`vtxc` from input + `descf` first-order correction (same result at convergence).

**Tangential from QEVL:** Fe nspin=2 at ecut=15 Ry collapses to NM under PseudoDojo; `Fe_dalcorso.upf` at higher cutoff may be needed to exercise magnetism. C diamond non-convergence deserves own proposal once Si offset clears.

## 2026-04-16 — Handoff + KBTF

**QE discrepancy baseline:**

| System | QE (eV) | Ours (eV) | ΔE (eV) |
|--------|---------|-----------|---------|
| Si (2 atoms, 15 Ry, 4×4×4) | -231.61 | -218.28 | 13.3 |
| Fe BCC (1 atom, 16 Ry, 4×4×4) | -3059.46 | -3104.83 | 45.4 |
| C diamond | — | — | Non-convergent |

Γ degeneracy broken ⇒ bug in V_local/V_NL form factors, not energy accounting. Error scales with Z. Formula audit: all 23 items verified correct — no formula bugs. Proposals SIMP + VERF spawned.

**KBTF classification:** `test_09` was a test bug (UPF D_ij is diagonal 6×6 post-QE diagonalisation of HGH h^l 3×3 blocks — published HGH values like h^0_11=2.95 Ry cannot be compared against UPF D[0,0]=11.13 Ry directly). test_07 trapezoidal artifact (SIMP); test_vloc bare-Coulomb (VERF).

**Low-priority open questions:** IBZ 10 vs QE 8 (→ SYKP); spin exchange at `xc.rs:249` non-standard weighted-average (correct, needs doc); `total_energy()` docstring omits V_local(G=0)·N_el (→ MADOC).

Refs: PZ PRB 23 5048; KB PRL 48 1425; NLCC PRB 26 1738; QE `vloc_mod.f90`, `simpsn.f90`, `setlocal.f90`.
