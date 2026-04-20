# Researcher Logbook

Entries: date, what was validated, discrepancies found (with numbers), references used. Physics findings only — not code quality or docs.

**Load-bearing cross-language convention (→ candidate CLAUDE.md promotion):** Rust's `f64::round()` is half-away-from-zero; numpy's `np.round` is half-to-even (banker's rounding). At exactly `d = ±0.5` they disagree in sign of the integer. Any Rust↔Python validation that maps real-valued distances/fractions to grid bins must either (a) avoid `.round()` entirely or (b) use an explicit wrap like `d - (d + 0.5).floor()` on both sides. VGCH Phase 1c (2026-04-19) lost half a day to this in a shell-average diagnostic — C diamond showed a bogus 0.2 e/Å³ asymmetry that evaporated on raw-sample diff. No production-code impact yet; diagnostic-only. Worth a §Conventions bullet if it bites again.

## 2026-04-20 — VGCH-2 Part C session-1: H-C1 / H-C3 CLEARED (PR #174)

Python Fermi-Dirac bisection over QE converged eigenvalues + weights
reproduces QE's reported E_F to:

| system    | |ΔE_F_py − E_F_QE| |
| ----------- | --------------------- |  |  |
| Cu FCC | 6 μeV |  |  |
| Fe BCC FM | 15 μeV (nspin=2) |  |  |
| C diamond | 10 μeV |  |  |
| MgO | 66 μeV |  |  |
| NaCl | 180 μeV |  |  |

pwdft-rs's `smearing::find_fermi_energy` matches QE's `efermig.f90` on:
(a) sign convention (`x = (ε−E_F)/σ; f = 1/(1+exp(x))` ≡ QE
`wgauss((E_F−ε)/σ, −99)`), (b) bracket initialization (`±10σ`), (c)
the bisection root agrees on `1e-14 eV` bracket-width tol vs QE's
`1e-10 e` count-tol. **H-C1 finder-algorithm: CLEARED.**

H-C3 smearing-function: Cu deck uses F-D + pwdft-rs uses F-D →
matched. Cross-smearing gap on Cu (F-D vs Gauss vs MP1): ≤ 0.10 eV.
**CLEARED on Cu.** On insulators the cross-smearing gap is 1.1-2.2 eV
by construction (band edge, no DOS in gap) — expected, not a bug.

**Convention landmine (→ CLAUDE.md candidate):** QE's `sumkg.f90`
weights already include the spin-degeneracy factor (`wk *= degspin`
at `setup.f90:673`; `Σ wk = 2` for nspin=1), and `wgauss` returns
values in [0,1] per state with NO additional spin factor applied.
pwdft-rs uses weights summing to 1 + explicit `spin_factor = 2/nspin`.
A naive Python reference using QE's `wk` + pwdft-rs's `spin_factor=2`
double-counts by 2× (original Cu run landed at E_F = 13.42 eV,
off by 5.79 eV from QE).

**Remaining Part C suspects, reordered:**

1. H-C4 (Cu/GaAs/MgO/Mg ρ_core(G) unpinned; Fe has a guard, heavies
   don't). Highest prior — fingerprint matches ΔE_xc = +8.87 eV at
   transplant.
2. H-C5 (V_NL `D_ij · β·β` contraction on Cu d-projectors). Middle
   prior — VNMT tests l=2 m-isolation on Si but not Cu double-sum.
3. H-C2 (n_bands margin). Lower — Cu has 26 eV headroom above E_F.

**What the 1.82 eV DOS-origin Cu gap means now:** With finder cleared,
the 1.82 eV Fermi mis-gauge at ρ_QE is a consequence of pwdft-rs and
QE producing **different eigenvalue distributions** at ρ_QE. Γ matches
to 0.26 ± 0.03 eV; non-Γ k-points are untested and must differ by
more (likely driven by V_NL contraction or NLCC V_xc on the
non-Γ Hamiltonian). Part C session-2's V_NL matrix-element test on
Cu will probe this directly.

**Artifacts:** `scripts/validate/vgch2c_fermi_reference.py` +
`.csv` (30 rows), `proposals/VGCH-MECH-mechanism-taxonomy.md` +
`VGCH-2-total-energy-assembly.md` updates.

**Gates:** Tier-1 357 pass / 0 fail / 49 ign. Clippy 18 / 24 baseline.
rustdoc clean.

## 2026-04-19 — VGCH-2 Part B: H3 CLEARED on Cu transplant iter-1 (PR #167)

Seed pwdft-rs with QE's converged Cu FCC ρ_QE, diagonalize one SCF iteration, compare per-term to QE.

**Key numbers (post-VGCH-SiEF-B1 #166, merged during session):**

- E_HF at ρ_QE = −4837.30 eV, QE total = −4853.64 eV → **+16.34 eV gap at SAME density**. E_HF is gauge-invariant so this survives the V_loc(G=0) re-gauge.
- Γ eigenvalues now +0.26 ± 0.03 eV offset from QE (pre-SiEF-B1: −7.47 eV uniform; SiEF-B1 closed that).
- pwdft E_F = 21.29 eV, QE E_F = 19.21 eV → Δ = +2.08 eV. **1.82 eV is DOS/occupation origin** (subtract the 0.26 eV eigenvalue offset).
- Δρ(in vs out) = 0.124 e/Å³ — iter-1 ρ_out ≠ ρ_QE.
- ρ_in / ρ_out both integrate to 19.00 electrons.

**H3 CLEARED.** pwdft and QE give different E_HF at the same density, so the 16.6 eV Cu residual is not a mixer-basin effect.

**New leading Part C suspect:** Fermi-finder / smearing on dense 3d DOS at E_F (1.82 eV DOS-origin mis-gauge even at matched eigenvalues). Prior suspects (NLCC Cu/GaAs/MgO, projector scaling) still in play.

**Load-bearing conventions:**

- QE `charge-density.dat` Fortran sequential-access binary, `mill_g(3, ngm_g)` column-major → Rust reshape as C-order `(ngm, 3)`. Verified via `rho(G=0)·Ω ≈ N_el`.
- pwdft `1/N`-forward FFT normalization matches QE `fwfft('Rho', ..., dfftp)`; no extra scaling needed.
- ρ(G) unit: QE e/Bohr³ → pwdft e/Å³ via `1/BOHR_TO_ANG³` = 6.7483…
- numpy `.npz` uses ZIP64 for large entries (comp_size = 0xFFFFFFFF, real in extra field). Avoided by emitting a flat VGCH2BIN bundle.

**Artifacts:** `scripts/validate/vgch2_parse_qe_density.py`, `qe_validation/cu_rho_qe.bin`, `src/scf/transplant.rs`, `tests/vgch_transplant_cu.rs`.

**Next:** Python Fermi-Dirac bisection on Cu iter-1 eigenvalues vs pwdft `smearing::find_fermi_energy`; log QE `verbosity='high'` per-iter bracket. If Fermi-finder closes E_F, re-run Cu transplant and expect Δρ to collapse. If not, escalate to NLCC Cu/GaAs/MgO pin tests and projector cross-check.

## 2026-04-19 — VGCH-SiEF: Si E_F 1.35 eV offset localized to V_loc(G=0) gauge (PR #164)

Per-band δ_n = ε_n^pwdft − ε_n^QE on Si diamond (4×4×4 Γ-centered,
ecut=15 Ry, LDA) at Γ: **mean −1.3515 eV, std 0.6 meV** — pure
rigid offset. ΔE_F (pwdft − QE) = −1.3495 eV matches mean band
shift to 2 meV. E_total green at 45 meV unchanged.

**Mechanism:** pwdft-rs zeroes `v_local_fft[0]` before H assembly
(`src/scf/context.rs:125`) and compensates total energy via
`with_g0_shift` (`src/scf/energy.rs:248`). **QE keeps V_loc(G=0) in
vltot** (`qe-7.5/PW/src/setlocal.f90:91-96`, lines 91-92: `v_of_0 =
DBLE(aux(1))` is read out but `aux(1)` stays in the IFFT input), so
QE's KS eigenvalues carry the DC offset and pwdft-rs's do not.
Predicted shift: `−N_atoms · V_loc(G=0)_per_atom`. Matches Si
1.343 predicted vs 1.3515 observed (9 meV margin). Matches Fe
(VGCH-1a: ~5.1–5.3 eV; predicted 5.174 eV) and C diamond (VGCH-1c:
−3.09 eV; predicted 3.092 eV).

**Fix candidate A** (scoped for Part B1, ~20 LOC, 4 files): drop
zeroing in context.rs:125, drop `e_local_g0_shift` from driver(s),
no-op or delete `with_g0_shift`, unignore `test_si_diamond_fermi_vs_qe`.
Candidates B (Ewald G=0) and C (Fermi solver convention) both ruled
out by inspection — Ewald never touches Hamiltonian diagonal; Fermi
solver matches QE's ef.f90 bisection. σ-independence of the shift
(std 0.6 meV at σ=0.136 eV) also rules out C.

**Per-system predictions** (independent of heavy-atom VGCH-2 work):
Cu 7.75, Ga+As 4.70, Na+Cl 1.67, Mg+O 3.35 eV — all would close
with Part B1. Total-energy arms unaffected.

**Artifacts** (all read-only, no src/ changes):

- `scripts/validate/si_ef_shift_trace.py` + `si_ef_shift.csv`
- `proposals/VGCH-heavy-atom-vloc-residual.md` § Light-atom E_F shift

## 2026-04-19 — GGAP Phase F-light: 6 remaining PBE tests wired; VQEF matrix populated (PR #161)

Wired Al/C/Cu/GaAs/NaCl/MgO PBE tests in `tests/qe_validation.rs`
mirroring the LDA arms (same cell/k-grid/mixer; new PseudoDojo
NC/PBE PP + `XcFunctional::Pbe`). Also tightened Si PBE from 100 meV
→ 20 meV (observed 12.4 meV). Matrix goes `1G/8Y/8R → 1G/15Y/0R`.

**Per-system PBE residuals (eV):**

| System | \|ΔE\| PBE | \|ΔE\| LDA (ref) | Ratio (LDA/PBE) |
|--------|-----------|------------------|-----------------|
| Si     | 0.012     | 0.033            | 2.7 (both GREEN-class) |
| Al     | **0.108** | 0.075            | **0.69** — PBE WORSE |
| C      | 0.322     | 1.45             | 4.5 |
| Fe     | 1.97      | 11.5             | 5.8 |
| Cu     | 10.06     | 16.2             | 1.6 |
| GaAs   | 17.30     | 33.6             | 1.9 |
| NaCl   | 4.86      | 7.7              | 1.6 |
| MgO    | **1.56**  | 10.1             | **6.5** — largest closer |

**Physics findings (important, log them here so VGCH-2A / next
light-atom investigator can pick up):**

1. **Al is functional-insensitive.** Al PBE is *worse* than Al LDA
   by 33 meV. This kills the hypothesis that Al's VGCH light-atom gap
   is an XC artifact. Root cause must be density basin / projector /
   symmetry. No proposal yet for this split.

2. **Heavy-atom VGCH-2 residual is partially functional-sensitive.**
   All five heavy/semicore systems (Fe, Cu, GaAs, NaCl, MgO) show
   1.6×–6.5× PBE improvement over LDA. This tells VGCH-2 Part A that
   the partial-cancellation signature (Δone-e vs ΔE_H opposite-sign)
   has a gradient-term-sensitive component. MgO (6.5×, Mg 2s/2p
   semicore) is the cleanest signal; worth looking first at the
   Mg pp semicore region in the VGCH-2 trace.

3. **C is partially functional-sensitive.** 4.5× PBE improvement on
   the 1.45 eV LDA gap. Combined with Al's functional-insensitivity,
   the "light-atom VGCH class" is at least two mechanisms.

**Gates:** Tier-1 0-fail; clippy 18/24 baseline; doc clean; all 7 PBE
Tier-2 tests pass under `--ignored` (392s wall on M3 Max).

## 2026-04-19 — VQEF-AL: Al ecut=24 QE regen, VGCH light-atom reclassification (PR #146)

Regenerated QE 7.5 Al FCC reference at basis-converged cutoff (ecut=24 Ry = PseudoDojo `.standard`, up from 15 Ry) to close the basis-truncation side of VQEF Al.

**QE itself moves 48 meV from ecut=15 → ecut=24**, confirming ecut=15 was under-converged on the QE side. New reference: E_total = −4.72724484 Ry = −64.317443 eV, E_F = 7.5876 eV, 6 SCF iterations, 229 PWs at Γ.

**pwdft-rs vs QE at matched ecut=24 (8×8×8):** |ΔE| = **74.9 meV**, essentially unchanged from pre-regen 83 meV at ecut=15. The pure-basis-truncation hypothesis (pwdft ecut sweep projected ~27 meV at ecut=24) was only partly right: both codes carry basis-set truncation at ecut=15, but pwdft-rs's basis convergence slope vs QE's is not a simple offset — aligning the basis does not close the gap.

**Reclassification:** Al moves from "SYKP/MPSH basis truncation" to **VGCH light-atom "different converged density"** class — same family as C diamond's 1.45 eV gap (PR #144), opposite-sign Δone-e / ΔE_H signature per the per-component audit. 75 meV is much closer to 0 than heavy-atom VGCH (7.7–33.6 eV), but structurally the same shape: pwdft-rs and QE converge to different Kohn-Sham densities.

No Fermi level comparison this PR (|ΔE_F| = 145 meV).

## 2026-04-19 — VQEF-QC: Si/Al/C quickchecks with new ECUT + MIXL + MPSH tools (PR #144)

Per-arm measurement for three light-atom YELLOW cells in the VQEF scoreboard using new per-PP ecut policy (ECUT, PR #136), mixer event logging (MIXL, PR #133), and Γ-centered MP (MPSH, PR #110).

**Si LDA split (GREEN energy, YELLOW Fermi).** Single combined `test_si_diamond_vs_qe` ignored; now two arms: `test_si_diamond_energy_vs_qe` passing at ΔE = 33.2 meV (< 40 meV tol), `test_si_diamond_fermi_vs_qe` still ignored with |ΔE_F| = 1.3495 eV. Band-to-band differences at Γ agree with QE to <10 meV as the old `#[ignore]` text promised; absolute eigenvalue shift is V_loc(G=0) convention and lives in the Fermi arm only. **First GREEN scoreboard cell** (was 0 GREEN / 8 YELLOW / 8 RED).

**Al ecut sweep at 8×8×8 Γ-centered.** Seven configs vs QE ecut=15 Ry ref (E_QE = −64.269 eV):

| ecut (Ry) | E_pwdft (eV) | ΔE (meV) | iters |
|---|---|---|---|
| 15 (baseline) | −64.186 | 83.1 | 10 |
| 20 | −64.227 | 42.7 | 10 |
| 24 | −64.243 | 26.9 | 10 |
| 30 | −64.261 | 8.8 | 10 |

Mixer variants at ecut=15 agree on E to 0.001 meV — residual is pure basis-set truncation on pwdft's side. Hypothesis at the time: regenerate QE ref at ecut=24+ closes it (VQEF-AL #146 later proved this partially wrong — real gap is 75 meV at matched ecut, VGCH light-atom class).

**C diamond SCF stall fixed (YELLOW stays YELLOW with root cause).** Nine configs at ecut=30 Ry, 4×4×4 Γ-centered vs QE's E = −324.407 eV. Plain Anderson stalls at Δρ≈1.7e-8 past 150 iters. **Broyden+Kerker converges in 12 iters** at Δρ=8.0e-11 → E_pwdft = −322.957 eV, ΔE = 1.45 eV. Pinned `MixingMode::Broyden { kerker: true }` on the test. Per-component at convergence (eV, ours − QE): Δone-e = +1.76, ΔE_H = −0.59, ΔE_xc = +0.29, ΔE_ewald ≈ 0. **Opposite-sign split across one-electron and Hartree matches the Cu/Fe VGCH signature.**

Seeded the hypothesis that became VGCH Phase 1b/1c (light-atom extension): C's 1.45 eV gap is in the same "different converged density" class as Cu/Fe, NOT mixer/basis/V_loc(G=0).

## 2026-04-19 — VGCH Phase 1c H2 cleared — SAD initial density bit-perfect vs QE; VGCH-2 spawned for energy assembly

Added `scripts/validate/vgch_sad_heavy.py` (Python reference for
`qe-7.5/PW/src/atomic_rho.f90` recipe: Simpson over log mesh,
per-species structure factor, IFFT, G=0 renormalization) and
`tests/vgch_sad_heavy.rs` (Tier-2) that pins pwdft-rs'
`generate_initial_density` on all 7 VGCH-class systems (C, Al, Fe,
Cu, GaAs, NaCl, MgO).

**Verdict: H2 CLEARED.** Raw-sample point-wise max |Δρ(r)|:

| system | max |Δρ| (e/Å³) | neg mass clamped (e) |
| --- | --- | --- |  |  |
| C diamond | 5.4e-11 | 0 |  |  |
| Al FCC | 5.5e-12 | 0 |  |  |
| Fe BCC | 4.0e-10 | 0 |  |  |
| Cu FCC | 4.0e-10 | 0 |  |  |
| GaAs | 1.2e-5 | 2.1e-5 |  |  |
| NaCl | 8.6e-11 | 0 |  |  |
| MgO | 4.8e-10 | 0 |  |  |

C diamond bit-perfect pre-clamp AND post-clamp — the 1.45 eV C
residual does NOT live in SAD. GaAs's 1.2e-5 outlier is the
clamp step (`initial_density.rs:144-148`) zeroing O(1e-5 e) of
Gibbs ringing near the As core that QE keeps
(`atomic_rho.f90:186-188`: "useless to set negative terms to zero,
they re-appear on FFT round-trip"). Total-energy impact O(1e-5 eV).

**Debug journey.** Initial shell-average run showed C/GaAs with
0.2 / 0.07 e/Å³ asymmetry between the two atoms. Raw-sample
point-wise diff was bit-perfect — the asymmetry lived entirely in
the shell-average function, not the ρ(r) arrays. Root cause:
Rust's `.round()` is half-away-from-zero; numpy's `np.round` is
half-to-even. At `d = ±0.5` fractional displacement (grid points
exactly 1/2 cell from the atom), the min-image wrap produced
opposite-sign cartesian displacements in Rust vs Python, giving
different distances and therefore different bins. Fix: both codes
now use `d - (d + 0.5).floor()` wrap, which is platform-independent.
No production-code impact — shell-average is a diagnostic only.

**H3 target (VGCH-2).** With H1 (β_l(q), landed PR #148) and H2
(SAD, this PR #156) both cleared, the 7-34 eV residuals must be
in:

- (H3a) Total-energy assembly — `src/scf/energy.rs` `with_g0_shift`
  - ρ_in/ρ_out pairing in Harris-Foulkes double counting. Phase 1a
  fingerprint: opposite-sign one-electron vs Hartree partial
  cancellation is consistent with an assembly-pairing bug. C
  diamond specifically shows every Γ eigenvalue offset by exactly
  `−2·V_loc(G=0) = −3.09 eV`; expected compensation via
  `e_local_g0_shift = v_local_g0·N_el` works for Si but not C/Cu/Fe.
- (H3b) SCF mixer basin — escalated only if H3a clears.

Spawned `proposals/VGCH-2-total-energy-assembly.md`. Part A:
`scripts/validate/vgch_energy_assembly_heavy.py` + `tests/vgch_energy_assembly_heavy.rs`
to trace per-term agreement vs QE. Part B (transplant experiment)
only if Part A clears.

**Gates.** 342 passed / 0 failed / 31 ignored (all Tier-1 baseline
preserved). Clippy 18/24 at baseline. Rustdoc clean. New
`build_sad_density_for_diagnostic` + `..._verbose` public helpers
in `src/scf/initial_density.rs` expose SAD output + pre-clamp /
post-clamp snapshots + clamp/renorm statistics for integration
tests. No production-path change.

## 2026-04-19 — VGCH Phase 1b H1 cleared — β_l(q) is bit-perfect vs QE

Added `scripts/validate/vgch_beta_l_heavy.py` (QE-convention Simpson
over log mesh, matches `qe-7.5/upflib/beta_mod.f90:111-116`
byte-for-byte) and `tests/vgch_beta_l_heavy.rs` (pins pwdft-rs'
`bessel_transform_projector` output against the Python CSV).

**Verdict: H1 CLEARED.** Across 11 elements (Si, C, Al, Fe, Cu, Ga,
As, Na, Cl, Mg, O) × all projectors × 10 q-values in [0, 7] Bohr⁻¹ =
590 rows, **max |Δ| = 3.17e-12 Bohr^{3/2}** (Fe l=2 d-projector at
q=6). 4 orders below the 1e-8 tolerance. C diamond (early-verdict
case, no semicore) bit-perfect at 2.14e-12 Bohr^{3/2}. Semicore
shells (Cu 3s/3p/3d, Fe 3s/3p, Mg 2s/2p) all indistinguishable from
Si's lighter 4-projector layout. Form factors are NOT the bug.

**What remains.** C diamond E_total still 1.45 eV off QE; every pwdft
Γ eigenvalue is −3.16 eV vs QE (consistent with v_local_g0 = 2·1.5458
= 3.09 eV / cell). Band-structure GAPS agree to <50 meV (22.121 vs
22.169 eV Γ_v→Γ_c for C), so the physics is the same — the shift is
pure V_loc(G=0) absolute-reference convention, cancelled inside
E_total by `e_local_g0_shift = v_local_g0·N_el`. Phase 1a already
verified that cancellation works for Si; it also works for the absolute
eigenvalue column but NOT apparently for E_total on C/Cu/Fe. Since
VGCMP Phases 1-4 proved the full PP→H pipeline is bit-correct on Si
(10⁻⁸ eV), this residual is either (a) in SCF mixing dynamics reaching
different fixed points across materials, or (b) a subtle E_total
accounting term I missed.

**Phase 1b H2 (SAD initial density) — not yet run.** Script template
and harness pattern from H1 are directly reusable (parse UPF
PP_RHOATOM, compare against pwdft's `generate_initial_density`
output on the FFT grid). PP_RHOATOM integrates to z_valence exactly
(C=4, Cu=19, Fe=16, Na=9, Cl=7) for every VGCH PP — so the PP-level
atomic density is well-defined; the open question is whether the
Bessel-transform + FFT assembly in pwdft's SAD matches QE's atomic
superposition to high precision.

**Phase 1b H3 (mixer basin) — diagnostic-only.** Hardest to close;
the most informative test is to transplant QE's iter-1 density into
pwdft and see which fixed point pwdft reaches. Out of scope for this
session; proposed as VGCH-2 follow-up.

**Post-H1 scoreboard:** no `#[ignore]` flipped — H1 was diagnostic-
only. 18 / 24 clippy warnings (unchanged baseline). No src/ changes.

## 2026-04-19 — VGCH Phase 1a — heavy-atom per-component diagnostic (no fix yet)

Added `tests/vgch_per_component_heavy.rs` (Cu 4×4×4 and Fe 8×8×8
nspin=1) and `scripts/validate/vgch_vloc_heavy.py`. Fix NOT shipped —
Phase 1a completed; Phase 1b scoped.

**Per-component residuals (eV, ours − QE):**

| system | one-electron | E_H | E_xc | E_ewald | E_total |
|---|---|---|---|---|---|
| Fe 8×8×8 | +10.99 | −0.70 | +0.93 | +0.006 | **+11.50** |
| Cu 4×4×4 | +37.03 | **−24.57** | +4.73 | +3e-5 | +17.31 |

**Key finding:** residual splits across one-electron and Hartree with
OPPOSITE SIGNS — signature of a different converged density, not a
form-factor bug. Cu one-e +37 eV partially canceled by E_H −24.6 eV
→ net +17 eV.

**V_local(G=0) ruled out.** Python ref `(4π/Ω) ∫ r²[V+Ze²/r] dr`
agrees with `v_local_of_g(0, Ω)` to all printed digits for all 11
heavy-atom PPs (Si/Fe/Cu/Ga/As/Na/Cl/Mg/O/Al/C). Fe Z=16 ref: 5.1736
eV, pwdft: 5.1736 eV.

**Γ eigenvalue shift.** Every Fe eigenvalue offset ~5.1–5.3 eV vs QE,
matching V_loc(G=0) convention (pwdft zeros G=0 in H, compensated by
N·V_loc(G=0) additive). Does NOT contribute to E_total residual.

**Phase 1b hypotheses (next session):**

1. Different SCF fixed point — mixer/initial-density issue on
   heavy-atom cells, not a PP bug. Needs ρ(G) shell-by-shell diff
   between pwdft and QE save files.
2. E_nonlocal d-projector scaling — Cu E_NL = −508.80 eV vs Fe
   +37.86 eV is a 13× swing. Need beta_q reference extended to Cu
   and Γ-point assembled-H cross-check for Cu.
3. Semicore ecut under-convergence. Test by ecut sweep Fe/Cu at
   ecut ∈ {15, 25, 40, 60} Ry on both codes.

Five VGCH `#[ignore]` strings unchanged (already cite VGCH/TBD). No
src/ changes; no regression risk.

## 2026-04-18 — VQEF roadmap (PR #82) — proposal only

Scoped the 8 × 2 (system × functional) QE validation matrix as a gating roadmap. Current state at `origin/main` 0520c8c = **0 GREEN / 8 YELLOW / 8 RED** out of 16 cells. All 7 `#[ignore]` markers in `tests/qe_validation.rs` attributed to one of:

- **MPSH (unfiled)** — Si/C/Al residuals 0.26 eV / non-convergence / 73 meV; MP shifted vs Γ-centered grid. 1 CE-day fix (shift param in `KPointSettings`). Recommended P1 in the critical path.
- **VGCMP Phase 5 (active)** — Fe/GaAs/Cu/NaCl/MgO residuals 9.5 / 33.6 / 16.2 / 7.7 / 10.1 eV. Note MgO's 10 eV attributed to Mg 2s/2p semicore PP (not high Z per se).
- **GGAP B+C+D (draft)** — all 8 PBE cells; no PBE code yet, dispatcher Phase A in flight.

**Fe FM decision: Path C** — accept NM under LDA (PseudoDojo NC/LDA drives collapse at any reasonable ecut), validate FM under PBE only with PseudoDojo NC/PBE Fe at ecut=50 Ry, target M = 2.22 μB ± 0.05. Both tests coexist.

**`pseudopotentials/nc/pbe/` verified present with all 8 target elements** (72 PPs total). No external downloads. `qe_validation/pseudo-pbe/` symlink is a VQEF deliverable.

Total VQEF-owned work: ~7 res-days + 1 CE-day (QE ref generation + test arms + VGC5 cross-checks). Grand total incl. upstream MPSH + VGCMP Phase 5 + GGAP B-D: ~20-25 working days, compressible to ~3 calendar weeks with 2 parallel agents.

**Next handoff:** EM approval → file MPSH as a separate CE proposal (SYKP §D2 has the shape already drafted).

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
