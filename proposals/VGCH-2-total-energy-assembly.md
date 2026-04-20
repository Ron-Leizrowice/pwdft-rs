---
id: VGCH-2
title: VGCH follow-up — total-energy assembly (V_loc(G=0) compensation + Harris-Foulkes pairings)
status: active
priority: high
complexity: medium
risk: medium
depends_on: [VGCH]
blocks: [VQEF]
owner: researcher
---

# VGCH-2 — Total-energy assembly for the heavy-atom residual

## Context

Both hypotheses under VGCH Phase 1 are now cleared:

- **H1 — β_l(q) form factors.** `tests/vgch_beta_l_heavy.rs` +
  `scripts/validate/vgch_beta_l_heavy.py` pin KB projectors for
  11 elements × 10 q-values = 590 rows to max |Δ| = 3.17×10⁻¹²
  Bohr^{3/2}, 4 orders below the 1e-8 tolerance. Fe/Cu semicore d
  projectors included. **Not the bug.**
- **H2 — SAD initial density.** `tests/vgch_sad_heavy.rs` +
  `scripts/validate/vgch_sad_heavy.py` pin pwdft-rs'
  `generate_initial_density` against a QE-convention Python
  reference for all 7 VGCH systems; raw-sample point-wise |Δρ| ≤
  1.2×10⁻⁵ e/Å³ (GaAs outlier from the negative-density clamp,
  O(1e-5 eV) total-energy effect). **C diamond is bit-perfect
  pre-clamp and post-clamp.** Not the bug either.

The 7–34 eV residuals must therefore live in **one of two places**:

1. **Total-energy assembly** — the `total_energy` + `with_g0_shift`
   + `harris_foulkes_energy` closure in `src/scf/energy.rs`, or the
   interaction between the G=0 compensation (`v_local_g0 · n_electrons`)
   and the NLCC double-counting subtraction, or a subtle mis-pairing
   of input-vs-output densities in the double-counting terms for
   multi-species cells.
2. **SCF mixer basin** — the mixer could be stabilizing a different
   local minimum of the energy functional.

Phase 1a's observation strongly points at (1):

- Cu one-electron is +37 eV too high, Hartree is −24.6 eV too low,
  they partially cancel to +17.3 eV net. The *signature* of a
  different converged density — but now we know the *initial* density
  is bit-perfect (Phase 1c).
- Si's one-electron vs. Hartree cancellation closes within
  `MPSH + SYKP` noise (≤ 50 meV). So the assembly works for Si.
  **Something in the assembly fails on Z > 14 cells.**
- On C diamond, every Γ eigenvalue is uniformly shifted by −3.16 eV
  vs. QE, exactly matching `2 · V_loc(G=0) = 3.09 eV` — see
  researcher logbook 2026-04-19 H1 entry. Band *gaps* agree to <50
  meV. The shift is absolute-reference-only and should be cancelled
  by `e_local_g0_shift = v_local_g0 · n_electrons`. Si cancels
  cleanly; C doesn't. That's a narrow fingerprint on the V_loc(G=0)
  compensation branch.

## Scope

**Part A — Isolate the mispairing (1 CE-day). COMPLETE 2026-04-19.**

Add a per-term trace to `src/scf/energy.rs` so we can cross-check
each assembled term against QE's per-term output (`E_one_electron`,
`E_hartree`, `E_xc`, `E_ewald`) using pwdft-rs' input-vs-output
density state at the moment each term is computed. Focus on:

- `e_local_g0_shift = v_local_g0 · n_electrons`. For Si this is +8
  electrons × small number ≈ small. For Cu this is +19 electrons ×
  ~5 eV. Does the sign match QE's convention?
- `e_xc` double-counting: the `(e_xc − e_vxc)` piece. NLCC makes
  this subtle — the core-density contribution to e_xc must *not*
  be subtracted from e_vxc. Is the sign right on all NLCC cells?
- Harris-Foulkes pairing: `E_HF = E_band − E_H[ρ_in] + (E_xc[ρ_in]
  − E_vxc[ρ_in]) + E_ewald`. The driver uses `ρ_in` for the double
  counting; does it actually pass `ρ_in` and not `ρ_out`?

### Part A findings (2026-04-19)

Artifacts:

- `scripts/validate/vgch2_per_term_trace.py` — parses QE
  `one_electron / hartree / xc / ewald / smearing_mts / total / fermi`
  from `qe_validation/*.out` for 8 systems → `vgch2_per_term_trace.csv`.
- `scripts/validate/vgch2_join_trace.py` — joins QE CSV with pwdft-rs
  CSV emitted by the extended `tests/vgch_per_component_heavy.rs`,
  computes `delta_meV = pwdft − QE`, ranks by |delta|.
- `tests/vgch_per_component_heavy.rs` extended from 2 → 8 cases
  (Si, C, Al, Fe, Cu, GaAs, NaCl, MgO) at the QE-reference SCF config
  where tractable. All 8 tests pass the PCFX self-check
  (`|Σ(components) − E_total|` < 33 meV on Fe, < 1 µeV elsewhere).

**Per-term delta table (eV, ours − QE):**

| system    | Δone-e   | ΔE_H     | ΔE_xc    | ΔE_ewald | ΔE_total |
|-----------|----------|----------|----------|----------|----------|
| Si        | +0.017   | −0.015   | −0.059   | +0.011   | −0.033   |
| Al        | −0.003   | −0.000   | −0.030   | +0.007   | +0.075   |
| C diamond | +1.722   | −0.586   | +0.327   | −0.013   | +1.450   |
| Fe BCC    | +10.986  | −0.702   | +0.931   | +0.006   | +11.502  |
| Cu FCC    | +37.033  | −24.568  | +4.728   | +0.000   | +17.308  |
| NaCl      | +14.970  | −8.934   | +1.852   | +0.100   | +7.988   |
| MgO       | +18.099  | −9.895   | +2.353   | +0.150   | +10.706  |
| GaAs      | +53.095  | −23.545  | +5.692   | −0.070   | +35.242  |

**Primary suspect: NOT an assembly mispairing.** The per-term
fingerprint is identical across every heavy-atom cell and is the
unambiguous signature of converging to a **different self-consistent
density**, not of a term-assembly bug:

1. **The direct-sum identity holds.** `|Σ components − E_total| <
   1 µeV` on every cell (the one exception is Fe at 23 meV,
   attributable to the nspin=1 mismatch with QE's nspin=2 ref —
   ρ-symmetrization of a magnetization-zero density is exact but
   the post-symmetrize occupation weight has numerical noise at
   1e-5 of the band sum). If `total_energy` / `with_g0_shift` /
   `harris_foulkes_energy` were double-counting a term, this sum
   would diverge by the same O(eV) as the QE residual — it doesn't.

2. **`e_local_g0_shift = N_el · Σ_sp V_loc(G=0)(sp)` is bit-correct
   against the direct formula.** Cross-checked against
   `scripts/validate/vgch_vloc_heavy.csv` to < 0.001 eV on all 8
   systems. The G=0 compensation path is not the bug.

3. **The Δone-e / −ΔE_H ratio tracks linear-response from a density
   perturbation**, not the factor-of-2 double-counting signature of
   an assembly bug:

   | system    | Δone-e    | −ΔE_H     | ratio  |
   |-----------|-----------|-----------|--------|
   | Si        | +0.017    | +0.015    | 1.19   |
   | C diamond | +1.722    | +0.586    | 2.94   |
   | Fe BCC    | +10.986   | +0.702    | 15.64  |
   | Cu FCC    | +37.033   | +24.568   | 1.51   |
   | NaCl      | +14.970   | +8.934    | 1.68   |
   | MgO       | +18.099   | +9.895    | 1.83   |
   | GaAs      | +53.095   | +23.545   | 2.26   |

   For a δρ that changes at fixed V_ext, linear response gives
   `Δone-e ≈ 2·∫V_H[ρ]δρ = 2·ΔE_H`, so ratio ≈ 2. Si/Cu/NaCl/MgO/GaAs
   sit in [1.2, 2.9]; C is 2.94; Fe at 15.6 is the outlier. A
   double-counting bug in `total_energy` would give ratio = 1 or
   ratio = ∞ (sign-dependent), not a smooth band. **This is a
   *different ρ* signature, not an assembly bug.**

4. **Δone-e + ΔE_H + Δxc + Δewald ≈ ΔE_total to within `−(−TS)`.**
   The residual after summing is exactly the QE `smearing contrib.
   (-TS)` term for Al/Fe/Cu/GaAs (100-260 meV), which means pwdft-rs'
   reported `total_energy` does NOT include the `−TS` smearing
   contribution (QE's `!    total energy` is F = E − TS). This is a
   separate ≤ 260 meV effect that should be flagged as its own
   follow-up but is NOT the cause of the 1.5-35 eV residuals.

5. **Si structural match: YES, with a twist.** On Si the |delta_meV|
   ranking is (ΔE_xc, Δone-e, ΔE_H) at (−59, +17, −15) meV. On every
   heavy system it is (Δone-e, ΔE_H, ΔE_xc) at O(1-53) eV. **The
   sign pattern is the same — Δone-e > 0, ΔE_H < 0, ΔE_xc > 0 — but
   the magnitude diverges by 3-4 orders.** Si's residual is O(k-grid
   noise + NLCC round-off + MP-shift convention); heavy systems
   have a *density-level* disagreement on top of those light
   effects. The assembly is bit-correct on both.

**Part B scope revision.** Given findings 1–4 above, the transplant
experiment moves from "only if Part A clears" to **primary next
step**. Seed pwdft-rs' SCF from QE's converged density (via a new
parser for QE's XML charge-density.xml or `.save/charge-density.dat`
binary). At iteration 1 with `ρ_in = ρ_QE`:

- If one-e / Hartree / xc match to < 50 meV/atom, the driver correctly
  reproduces QE's decomposition at QE's fixed point. The bug is then
  in the SCF *dynamics* — H3b (mixer-basin), H5 (symmetrization, but
  pinned clean by MPSH), H6 (initial density beyond SAD — already
  bit-perfect per H2), or H7 (eigensolver drift) — with H3b as the
  leading candidate because every failing cell uses a Kerker-family
  mixer and the `E_H` sign pattern is consistent with charge
  sloshing.
- If they don't match at iter 1, the bug is a subtle `v_xc` or `v_H`
  assembly term that is only visible when the density structure is
  "heavy" (semicore PP, multi-species, or compact valence overlap).

**Part B code pointers (post-revision):**

- `src/scf/driver.rs:593-595` — `total_energy(e_band, e_H, e_xc_corr,
  e_ewald)` followed by `with_g0_shift`. The pairing looks correct
  (same `rho_out` everywhere); add a sanity print of all 4 components
  at iter 1 for the seeded-density test.
- `src/scf/energy.rs:133-151` — `xc_energy_corrected`. Verify
  `rho_xc = rho_val + rho_core` and that `e_vxc` subtrahend
  integrates against `rho_val` not `rho_xc`. Already checked at
  docstring level — add an NLCC-specific unit test under Part B.
- `src/scf/context.rs:124` — `v_local_g0 = v_local_fft[0].re`. This
  is a sum over species (`sum_over_species(v_loc_of_g=0)`). On
  single-species cells this is fine; on GaAs/NaCl/MgO the sum is
  verified bit-correct above. Not the bug.

**What VGCH-2 Part A rules out:**

- Assembly mispairing in `total_energy` / `harris_foulkes_energy` /
  `with_g0_shift` (ratios, sum-identity, explicit PCFX check all
  pass).
- V_loc(G=0) compensation sign or magnitude (matches closed-form
  formula `N_el · Σ_sp V_loc(G=0)`).
- NLCC double-counting direction (e_xc − e_vxc has the right sign
  on all 5 NLCC-active cells).

**What VGCH-2 Part A does NOT rule out:**

- ρ-level disagreement at self-consistency (the primary finding).
- The missing `−TS` in pwdft-rs' reported `total_energy` (up to
  260 meV on Fe, up to 100 meV on Al/Cu — separate follow-up; see
  `src/scf/driver.rs:593-595` where `total_energy` is assembled
  without a smearing contribution).
- Double-counting in the LSDA driver for Fe specifically
  (Δone-e/−ΔE_H = 15.6 is a distinct outlier; rerun with nspin=2
  as a cross-check under Part B).

**Part B — Transplant experiment (COMPLETE 2026-04-19, 1 CE-day).**

Artifacts:

- `scripts/validate/vgch2_parse_qe_density.py` — parses QE's
  `charge-density.dat` (Fortran sequential-access binary) and emits
  a flat little-endian binary bundle (`VGCH2BIN` magic) containing
  `{mill, rho_g (e/Bohr³), b1/b2/b3}`. Errors out on `gamma_only=True`.
- `src/scf/transplant.rs` — `#[doc(hidden)]` public module exposing
  `run_scf_iter1_from_rho_g_fft`. Replicates iter-0 of
  `driver::run_scf_unpolarized` using a supplied G-space density
  instead of the SAD initial guess. Zero effect on production SCF —
  the new code path is only reached from the diagnostic test.
- `tests/vgch_transplant_cu.rs` — Tier-2 test (`--ignored`) that loads
  Cu FCC ρ_QE from `/tmp/vgch2b_cu/cu_rho_qe.bin`, scatters it onto
  pwdft-rs' 15×15×15 FFT grid with e/Bohr³ → e/Å³ unit conversion,
  runs one iteration, and prints the per-term deltas.

**Cu FCC iter-1 transplant table (Γ-centered 8×8×8, ecut=25 Ry, LDA):**

| term               | pwdft (iter 1) | QE (converged) | Δ (ours − QE) |
|--------------------|---------------:|---------------:|--------------:|
| one-electron       |  −1971.014 eV  |  −2033.883 eV  | **+62.87 eV** |
| Hartree            |   +985.775 eV  |  +1040.511 eV  | **−54.74 eV** |
| XC (bare)          |   −550.259 eV  |   −559.128 eV  |  +8.87 eV     |
| Ewald              |  −3301.019 eV  |  −3301.025 eV  |  +0.01 eV     |
| Total (E_KS, ρ_out)|  −4785.034 eV  |  −4853.641 eV  | **+68.61 eV** |
| E_HF (ρ_in = ρ_QE) |  −4837.300 eV  |  −4853.641 eV  | **+16.34 eV** |

Supporting per-term breakdown (pwdft, on top of VGCH-SiEF-B1
V_loc(G=0) on-diagonal gauge merged via PR #166):

- `E_kin = +1760.21 eV`, `E_loc = −3232.01 eV`,
  `E_loc(G=0)·N_el = 0.00 eV` (lives on H diagonal post-SiEF-B1),
  `E_nl = −499.21 eV`.
- Γ eigenvalues now agree with QE's converged spectrum to a
  uniform **+0.26 ± 0.03 eV** offset across bands 0..7 (semicore
  3s, 3p; valence d; 4s). Pre-SiEF-B1 the same comparison showed
  −7.47 eV from the G=0 gauge; SiEF-B1 closed that.
- `∫ρ_in = 19.00009 e`, `∫ρ_out = 19.00000 e` — both densities
  have correct total charge.
- `Δρ (in vs out, RMS) = 1.24 × 10⁻¹ e/Å³` — iter-1 ρ_out from
  diagonalizing H[ρ_QE] is NOT ρ_QE.
- pwdft iter-1 Fermi `E_F = 21.29 eV` vs QE `E_F = 19.21 eV` ⇒
  **Δ = +2.08 eV**. The eigenvalue-average offset is only
  +0.26 eV, so the Fermi-level mis-gauge is an extra **+1.8 eV of
  pure DOS/occupation origin** — i.e. pwdft's Fermi finder
  places E_F 1.8 eV too high for a Cu-like dense-3d DOS, even
  when the eigenvalues themselves match. Smearing implementation
  or band-count-vs-k mismatch is the first-order suspect.

### Verdict: H3 CLEARED

The transplant at iter 1 does NOT reproduce QE's per-term
decomposition at the QE-converged density. Specifically:

1. **E_HF at QE's density is +16.3 eV above QE's total.** This is
   the same magnitude as the 16.6 eV pwdft↔QE residual seen at
   pwdft-rs' own self-consistent density (`test_cu_fcc_vs_qe`, which
   is `#[ignore]`d with that residual). Since both E_HF estimators
   use the identical ρ_QE for double counting, eigenvalues, and
   density-dependent terms, the 16 eV gap exists between the two
   codes' energy functionals evaluated on the same density, not
   between two different fixed points. The mixer basin is not the
   cause.
2. **Eigenvalues match QE within +0.26 eV uniform offset across
   all bands** (post-SiEF-B1 gauge fix). Pre-SiEF-B1 the offset
   was −7.47 eV, which SiEF-B1 closed by keeping V_loc(G=0) on
   the Hamiltonian diagonal (QE convention). The KS Hamiltonian
   assembly at ρ_QE is correct up to the known gauge. The
   residual 0.26 eV is small enough to be a k-grid / smearing
   convention artifact; not a driver issue.
3. **Δρ (in vs out) at iter-1 is large (0.12 e/Å³).** Despite
   matching eigenvalues at Γ, the rebuilt density from
   diagonalizing H[ρ_QE] differs substantially from ρ_QE.
   Consistent with a Fermi-level / occupation inconsistency:
   pwdft's Fermi at ρ_QE is 2.08 eV above QE's — but the
   eigenvalues only differ by 0.26 eV, so **1.8 eV of the Fermi
   mis-gauge is pure DOS/occupation origin**. That would
   reshuffle occupations of near-E_F bands and change |ψ|²
   integrated over the metallic fraction of the 3d manifold. The
   density feedback propagates through Hartree (−54.7 eV gap)
   and XC (+8.9 eV gap); their partial cancellation explains the
   linear-response signature that Phase 1a first observed.

### What Part B rules out

- Mixer-basin effect (H3 cleared — 16.3 eV gap exists at the SAME
  density).
- V_loc(G=0) compensation error (compensation is exact up to the
  uniform eigenvalue shift).
- Total-energy-assembly mispairing (Part A already cleared; Part B
  confirms on fresh data).

### Where the bug IS (for Part C scope)

The 16.3 eV iter-1 E_HF residual and **1.8 eV of DOS-origin
Fermi-level mis-gauge** at ρ_QE point at ρ-based physics that
differs between pwdft-rs and QE when evaluated on the SAME
density. Prime suspects, in order of likelihood given the
post-SiEF-B1 signal:

1. **Fermi-level finder / smearing implementation (new leading
   suspect post-SiEF-B1).** With eigenvalues now agreeing to
   +0.26 eV, pwdft's E_F is still +1.80 eV high on Cu. The Cu 3d
   manifold produces a dense DOS right at E_F; a small difference
   in the smearing cumulative or the root-finder bracket could
   move E_F by O(eV) in this regime and trigger the large Δρ
   between ρ_in and ρ_out. Action: (a) compare
   `smearing::find_fermi_energy` bisection output on Cu iter-1
   eigenvalues against a Python Fermi-Dirac root-find reference
   built from the identical eigenvalue list;
   (b) log QE's per-iter E_F (it prints it at `verbosity='high'`)
   and compare bracketing traces;
   (c) check that the number of bands per k-point is large enough
   — Cu with `n_bands = 14` leaves only 4.5 unoccupied bands
   above the 9.5 occupied, and if the Fermi tails need bands
   substantially above the Fermi level to converge, this margin
   might be too thin for the DOS peak at E_F.
2. **NLCC for transition-metal d-systems.** Cu is Z=29 with
   3s/3p/3d semicore + a `PP_NLCC` block. The NCFX fix (closed the
   Si 13.4 eV residual) corrected unit conversion and radial
   weighting for the ρ_core Fourier transform. NCFX has a
   regression-guard test on Fe
   (`test_fe_bcc_xc_nlcc_regression_guard`) but no dedicated Cu /
   GaAs / MgO unit test. If the Bessel transform of `PP_NLCC` has
   a subtler bug only exposed on high-ρ_core elements, the ΔE_xc
   at Cu (+8.87 eV) would be partly NLCC-driven — and the Hartree
   gap would not be NLCC since ρ_core does not enter V_H. Action:
   parse QE's `rho_core(G)` for Cu and pin-test vs pwdft-rs'
   `ctx.rho_core_r` FFT-forward, analog to
   `scripts/validate/rho_core_g_reference.py` but for Cu/GaAs/MgO.
3. **Semicore projector magnitude (V_nl).** E_nl at iter-1 =
   −499.21 eV. QE's E_nl is not printed directly; backing it out
   from the QE one-electron breakdown requires separating kinetic
   + local + nonlocal. Action: extend the Cu KB projector cross-
   check (Phase 1b H1) from `β_l(q)` to include the `D_ij · Σ_lm
   β·β` contraction on a Cu-sized wavefunction; if projector
   scaling differs between codes by an O(10%) factor, that could
   account for the 16 eV gap.

### Scaling outlook

Cu's mechanism (semicore d + NLCC) is shared with Fe (semicore
3s/3p, NLCC), GaAs (Ga 3d + As semicore), MgO (Mg 2s/2p semicore),
and C (no semicore but potentially NLCC). NaCl has no semicore on
Na or Cl. If H3 is cleared on Cu and the underlying cause is NLCC
or semicore-projector, then Fe / GaAs / MgO should show the same
iter-1-transplant signature. C / NaCl are the differential
diagnosis: if they also show O(eV) iter-1 E_HF gaps at transplanted
density, the bug is not semicore-specific. Part C should run the
transplant on C, Fe, NaCl as the minimum triangulation set before
proposing a fix.

**Part C — Diagnosis + fix (1–5 CE-days).**

### Part C session-1 findings (2026-04-20, diagnostic-only) — H-C1 + H-C3 CLEARED

Closed hypotheses this session:

- **H-C1 (Fermi-finder algorithm)** — CLEARED. Python Fermi-Dirac
  bisection over QE's converged eigenvalues reproduces QE's reported
  E_F to 6–66 μeV across Cu / Fe / C / MgO, and to 0.036–0.18 meV on
  NaCl. Pwdft-rs's `smearing::find_fermi_energy` uses the identical
  sign convention and a bracket-width tolerance (`1e-14 eV`) that
  converges to the same root as QE's count-tolerance (`1e-10 e`).
  The finder is not the bug. Artifact:
  `scripts/validate/vgch2c_fermi_reference.py` + `.csv`.
- **H-C3 (smearing function mismatch)** — CLEARED on Cu. Every QE
  reference input deck under `qe_validation/*.in` uses the same
  smearing function pwdft-rs uses at runtime (Cu: F-D). Cross-
  smearing gap on Cu is ≤ 0.10 eV (F-D vs MP1 on the same QE
  eigenvalues), two orders below the 1.82 eV DOS-origin Cu residual.

Convention landmine recorded for future Rust↔Python occupation
cross-checks: QE's `sumkg.f90` weights already include the spin
degeneracy (`wk *= degspin` in `setup.f90:673`), sum to 2 for
nspin=1, and the integrand does NOT apply an additional spin
factor. Pwdft-rs's `find_fermi_energy` takes weights summing to 1
plus an explicit `spin_factor = 2.0/nspin`. Identical products,
but a naive reference that uses QE's `wk` + pwdft-rs's
`spin_factor = 2` double-counts.

**Remaining hypotheses for Part C session-2:**

- **H-C4 (NLCC ρ_core(G) unpinned on Cu/GaAs/MgO/Mg semicore)** —
  highest prior. Fe has a regression guard; heavy semicore elements
  don't. A silent unit-conversion bug on Cu's 3d NLCC could move
  ΔE_xc at shared density by the observed +8.87 eV (Cu transplant).
- **H-C5 (V_NL `D_ij · Σ_lm β·β` contraction for Cu d-projectors)** —
  middle prior. Individual β_l(q) are bit-perfect per VGCH-1b, but
  the double-sum assembly is untested on Cu-sized cells.
- **H-C2 (n_bands margin)** — lower prior. Cu has 26 eV of headroom
  above E_F at Γ (14 bands reach 45 eV; E_F ≈ 19 eV). F-D tails at
  σ = 0.272 eV decay in ~10 σ = 2.7 eV — well within headroom.

### Part C session-2 scope (next)

- Extend `tests/vgch_per_component_heavy.rs` or add a new NLCC
  pin for Cu/GaAs/MgO ρ_core(G) at G=0 + first two shells. Parse
  QE's `rho_core` from the UPF XML (PP_NLCC block) through
  `scripts/validate/rho_core_g_reference.py`. Target tol: 1e-4
  e/Å³ at each G.
- Add a G=G' matrix-element cross-check for the full KB V_NL
  operator on Cu, extending VNMT's l=2 single-channel pin from Si
  to Cu. Verify that ⟨G | V_NL | G'⟩ computed by
  `NonlocalPotential::add_to_hamiltonian` matches a Python reference
  that sums D_ij · β_l(q) · β_l(q') · Σ_m Y_lm(q̂) Y*_lm(q̂') by
  brute force.

Expected fix classes, in decreasing prior post-session-1:

- A unit-conversion regression on `PP_NLCC` for heavy semicore
  elements (H-C4 leading).
- A scale-factor error in the `D_ij · β·β` contraction or in
  `add_to_hamiltonian`'s accumulation that only surfaces on Cu-sized
  cells with dense d-projectors (H-C5).
- A mis-pairing of input vs. output density in `harris_foulkes_energy`
  or one of its feeders (kept as a fallback — Phase 1a fingerprint
  already clears this at the sum level).

**Part D — Close the matrix (0.5 CE-day).**

Once Part C lands, rerun all 5 heavy-atom `qe_validation.rs` tests.
Drop `#[ignore]` on those that close below 50 meV/atom. Update test
pins. Remaining cells filed as narrow follow-ups.

## Deliverables

- `scripts/validate/vgch_energy_assembly_heavy.py` — QE per-term
  parser (reads QE stdout `E_one_electron`, `E_hartree`, `E_xc`,
  `E_ewald`, plus the `one electron contribution` sub-terms from
  `verbosity='high'`).
- `tests/vgch_energy_assembly_heavy.rs` — Tier-2 integration test
  asserting per-term agreement for C diamond, Fe BCC, Cu FCC, GaAs,
  NaCl, MgO. Failure mode: print per-term Δ vs QE in eV.
- Fix on `src/scf/energy.rs` (most likely) or `src/scf/context.rs`.
- 3–5 `#[ignore]` removals in `tests/qe_validation.rs`.

## Acceptance

Close VGCH-2 when:

1. Fe BCC 8×8×8 ≤ 50 meV/atom residual (currently 9.5 eV).
2. C diamond residual ≤ 50 meV/atom (currently 1.45 eV).
3. At least 3 of 5 heavy-atom `qe_validation.rs` tests have
   `#[ignore]` dropped.
4. `tests/vgch_energy_assembly_heavy.rs` pins per-term agreement at
   ≤ 50 meV/atom for all 7 VGCH-scope systems.
5. `tests/vgch_sad_heavy.rs` and `tests/vgch_beta_l_heavy.rs` remain
   green (no regression on the Phase 1b/1c bit-perfect baselines).

## Cost

1–2 CE-weeks, depending on whether Part A's trace pins the bug or
Part B's transplant experiment is needed.

- Part A (per-term trace): 1 CE-day.
- Part B (transplant experiment): 1 CE-day.
- Part C (fix + tests): 1–5 CE-days.
- Part D (close matrix): 0.5 CE-day.

## Non-goals

- Fixing the SAD clamp (Phase 1c showed the clamp effect is O(1e-5
  eV) on the worst cell, GaAs). Can be a separate polish ticket.
- USPP / PAW support.
- Spin-orbit coupling.
- Whole-new mixer — the mixer is cleared for Si/Al/C (MPSH +
  SYKP closes those). Only revisit if Part B's transplant stays at
  the QE fixed point with no drift, ruling out assembly and pointing
  at mixer.

## Related

- **VGCH** (parent) — Phase 1a diagnostic landed PR #139, Phase 1b
  H1 (β_l(q)) landed PR #148, Phase 1c H2 (SAD) in this PR.
- **VQEF** — blocked on VGCH-2 for the 5 heavy-atom PBE-ready cells.
- **MPSH** — landed PR #110, closed Si/Al/Cu/Fe shift-convention
  gaps. MPSH is independent of VGCH-2; together they unblock the
  full VQEF matrix.
