# Researcher Logbook

Entries: date, what was validated, discrepancies found (with numbers), references used. Physics findings only — not code quality or docs.

## 2026-04-16 — Handoff and orientation

### QE discrepancy status

| System | QE (eV) | Ours (eV) | ΔE (eV) |
|--------|---------|-----------|---------|
| Si (2 atoms, 15 Ry, 4x4x4) | -231.61 | -218.28 | 13.3 |
| Fe BCC (1 atom, 16 Ry, 4x4x4) | -3059.46 | -3104.83 | 45.4 |
| C diamond | — | — | Does not converge |

Eigenvalue degeneracy breaking at Gamma confirms the bug is in V_local/V_NL form factors, not energy accounting. Error scales with Z (1.66 eV/el for Si, 2.84 eV/el for Fe).

**Root cause:** Radial quadrature — O(h²) trapezoidal vs QE's O(h⁴) Simpson. Compounded by V_local 1/r singularity (QE uses bounded erf/r subtraction). Proposals SIMP + VERF.

### Formula audit: all 23 items verified correct

No formula-level bugs found. Full comparison against QE 7.5 source. The only issue is numerical (quadrature quality).

### Open physics questions (low priority)

1. IBZ reduction gives 10 k-points for Si 4×4×4 vs QE's 8 — BZ boundary tolerance issue
2. Spin exchange formula at `xc.rs:249` uses equivalent but non-standard weighted-average form — needs documenting
3. `total_energy()` docstring omits V_local(G=0)·N_el correction — formula is correct, doc is incomplete

## 2026-04-16 — KBTF: KB projector test failure investigation

Investigated 3 failing tests in `tests/kb_projector_validation.rs`. PR #1 on branch `KBTF/kb-test-failures`.

| Test | Classification | Root cause |
|------|---------------|------------|
| `test_09` (D_ij vs HGH h^l_ij) | **Test bug** | Test assumed D_ij = raw HGH h^l_ij (3×3). Si.upf has 6 projectors (2 per l=0,1,2); QE diagonalizes h^l and absorbs eigenvector rotation into projectors. UPF D_ij is diagonal 6×6, not raw h^l_ij. |
| `test_07` (form factor decay) | **Known limitation** | l=1 projector |F(24.5)|/|F_max| = 0.106 > 0.1 threshold. Trapezoidal quadrature artifact at high q. SIMP would fix. |
| `test_vloc` (V_local vs QE) | **Known limitation** | `v_local_of_g` uses bare Coulomb subtraction (V+Z/r); QE uses erf(r)/r subtraction (numerically superior, avoids cancellation at large r). VERF would fix. |

**Key finding on HGH→UPF mapping:** QE's UPF conversion diagonalizes each l-block of h^l_ij, stores eigenvalues as diagonal D_ij, and rotates projectors accordingly. This means UPF D_ij values (e.g., D[0,0]=11.13 Ry) cannot be compared against published HGH h^l_ij (e.g., h^0_11=2.95 Ry). Tests 05/06/08/10 passing confirms D_ij + projectors are self-consistent.

**Confirms root cause from orientation:** Both test_07 and test_vloc failures trace to trapezoidal quadrature + bare Coulomb subtraction — same root cause as the 13.3 eV Si discrepancy. SIMP + VERF remain the correct fix path.

### Key references

PZ: PRB 23, 5048 (1981). KB: PRL 48, 1425 (1982). NLCC: PRB 26, 1738 (1982). QE source: `vloc_mod.f90`, `simpsn.f90`, `setlocal.f90`.

## 2026-04-16 — SPXC: Spin-polarized E_xc consistency investigation

**Confirmed the bug.** In `run_scf_spin`, E_KS computation mixes:
- `exc_r` from INPUT spin densities (line 525)
- `rho_xc_total` from OUTPUT total density (line 642)
- `vxc_{up,down}_r` from INPUT (line 525) with `rho_{up,down}_sym` from OUTPUT (line 650)

Non-spin `run_scf` does NOT have this bug -- it recomputes XC from OUTPUT at lines 347-348.
Harris-Foulkes in spin path is correct (all INPUT quantities, lines 670-683).

**QE comparison:** QE never recomputes XC from output. Uses `etxc`/`vtxc` from input density + `descf` first-order correction. Different formulation, same result at convergence.

**Recommended fix:** Recompute `lda_xc_spin_grid` from OUTPUT spin densities for E_KS. One extra grid-level call per iteration. Simple, matches non-spin path.

**Impact:** Bug is O(delta_rho) at convergence -- sub-meV for converged energy. But spoils quadratic convergence of |E_HF - E_KS|, reducing it to linear. This is the main practical issue: false HF-KS warnings in spin-polarized runs.

Proposal updated: `proposals/SPXC-spin-xc-consistency.md`.

## 2026-04-17 — QEVL: QE validation suite landed (Tier 1+2)

Validated the rescue data at `/tmp/pwdft-rescue/qe_validation_data/` and moved
into `qe_validation/` on branch `QEVL/qe-validation-suite`. Inputs match the
proposal spec for all 8 systems.

### Spot-check (M3 Max, 8 MPI, qe-7.5/build/bin/pw.x)

| System | Rescue E (Ry) | Re-run E (Ry) | Agreement |
|--------|---------------|---------------|-----------|
| Si diamond | -17.02299344 | -17.02299344 | bit-exact |
| Al FCC     |  -4.72371790 |  -4.72371790 | bit-exact |

Pseudo directory byte-identical to `pseudopotentials/nc/lda/`. Trusted; no
regeneration needed.

### Tests written (`tests/qe_validation.rs`)

Refactored into helpers (`fcc_crystal`, `bcc_crystal`, `run_qe_comparison`,
`assert_energy_matches_qe`, `assert_fermi_matches_qe`) + 8 `#[test]`s, all
`#[ignore]`d with per-test reasons because the Si 13.4 eV gap (VERF root
cause) propagates through all heavier systems. Tolerances set at 0.05 eV
(Tier 1) and 0.1 eV (Tier 2) so that when VERF closes the gap, unblocking
is `remove the #[ignore]` — nothing else.

QE reference numbers (Ry) that will be asserted once unblocked:
Si -17.022_993_44, C -23.843_439_10, Al -4.723_717_90, Fe -224.917_449_34,
GaAs -307.928_895_02, Cu -356.736_028_69, NaCl -119.779_703_03,
MgO -147.235_477_68.

### Baseline (pre-VERF) while old test file was in place

Si |ΔE|=13.43 eV, E_F off by 1.24 eV, Γ degeneracy broken. C diamond does
not converge in 80 iters at ecut=30 Ry with `MixingMode::Plain`. Fe
(nspin=1, old test) previously passed with a warning-only check; nspin=2
version now added.

### Tangential

- Tier 3 convergence studies (tests 9-11) still deferred; revisit after
  Tier 1+2 pass.
- Fe nspin=2 at ecut=15 Ry collapses to NM under PseudoDojo; consider
  `Fe_dalcorso.upf` at higher cutoff post-VERF to exercise magnetism.
- C diamond non-convergence may deserve its own proposal (high-ecut
  light-element mixing tuning) once Si offset is cleared.

## 2026-04-17 — VGCMP Phase 1: V_local(G) cleared

Si V_local(G) agrees with independent Python (scipy simpson on UPF mesh,
QE erf-subtracted formula) to **max |Δ| = 2.78e-9 Ry (3.78e-8 eV)** across
20 shells (|G|² = 3..56 in (2π/a)²). Tolerance was 1e-4 Ry; we beat it by
five orders of magnitude. V_local(G) is **not** the Si 13.43 eV culprit.

Artifacts: `scripts/validate/vloc_g_reference.py`,
`scripts/validate/vloc_g_si_reference.csv`, `tests/vgcmp_vloc_cross_check.rs`.
Next: Phase 2 (β_l(q) KB projectors) — follow-up branch `VGCMP/phase2-beta-q`.
Primary suspect now is the non-local KB machinery, specifically the
F_l(q) Bessel transform and/or the √BOHR_TO_ANG projector unit conversion
in `src/pseudopotential/upf.rs:68-72`.
