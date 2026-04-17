---
id: NLCC
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# NLCC: Nonlinear Core Correction Audit

## Problem

Nonlinear core correction (NLCC; Louie, Froyen, Cohen, PRB 26, 1738 (1982)) is
required for *correct* LDA/GGA total energies and magnetic moments on elements
where the valence and core densities overlap — typically transition metals
(Fe, Cu, Ni) and alkali metals with core $ns$ + valence $n(s+1)$ configurations.
Without NLCC, $E_{xc}[\rho_{\text{val}}]$ misses the nonlinear coupling
$E_{xc}[\rho_{\text{val}} + \rho_{\text{core}}] - E_{xc}[\rho_{\text{val}}]$,
which can be tens of meV to eV for 3d metals.

**Current status:** pwdft-rs has NLCC infrastructure at
`src/pseudopotential/upf.rs:94-111` (parser), `src/pseudopotential/mod.rs:36-89`
(data + `has_nlcc()`), `src/scf/potentials.rs:61-113` (grid placement), and
`src/scf/energy.rs:53-73` + `src/scf/energy.rs:153-165` (XC consumption). The
spin-polarized path at `src/scf/mod.rs:512-526,664-665` splits $\rho_{\text{core}}/2$
per channel. All of this was audited in this proposal and found correct —
see §"Audit findings" below.

**Gap:** the Si test suite (which dominates coverage) has no NLCC, so the code
path is **entirely untested against QE**. The existing tests for NLCC — if any
— only check that the `rho_core_r` array is non-empty and finite (cursory).
No end-to-end validation has been done against a QE reference for an
NLCC-requiring element.

This proposal is a **documentation + test-coverage pass**, not a bug fix.
The audit below justifies not needing a code change.

## Background

Louie-Froyen-Cohen NLCC prescription (PRB 26, 1738 (1982), Eq. 5):

$$
E_{xc}^{\text{NLCC}} = \int \epsilon_{xc}[\rho_{\text{val}}(r) + \rho_{\text{core}}(r)]\,
    \big(\rho_{\text{val}}(r) + \rho_{\text{core}}(r)\big) \, d^3r
$$

but the XC *potential* is constructed from the same total density and acts
**only on valence orbitals** (core is frozen). The double-counting correction
in the KS total energy is

$$
E_{\text{dc}} = \int \rho_{\text{val}}(r)\, v_{xc}[\rho_{\text{val}} + \rho_{\text{core}}](r)\, d^3r
$$

so the final contribution is

$$
E_{xc} - E_{\text{dc}} = E_{xc}[\rho_{\text{val}} + \rho_{\text{core}}]
    - \int \rho_{\text{val}}\, v_{xc}[\rho_{\text{val}} + \rho_{\text{core}}]\, d^3r.
$$

Key invariants the implementation must preserve:
1. $\rho_{\text{core}}$ is **not** added to the Hartree source (core is
   assumed not to contribute to $V_H$ at the NLCC approximation level).
2. $\rho_{\text{core}}$ is **not** added to the valence electron count
   (normalization).
3. In LSDA, $\rho_{\text{core}}$ is split equally between channels
   because the core is spin-unpolarized.
4. $\rho_{\text{core}}$ enters only $E_{xc}$ and $v_{xc}$.

QE reference: `qe-7.5/PW/src/v_of_rho.f90:440-620` (subroutine `v_xc`), which
adds `rho_core(ir)` to `rho%of_r(ir,1)` before calling `xc_lda`, and subtracts
it before return (lines 511, 523, 540, 579).

## Audit findings

### Code paths inspected

| File / line | Role | Verdict |
|---|---|---|
| `src/pseudopotential/upf.rs:94-111` | UPF NLCC parser | Correct |
| `src/pseudopotential/mod.rs:36-89` | `PseudopotentialData::core_charge`, `has_nlcc()` | Correct |
| `src/scf/potentials.rs:58-113` | `compute_core_density` | Correct |
| `src/scf/context.rs:105-110,137` | Core density stored on ScfContext | Correct |
| `src/scf/energy.rs:53-73` | `xc_energy_corrected` | Correct |
| `src/scf/energy.rs:153-165` | `add_core_density` | Correct |
| `src/scf/mod.rs:275-284,349-350` | Non-spin XC uses $\rho_{\text{val}} + \rho_{\text{core}}$ | Correct |
| `src/scf/mod.rs:512,525-527,664-665` | Spin XC uses $\rho_{\text{val}\sigma} + \rho_{\text{core}}/2$ | Correct |
| `src/scf/initial_density.rs:70-85,149-184` | Initial density **does not** include $\rho_{\text{core}}$ | Correct |
| `src/scf/density.rs:93` | Density normalized to $N_{\text{val}}$, not $N_{\text{val}} + N_{\text{core}}$ | Correct |

### Specific invariant checks

**(1) Hartree does not see $\rho_{\text{core}}$:**
`src/scf/mod.rs:269` (non-spin) passes `rho_g` (built from `rho_r`, the
*valence* density) to `hartree_on_fft_grid`. `src/scf/mod.rs:522` (spin)
builds Hartree from `rho_total_r = rho_up_r + rho_down_r` which is also
valence-only. Confirmed correct.

**(2) Electron count excludes core:**
`src/scf/density.rs:93` normalizes to `n_electrons = grid.n_electrons`.
Caller supplies `n_electrons = sum_atoms(pp.z_valence)` at
`src/scf/context.rs` — $Z_{\text{valence}}$, not $Z_{\text{total}}$. Confirmed.

**(3) Spin split is $\rho_{\text{core}}/2$:**
`src/scf/mod.rs:512` `rho_core_half = rho_core / 2`. Confirmed.

**(4) XC-only consumption:**
`add_core_density` (`src/scf/energy.rs:155`) is called only at XC sites
(`src/scf/mod.rs:276,349,525,526,654,664,665,693`), never for V_H or
density-normalization. Confirmed by grep.

**(5) NLCC in both E_KS and E_HF:**
Non-spin: `rho_for_xc` (INPUT) for E_HF, `rho_new_for_xc` (OUTPUT) for E_KS —
both use `add_core_density`. Spin: `rho_up_xc` + `rho_down_xc` (INPUT) for
E_HF, `rho_up_xc_out` + `rho_down_xc_out` (OUTPUT) for E_KS — both use
`add_core_density(.., rho_core_half)`. Confirmed consistent with SPXC fix.

### Matches QE convention

`qe-7.5/PW/src/v_of_rho.f90:511` does `rho%of_r(ir,1) = rho%of_r(ir,1) + rho_core(ir)`
before the XC call, then line 523/540/579 subtracts it back. pwdft-rs's
`add_core_density` is the pure-function equivalent (returns a new vector
instead of mutating in place). No algorithmic difference.

### Conclusion

**No bug found.** The NLCC path is mathematically correct end-to-end. The
risk is purely on the test-coverage side: there is no regression against
QE for an NLCC-requiring element, so a future refactor could silently
break NLCC and no test would catch it.

## Proposed work (test coverage and docs)

### Part A — unit tests for invariants

Add `src/scf/energy.rs` `#[cfg(test)] mod tests` cases:

1. **`test_add_core_density_empty_is_identity`.**
   `add_core_density(&rho_val, &[])` returns `rho_val` unchanged.
2. **`test_add_core_density_sums_elementwise`.**
   For $\rho_{\text{val}} = [1, 2, 3]$, $\rho_{\text{core}} = [0.1, 0.2, 0.3]$:
   output is $[1.1, 2.2, 3.3]$.
3. **`test_add_core_density_clamps_nonnegative`.**
   For $\rho_{\text{val}} = [-0.5, -0.1, 0.1]$, $\rho_{\text{core}} = [0.3, 0.0, 0.0]$:
   output is $[0.0, 0.0, 0.1]$ — the clamp at `src/scf/energy.rs:162` fires.
4. **`test_xc_energy_corrected_with_nlcc`.**
   Construct a toy case where $\rho_{\text{val}} \ne \rho_{\text{xc}}$ and
   verify $E_{xc}$ uses $\rho_{\text{xc}}$ (val+core) while the double-counting
   subtracts $\int \rho_{\text{val}} v_{xc}$ (not $\int \rho_{\text{xc}} v_{xc}$).

### Part B — integration test against QE for an NLCC element

Add `tests/nlcc_fe_validation.rs` (or extend `tests/qe_validation.rs`):

1. Select a pseudopotential with NLCC. Options:
   - `pseudopotentials/nc/lda/Fe.upf` if it has `core_correction="T"`; or
   - `Fe_dalcorso.upf` (researcher logbook 2026-04-17 mentions this as a
     higher-cutoff Fe PP with better magnetism behavior). Needs to be
     added to `pseudopotentials/nc/lda/` if not already there.
2. Run an SCF on Fe BCC at ecut high enough to converge the relevant
   orbitals (logbook suggests ≥ 30 Ry for Fe). `qe_validation/fe_bcc.in`
   already exists.
3. Compare total energy to QE reference with a 0.1 eV tolerance (same as
   existing Tier-2 qe_validation tests).
4. Add a companion test that **disables NLCC** in the same PP (manually
   zero the `core_charge` field after loading) and asserts the energy is
   *different* by > 0.1 eV. This proves the test actually exercises the
   NLCC pathway.
5. This test will be `#[ignore]`d initially (behind VERF and any remaining
   per-component discrepancies), exactly like the existing
   `tests/qe_validation.rs` tests. The cost is zero until unblocked.

### Part C — documentation pass

1. Extend the module docstring in `src/scf/energy.rs` to include the NLCC
   formula and cite Louie-Froyen-Cohen PRB 26, 1738 (1982).
2. Add a note to `src/scf/potentials.rs` `compute_core_density` docstring
   explaining the (1/Ω) normalization convention and citing the QE
   subroutine name (`v_of_rho.f90:v_xc`).
3. Update the main SCF module docstring (`src/scf/mod.rs`) to mention NLCC
   as one of the supported features.
4. Add a `CLAUDE.md` / README mention: "LDA + NLCC is supported and used
   automatically when the UPF file has `core_correction="T"`."

## Risk assessment

- **Zero behavioral risk for Parts A and C** — test-only and documentation-only.
- **Part B** may reveal a real discrepancy when run against QE. If so, the
  audit above is wrong and this proposal morphs into a bug-fix proposal
  with evidence; that is the desired outcome of a real test.
- **Part B cost depends on available Fe PP**. If `Fe.upf` in the current
  pseudopotential set has `core_correction="F"`, we need to either
  (a) regenerate it from PseudoDojo with NLCC enabled, or (b) download
  `Fe_dalcorso.upf` from the QE site. Option (b) is ~1 MB and reproducible.

## Verification plan

The verification is the work itself — Parts A, B, C *are* the test coverage.
Additionally:

- `cargo test add_core_density` should show 3 passing tests after Part A.
- `cargo test --ignored nlcc_fe` should show 1 passing test when unblocked.
- `cargo doc --open` on `src::scf::energy` should surface the LFC reference.

## Estimated effort

4 hours total:
- Part A (unit tests): 1 h.
- Part B (integration test + Fe PP setup): 2 h.
- Part C (docs): 1 h.

If Part B uncovers a real bug, budget +4 h for investigation (triangulate
with the NLCC-off control test, then close the loop via a new proposal).

## Success criteria

1. Four new unit tests for `add_core_density` and NLCC-aware
   `xc_energy_corrected` are green.
2. One integration test exercises NLCC end-to-end against QE reference
   data (ignored until VERF blockers lift; passes when unignored).
3. Module docstrings in `src/scf/energy.rs` and `src/scf/potentials.rs`
   cite Louie-Froyen-Cohen 1982.
4. The audit findings in this proposal are linked from the `src/scf/mod.rs`
   module docstring so future readers know NLCC has been validated.
