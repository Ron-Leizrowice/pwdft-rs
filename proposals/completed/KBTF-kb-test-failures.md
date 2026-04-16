---
id: KBTF
status: completed
priority: critical
complexity: small
risk: low
depends_on: []
blocks: [SIMP]
---

# KBTF: Investigate 3 Failing KB Projector Validation Tests

## Problem

Three tests in `tests/kb_projector_validation.rs` are failing. These must be understood and resolved before SIMP (Simpson quadrature) work begins, because SIMP modifies the same radial integral machinery and we need a clean baseline to distinguish pre-existing bugs from regressions.

## Failing Tests

### 1. `test_09_hgh_parameter_crosscheck` — D_ij mismatch

**Symptom:** D[0,0] parsed = 151.46 eV, expected = 40.18 eV (~3.77x ratio).

The test hard-codes HGH h^l_ij reference values from PRB 58, 3641 (1998) Table I and compares against the D_ij parsed from `Si.upf`. The 3.77x factor is suspiciously close to structural constants that appear in the HGH→UPF D_ij transformation.

**Investigation needed:**
- Determine whether the UPF D_ij includes normalization factors not present in the raw HGH h^l_ij (e.g., overlap integrals of the projector Gaussians). QE's `upflib/upf_to_internal.f90` may apply such factors.
- Check if this is a **test bug** (wrong reference values) or a **parser bug** (wrong conversion).
- The non-local potential tests (05, 06, 08, 10) all pass, which suggests the parsed D_ij and projectors are internally consistent — the test's expected values may simply be wrong.

### 2. `test_07_form_factor_behavior` — High-q decay failure

**Symptom:** Projector 2 (l=1) form factor |F(q_max)|/|F_max| = 0.106, threshold is 0.1.

The test asserts that F(q) decays to <10% of its peak at q = 24.5 Å⁻¹. Projector 2 (l=1) barely fails this threshold. The l=0 projectors also show elevated ratios (0.032) but pass.

**Investigation needed:**
- Determine if this is a quadrature issue (trapezoidal integration at high q produces oscillatory artifacts) that SIMP would fix.
- Check if the decay threshold of 0.1 is too aggressive for the l=1 projector on this grid.
- If the test threshold is wrong, adjust it. If it's a real quadrature issue, document it as motivation for SIMP.

### 3. `test_vloc_comparison_with_qe` — V_local(G) mismatch with QE

**Symptom:**
- V_local(G=0): ours = +1.343 eV, QE = -1.003 eV (diff = 2.35 eV, wrong sign)
- |V_local(G=(1,0,0))|: ours = 5.468 eV, QE = 6.968 eV (diff = 1.50 eV)

**Investigation needed:**
- This is likely the same root cause as QEDX — trapezoidal quadrature + lack of erf subtraction in V_local Fourier transform.
- Determine if the QE reference values in the test are correct (extracted from pp.x Cube file FFT).
- Check whether this test should be expected to fail until SIMP+VERF land, and if so, mark it `#[ignore]` with a comment referencing those proposals.

## Deliverables

1. For each failing test, classify as: **test bug** (fix the test), **code bug** (fix the code), or **known limitation** (mark `#[ignore]` with proposal reference).
2. Write up findings in the Researcher logbook.
3. If any test reveals a genuine code bug, file a separate proposal or update an existing one.

## Approach

1. Read `qe-7.5/upflib/` sources to understand how HGH parameters map to UPF D_ij values.
2. Read our `src/pseudopotential/upf.rs` parser to trace the D_ij conversion chain.
3. For test_07, re-run with Simpson quadrature (already implemented in the test file) to see if decay improves.
4. For test_vloc, compare our `v_local_of_g()` against QE's `vloc_of_g()` in `qe-7.5/PW/src/`.

## Estimated Effort

~1-2 hours of investigation. No production code changes expected (test fixes only).
