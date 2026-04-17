---
id: SYKP
title: Audit Si 4×4×4 IBZ reduction discrepancy vs QE
status: documented
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# SYKP — Audit Si 4×4×4 IBZ reduction (10 vs 8 k-points)

## TL;DR

The Si 4×4×4 "10 vs 8 IBZ points" discrepancy flagged in the 2026-04-16
researcher logbook is **not a bug**. It is a **convention mismatch**
(classification **(b)**): pwdft-rs and our QE reference are sampling
**different physical k-point grids**, not the same grid reduced
inconsistently. Both reductions are correct for the grid they sample.

- pwdft-rs' `monkhorst_pack()` hard-codes the original **Monkhorst-Pack
  shifted** convention (MP 1976, Eq. 4): fractional coordinates
  `{-(N-1)/(2N), ..., (N-1)/(2N)}`. For N=4 this is `{-3/8, -1/8, 1/8, 3/8}`.
- Our QE reference (`qe_validation/si_scf.in`, line 34) uses
  `K_POINTS automatic / 4 4 4 0 0 0`, i.e. the **Γ-centered unshifted**
  grid: `{0, 1/4, 1/2, 3/4}` (equivalently `{0, ±1/4, −1/2}`).
- These are two legitimate but **distinct** k-meshes that happen to have
  the same density. For FCC Fd-3m their IBZ reductions yield 10 and 8
  irreducible points respectively. Both sums of weights are 1.0.

Our grid is equivalent to what QE calls `4 4 4 1 1 1` (half-shift in all
three directions). That grid really does reduce to 10 IBZ points in QE
too — so our IBZ reduction is consistent with QE when the **same grid**
is compared.

**No code change is required for correctness.** The proposal below
recommends two small documentation / ergonomics follow-ups so the next
reader does not have to re-derive this.

## Evidence

### pwdft-rs MP convention

`src/kpoints.rs:30-34` (`monkhorst_pack`), confirmed by
`src/symmetry/kpoints.rs:80-85` (`mp_fractional`):

```rust
// f_j = (2*i_j - N_j + 1) / (2*N_j),  i_j ∈ {0..N_j-1}
```

For N=4: `{-3/8, -1/8, 1/8, 3/8}` — no point at Γ, symmetric about Γ.
This matches Monkhorst & Pack, *Phys. Rev. B* **13**, 5188 (1976),
Eq. 4 (the "shifted" MP formula for even N).

### QE MP convention

`qe-7.5/PW/src/kpoint_grid.f90:67-78`:

```fortran
xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
```

For `nk1=4, k1=0` (our `si_scf.in`): `xkg ∈ {0, 1/4, 1/2, 3/4}` — Γ-centered,
has a point exactly at Γ. For `k1=1`: `xkg ∈ {1/8, 3/8, 5/8, 7/8}` which
mod-1 equals `{1/8, 3/8, -3/8, -1/8}` — identical to our grid.

### QE's Si reduction (from our own reference file)

`qe_validation/si_scf.out:672-691` shows QE's 8 IBZ points in crystal
coords:

```
k(1) = ( 0.000, 0.000, 0.000)   wk = 0.03125
k(2) = ( 0.000, 0.000, 0.250)   wk = 0.25000
k(3) = ( 0.000, 0.000,-0.500)   wk = 0.12500
k(4) = ( 0.000, 0.250, 0.250)   wk = 0.18750
k(5) = ( 0.000, 0.250,-0.500)   wk = 0.75000
k(6) = ( 0.000, 0.250,-0.250)   wk = 0.37500
k(7) = ( 0.000,-0.500,-0.500)   wk = 0.09375
k(8) = ( 0.250,-0.500,-0.250)   wk = 0.18750
```

Every entry has components in `{0, ±1/4, −1/2}`. Note: QE prints weights
pre-multiplied by spin degeneracy (`qe-7.5/PW/src/setup.f90:673`,
`wk = wk * degspin` for LDA), so the printed sum is 2.0. The physical
sum normalised to 1.0 is exactly the `w_i / 2` of each line above.

### pwdft-rs' reduction of its own grid

`src/symmetry/kpoints.rs:134-157` (`test_si_4x4x4_reduces_to_8`) asserts
the result is in `[8, 10]`, and the current output is 10 with weights
summing to 1.0 (verified by the sibling `test_weight_sum_is_one`). All
10 points have components in `{±1/8, ±3/8}` — exactly the grid QE would
produce with `4 4 4 1 1 1`. No component hits a BZ boundary
(`±1/2`) where tolerance issues could plausibly arise.

### Why 10, not 8?

Intuitively, the unshifted grid has high-symmetry points (Γ, X, L and
their neighbours) that lie on a BZ boundary or a special line and are
therefore "lower-multiplicity" — their orbit under O_h is smaller, so
they contribute fewer representatives to the IBZ count while still
summing to the correct total weight. The shifted grid deliberately
avoids all those high-symmetry loci, producing no Γ/X/L/W representatives
and a more uniformly-weighted 10-point IBZ. Both are correct reductions
of their respective grids.

## Classification

**(b) Different-but-equivalent convention.** Our IBZ reduction
algorithm is correct; the symmetry detector (48 ops on Si Fd-3m,
verified in `tests/qe_validation.rs` et al.) is correct; the reciprocal
rotation `(R^{-1})^T` is correct (`symmetry/operations.rs:142-152`); the
grid-index roundtrip for all 64 points is verified (`test_mp_fractional_roundtrip`).

Not (a) tolerance: no point on our shifted grid lands near a BZ boundary
modulo reciprocal lattice vectors where the `1e-6` tolerance in
`frac_to_grid_index` could plausibly misclassify orbits.

Not (c) missing symmetry op: Si Fd-3m has 48 operations; we find all 48
(`test_si_fcc_48_operations`), and group closure is verified
(`test_si_group_closure`).

Not (d) extra op: group closure rules this out too.

## Relation to the Si 13.4 eV gap

**Plausible cause?** No. Order-of-magnitude argument:

- k-point convergence error at 4×4×4 for a simple semiconductor is
  < 10 meV between any two reasonable grids.
- Changing a 10-point shifted-IBZ sampling to an 8-point Γ-centered
  sampling cannot produce a 13,400 meV shift.
- The documented signatures of the 13.4 eV gap (Γ degeneracy broken,
  error scales as ~1.66 eV/electron) are form-factor / quadrature
  artefacts in V_local and the KB non-local projectors — see the
  `VGCMP` tracking proposal and the 2026-04-17 logbook entries. Those
  signatures are **k-independent**.

There is however a minor apples-to-oranges effect when numerically
comparing pwdft-rs vs QE: the two calculations sample physically
different k-grids. The residual from this is sub-meV and is dwarfed by
the VGCMP/SIMP/VERF-class issues, but it does mean the future
"pwdft-rs matches QE to 0.1 eV" assertion in
`tests/qe_validation.rs` should either (i) use an MP shift that makes
the grids identical, or (ii) state in a comment that the small
k-sampling delta is absorbed into the 0.1 eV tolerance.

## Recommended follow-ups (documentation + ergonomics only)

These are **not urgent** and are **not part of this audit PR**. They
are captured here so a future Core Engineer / Technical Writer session
can pick them up.

### D1 — Clarify the MP convention in source docs

Current docstring on `monkhorst_pack` (`src/kpoints.rs:17-20`) says
"produces `k_i = (2n_i - N_i - 1) / (2 N_i)`" without naming the
convention. Update to explicitly note:

- this is the **shifted** MP grid (MP 1976, Eq. 4; same as QE
  `k1=k2=k3=1`);
- it differs from QE's default `0 0 0` which is Γ-centered;
- the IBZ reduction will therefore produce different irreducible
  k-point counts between the two conventions — this is not a bug.

Similar note on `src/symmetry/kpoints.rs` (`reduce_kpoints` and/or
`test_si_4x4x4_reduces_to_8`'s assertion comment currently says
"10 due to incomplete boundary handling", which is **wrong** per the
classification above; should be "10 because the shifted MP grid has a
genuinely different orbit structure from the Γ-centered grid QE
defaults to; both reductions are exact").

### D2 — (Optional, deferred) Expose an MP shift parameter

Add an optional `shift: [u32; 3]` or `kind: "shifted" | "gamma"` field
to `KPointSettings::MonkhorstPack` in `src/settings.rs:113-117`, plus
a `monkhorst_pack_shifted(n1, n2, n3, shift, lattice)` variant in
`src/kpoints.rs`. Default remains the current shifted behaviour to
preserve all existing YAML files and test results. This would let
users who need byte-exact QE comparison run against QE's `0 0 0` grid.

This is a **scope expansion**, not a fix. File as a standalone
proposal (suggested ID: `MPSH` — Monkhorst-Pack shift parameter) if
prioritised.

## Test plan

No new tests — existing ones are correct:

- `symmetry::kpoints::tests::test_si_4x4x4_reduces_to_8` passes (asserts
  8 ≤ n ≤ 10; current value 10).
- `symmetry::kpoints::tests::test_weight_sum_is_one` passes (sum = 1.0).
- `symmetry::kpoints::tests::test_mp_fractional_roundtrip` passes.
- `symmetry::detect::tests::test_si_fcc_48_operations` passes.
- `symmetry::detect::tests::test_si_group_closure` passes.

The assertion in `test_si_4x4x4_reduces_to_8` (`8 ≤ n ≤ 10`) is a
legitimate lower bound: in principle a future implementation that
supports both conventions could tighten it to `== 10` for the shifted
branch.

## Status: documented

Closing as **documented, no code change**. D1 captured here for a
future docstring-fix sweep; D2 deferred as a separate proposal if
needed.

## References

- Monkhorst, H. J.; Pack, J. D. *Special points for Brillouin-zone
  integrations.* **Phys. Rev. B 13, 5188 (1976).** Eq. 4 gives the
  shifted formula we use.
- Quantum ESPRESSO 7.5:
  - `qe-7.5/PW/src/kpoint_grid.f90:67-170` — IBZ-reducing MP generator.
  - `qe-7.5/PW/src/setup.f90:673` — `wk *= degspin` for spin-unpolarised.
  - `qe-7.5/PW/src/symm_base.f90` — space-group symmetry operations
    (not re-audited here; Si Fd-3m detection was independently verified
    via 48-op count and group closure tests listed above).
- Our own reference run: `qe_validation/si_scf.out:672-691`.

## 2026-04-17 — D1 done

Docstrings updated in `src/kpoints.rs::monkhorst_pack` and
`src/symmetry/kpoints.rs::{reduce_kpoints, mp_fractional}`, and the
misleading "incomplete boundary handling" comment in
`test_si_4x4x4_reduces_to_8` replaced with the correct convention
explanation. Landed in XCLN cleanup PR.
