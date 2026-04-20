---
id: SPXC
status: active
priority: medium
complexity: small
risk: medium
depends_on: []
blocks: []
---

# SPXC: Fix Spin-Polarized E_xc Density Consistency

## Problem

In `run_scf_spin` (src/scf/mod.rs), the Kohn-Sham energy E_KS mixes input and
output quantities in the XC energy computation. Specifically:

- `exc_r` (XC energy density per electron) is computed at line 525 from the
  **INPUT** spin densities `rho_up_r`, `rho_down_r` plus half-core.
- `rho_xc_total` (total density for E_xc integration) is computed at line 642
  from the **OUTPUT** total density `rho_total_new` plus full core.
- `vxc_up_r`, `vxc_down_r` (XC potentials) are from the **INPUT** spin densities.
- `rho_up_sym`, `rho_down_sym` (for E_vxc integration, line 650) are **OUTPUT** densities.

The E_KS XC correction is:

    E_xc_corrected = sum(rho_xc_total * exc_r) * dvol - sum(rho_out_spin * vxc_in_spin) * dvol

This uses OUTPUT total density with INPUT exc_r and OUTPUT spin densities with
INPUT vxc potentials. The exc_r depends on spin polarization
zeta = (rho_up - rho_down)/(rho_up + rho_down), so using it with a different
total density is not physically meaningful.

The non-spin `run_scf` path does NOT have this bug: at lines 347-348 it
recomputes `(exc_r, vxc_r_energy) = xc::lda_xc_grid(&rho_new_for_xc)` from
the OUTPUT density before computing E_KS.

The Harris-Foulkes energy computation (lines 670-683) is correct: all quantities
(rho_xc_total_in, exc_r, vxc, rho_in) are consistently from the INPUT density.

## QE Comparison

QE (PW/src/electrons.f90) uses a different but equivalent formulation:

1. `v_of_rho(rhoin, ...)` computes `etxc`, `vtxc`, and `V` all from the same
   MIXED/INPUT density `rhoin` (lines 937-938).
2. `deband = delta_e()` computes `-int(rho_out * V[rho_in])` (line 849).
3. Total energy (line 1133):
   `etot = eband + (etxc - etxcc) + ewld + ehart + deband + demet + descf`
4. `descf = delta_escf()` is a first-order correction
   `-int((rhoin - rho_out) * V[rho_in])` that accounts for the difference
   between using `rhoin` vs `rho_out`.

Key point: **QE never recomputes exc from the output density.** Instead, it uses
`etxc` from the input/mixed density and adds `descf` as a first-order
variational correction. At self-consistency, `descf -> 0` because
`rhoin = rho_out`.

QE's E_xc integral in `v_xc()` (v_of_rho.f90, line 539 for nspin=2) is:
`etxc += e2 * (ex(ir) + ec(ir)) * rho%of_r(ir,1)` where `rho%of_r(ir,1)` is
the same (input) density used to compute `ex` and `ec`. Fully consistent.

## Recommended Fix

**Recompute XC from OUTPUT spin densities for E_KS**, matching the non-spin path.

After line 641 (`density_r_to_g`), add:

```rust
// Recompute XC from OUTPUT spin densities for E_KS energy
let rho_up_xc_out = add_core_density(&rho_up_sym, &rho_core_half);
let rho_down_xc_out = add_core_density(&rho_down_sym, &rho_core_half);
let (exc_r_out, vxc_up_r_out, vxc_down_r_out) =
    xc::lda_xc_spin_grid(&rho_up_xc_out, &rho_down_xc_out);
```

Then use `exc_r_out`, `vxc_up_r_out`, `vxc_down_r_out` with
`rho_xc_total`/`rho_up_sym`/`rho_down_sym` for E_KS, keeping the existing
INPUT-based quantities for E_HF.

This adds one extra `lda_xc_spin_grid` call per iteration (grid-level, no FFT),
costing about the same as the existing call. Acceptable for correctness.

Alternative (QE-like): keep everything from the input density and add a
`descf`-like correction term. This avoids the extra XC evaluation but adds
code complexity. The simple "recompute from output" approach is preferred
since the non-spin path already does this.

## Expected Numerical Impact

**Near convergence (delta_rho ~ 1e-6):** The error is O(delta_rho), roughly:

    delta_E_xc ~ integral((rho_out - rho_in) * exc[rho_in]) dr ~ O(delta_rho * omega)

For Fe BCC with omega ~ 22 Ang^3 and delta_rho = 1e-6, this is sub-meV.
The final converged energy is not significantly affected.

**Effect on |E_HF - E_KS| convergence:**

The Harris-Foulkes gap |E_HF - E_KS| should shrink as O(delta_rho^2) at
self-consistency (this is the defining property of HF). With the bug, E_KS has
an O(delta_rho) error that does not cancel, so the gap converges only as
O(delta_rho) -- linearly, not quadratically. This spoils the HF energy's
value as a convergence quality indicator for spin-polarized calculations.

This is the primary practical impact: the HRFK convergence monitoring will
log anomalously large |E_HF - E_KS| values for spin-polarized systems,
potentially triggering false warnings.

**Summary:** The bug does not corrupt converged energies but degrades the
quadratic convergence of |E_HF - E_KS|, which is the whole point of computing
Harris-Foulkes. Fix is small and should be done.

## Implementation Plan

1. After computing `rho_total_new` and `rho_total_new_g`, recompute spin XC
   from OUTPUT spin densities (rho_up_sym, rho_down_sym) with half-core.
2. Use the OUTPUT XC quantities (exc_r_out, vxc_up_r_out, vxc_down_r_out)
   for the E_KS computation.
3. Keep the existing INPUT quantities (exc_r, vxc_up_r, vxc_down_r) for E_HF.
4. Verify: run Fe BCC spin-polarized and confirm |E_HF - E_KS| converges
   quadratically (as O(delta_rho^2)) instead of linearly.

## References

- QE 7.5: PW/src/electrons.f90 lines 849, 937-938, 1133
- QE 7.5: PW/src/v_of_rho.f90 lines 529-548 (spin v_xc)
- Harris, Phys. Rev. B 31, 1770 (1985) -- HF functional stationarity
- Foulkes & Haydock, Phys. Rev. B 39, 12520 (1989)

## Origin

Discovered by Core Engineer during HRFK implementation (2026-04-16).
Logged in `.claude/logbooks/core-engineer.md`.

## 2026-04-17 — Attempt 1: Partial

A Core Engineer agent drafted the fix (archived at `/tmp/pwdft-rescue/SPXC-attempt.diff`, ~20-line change in `src/scf/mod.rs`). The implementation matches the "Recommended Fix" section above exactly: recomputes `(exc_r_out, vxc_up_r_out, vxc_down_r_out)` from the OUTPUT spin densities after `density_r_to_g`, uses the OUTPUT triple in the E_KS integration, keeps INPUT-derived quantities for E_HF.

A new regression test was started at `tests/spin_polarization.rs` (archived at `/tmp/pwdft-rescue/SPXC-tests.diff`, 52-55 lines). The test exercises BCC Fe at nspin=2 and asserts |E_HF - E_KS| decreases quadratically with decreasing `conv_threshold`.

### Issue encountered

With `conv_threshold = 1e-7` the Fe BCC test failed to converge in 120 iterations. The existing `test_fe_ferromagnetic_fixed_moment` uses `conv_threshold = 1e-6` and does converge. It is unclear whether this is:

- A genuine convergence issue introduced by the fix (unlikely — the change only affects energy accounting post-density-update)
- An artifact of the tighter threshold (Anderson mixing's default history depth may be insufficient)
- Test-scaffolding noise (the agent reported SCF bailing out at ~12 ms per iteration, which suggests an early exit not a real iteration)

### Recommended Next Step

1. Re-implement the fix (the ~20-line change is straightforward and was correct on inspection).
2. Write the regression test with `conv_threshold = 1e-6` matching the existing Fe test; verify |E_HF - E_KS| improves at least one order of magnitude vs baseline at that threshold.
3. Defer tight-threshold (1e-7) convergence investigation to a separate proposal if mixing history depth turns out to be the issue.
