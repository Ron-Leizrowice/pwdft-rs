---
id: SPNC
status: active
priority: high
complexity: small
risk: low
depends_on: [SPXC]
blocks: []
---

# SPNC: Use Per-Spin Density Diff for nspin=2 Convergence

## Problem

After SPXC landed (PR #12, completed 2026-04-17), the Harris-Foulkes gap on Fe
BCC fixed-mag dropped from |E_HF - E_KS| ≈ 22.2 eV to ≈ 13.0 eV. SPXC fixed
the input/output inconsistency in the XC energy terms but left a residual
~13 eV gap that is not explained by O(delta_rho²).

The SPXC implementor (2026-04-17 logbook entry) identified the cause:
`run_scf_spin` in `src/scf/mod.rs` declares convergence based **only on the
total density difference** `||rho_total_new - rho_total_r||`, computed at
line 639:

```rust
let rho_total_new: Vec<f64> = rho_up_sym.iter().zip(rho_down_sym.iter())
    .map(|(&u, &d)| u + d).collect();
let delta = density_diff(&rho_total_r, &rho_total_new, ctx.omega, ctx.n_grid);
```

The spin polarization `zeta = (rho_up - rho_down)/(rho_up + rho_down)` can
differ substantially between input and output densities even when the total
densities match, because a simultaneous `+epsilon` in `rho_up` and
`-epsilon` in `rho_down` leaves `rho_total` unchanged but keeps `zeta_in ≠
zeta_out`.

Consequently, the spin-polarized XC energy density `exc[rho, zeta]` —
which is what QE and our SPXC-fixed code both integrate — remains mismatched
between input- and output-density evaluations. The Harris-Foulkes gap
converges only as O(delta_zeta), linearly, rather than O(delta_rho²),
quadratically. At `conv_threshold = 1e-6` this leaves ≈13 eV of residual
|HF-KS| on Fe.

## QE Comparison

QE's `scf_run` (`PW/src/electrons.f90`) computes the density error per spin
channel and mixes per-channel. The mixing driver `mix_rho` uses the full
2-component density (up + down for nspin=2) to compute the residual norm.
See `PW/src/mix_rho.f90:234-260`. In particular the `dr2` (density-squared
residual) accumulator is summed over ALL nspin components, which for
nspin=2 means the total density diff plus the magnetization diff — strictly
stronger than total-only.

Our convergence check should match this in spirit: terminate only when
both spin channels have stopped moving.

## Proposed Metric

Replace the total-density-only delta with a per-spin max:

```rust
let delta_up = density_diff(&rho_up_r, &rho_up_sym, ctx.omega, ctx.n_grid);
let delta_down = density_diff(&rho_down_r, &rho_down_sym, ctx.omega, ctx.n_grid);
let delta = delta_up.max(delta_down);
```

### Rationale for per-spin max (not L2 sum, not total)

- **Strictly stronger than total.** If both channels converge individually,
  the total also converges (triangle inequality). The converse is false —
  a spin-flip fluctuation `+epsilon/-epsilon` is invisible to the total but
  present in the max.
- **max() matches our existing single-scalar threshold.** We already have
  one knob `conv_threshold`. Replacing the scalar with `max(|delta_up|,
  |delta_down|)` keeps the knob semantics identical — users calibrating
  `conv_threshold = 1e-6` continue to get "no channel has moved more than
  1e-6 in L2 norm".
- **Not the L2 sum** (sqrt(delta_up² + delta_down²)) because that would
  shift the effective threshold by up to √2, breaking backwards
  compatibility for well-tuned tests. max() preserves calibration.
- **Not the total-plus-magnetization formulation** that QE uses, because
  that requires picking a relative weight between rho_total residual and
  magnetization residual. The per-spin max implicitly weights them equally
  and avoids introducing a new parameter.

## Verification Plan

1. **Unit correctness:** `cargo test` must remain green with the new
   metric. Existing tests all use nspin=1 or nspin=2 with reasonably
   symmetric starting conditions; convergence should still be reachable.

2. **Regression test tightening (`tests/spin_polarization.rs::
   test_fe_spin_xc_consistency_regression`):** with both SPXC and SPNC in
   place, |E_HF - E_KS| on Fe BCC fixed-mag (4x4x4 k, 15 Ry, starting_mag=
   0.5, conv=1e-6) should collapse to the true O(delta_rho²) floor. The
   SPXC-era baseline is 13.03 eV. Empirically measure the post-SPNC value
   and set the assertion threshold to 2× the empirical number (with
   sensible headroom). Expected post-SPNC: O(1e-3) to O(1e-1) eV.

3. **Iteration count impact:** expect +0% to +30% iterations on Fe
   fixed-mag because the converge criterion is strictly stricter. If the
   iteration count doubles, that is a signal we have a separate
   convergence issue elsewhere and should revisit.

4. **Fe ferromagnetic (free-moment) test:** `test_fe_ferromagnetic_fixed_moment`
   currently does not converge in 100 iterations and is gated on `Ok`.
   SPNC will not change this — noted for tracking only.

## Implementation Plan

1. In `src/scf/mod.rs::run_scf_spin`, replace the single `density_diff`
   call at line 639 with per-spin max computation. Leave
   `rho_total_new` construction in place (still needed for
   `rho_total_new_g`, Hartree recompute, and energy integrals).
2. Add a short comment citing SPNC proposal to explain the per-spin
   rationale.
3. Update `tests/spin_polarization.rs::test_fe_spin_xc_consistency_regression`
   to assert the new tighter threshold. Update the comment block to
   reference both SPXC (residual artifact diagnosis) and SPNC (resolution).
4. Run full test suite; log iteration-count and |HF-KS| deltas.

## Expected Numerical Impact

**Fe BCC fixed-mag (SPXC baseline):** |HF-KS| = 13.03 eV, 244 iters.

### Empirical Result

Running the same Fe BCC fixed-mag=2 setup with per-spin criterion reveals
that the previous "convergence" was a false positive. Both spin channels
enter a steady-state limit cycle at `Δρ_up = Δρ_down ≈ 0.254` from iter
~5 onwards, with `Δρ_total ≈ 1e-7` (a +ε/−ε cancellation between channels).
The SCF now correctly reports `ConvergenceFailure` after 300 iterations
rather than declaring victory at a spurious state with 13 eV of |HF-KS|
spurious gap.

**Root cause:** the fixed-magnetization=2 constraint forces a state that is
not a stable SCF fixed point for this PP (the LDA ground state of Fe BCC
on `pseudopotentials/nc/lda/Fe.upf` is non-magnetic). The per-channel
Anderson mixers, operating independently, cannot co-ordinate the transfer
of charge density between channels that would be needed to close the gap.
Fixing this requires a better mixer or a different PP — out of SPNC scope.

**Verification on a well-behaved system (Si nspin=2, no constraint,
starting_mag=0.2):** SCF converges in 23 iterations with
**|HF-KS| = 7.19e-7 eV** at conv_threshold=1e-6. Magnetization relaxes to
0 as expected (Si is non-magnetic). This demonstrates true O(Δρ²)
quadratic convergence of the Harris-Foulkes gap with both SPXC and SPNC
in place — exactly the convergence floor we expected.

**Physics:** none changed. This is a pure convergence-criterion tightening.
Converged values (when reached) are identical under either metric — only
the point at which we stop shifts. The new metric correctly refuses to
claim convergence for states that were never actually converged.

## References

- QE 7.5: `PW/src/mix_rho.f90` lines 234-260 (multi-component residual)
- `.claude/logbooks/core-engineer.md` 2026-04-17 SPXC entry
- `proposals/completed/SPXC-spin-xc-consistency.md`

## Origin

Identified by the SPXC implementing agent during post-implementation
analysis (2026-04-17). The residual 13 eV |HF-KS| gap on Fe fixed-mag
after SPXC could not be reduced by tighter `conv_threshold` because the
convergence check was blind to spin channel drift.
