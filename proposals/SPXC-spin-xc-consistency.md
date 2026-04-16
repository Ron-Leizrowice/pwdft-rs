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

Found during HRFK implementation: in the spin-polarized SCF loop (`run_scf_spin`), the E_xc integral uses `exc_r` computed from **input** spin densities but integrates it against `rho_xc_total` which is the **output** total density + core charge. This mixes input and output quantities in the Kohn-Sham energy, breaking the variational property.

The non-spin path (`run_scf`) does not have this issue.

## Investigation Needed

1. Trace `exc_r` and `rho_xc_total` through the spin SCF loop to confirm the inconsistency
2. Determine whether `exc_r` should be recomputed from output densities, or `rho_xc_total` should use input densities for E_KS
3. Check how QE handles this in the spin-polarized case (`qe-7.5/PW/src/v_of_rho.f90`)
4. Quantify the numerical impact (may be negligible near convergence but matters for E_KS vs E_HF comparison)

## Origin

Discovered by Core Engineer during HRFK implementation (2026-04-16). Logged in `.claude/logbooks/core-engineer.md`.
