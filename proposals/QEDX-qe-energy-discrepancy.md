---
id: QEDX
status: active
priority: critical
complexity: large
risk: high
depends_on: [SIMP, VERF]
blocks: [QEVL]
---

# QEDX: Systematic Energy Discrepancy vs QE

## Problem

pwdft-rs SCF total energies differ from QE 7.5 by 13–45 eV depending on the material, using the same pseudopotentials (PseudoDojo ONCV LDA), same ecut (15 Ry), same k-grid (4×4×4), and same smearing (Fermi-Dirac σ=0.01 Ry). Eigenvalue degeneracies that should be exact by symmetry are broken.

### Reproduction

```bash
# Run QE references (results already in qe-7.5/runs/{si,c,fe}_pseudodojo/)
cargo test --release --test qe_validation -- --nocapture
```

### Measured discrepancies

| System | Atoms | Z_val | QE energy (eV) | Our energy (eV) | ΔE (eV) | ΔE/el (eV) | Status |
|--------|-------|-------|-----------------|------------------|---------|-------------|--------|
| Si diamond | 2 | 4 | -231.61 | -218.28 | 13.3 | 1.66 | Converges (11 iter) |
| C diamond | 2 | 4 | -312.51 | — | — | — | Does not converge |
| Fe BCC | 1 | 16 | -3059.46 | -3104.83 | 45.4 | 2.84 | Converges (16 iter) |

### Eigenvalue degeneracy breaking

In Si (FCC) at Gamma, bands 2-4 should be triply degenerate. QE gives 6.080, 6.080, 6.080 eV. We give -1.04, 3.75, 3.76 eV — completely wrong pattern and magnitude.

In Fe (BCC) at Gamma, bands 2-4 (3p semicore) should be triply degenerate. QE gives -46.50, -46.50, -46.50 eV. We give -44.86, -44.86, -38.56 eV — one state is split by 6 eV.

## Research: What's validated, what's suspect

### Root cause identified (April 2026 audit)

A comprehensive line-by-line comparison against QE 7.5 source code verified all 18 core formulas, unit conversions, and conventions as correct. The discrepancy is traced to **radial quadrature quality** — see Proposals 38 and 39 for the fix:

1. **Simple sum vs Simpson's rule (Proposal 38):** All radial integrals use `integral += f * dr` (O(h^2)), while QE uses Simpson's rule (O(h^4), ~100x more accurate). This affects V_local, beta projectors, core density, and atomic density.

2. **V_local Coulomb subtraction singularity (Proposal 39):** Our V_local integrand has a divergent `Ze^2/r` term near r=0. QE uses erf subtraction to keep the integrand bounded. Combined with lower-order quadrature, this produces G-dependent errors that break eigenvalue degeneracies.

The per-electron error scaling with Z (Si 1.66 eV/el vs Fe 2.84 eV/el) is consistent: the near-origin singularity grows with Z_valence.

### Validated components (expanded in April 2026 audit)

| Component | Test | Status |
|-----------|------|--------|
| **Ewald energy** | Fe: ours = -2337.167 eV, QE = -2337.173 eV | **Match (0.006 eV)** |
| **Basis size** | Fe Gamma: 79 PWs, matches QE exactly | **Match** |
| **V_NL Hermiticity** | Max violation = 2.2e-16 | **Correct** |
| **Free-electron eigenvalues** | Band 0 = 0, band 1 = 36.52 eV (BCC Fe) | **Correct** |
| **PZ correlation constants** | Verified against Table I of PZ 1981 | **Correct** |
| **Smearing/occupations** | Extensive unit tests, nspin=1/2 consistency | **Correct** |
| **All unit conversions** | UPF parser: r, V, beta, D_ij, rho_atom, core_charge | **Correct** |
| **XC NLCC double-counting** | E_xc uses rho_total, E_vxc uses rho_val | **Correct (matches QE)** |
| **KB form factors** | `∫ r·β(r) j_l(qr) r dr` matches QE beta_mod.f90:112-113 | **Correct** |
| **Hamiltonian assembly** | T + V_eff(G-G') via Miller index, FFT convention | **Correct** |
| **V_local(G=0) convention** | Excluded from H, added to E_total | **Correct (matches QE)** |
| **Slater exchange + PZ correlation** | All parameters, derivatives, unit chain | **Correct** |
| **Spin-polarized LSDA** | f(zeta), exchange, correlation interpolation | **Correct** |

### Original suspect components (pre-audit)

**1. V_local Fourier transform — CONFIRMED: QUADRATURE QUALITY ISSUE**

The formula in `v_local_of_g()` is correct, and the Coulomb subtraction approach is mathematically valid. The issue is purely numerical: O(h^2) quadrature with a near-singular integrand vs QE's O(h^4) Simpson with a smooth integrand.

**2. XC double-counting with NLCC — CLEARED**

Verified correct by reading both our code and QE's `v_of_rho.f90`.

**3. Projector unit conversion — CLEARED**

Confirmed matching QE: `upf%beta * besr * r` (QE beta_mod.f90:113) = `rb * jl * r` (our nonlocal.rs:235). The `4π/√Ω` prefactor and `D_ij` conventions are consistent.

**4. FFT normalization or grid mismatch — CLEARED**

FFT convention is consistent. Hamiltonian assembly via Miller index lookup is correct. The broken degeneracies are explained by G-dependent quadrature errors, not grid asymmetry.

**5. V_local(G=0) convention — CLEARED**

Matches QE exactly.

## Implementation: Fix Plan

The root cause is identified. Fix via Proposals 38 and 39:

1. **Proposal 38: Simpson's rule** — Replace O(h^2) sum with O(h^4) Simpson in all 4 radial integral sites. Highest impact, lowest risk. ~1-2 hours.
2. **Proposal 39: erf subtraction** — Adopt QE's erf/r decomposition for V_local(G!=0). Only needed if Proposal 38 alone is insufficient. ~1-2 hours.

### Remaining diagnostics (if Proposals 38+39 are insufficient)

These steps from the original investigation remain valid as fallback diagnostics:

- **Binary search (D_ij=0):** Disable V_NL and compare eigenvalues to isolate V_local vs V_NL.
- **V_local(G) point-by-point comparison:** Extract QE V_local(G) via pp.x and compare.
- **Per-component energy printout:** Log E_band, E_H, E_xc, E_vxc, E_ewald individually.

## Verification

1. Si total energy within 0.1 eV of QE (-231.61 eV)
2. C diamond converges and matches QE within 0.1 eV
3. Fe total energy within 0.1 eV of QE (-3059.46 eV)
4. All eigenvalue degeneracies exact to 1e-4 eV at Gamma
5. `cargo test --release --test qe_validation` passes with all assertions < 0.5 eV

## Estimated Effort

2-4 hours total via Proposals 38 and 39. The root cause is understood; implementation is mechanical.
