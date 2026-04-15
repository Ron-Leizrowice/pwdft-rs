# Proposal 25: Investigate Fe BCC 210 eV Energy Discrepancy

## Problem

BCC Fe (non-magnetic, NC PP, 15 Ry cutoff) produces a total energy of -391.16 eV, while QE gives -600.93 eV — a 210 eV discrepancy. The Gamma eigenvalues are qualitatively wrong:

| Band | QE (eV) | pwdft-rs (eV) |
|------|---------|---------------|
| 1 | 5.16 | 18.09 |
| 2 | 26.26 | 18.09 |
| 3 | 26.26 | 19.49 |
| 4 | 27.15 | 31.49 |
| 5 | 27.15 | 31.49 |
| 6 | 27.15 | 33.22 |

## What works

- PP parsing is correct: D_ij = [0.152, 0, 0, -135.85] eV, projectors at l=1 (p) and l=2 (d)
- Spherical Bessel j_l and Legendre P_l work for arbitrary l (verified in unit tests)
- Si (l=0,1 projectors) works to within multi-minimum tolerance

## Suspected causes (in order of likelihood)

1. **V_local Bessel transform for Fe PP**: The V_local radial function may have different characteristics than Si's HGH PP. The Coulomb subtraction in `v_local_of_g` assumes a specific form — verify the integral converges for Fe.

2. **Structure factor convention for BCC**: BCC has 1 atom at origin. The structure factor S(G) = exp(-iG·τ) = 1 for all G when τ=0. This is trivial but verify it's not introducing a factor-of-2 error (the primitive BCC cell has 1 atom, not 2).

3. **Non-local potential angular terms**: With l=2 (d-projectors), the angular factor is (2l+1)/(4π) × P_2(cos θ) = 5/(4π) × (3cos²θ-1)/2. The P_2 recurrence was just added — verify against explicit formula.

4. **Ewald energy for BCC**: Single-atom primitive cell has trivial Ewald sum. Verify the Z_val=8 for Fe is used correctly.

5. **Basis set for BCC**: The G-vector enumeration for the BCC reciprocal lattice (FCC in reciprocal space) may have issues at the BZ boundary.

## Investigation plan

1. Compare V_local(G=0) between QE and pwdft-rs for Fe
2. Compare eigenvalues of kinetic-only Hamiltonian (should match exactly)
3. Add V_local, compare eigenvalues
4. Add V_NL, compare eigenvalues — isolate which component introduces the error
5. Compare Ewald energy
