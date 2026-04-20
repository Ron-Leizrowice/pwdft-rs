# Radial Integration

## Current Implementation

All radial integrals use a simple sum:

```rust
for i in 0..n {
    integral += f(r_i) * rab[i];
}
```

This is O(h^2) trapezoidal-like accuracy. Four sites use this pattern:

| Location | Integral |
|----------|----------|
| `src/pseudopotential/mod.rs:104-132` | V_local(G) form factor |
| `src/potential/nonlocal.rs:224-238` | Beta projector form factors |
| `src/scf/potentials.rs:83-97` | NLCC core density Bessel transform |
| `src/scf/initial_density.rs:162-173` | Atomic density Bessel transform |

## QE Convention: Simpson's Rule

Quantum ESPRESSO uses Simpson's rule throughout (`qe-7.5/upflib/simpsn.f90`),
which is O(h^4) — roughly 100-10000x more accurate on the same grid.

QE Simpson weights: `1/3, 4/3, 2/3, 4/3, ..., 4/3, 1/3` for odd mesh,
with a boundary correction for even mesh.

## V_local Singularity

Our V_local integrand subtracts the full Coulomb tail:

```text
v_short(r) = V_local(r) + Z*e^2/r
```

Since pseudopotentials are smooth at the origin (V_local(0) is finite),
this sum diverges as `Z*e^2/r` at r -> 0. The integral converges (damped by
r^2 * sinc(Gr)), but the integrand has a very large spike near r=0.

QE (`qe-7.5/upflib/vloc_mod.f90:138`) uses erf subtraction instead:

```text
v_short(r) = V_local(r) + Z*e^2*erf(r)/r
```

Since `erf(r)/r -> 2/sqrt(pi)` as r -> 0, the integrand is bounded everywhere.
QE then adds back `4*pi*Z*e^2*exp(-G^2/4)/(Omega*G^2)` analytically (the FT
of `erfc(r)/r`).

Both approaches give the same V_local(G) mathematically. The difference is
purely numerical: our integrand is singular, QE's is smooth.

## Impact on Results

The combination of lower-order quadrature and a singular integrand is the
identified root cause of the 13-45 eV energy discrepancy vs QE documented
in Proposal 30. The error is G-dependent (different |G| sample different
parts of the radial integrand), which explains the broken eigenvalue
degeneracies at high-symmetry k-points.

The per-electron error scales with Z_valence:

- Si (Z=4): 1.66 eV/electron
- Fe (Z=16): 2.84 eV/electron

This is consistent with the near-origin singularity growing with Z.

## Fix Plan

1. **Proposal 38:** Implement Simpson's rule for all 4 radial integral sites.
   Port QE's `simpsn.f90` (~20 lines of logic). Highest impact.

2. **Proposal 39:** Adopt erf subtraction for V_local (G != 0 branch only).
   Makes the integrand bounded everywhere. Only needed if Proposal 38
   alone is insufficient.
