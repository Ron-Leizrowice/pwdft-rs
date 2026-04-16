# Documentation Index

Reference documents for the physics and mathematics in pwdft-rs.
Each document covers one component: the formula, our convention, code location,
audit status, and known issues.

## Core Physics

- [Total Energy](total-energy.md) — Kohn-Sham total energy, band energy double-counting, NLCC modification
- [Potentials](potentials.md) — Local pseudopotential (Bessel transform), Hartree, XC (Slater exchange + PZ correlation)
- [Non-Local Potential](nonlocal.md) — Kleinman-Bylander separable form, projector form factors, angular sums
- [Ewald Summation](ewald.md) — Ion-ion energy: reciprocal, real-space, self-energy, background terms
- [Electron Density](density.md) — Density from wavefunctions, superposition of atomic densities (SAD), normalization
- [SCF Convergence](scf-convergence.md) — Anderson/Pulay DIIS mixing, Kerker preconditioning, convergence criteria
- [Smearing and Entropy](smearing.md) — Fermi-Dirac, Gaussian, Methfessel-Paxton, Cold (Marzari-Vanderbilt); entropy; sigma to 0

## Infrastructure

- [Basis Set and FFT](basis-and-fft.md) — Plane-wave basis, G-vectors, energy cutoff, FFT convention and normalization
- [Symmetry](symmetry.md) — Space group detection, k-point IBZ reduction, density symmetrization
- [Units and Constants](units.md) — eV/Angstrom convention, physical constants, UPF unit conversion chain
- [Radial Integration](radial-integration.md) — Quadrature methods: current trapezoidal vs QE Simpson; known accuracy issue

## Reference

- [Pitfalls](pitfalls.md) — Lessons learned: NLCC omission, V_local(G=0), FFT normalization, radial quadrature
- [Key References](references.md) — Canonical papers and textbooks cited throughout

## How to Use These Docs

When auditing or modifying a component, read the relevant doc first to understand:
1. The mathematical formula being implemented
2. Our convention choices (units, normalization, sign)
3. The code location (file and line range)
4. Audit status: what has been verified against QE and what hasn't
5. Known issues or pitfalls specific to that component
