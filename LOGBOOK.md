# pwdft-rs Logbook

## 2026-04-16: Project State Assessment

### Overview

pwdft-rs is a plane-wave DFT solver in Rust targeting Apple Silicon (Metal GPU). The code implements norm-conserving pseudopotential LDA/LSDA with Kleinman-Bylander non-local projectors, validated against Quantum ESPRESSO 7.5.

**Codebase:** ~9,000 lines of Rust (src/), ~3,300 lines of tests. 159 unit/integration tests, zero clippy warnings (curated pedantic lints enabled). 19 completed proposals, 24 open.

### What works

| Feature | Status | Notes |
|---------|--------|-------|
| Plane-wave basis | Complete | Automatic cutoff, Miller indices, FFT mapping |
| UPF v2 parser | Complete | PseudoDojo ONCV LDA/PBE (74+72 elements) |
| Local pseudopotential | Complete | Spherical Bessel transform on PP radial grid |
| KB non-local potential | Complete | Arbitrary l, precomputed form factors, Hermitian |
| Ewald ion-ion | Complete | Matches QE to 0.006 eV (validated on Fe BCC) |
| LDA XC (Perdew-Zunger) | Complete | Slater exchange + PZ correlation, NLCC support |
| LSDA spin-polarized XC | Complete | f(zeta) interpolation, per-spin potentials |
| Hartree potential | Complete | Reciprocal space, G=0 excluded |
| SCF loop (nspin=1,2) | Complete | Anderson/Pulay mixing, Kerker preconditioning |
| Fermi-Dirac/Gaussian/MP/cold smearing | Complete | With entropy and free energy |
| Density symmetrization | Complete | Space group detection, point group ops |
| Free-electron band structure | Complete | Validated Si, C, Fe against analytic |
| faer eigensolver | Complete | Pure Rust, 1.8-2.8x faster than LAPACK |
| ndrustfft 3D FFT | Complete | Zero unsafe code |
| GPU acceleration (wgpu/Metal) | Partial | Hartree, XC, V_eff assembly in f32 |
| K-point symmetry reduction | Implemented | **Not wired into SCF path** |

### Known issues

#### Critical: Energy discrepancy vs QE (Proposal 30)

Total energy differs from QE by 13-45 eV using identical pseudopotentials and parameters. Root cause identified as **radial quadrature quality** (see Proposals 38-39):

| System | Z_val | QE (eV) | pwdft-rs (eV) | |Delta| (eV) | |Delta|/el |
|--------|-------|---------|---------|----------|-----------|
| Si diamond | 4 | -231.61 | -218.28 | 13.3 | 1.66 |
| C diamond | 4 | -312.51 | (no conv) | — | — |
| Fe BCC | 16 | -3059.46 | -3104.83 | 45.4 | 2.84 |

Eigenvalue degeneracies are broken (e.g., Si 3-fold at Gamma split by ~4 eV). This is consistent with G-dependent quadrature errors in V_local.

**Fix path:** Proposal 38 (Simpson's rule for radial integrals) + Proposal 39 (erf subtraction for V_local singularity). Estimated 2-4 hours.

#### Medium: K-point reduction not used in SCF

`symmetry/kpoints.rs` has `reduce_kpoints()` that reduces Si 4x4x4 from 64 to 8 IBZ k-points. But the SCF path always uses the full unreduced grid. This is an 8x performance penalty for cubic systems.

#### Low: C diamond does not converge

C diamond fails to converge in 80 iterations with plain mixing at 15 Ry. Likely related to the V_local quadrature issue (broken degeneracies → oscillating density).

### Performance benchmarks

**Hardware:** Apple M3 Max, 36 GB. **Date:** 2026-04-16.

All runs: ecut=15 Ry, 4x4x4 MP grid, FD smearing sigma=0.01 Ry, PseudoDojo NC LDA PPs.

| System | PWs | k-pts | Iters | pwdft-rs (wall) | QE 7.5 (wall) | QE config |
|--------|-----|-------|-------|----------------|---------------|-----------|
| Fe BCC (1 atom) | 79 | 64 vs 8 | 16 vs 8 | **0.14s** | **0.11s** | 12 MPI |
| Si diamond (2 atoms) | 283 | 64 vs 8 | 11 vs 7 | **1.4s** | **5.8s** | 12 MPI |
| C diamond (2 atoms) | 65 | 64 vs 8 | 80 (no conv) | — | **0.09s** | 12 MPI |

Notes:
- **pwdft-rs uses 64 unreduced k-points** while QE uses 8 symmetry-reduced. This makes Si appear fast but it's doing 8x more work per iteration. With k-point reduction, Si would be ~0.2s.
- Fe is fast due to tiny basis (79 PWs) — eigensolve is trivial. The 16 vs 8 iteration count is due to the quadrature-induced potential errors.
- QE timings include MPI startup overhead (~0.05s), which dominates for small systems.
- Our code uses rayon for k-point parallelism (all cores) vs QE's MPI (12 ranks).

### Architecture

```
main.rs → input.rs → scf/mod.rs (SCF loop)
                      ├── scf/context.rs (immutable precomputed state)
                      ├── scf/energy.rs (band, Hartree, XC, total)
                      ├── scf/density.rs (wavefunction → charge density)
                      ├── scf/mixing.rs (Anderson/Pulay + Kerker)
                      ├── scf/smearing.rs (occupations + entropy)
                      ├── scf/potentials.rs (V_local, core density, H assembly)
                      └── scf/grid.rs (FFT grid management)

potential/xc.rs          LDA/LSDA exchange-correlation
potential/nonlocal.rs    KB projectors and form factors
potential/hartree.rs     Coulomb potential
ewald.rs                 Ion-ion energy
eigensolver/dense.rs     faer Hermitian eigendecomposition
fft.rs                   ndrustfft 3D FFT (zero unsafe)
pseudopotential/upf.rs   UPF v2 parser
symmetry/                Space group detection + k-point reduction
gpu/                     wgpu Metal compute shaders
```

### Pseudopotential library

```
pseudopotentials/
  nc/lda/  — 74 PseudoDojo ONCV elements (H-Zr, La, Hf-Bi)
  nc/pbe/  — 72 PseudoDojo ONCV elements
  uspp/pbe/ — 41 SSSP efficiency ultrasoft (for future USPP support)
  paw/pbe/  — 12 SSSP efficiency PAW (for future PAW support)
```

### Next priorities

1. **Fix radial quadrature** (Proposals 38-39) — this is blocking all QE validation
2. **Wire k-point reduction into SCF** — free 8x speedup for symmetric systems
3. **Indicatif progress bars** (Proposal 05) — UX improvement
4. **Input validation** (Proposal 19) — catch bad inputs early
5. **XC GPU shader audit** (Proposal 23) — CPU/GPU consistency
