# Proposal Index

Proposals use 4-letter IDs (e.g., `SIMP`) to avoid numbering conflicts when multiple agents work concurrently. Each proposal file is named `XXXX-slug.md` with YAML frontmatter containing metadata.

## Active

### Critical — Physics Correctness

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| SIMP | Simpson's Rule for Radial Integrals | medium | medium | — | VERF, QEVL |
| VERF | V_local erf Coulomb Subtraction | small | medium | SIMP | QEVL |
| QEDX | Systematic Energy Discrepancy vs QE | large | high | SIMP, VERF | QEVL |

### High — Foundation & Code Quality

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| QEVL | QE Validation Test Suite | medium | medium | SIMP | — |

### Medium — Enhancements & Performance

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| SDED | Deduplicate Settings Enums (MixingModeType + OccupationType) | small | low | — | — |
| HRFK | Harris-Foulkes Energy | small | low | — | — |
| FFTB | FFT Buffer Reuse | small | low | — | — |
| XCPR | XC and Spin Diagonalization Parallelization | small | low | — | — |
| ERRH | Error Handling Cleanup (16 production unwrap/panic/expect) | medium | medium | — | — |
| BROY | Broyden Mixing and Adaptive Beta | medium | medium | — | — |
| CFGN | Expose Hardcoded Numerics as Settings | large | medium | SIMP | — |

### Low / Deferred

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| WFRX | Wavefunction Reuse Between SCF Iterations | medium | low | — | DVSN |
| DVSN | Iterative Eigensolver (Davidson / LOBPCG) | large | medium | — | SPRS |
| HD5I | HDF5 Restart and Structured Output | large | medium | — | — |
| SPRS | Sparse Matrix Support | large | medium | DVSN | — |
| CUCL | CubeCL GPU Kernels | large | high | — | — |

## Reference Documents

- `psuedopotentials-sota.md` — Survey of pseudopotential formalisms, libraries, and formats

## Completed

| ID | Title |
|----|-------|
| DLTB | Fix Convergence Failure Delta Bug |
| CBRT | Replace powf(1/3) with cbrt() |
| CNST | Consolidate Physical Constants |
| FAER | faer eigensolver (legacy #01) |
| FFTW | FFTW/ndrustfft integration (legacy #02) |
| NDAR | ndarray grid ops (legacy #03) |
| SPFN | Special functions / puruspe (legacy #04) |
| PROG | Indicatif progress bars (legacy #05) |
| KRKR | Kerker preconditioning (legacy #09) |
| ECON | Energy convergence tracking (legacy #10) |
| ENTR | Entropy and free energy (legacy #11) |
| SMER | Smearing schemes (legacy #12) |
| CNLP | Cache non-local potential (legacy #17) |
| PSP8 | PSP8 D_ij fix (legacy #18) |
| IVAL | Input validation (legacy #19) |
| ELEM | Elements crate (legacy #21) |
| NUMC | Numeric constants cleanup (legacy #22) |
| XCGA | XC GPU shader audit (legacy #23) |
| SDWF | Symmetry density wrapping fix (legacy #24) |
| EWSF | Ewald structure factor cleanup (legacy #25) |
| FEDB | Fe d-electron bug (legacy #25) |
| MDOC | Math documentation (legacy #26) |
| CLIP | Pedantic clippy lints (legacy #27) |
| GITC | Git repo cleanup (legacy #28) |
| MSTR | Math docstring gaps (legacy #29) |
| YAML | YAML input migration (legacy #32) |
| KBTF | Investigate KB Projector Test Failures |
| CLEN | Minor Code Quality Cleanups |
| DDUP | SCF Code Deduplication |

## Notes

- **QEDX** is a tracking proposal — actual fixes are SIMP + VERF. Archive when those land.
- **ERRH** counts refreshed: 11 unwrap + 2 panic + 3 expect in production code (was 25).
- **SDED** scope reduced: SmearingType dedup done, MixingModeType + OccupationType remain.
- **CFGN** dependency on CNST satisfied (archived). DDUP done, SIMP still pending.
- **HD5I** references deleted `src/input.rs` — update to YAML Settings when implementing.
- **SIMP** is now unblocked (KBTF completed). Critical path item.
