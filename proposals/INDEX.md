# Proposal Index

Proposals use 4-letter IDs (e.g., `SIMP`) to avoid numbering conflicts when multiple agents work concurrently. Each proposal file is named `XXXX-slug.md` with YAML frontmatter containing metadata.

## Active

### Critical — Physics Correctness

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| VGCMP | V_local(G) & KB projector cross-check vs QE (Phase 1 done; Phase 2-4 pending) | medium | low | VERF | QEVL |

**2026-04-17:** VERF landed (cosmetic but adopted). **VGCMP Phase 1 result: V_local(G) is CORRECT** (max |Δ| = 2.78e-9 Ry vs independent Python Simpson ref, across Si's first 20 G-shells). **The 13.4 eV Si gap is NOT in V_local(G).** Investigation now focuses on KB non-local: Phase 2 will cross-check β_l(q) form factors. Primary suspects: `NonlocalPotential::F_l(q)`, UPF projector √(BOHR_TO_ANG) unit conversion, D_ij diagonalization. QEDX archived as superseded by VGCMP.

### High — Foundation & Code Quality

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| SPNC | Per-Spin Density Diff for nspin=2 Convergence | small | low | SPXC | — |
| CCMX | Coupled-Channel Mixer for nspin=2 (mix (ρ_total, m)) | medium | medium | SPNC | — |
| QEVL | QE Validation Test Suite | medium | medium | — | — |
| TAUD | Test Suite Quality Audit (silent-pass, tolerance gaps) | medium | low | — | — |

### Medium — Enhancements & Performance

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| FFTB | FFT Buffer Reuse | small | low | — | — |
| XCPR | XC and Spin Diagonalization Parallelization | small | low | — | — |
| CFGN | Expose Hardcoded Numerics as Settings | large | medium | — | — |

### Medium — Code Quality Lints (follow-up to CLIP)

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| QLNT | Tier-1 Quality Lints (follow-up to CLIP) | small | low | — | — |
| FMAD | Fused Multiply-Add via `suboptimal_flops` | medium | medium | — | — |

### Low / Deferred

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| CLSS | `cast_lossless` + Doc Hygiene | small | low | — | — |
| MUST | `must_use_candidate` Annotations | medium | low | — | — |
| CAST | Numeric Cast Safety Audit | medium | medium | — | — |
| WFRX | Wavefunction Reuse Between SCF Iterations | medium | low | — | DVSN |
| DVSN | Iterative Eigensolver (Davidson / LOBPCG) | large | medium | — | SPRS |
| HD5I | HDF5 Restart and Structured Output | large | medium | — | — |
| SPRS | Sparse Matrix Support | large | medium | DVSN | — |
| CUCL | CubeCL GPU Kernels | large | high | — | — |
| SYKP | Audit Si 4×4×4 IBZ reduction (documented: convention mismatch, no code change) | small | low | — | — |

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
| SIMP | Simpson's Rule for Radial Integrals |
| HRFK | Harris-Foulkes Energy |
| SDED | Deduplicate Settings Enums |
| BROY | Broyden Mixing (core algorithm; adaptive beta deferred) |
| ERRH | Error Handling Cleanup |
| VERF | V_local erf Coulomb Subtraction (landed as QE convention; regression test only) |
| SPXC | Fix Spin-Polarized E_xc Density Consistency (|HF-KS| 22→13 eV on Fe; residual from nspin=2 convergence — see follow-up) |
| QEDX | Systematic Energy Discrepancy vs QE (superseded by VGCMP) |

## Notes

- **QEDX** (archived 2026-04-17): superseded by VGCMP. SIMP closed the Fe gap; VERF cosmetic; VGCMP Phases 1+2 ruled out V_local(G) and β_l(q). Remaining Si gap continues under VGCMP Phase 3+4.
- **VERF** (done): landed the QE erf-subtraction convention. Numerically equivalent to bare-Coulomb on our current log mesh — the change is adopted for alignment with QE and robustness for future high-Z PPs, not as a fix. Added `tests/vloc_erf_consistency.rs` to pin the equivalence.
- **VGCMP** (new): active investigation into the Si 13.4 eV gap via form-factor-level cross-check against QE (V_local(G), β_l(q), D_ij).
- **SIMP** result: Fe BCC validation passes. Si still 13.4 eV off — VGCMP will determine where.
- **SPXC** found during HRFK: spin-polarized E_KS mixes input exc_r with output rho_xc_total.
- **CCMX** (new): opened during SPNC post-mortem. Independent Anderson mixers on `(ρ↑, ρ↓)` can't converge Fe fixed-mag=2 (limit cycle). Fix is to mix `(ρ_total, m)` instead, matching QE's `rhoz_or_updw` basis change (`qe-7.5/PW/src/sum_band.f90:307`).
- **CFGN** all dependencies satisfied (DDUP + SIMP done).
- **BROY** landed core algorithm only; adaptive-beta and periodic Pulay deferred to follow-ups.
- **HD5I** references deleted `src/input.rs` — update to YAML Settings when implementing.
