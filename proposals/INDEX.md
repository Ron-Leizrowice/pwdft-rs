# Proposal Index

Proposals use 4-letter IDs (e.g., `SIMP`) to avoid numbering conflicts when multiple agents work concurrently. Each proposal file is named `XXXX-slug.md` with YAML frontmatter containing metadata.

## Active

### Critical — Physics Correctness

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| VGC5 | Per-Component Energy Accounting (Si vs QE) — VGCMP Phase 5 | medium | low | — | — |

**2026-04-17:** VGCMP Phases 1+2+3+4 all done (PR #29). The entire pseudopotential → Hamiltonian assembly pipeline is bit-correct vs QE: V_local(G), β_l(q), D_ij, and assembled diagonal H[G,G] all clear to machine precision. **The 13.4 eV Si gap is OUTSIDE the matrix assembly.** VGC5 (Phase 5) will tabulate per-component energies side-by-side. Prime suspect: the V_local(G=0) compensating background shift in `total_energy()` (`src/scf/context.rs:93-94` zeroes `v_local_fft[0]` and stashes it separately; may not be added back). Geometry-dependent — explains why Fe (matches to 0.02 eV) and Si (off by 13.4 eV) diverge.

### High — Foundation & Code Quality

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| CCMX | Coupled-Channel Mixer for nspin=2 (mix (ρ_total, m)) | medium | medium | — | — |
| ITEV | Iterative Eigensolver via `faer::partial_self_adjoint_eigen` (supersedes DVSN) | medium | medium | — | — |
| VLQR | V_local QE Reference Data Re-extraction (TAUD PR D follow-up) | small | low | — | — |

### Medium — Enhancements & Performance

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| XCPR | XC Parallelization — Step 3 (spin-channel `rayon::join`) remaining | small | low | — | — |
| CFGN | Expose Hardcoded Numerics as Settings | large | medium | — | — |
| MXBA | Adaptive Mixing Beta (BROY follow-up — Eyert 1996 residual-monitor rule) | medium | medium | — | — |
| PRPL | Periodic Pulay Mixing (BROY follow-up — Banerjee et al. JCTC 2016) | small | low | — | — |
| NLCC | Nonlinear Core Correction Audit (docs + test coverage; no bug found) | small | low | — | — |

### Medium — Code Quality & Refactor

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| SOPT | Drop `Option<&SymmetryInfo>` from `ScfContext` | small | low | — | — |
| DBGC | `ScfResult` Debug derive + GPU-test reference-value const (TAUD nits) | small | low | — | — |

### Low / Deferred

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| CLSS | `cast_lossless` + Doc Hygiene | small | low | — | — |
| CAST | Numeric Cast Safety Audit | medium | medium | — | — |
| WFRX | Wavefunction Reuse Between SCF Iterations (re-scope as ITEV warm-start) | medium | low | — | — |
| DVSN | Iterative Eigensolver (Davidson / LOBPCG) — SUPERSEDED BY ITEV | large | medium | — | — |
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
| SIMP | Simpson's Rule for Radial Integrals |
| HRFK | Harris-Foulkes Energy |
| SDED | Deduplicate Settings Enums |
| BROY | Broyden Mixing (core algorithm; adaptive beta deferred) |
| ERRH | Error Handling Cleanup |
| VERF | V_local erf Coulomb Subtraction (landed as QE convention; regression test only) |
| SPXC | Fix Spin-Polarized E_xc Density Consistency (|HF-KS| 22→13 eV on Fe) |
| QEDX | Systematic Energy Discrepancy vs QE (superseded by VGCMP) |
| QEVL | QE Validation Test Suite (Tier 1+2) |
| MUST | `must_use_candidate` Annotations (70 of 80 sites) |
| SPNC | Per-Spin Density Diff for nspin=2 Convergence |
| SYKP | Audit Si 4×4×4 IBZ reduction (convention mismatch — docstrings clarified, no code change) |
| FFTB | FFT Buffer Reuse (23–33% speedup on `fft/scf_iter_20x_*`) |
| FMAD | Fused Multiply-Add via `suboptimal_flops` (-3 to -4% on `lda_xc_grid_*`) |
| QLNT | Tier-1 Quality Lints (11 new lints, 15 fixes; GPU spillover → QLN2) |
| TAUD | Test Suite Quality Audit (5 PRs landed; PR D uncovered VLQR follow-up) |
| QLN2 | GPU + benches lint follow-up (closes the CLIP `--features gpu` gap) |
| CIGP | Document `--features gpu` in clippy CI gate (CLAUDE.md + agent defs) |

## Notes

- **VGCMP** (Phases 1-4 done): the entire PP→H assembly pipeline is bit-correct vs QE. The 13.4 eV Si gap is OUTSIDE matrix assembly. **VGC5** is the next-step Phase 5 (per-component energy accounting). Prime suspect: V_local(G=0) compensating shift missing in `total_energy()`.
- **CCMX** (active): Independent Anderson mixers on `(ρ↑, ρ↓)` can't converge Fe fixed-mag=2 (limit cycle). Fix is to mix `(ρ_total, m)` instead, matching QE's `rhoz_or_updw` basis change.
- **TAUD** (done): all 5 PRs landed. PR D uncovered a sign-flipped V_local in the test's QE Cube reference (NOT in our Rust code — VGCMP Phase 1 already proved Rust correct). Captured as VLQR. Test re-`#[ignore]`'d with diagnostic numbers in the reason string.
- **XCPR** (active): Step 1+2 (XC grid parallelization) merged. Step 3 (spin-channel `rayon::join`) remains.
- **CFGN** all dependencies satisfied (DDUP + SIMP done).
- **BROY** landed core algorithm only; adaptive-beta (MXBA) and periodic Pulay (PRPL) follow-ups are now open proposals.
- **NLCC** audit (2026-04-17): all code paths verified correct against QE `v_of_rho.f90`. Hartree excludes core, electron count excludes core, LSDA splits core/2 per spin, XC uses val+core with val-only double-counting. No bug. Proposal is documentation + integration test against an NLCC element (Fe).
- **HD5I** references deleted `src/input.rs` — update to YAML Settings when implementing.
- **SOPT** (new 2026-04-17): refactor `ScfContext.symmetry: Option<&SymmetryInfo>` to always-present (identity fallback). Every material has identity group; Option encodes a setting not a structural fact.
- **ITEV** (new 2026-04-17, Performance Engineer): Post-FFTB/FMAD profiling shows eigensolver at 85-90% of SCF user CPU (n_pw=259: 56 ms/call; n_pw=725: 836 ms/call). `faer 0.24` ships `matrix_free::eigen::partial_self_adjoint_eigen` (implicitly-restarted Arnoldi, matrix-free via `LinOp`, warm-start via `v0`). Supersedes DVSN's hand-rolled Davidson plan. Projected 2.5-4× SCF wall-time speedup at production sizes. WFRX becomes the warm-start knob.
