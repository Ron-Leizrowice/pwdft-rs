# Proposal Index

Proposals use 4-letter IDs (e.g., `SIMP`) to avoid numbering conflicts when multiple agents work concurrently. Each proposal file is named `XXXX-slug.md` with YAML frontmatter containing metadata.

## Active

### Critical — Physics Correctness

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| PCFX | Density symmetrization in G-space (non-symmorphic τ fix) | medium | medium | — | — |

**2026-04-18:** NCFX landed (closed the 13.4 eV Si gap to 0.26 eV — 52× reduction). Remaining critical item is PCFX: real-space density symmetrization uses `nint`-based rotation on 18³ grid and can't represent Fd-3m's τ=(1/4,1/4,1/4); fix is G-space phase-factor symmetrization (QE-style). Residual on Si traced to PCFX + MP-shifted grid (SYKP territory).

### High — Foundation & Code Quality

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| CCMX | Coupled-Channel Mixer for nspin=2 (mix (ρ_total, m)) | medium | medium | — | — |
| ITEV | Iterative Eigensolver via `faer::partial_self_adjoint_eigen` (supersedes DVSN) | medium | medium | — | — |

### Medium — Enhancements & Performance

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| CFGN | Expose Hardcoded Numerics as Settings | large | medium | — | — |
| MXBA | Adaptive Mixing Beta (BROY follow-up — Eyert 1996 residual-monitor rule) | medium | medium | — | — |
| PRPL | Periodic Pulay Mixing (BROY follow-up — Banerjee et al. JCTC 2016) | small | low | — | — |
| NLCC | Nonlinear Core Correction Audit (docs + test coverage; no bug found) | small | low | — | — |

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
| SOPT | Drop `Option<&SymmetryInfo>` from `ScfContext` (incl. P1+TR regression test) |
| VLQR | Retire Cube-based `test_vloc_comparison_with_qe` (superseded by VGCMP Phase 1) |
| XCPR | XC + Spin Diagonalization Parallelization (3 steps; Step 3 = 1.12–1.24×, ceiling at 2× post-ITEV) |
| VGC5 | Per-Component Energy Accounting (VGCMP Phase 5) — localized 13.4 eV Si gap to E_xc / NLCC |
| DBGC | `ScfResult` Debug derive + GPU reference consts (TAUD nits) |
| PCRS | Per-Component Energy Residual Investigation (traced 1.2 eV to non-symmorphic τ symmetrization → PCFX) |
| NCFX | NLCC Core-Density Unit and Radial-Weight Fix (closes the 13.4 eV Si gap: 13.43 → 0.26 eV) |

## Notes

- **VGCMP + VGC5 + NCFX** (all done): the entire PP→H assembly pipeline is bit-correct vs QE, and NCFX (NLCC unit conversion + missing r²·4π in Bessel FT) closed the 13.4 eV Si gap to 0.26 eV. Residual attributed to Monkhorst-Pack shifted-vs-Γ-centered grid convention (SYKP). Prime suspect (V_local(G=0) shift) ruled out — already present.
- **CCMX** (active): Independent Anderson mixers on `(ρ↑, ρ↓)` can't converge Fe fixed-mag=2 (limit cycle). Fix is to mix `(ρ_total, m)` instead, matching QE's `rhoz_or_updw` basis change.
- **TAUD** (done): all 5 PRs landed. PR D uncovered a sign-flipped V_local in the test's QE Cube reference (NOT in our Rust code — VGCMP Phase 1 already proved Rust correct). Captured as VLQR. Test re-`#[ignore]`'d with diagnostic numbers in the reason string.
- **XCPR** (done 2026-04-17): Steps 1+2 (XC grid) + Step 3 (spin-channel `rayon::join`) all landed. Step 3 speedup 1.12–1.24× (faer's internal gemm already saturates 8 cores during eigensolve); ceiling ~2× after ITEV drops per-k eigensolve cost.
- **CFGN** all dependencies satisfied (DDUP + SIMP done).
- **BROY** landed core algorithm only; adaptive-beta (MXBA) and periodic Pulay (PRPL) follow-ups are now open proposals.
- **NLCC** audit (2026-04-17): code-path audit against QE `v_of_rho.f90` verified Hartree/electron-count/LSDA-split/XC-double-counting invariants. **NCFX** (now landed) fixed the underlying storage-convention bug (`PP_NLCC` bare ρ_core in e/Bohr³, not 4πr²·ρ in e/Bohr) + missing r²·4π in the Bessel FT. The NLCC audit's Part A/B/C (unit tests, Fe integration, docs) is now unblocked.
- **HD5I** references deleted `src/input.rs` — update to YAML Settings when implementing.
- **SOPT** (done 2026-04-17): refactored `ScfContext.symmetry: Option<&SymmetryInfo>` to always-present with identity-only fallback. Review caught a P1+TR regression in the initial `is_trivial()` short-circuit in main.rs; fixed by dropping the branch and tightening the predicate. Regression test pins the behavior.
- **ITEV** (new 2026-04-17, Performance Engineer): Post-FFTB/FMAD profiling shows eigensolver at 85-90% of SCF user CPU (n_pw=259: 56 ms/call; n_pw=725: 836 ms/call). `faer 0.24` ships `matrix_free::eigen::partial_self_adjoint_eigen` (implicitly-restarted Arnoldi, matrix-free via `LinOp`, warm-start via `v0`). Supersedes DVSN's hand-rolled Davidson plan. Projected 2.5-4× SCF wall-time speedup at production sizes. WFRX becomes the warm-start knob.
- **PCRS → PCFX** (2026-04-17, Researcher): per-component identity residual of 1.204 eV on Si is NOT SCF noise (plateau across 4 orders of conv_threshold). Root cause: `symmetrize_density` applies Fd-3m's fractional translation τ=(1/4,1/4,1/4) via `nint` on an 18³ grid, and 18 is not divisible by 4 — every application of the glide bleeds ρ into the wrong grid point. Turning symmetry OFF drops the residual to 1e-8 eV. Fix = symmetrize ρ(G) in reciprocal space (QE-style phase factors); see **PCFX**.
