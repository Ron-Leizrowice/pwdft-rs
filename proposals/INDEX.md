# Proposal Index

Proposals use 4-letter IDs (e.g., `SIMP`) to avoid numbering conflicts when multiple agents work concurrently. Each proposal file is named `XXXX-slug.md` with YAML frontmatter containing metadata.

## Active

**2026-04-18:** Critical-Physics slate is **empty**. NCFX closed the 13.4 eV Si gap to 0.26 eV; PCFX (this release) closed the 1.2 eV per-component residual to 3.5e-11 eV by moving density symmetrization to G-space. The remaining ~23 meV Si residual vs QE is attributed to Monkhorst-Pack shifted-vs-Γ-centered grid convention (SYKP territory), below proposal-priority threshold.

**2026-04-19 (GRUM grooming pass):** VNLM closed (PR #49, 2.9–4.5× V_NL speedup); MODR closed (all phases A–D landed as PRs #46/#50/#48/#47). ITEV moved to "Deferred — Blocked on upstream" pending a faer 0.24 `iterate_lanczos` reorthogonalization bug fix. WFRX elevated from Low to High under a new "High — Performance" subsection: technique 1 (subspace diag) is independent of ITEV and delivers 20–30% SCF speedup on the current dense eigensolver.

**2026-04-19 (GRM2 grooming pass — 22 PRs #80–#101 merged today):** WFRX Technique 1 landed (PR #99, opt-in `scf.subspace_diag`, 7% at n_pw=725; Technique 2 stays deferred on ITEV). Promoted WFRX to Completed. Added MLFX/QELK/UNTS/DOCX/CLNP to Completed as small/reactive landings (no proposal files). GGAP Phase A landed (PR #85) — title annotated with phase state. ALOC Finding F-5 landed (PR #100, alloc traffic 16.8 GB → 0 per SCF at production sizes); F-7 and F-12 remain. TRV2 F1 (PR #98) and F3 (PR #96) landed; F2 (CCMX extraction) deferred on WFRX/driver refactor; 10 Category 2–5 findings remain. ERR2 P0 landed (PR #86); ERR2-AX (operations.rs annotations) and P1 (InvalidInput split) remain. MAUD-AC still in flight — title left as-is this pass. Machine-lock enforcement is now owner-scoped end-to-end (MLFX). Next strategic item: GGAP Phase B (PBE semilocal + gradient FFT helper) once MPSH drivers land; that unblocks 7 PBE validation cells in VQEF. Path forward — Validation: MPSH drivers (in flight) → 3 LDA cells; VGCH Phase 1 (Fe ecut sweep) after Phase 0 landed via CLNP. Perf: WFRX Technique 1 done, ALOC F-7/F-12 + GOPT PR-B next. "High — Foundation & Code Quality" subsection retained as a header slot but empty (MODR phases A–D all landed).

### High — Foundation & Code Quality

_No active entries (MODR's 4 phases all landed; see Completed)._

### High — Performance

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| TPRF | `[profile.test] opt-level=3` — shrink full gate from ~11 min to seconds (in flight — this PR) | trivial | low | — | TSPL |
| TSPL | Bifurcate test suite — fast Tier-1 default + heavy Tier-2 opt-in via `#[ignore]` (follow-up to TPRF) | small | low | TPRF | — |

### High — Validation

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| MPSH | Monkhorst-Pack shift alignment with QE convention (closes 4 light-atom `#[ignore]`s) | small-medium | low | — | VQEF |
| VGCH | Heavy-atom V_local(G) residual — post-VGCMP continuation (closes 5 Z>14 `#[ignore]`s; Phase 0 landed as CLNP PR #94; Phase 1 Fe ecut sweep next) | medium-large | medium | — | VQEF |
| VQEF | Full LDA+PBE QE validation matrix (8 systems × 2 functionals) — roadmap | medium | low | VGCMP, GGAP, QELK | — |

### Medium — Enhancements & Performance

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| CFGN | Expose Hardcoded Numerics as Settings | large | medium | — | — |
| MADOC | Mathematical documentation push (phased; MADOC-A first) | large | low | — | DLNT |
| GGAP | GGA/PBE exchange-correlation functional (Phase A dispatcher LANDED PR #85; Phases B–F ~7–10 CE-days; B unblocks 7 VQEF PBE cells) | large | medium | — | HYBR |
| HYBR | Hybrid functional (PBE0, HSE06) with ACE compression (phased 0–6; ~7–11 CE-weeks) | large | high | GGAP | — |
| ERR2 | Panic-free production (clippy::unwrap_used + structured InvalidInput split; P0 landed PR #86; ERR2-AX operations.rs annotations in flight; P1 InvalidInput split remains) | medium | low | — | — |
| GOPT | GPU kernel + wgpu host path optimization audit (scoping; 3 major + 5 modest + 4 micro findings; PR-B in flight; PRs A/C/D/E to follow) | medium | low-medium | — | — |
| TRV2 | Fresh test-suite review — post-PCFX/CCMX/NCFX/GGAP coverage pass (F1+F3 landed PRs #98/#96; F2 CCMX-extraction deferred on WFRX/driver refactor; 10 Categories 2–5 findings remain) | medium | low | — | — |
| MAUD | Mathematical accuracy audit of core physics modules (post-MADOC-A cold read; 1 docstring A + 9 C findings; MAUD-AC in flight will address top 2) | small | low | MADOC | — |
| ALOC | Per-iteration allocation audit for SCF hot loop (F-5 landed PR #100: alloc traffic 16.8 GB → 0 per SCF at production sizes; F-7 in flight; F-12 remains) | medium | low | — | — |

### Deferred — Blocked on upstream

| ID | Title | Complexity | Risk | Depends On | Blocks | Blocked On |
|----|-------|-----------|------|------------|--------|------------|
| ITEV | Iterative Eigensolver via `faer::partial_self_adjoint_eigen` (supersedes DVSN) | medium | medium | — | — | faer 0.24 `iterate_lanczos` reorthogonalization bug — revisit when upstream releases 0.25+ |

### Low / Deferred

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| CLSS | `cast_lossless` + Doc Hygiene | small | low | — | — |
| HKIN | Drop unused `Option<&dyn Fn>` V_eff param from `build_hamiltonian` + `compute_band_structure` (zero `Some` call sites) | trivial | low | — | — |
| DVSN | Iterative Eigensolver (Davidson / LOBPCG) — SUPERSEDED BY ITEV | large | medium | — | — |
| HD5I | HDF5 Restart and Structured Output | large | medium | — | — |
| SPRS | Sparse Matrix Support | large | medium | DVSN | — |
| CUCL | CubeCL GPU Kernels | large | high | — | — |
| FLUP | Follow-up backlog — 7 unpromoted items seeded from today's code reviews | small | low | — | — |

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
| PRPL | Periodic Pulay mixing (BROY Step 4; landed as PR #39 2026-04-18) |
| VGCMP | V_local(G) + KB Projector Cross-Check vs QE 7.5 (Phases 1–4; root cause moved to NCFX/NLCC chain, heavy-atom follow-up tracked in VGCH) |
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
| NLCC | NLCC Audit: Part A/B/C — ρ_core(G) unit tests, Fe E_xc regression guard, LFC docs |
| CCMX | Coupled-Channel Mixer for nspin=2 (Fe 4×4×4 free-mag: limit cycle → 14 iters, \|HF-KS\|=1e-4 eV) |
| PCFX | Density Symmetrization in G-Space (Non-Symmorphic τ Fix): per-component residual 1.2 eV → 3.5e-11 eV |
| CAST | Numeric Cast Safety Audit (3 correctness lints enabled; ~148 sites triaged; 4 assertion-guarded rewrites) |
| MXBA | Adaptive Mixing Beta (Eyert 1996 residual-norm monitor; opt-in, default off — hurts Fe CCMX) |
| TACC | Test Accuracy + Relevance Audit (findings #2–#6 landed via TACC-I; finding #1 fixed in TACC-II) |
| TYPE | Numeric-type audit Phase A — i32→i8 SpaceGroupOp rotations + i32→i16 Miller + dead `index_map` (5% on `symmetrize_density_g`) |
| VNLM | V_NL Hamiltonian assembly via single GEMM (n_pw=725: 27 ms → 5.2 ms, 5.2×; 4.5× per-k over 15 SCF iters) |
| MODR | Modular refactor — split god-modules (all 4 phases A–D landed: PRs #46/#50/#48/#47) |
| WFRX | Wavefunction reuse — Technique 1 (subspace diag) landed PR #99 opt-in, 7% at n_pw=725; Technique 2 deferred on ITEV upstream fix |
| MLFX | Machine-lock hardening — owner-scoped acquire, atomic mkdir, wait mode, PID liveness (4 bugs; PR #101) |
| QELK | Machine-lock policy — require lock for all QE runs (pw.x/ph.x/pp.x/etc.; CLAUDE.md policy update; PR #81) |
| UNTS | CLAUDE.md units correction — internal units are eV/Å (not Ry/Bohr); Ry/Bohr only at UPF boundary (PR #88) |
| DOCX | Fix `cargo doc --no-deps -- -D warnings` recipe — rejected by current cargo; use `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` (PR #95) |
| CLNP | Cleanup bundle — VGCH Phase 0 (ignore-string relabel on 5 heavy-atom tests) + ERR2 `BUG:` prefix on broyden expect (PR #94) |

## Notes

- **VGCMP + VGC5 + NCFX** (all done): the entire PP→H assembly pipeline is bit-correct vs QE, and NCFX (NLCC unit conversion + missing r²·4π in Bessel FT) closed the 13.4 eV Si gap to 0.26 eV. Residual attributed to Monkhorst-Pack shifted-vs-Γ-centered grid convention (SYKP). Prime suspect (V_local(G=0) shift) ruled out — already present.
- **CCMX** (done 2026-04-18): Mixer moved to QE's `(ρ_total, m) = (ρ↑+ρ↓, ρ↑−ρ↓)` basis (see `rhoz_or_updw` in `qe-7.5/PW/src/scf_mod.f90`). Kerker disabled on the magnetization channel (it's a charge-response, not a spin-response). Fe BCC 4×4×4 free-mag regression: 14 iters, |HF-KS| = 1.06e-4 eV, M = 0 μB (pre-CCMX: 200+ iters, Δρ limit cycle at 0.254, |HF-KS| ≈ 13 eV). `qe_validation.rs::test_fe_bcc_fm_vs_qe` now converges at 8×8×8 but remains `#[ignore]` pending VGCMP heavy-atom V_local fix (9.5 eV energy gap). Fixed-mag=2 still fails (PP prefers M=0; not a mixer problem).
- **TAUD** (done): all 5 PRs landed. PR D uncovered a sign-flipped V_local in the test's QE Cube reference (NOT in our Rust code — VGCMP Phase 1 already proved Rust correct). Captured as VLQR. Test re-`#[ignore]`'d with diagnostic numbers in the reason string.
- **XCPR** (done 2026-04-17): Steps 1+2 (XC grid) + Step 3 (spin-channel `rayon::join`) all landed. Step 3 speedup 1.12–1.24× (faer's internal gemm already saturates 8 cores during eigensolve); ceiling ~2× after ITEV drops per-k eigensolve cost.
- **CFGN** all dependencies satisfied (DDUP + SIMP done).
- **BROY** landed core algorithm; periodic Pulay (PRPL) landed 2026-04-18 as PR #39; adaptive-beta (MXBA) landed 2026-04-18 as opt-in (default off — see `proposals/completed/MXBA-adaptive-mixing-beta.md` completion note for the Fe CCMX failure mode).
- **NLCC** audit (2026-04-18, done): Part A/B/C landed — 4 `ρ_core(G)` unit tests for Si/Fe in `src/pseudopotential/upf.rs` (Python/SciPy reference at `scripts/validate/rho_core_g_reference.py`), `test_fe_bcc_xc_nlcc_regression_guard` E_xc defensive guard in `tests/qe_validation.rs` (|Δ_xc| = 0.692 eV vs QE, pre-NCFX baseline 48.85 eV), and LFC references + NLCC invariants added to module docstrings for `scf::energy`, `scf::mod::ScfParams`, `PseudopotentialData::core_charge`, and CLAUDE.md SCF loop section.
- **HD5I** references deleted `src/input.rs` — update to YAML Settings when implementing.
- **SOPT** (done 2026-04-17): refactored `ScfContext.symmetry: Option<&SymmetryInfo>` to always-present with identity-only fallback. Review caught a P1+TR regression in the initial `is_trivial()` short-circuit in main.rs; fixed by dropping the branch and tightening the predicate. Regression test pins the behavior.
- **ITEV** (new 2026-04-17, Performance Engineer): Post-FFTB/FMAD profiling shows eigensolver at 85-90% of SCF user CPU (n_pw=259: 56 ms/call; n_pw=725: 836 ms/call). `faer 0.24` ships `matrix_free::eigen::partial_self_adjoint_eigen` (implicitly-restarted Arnoldi, matrix-free via `LinOp`, warm-start via `v0`). Supersedes DVSN's hand-rolled Davidson plan. Projected 2.5-4× SCF wall-time speedup at production sizes. WFRX becomes the warm-start knob.
- **PCRS → PCFX** (2026-04-17, Researcher): per-component identity residual of 1.204 eV on Si is NOT SCF noise (plateau across 4 orders of conv_threshold). Root cause: `symmetrize_density` applies Fd-3m's fractional translation τ=(1/4,1/4,1/4) via `nint` on an 18³ grid, and 18 is not divisible by 4 — every application of the glide bleeds ρ into the wrong grid point.
- **PCFX** (done 2026-04-18, Core Engineer): `symmetrize_density_g` in `src/symmetry/density.rs` projects ρ onto the S-invariant subspace via `ρ_sym(G) = (1/N) Σ_S exp(-i G·τ) ρ(R^T G)` — exact for any fractional translation on band-limited input. Si per-component identity residual: 1.204 eV → 3.5e-11 eV (5 orders of magnitude improvement below the 1e-5 eV target). Si E_total shifted by 23 meV (PCFX correction, matches proposal's ~17 meV estimate). Fe BCC (symmorphic Im-3m, τ=0) unchanged. Convention details + band-limitation requirement documented inline.
- **WFRX** (done 2026-04-19, Performance Engineer): Technique 1 (subspace diagonalization on the current dense eigensolver) landed as PR #99 opt-in behind `scf.subspace_diag` (default off). Measured ~7% SCF wall-time improvement at n_pw=725 on the existing benches. Technique 2 (warm-start iterative) is ~30 lines once ITEV unblocks — the cache plumbing landed with Technique 1. See `proposals/completed/WFRX-wavefunction-reuse.md` completion note.
- **MLFX** (done 2026-04-19, Core Engineer): machine-lock hardening — four bugs fixed in one PR (#101). Owner worktree recorded and enforced by the Bash PreToolUse hook (cross-worktree cargo denied); atomic `mkdir` acquire (no racing acquires both winning); `--wait` mode with configurable timeout (long benches no longer get stolen at 30 min); PID-liveness check plus a 3-hour time cap for PID reuse on torn-down worktrees. 17-case shell-test suite added at `.claude/bin/tests/machine-lock.test.sh`.
- **QELK** (done 2026-04-19): policy update — every QE run (`pw.x`/`ph.x`/`pp.x`/`bands.x`/`projwfc.x`/etc.) now requires the machine lock because QE saturates CPU and corrupts bench wall-time measurements. CLAUDE.md § Machine Coordination spells out the scope (it's about benchmark integrity, not test isolation).
- **UNTS** (done 2026-04-19): CLAUDE.md claimed Ry/Bohr internal units; in reality the engine runs in eV/Å/e-Å⁻³ from the UPF boundary onward. `src/pseudopotential/upf/convert.rs` is the only conversion site. CLAUDE.md § Conventions now states this explicitly.
- **DOCX** (done 2026-04-19): `cargo doc --no-deps -- -D warnings` is rejected by current cargo (the `-D warnings` flag has to travel through `RUSTDOCFLAGS`). The recipe in CLAUDE.md § Code Quality is now `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps`.
- **CLNP** (done 2026-04-19): small cleanup bundle — VGCH Phase 0 relabeled five heavy-atom `#[ignore]` strings in `tests/qe_validation.rs` to reference VGCH instead of the retired VGCMP ticket, and ERR2 added a `BUG:` prefix to the `broyden.rs` `.expect(...)` so the panic text is greppable in production-panic triage. Tiny but worth its own landing to keep the commit atomic.
