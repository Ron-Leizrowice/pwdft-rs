# Proposal Index

Proposals use 4-letter IDs (e.g., `SIMP`) to avoid numbering conflicts when multiple agents work concurrently. Each proposal file is named `XXXX-slug.md` with YAML frontmatter containing metadata.

## Active

**2026-04-19 (GRM11 end-of-day consolidation — 59 PRs merged today).** The day produced a full run of the GGA/PBE functional (Phases A+A.1+B+C+D+F-light all landed — Si PBE 12 meV GREEN, Al PBE 8 meV GREEN, Fe PBE retains M=2.16 μB at 1.97 eV VGCH-class), closed the VQEF matrix from `0 G / 8 Y / 8 R` → `4 G / 12 Y / 0 R` on E_total + `3 G` on E_F (Si/Al/C), and fractured the "VGCH heavy-atom residual" into three distinct mechanism classes (VGCH-MECH #168). Five root-cause hypotheses ruled out today: β_q projectors (VGCH-1a/b), SAD initial density (VGCH-1c), V_loc(G=0) gauge (SiEF-B1 #166), smearing entropy reporting (TSEN #162), mixer basin (VGCH-2B #167). Next-day entry point: **VGCH-2 Part C** on Cu Fermi-finder + Class B Fe LDA diagnostic + Class C C transplant in parallel. Full PR list at bottom Completed table. Archived today: **MIXA**, **TSPL**, **LOGH**, **MOAD**, **VGCH** (parent) — all landed.

**2026-04-18:** Critical-Physics slate is **empty**. NCFX closed the 13.4 eV Si gap to 0.26 eV; PCFX (this release) closed the 1.2 eV per-component residual to 3.5e-11 eV by moving density symmetrization to G-space. The remaining ~23 meV Si residual vs QE is attributed to Monkhorst-Pack shifted-vs-Γ-centered grid convention (SYKP territory), below proposal-priority threshold.

**2026-04-19 (GRUM grooming pass):** VNLM closed (PR #49, 2.9–4.5× V_NL speedup); MODR closed (all phases A–D landed as PRs #46/#50/#48/#47). ITEV moved to "Deferred — Blocked on upstream" pending a faer 0.24 `iterate_lanczos` reorthogonalization bug fix. WFRX elevated from Low to High under a new "High — Performance" subsection: technique 1 (subspace diag) is independent of ITEV and delivers 20–30% SCF speedup on the current dense eigensolver.

**2026-04-19 (Stack decision):** Engineering stack codified — observability stays on `log` + `env_logger` + `indicatif`; profiling adopts `samply` (PROF); benchmarks stay on `criterion`. The `tracing` ecosystem was considered and dropped; reopen triggers documented in PROF § "When to reconsider tracing". MIXL/LOGH/MOAD/DEAD/XCTH/PROF land independently — each is the simplest tool for its job rather than a piece of a unified observability framework.

**2026-04-19 (GRM5 grooming pass — stale-scope review):** Archived **DVSN** (superseded by ITEV — faer's native partial solver), **SPRS** (sparsity doesn't fit plane-wave basis — Hamiltonian is structurally dense because V_eff is a G-space convolution), **HD5I** (3 days deferred with no user demand; re-open fresh when MD/geometry-optimization lands). Re-scoped **CFGN** post-CFGN1 (#114): priority medium→low, complexity large→medium, 10 knobs left (Fermi search tuning, iterative eigensolver tuning, Ewald cutoffs, numerics floors); each should land individually when a user asks. Nine PRs merged today: MPSH #110, ERR2-AX #111, CLAU #113, CFGN1 #114, TPRF #115, PR #112 (observability), GRM4 #116, TPRFB #117, DCLN #118. New `Archived` section added between Deferred and Completed.

**2026-04-19 (GRM4 grooming pass — post-MY_THOUGHTS.md review):** Added 5 new proposals drawn from user review: **ESPL** (split ElectronSettings + max_iter 100→50), **ECUT** (per-PP recommended ecutwfc from PseudoDojo table, drop hardcoded 204.09 eV), **DCLN** (strip 57 proposal-ID + 51 QE tokens from public rustdoc), **ELMN** (trim atoms.rs to a `pub use`), **TYPB** (revert premature i16 Miller narrowing, `fft_grid_size` → u32). DCLN blocks MOAD (which writes 14 new module headers — land DCLN first so those headers are clean). TPRF landed as PR #115 (test profile opt-level=3; cargo test 11 min → 95 s, 7×). PR #112 landed as squash #116 (PROF/DEAD/XCTH/LOGH/MIXL/MOAD + TRCE delete). 6 total PRs merged this pass: MPSH #110, ERR2-AX #111, CLAU #113, CFGN1 #114, TPRF #115, PR #112. MPSH empirically refuted the shift-convention prior for C/Al/Fe residuals — only Si E_total closed; C stall and Al 83 meV gap are NOT shift-related and need separate proposals (C mixer/ecut, Al ecut/Kerker). VGCH V_loc(G=0) eigenvalue zero-reference (1.35 eV Si E_F shift) is a distinct issue from VGCH's heavy-atom residual — may warrant its own proposal rather than folding into VGCH Phase 1.

**2026-04-19 (GRM2 grooming pass — 22 PRs #80–#101 merged today):** WFRX Technique 1 landed (PR #99, opt-in `scf.subspace_diag`, 7% at n_pw=725; Technique 2 stays deferred on ITEV). Promoted WFRX to Completed. Added MLFX/QELK/UNTS/DOCX/CLNP to Completed as small/reactive landings (no proposal files). GGAP Phase A landed (PR #85) — title annotated with phase state. ALOC Finding F-5 landed (PR #100, alloc traffic 16.8 GB → 0 per SCF at production sizes); F-7 and F-12 remain. TRV2 F1 (PR #98) and F3 (PR #96) landed; F2 (CCMX extraction) deferred on WFRX/driver refactor; 10 Category 2–5 findings remain. ERR2 P0 landed (PR #86); ERR2-AX (operations.rs annotations) and P1 (InvalidInput split) remain. MAUD-AC still in flight — title left as-is this pass. Machine-lock enforcement is now owner-scoped end-to-end (MLFX). Next strategic item: GGAP Phase B (PBE semilocal + gradient FFT helper) once MPSH drivers land; that unblocks 7 PBE validation cells in VQEF. Path forward — Validation: MPSH drivers (in flight) → 3 LDA cells; VGCH Phase 1 (Fe ecut sweep) after Phase 0 landed via CLNP. Perf: WFRX Technique 1 done, ALOC F-7/F-12 + GOPT PR-B next. "High — Foundation & Code Quality" subsection retained as a header slot but empty (MODR phases A–D all landed).

### High — Validation (VQEF track)

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| RWHK | Reward-hacking audit (2026-04-19). **1 Critical** (C1: `test_gpu_vs_cpu_scf_direct_comparison` compares GPU-vs-GPU, not GPU-vs-CPU — f32/f64 consistency signal structurally absent), **4 Major** (E2 ITEV heuristic narrow coverage; F3 NLCC convention wrong-sign undefended; H1 PBE-actually-invoked soft assertion gap; A5 five `#[ignore]`d LDA tolerances don't match their disclosed residuals), **6 Minor** | small | low | — | VQEF |
| VGCH-2 | Heavy-atom residual hunt. Part A (#160) per-term trace; Part B (#167) shared-density transplant on Cu **refuted H3 mixer-basin** (+16.34 eV gap at ρ_QE — functional disagrees at same density). Part C (next) = Fermi-finder Python reference, `n_bands` margin audit, smearing function alignment, Cu/GaAs/MgO ρ_core pinning, `D_ij·β·β` contraction cross-check | medium | medium | — | VQEF |
| VGCH-MECH | Mechanism taxonomy — splits remaining 12 YELLOW cells into **Class A** (8 cells, energy-functional-at-shared-density, continued under VGCH-2 Part C), **Class B** (Fe LDA Hamiltonian-side outlier per BSUM ratio 0.99× — PZ-vs-PW92 / spin-path hypothesis, 1 cell, 1 CE-day diagnostic), **Class C** (C diamond — VGCH-2E reclassified to Class A at light-atom magnitude, ratio 2.25 at shared ρ). Classes B+C can run parallel to VGCH-2 Part C | medium | medium | VGCH-2, BSUM | VQEF |
| VGCH-2D | Fe LDA Class B diagnostic — discriminates PZ-vs-PW92 correlation vs. spin-driver (CCMX) bias via 2×2 matrix of (PZ/PW92) × (nspin=1/2) iter-1 transplants from QE's M-collapsed ρ. Key ground truth: every `qe_validation/*_scf.out` LDA deck reports `SLA+PW` (PW92), while pwdft-rs `XcEvaluator::Pz` uses PZ-81 — global functional mismatch on LDA, amplified on transition metals. Diagnostic-only, ~1 CE-day | small | low | VGCH-MECH, VGCH-2 | VQEF |
| VGCH-2E | C diamond Class C transplant — at shared ρ_QE, pwdft-rs E_total is +2.28 eV off QE (larger than SCF 1.45 eV); ΔE_1e=+1.87, ΔE_H=−0.83, ratio 2.25 — Class A fingerprint at light-atom magnitude. E_F matches to 11 meV (Fermi-finder NOT the bug); mixer-basin hypothesis refuted. Primary suspect: NLCC ρ_core(G) Bessel transform on C (C LDA PP has `core_correction="T"`). Follow-up: pin ρ_core(G) + NLCC-off ablation | small | low | VGCH-MECH, VGCH-2B | VQEF |
| VQEF | Full LDA+PBE QE validation matrix (8 systems × 2 functionals). **Scoreboard (EOD 2026-04-19): 4 GREEN / 12 YELLOW / 0 RED on E_total + 3 GREEN on E_F (Si / Al / C).** **Caveat per RWHK C1:** GPU-vs-CPU consistency signal is structurally absent — the `test_gpu_vs_cpu_scf_direct_comparison` test runs GPU on both branches. Closure path sequenced via VGCH-MECH | medium | low | RWHK, VGCH-2, VGCH-MECH, GGAP-E | — |

### Medium — Physics & Performance

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| GGAP | GGA/PBE functional — **Phases A+A.1+B+C+D+F-light all landed** (#85, #155, #145, #151, #158, #161). Si PBE 12 meV GREEN, Al PBE 8 meV GREEN, Fe PBE retains M=2.16 μB (VGCH-class residual). **Only Phase E remains** (GPU PBE shader, 2–3 CE-days, deferred until CPU path validates fully) | medium | medium | — | HYBR |
| HYBR | Hybrid functional (PBE0, HSE06) with ACE compression (phased 0–6; ~7–11 CE-weeks) | large | high | GGAP-E | — |
| ITEV | Iterative Eigensolver — faer Lanczos fix vendored (#129); both correctness defects closed (ITEV2 #140, Si ecut=100 Dense↔Iterative | ΔE | =4.52e-12 eV). **Only Phase-5 step-4 remains**: end-to-end SCF wall-time bench with WFRX active → decide default flip | small | low |
| GOPT | GPU audit — PR-B #106 + PR-A #138 landed (128³ hartree −34%, v_eff −40%); F5/F6/F8 empirically tested, no measurable gain; remaining lever is F1+F2 chain fusion (blocked on `src/scf/driver.rs` quiescence) | medium | low-medium | — | — |
| ALOC | Per-iteration allocation audit — F-5 landed (16.8 GB → 0 per SCF); F-7 (psi_g in band loop, 200-1000 µs/iter) + F-12 (FFT3D twiddle rebuild, 100-500 µs/iter) remain | medium | low | — | — |
| ERR2 | Panic-free production — P0 #86, ERR2-AX #111, P1.a #152 (variants), P1.b #159 (CRYSTAL 3), P1.c #163 (PARAM 9) all landed. **Only P1.d remains** (ELEMENT + UPF, 2 sites, mechanical ~30 min) | small | low | — | — |
| MADOC | Mathematical documentation push — A + B + C phases landed (MADOC-A series, MADOC-B #121, CLSS #122 doc hygiene). **Remaining phases D+E** cover unaddressed modules | medium | low | — | — |
| MAUD | Mathematical accuracy audit — MAUD-AC addressed top 2 findings; 7 C-level docstring accuracy items remain | small | low | MADOC | — |
| TRV2 | Fresh test-suite review — F1+F3 landed (#98/#96); F2 CCMX-extraction deferred on WFRX/driver refactor; 10 Categories 2–5 findings remain | medium | low | — | — |
| CFGN | Numerics knobs — CFGN1 #114 + DFLT #137 + G2ZT #153 landed; 10 knobs remain across Fermi-search / iterative-eigensolver / Ewald / floors; each lands as a separate small PR on demand | medium | low | — | — |
| ESPL | Split `ElectronSettings` — system physics vs convergence knobs; drop `scf.max_iter` default 100 → 50 | small | low | — | — |

### Low / Deferred

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| URES | Replace 82 `#[must_use]` annotations with `unused_results` lint — unconditional coverage without per-function annotation | small | low | — | — |
| CUCL | CubeCL GPU Kernels (deferred — explicit trigger conditions in proposal) | large | high | — | — |
| FLUP | Follow-up backlog — remaining unpromoted items: MXB2 (Fe CCMX retune), EIGV/EIGW (bench-noise triage), FLP3 (NLCC 60-element regression parametric expansion), TYPE-AX (folded into TYPB) | small | low | — | — |
| LTOB | Benchmark `lto = "fat"` vs `"thin"` on the SCF hot path; adopt iff ≥2% production SCF win and link time <5 min | small | low | — | — |

## Reference Documents

- `psuedopotentials-sota.md` — Survey of pseudopotential formalisms, libraries, and formats

## Archived

Proposals moved out of the active/deferred backlog after a stale-scope review. File bodies kept in `proposals/completed/` for historical reference; frontmatter `status: archived` + `archived_reason` documents why. Re-open as a fresh proposal if the premise changes.

| ID | Title | Archived | Reason |
|----|-------|----------|--------|
| DVSN | Iterative Eigensolver (Davidson / LOBPCG) | 2026-04-19 | Superseded by ITEV (faer's native `partial_self_adjoint_eigen`); custom Davidson/LOBPCG is wasted work |
| SPRS | Sparse Matrix Support | 2026-04-19 | Sparsity assumption wrong for plane-wave DFT — V_eff convolution fills the Hamiltonian; sparse storage only helps real-space-grid or atomic-orbital bases (neither is this project's direction) |
| HD5I | HDF5 Restart and Structured Output | 2026-04-19 | 3 days deferred with no user demand and no blocking dependency; re-open fresh when MD or geometry-optimization lands |

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
| SPXC | Fix Spin-Polarized E_xc Density Consistency ( |
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
| MPSH | Monkhorst-Pack Γ-centered default matches QE convention (PR #110; closes Si E_total convention gap; C/Al/Fe residuals empirically refuted as shift-related) |
| MADOC-B | Mathematical completeness for physics-kernel docstrings (9 P/A → E upgrades in `potential/xc.rs`; PR #121) |
| DCLN | Strip proposal-ID tokens + QE references from public rustdoc (physics-first prose; PR #118) |
| TPRF | `[profile.test] opt-level=3` + `TSPL` proposal seeded (cargo test 11 min → 95 s, 7×; PR #115) |
| TPRFB | Share dep artifacts between dev and test profiles (eliminate 30–60 s first-compile hit; PR #117) |
| CFGN1 | Expose initial-density `gaussian_sigma` as YAML-configurable input (first CFGN knob; PR #114) |
| CLAU | Consolidate 33 per-module ERR2 test allows → crate-level `cfg_attr(test)` (PR #113) |
| ERR2-AX | Annotate TYPE-A narrowing expect sites with `#[expect(reason=...)]` (PR #111) |
| HDWR | Fix Apple M2 → M3 Max in docs + auto-detect MPI rank counts (no hardcode; PR #108) |
| CLSS | `cast_lossless` + `missing_errors_doc` + `missing_panics_doc` enabled; 85-site migration (PR #122) |
| HKIN | Drop unused `Option<&dyn Fn>` V_eff param from `build_hamiltonian`/`compute_band_structure` (−26 LOC; PR #123) |
| XCTH | Remove `XC_PARALLEL_THRESHOLD`; always rayon in `lda_xc_grid` / `lda_xc_spin_grid` (consistency with rest of engine; PR #124) |
| ELMN | Trim `src/atoms.rs` to `pub use mendeleev::Element;` (−58 LOC; 4 call sites use `Element::iter().find(...)`; PR #126) |
| DEAD | DRSD drop `symmetry/density/real_space.rs` (deprecated real-space symmetrizer) + SMRT delete misleading test + DHPC narrow `diagonalize_hermitian` to `pub(crate)` (−239 LOC; PR #127) |
| MADOC-C | MADOC phase C — `CLSS`: `cast_lossless` lint + `# Errors` / `# Panics` doc hygiene gate (PR #122) |
| CUCL | Explicit trigger conditions for un-deferring CubeCL (GPU non-local / GGAP Phase E / wgpu regression / CubeCL 1.0) (PR #120) |
| GRM8 | Vendor faer v0.24.0 + MAX_REORTH `iterate_lanczos` fix; wire `[patch.crates-io]` at `./faer/` (PR #129) |
| GRM9 | Archive merged proposal files (HKIN/XCTH/ELMN/DEAD/MPSH/DCLN); drop ITEV "not vendored yet" caveat; move ITEV out of Deferred (PR #130) |
| PROF | Adopt `samply` as canonical profiler; codify observability vs. profiling vs. benchmarking split in CLAUDE.md + perf-engineer agent + logbook (PR #131) |
| LOGH (MELG) | `eprintln!` → `log::info!` migration for `src/main.rs` SCF summary (PR #132; PCEP phase for `src/pseudopotential/**` tracked as LOGH-2) |
| MIXL | Mixer init `info!` + auto-q_TF + DIIS truncation `debug!` + adaptive-β trigger `debug!` (PR #133; +217 LOC) |
| TYPB | Integer type cleanup — i16 Miller → i32, `fft_grid_size` → u32, 4 i8-rotation `expect` sites annotated, 3 cast suppressions removed (PR #134; −78 LOC, bench −0.1% noise) |
| MOAD | Module-orientation `//!` docstrings on 11 src files — crate-root cargo doc landing page populated (PR #135; +158 LOC) |
| ECUT | Per-PP recommended `ecutwfc` from PseudoDojo `.standard` table; `BasisSettings::ecutwfc` → `Option<f64>`; 85-element lookup table with 100 eV safety floor (PR #136) |
| DFLT | Hoist two `1e-15` density threshold literals into named `const`s with physics docstrings (FLUP entry; PR #137) |
| GOPT-A | BufferPool scratch `Vec<f32>` for CPU→GPU uploads; 128³ hartree −33.9%, v_eff −39.9%; 64³ hartree −28.6%, v_eff −20.4% (PR #138) |
| VGCH Phase 1a | Heavy-atom per-component diagnostic — rules out V_local(G=0) Z-scaling and Ewald; residual lives in one-electron/Hartree partial cancellation, not form factors. Phase 1b needs β_q / initial density / mixer basin investigation (PR #139; diagnostic-only, no src/ changes) |
| ITEV2 | Close ITEV defects 1+2 — adaptive `krylov_max_dim` scales with basis (Lehoucq & Sorensen §3.2); WFRX warm-start wired on Iterative path. Si ecut=100 Dense↔Iterative |
| GRM6 | Consolidate TYPB + TYPE-AX into one integer-type cleanup; archive merged CLSS + HKIN (PR #125) |
| GRM7 | Fold ITEVF findings into ITEV; fix CLAUDE.md perf claim + DEAD drive-by (PR #128) |
| GRM8 | Vendor faer v0.24.0 + MAX_REORTH `iterate_lanczos` fix; wire `[patch.crates-io]` at `./faer/` (PR #129) |
| GRM9 | Archive merged proposal files + move ITEV out of Deferred (PR #130) |
| GRM10 | Archive MIXL/TYPB/ECUT/PROF; refresh ITEV/VGCH/GOPT status (PR #141) |
| LGB2 | Consolidate logbooks — append 8 pending entries from 46-PR wave (PR #157) |
| LOGH-2 | Delete 13 redundant `eprintln!` in `src/pseudopotential/**` (PCEP; PR #142) |
| MOAD-2 | Finish `//!` module headers + fix stale `EigensolverKind` claim (PR #143) |
| VQEF-QC | Si/Al/C quickchecks — Si E_total flips first GREEN; Al reclassified basis-truncation → VGCH-light; C mixer-stall pinned via Broyden+Kerker (PR #144) |
| GGAP Phase B | PBE exchange (non-spin) — `pbex` with iflag=1 (PR #145) |
| VQEF-AL | Regenerate QE Al ref at ecut=24 Ry — reclassified to VGCH light-atom "different converged density" class (PR #146) |
| ERR2 P1 scoping | InvalidInput variant split plan — 15 sites → 4 variants across P1.a–d (PR #147) |
| VGCH Phase 1b | β_l(q) projector form factors bit-perfect vs QE (H1 CLEARED; PR #148) |
| MIXA | Pin Plain-Anderson-stall-on-wide-gap-insulators as negative regression (C diamond, PR #149) |
| TSPL | Bifurcate test suite — Tier-1 fast default (12s) + Tier-2 `#[ignore]` opt-in (58s); 8× speedup on default `cargo test` (PR #150) |
| GGAP Phase C | PBE correlation + PW92 helper; wire end-to-end (scaffolded, needs A.1 to run; PR #151) |
| ERR2 P1.a | Add InvalidCrystal / InvalidParam / UnknownElement / InvalidPseudopotential variants (PR #152) |
| FLUP G2ZT | Hoist `1e-12` |
| GGAP F-pre | QE PBE reference data for all 8 VQEF systems (PR #154) |
| GGAP Phase A.1 | Driver-side density gradient FFT + semilocal V_xc assembly — **Si PBE end-to-end 12.4 meV** (PR #155) |
| VGCH Phase 1c | SAD initial density bit-correct vs QE (H2 CLEARED; PR #156) |
| GGAP Phase D | Spin-polarized PBE (`pbex` spin-scaling + `pbec_spin`); **Fe BCC FM PBE retains ferromagnetism** (M=2.16 μB vs QE 2.34; residual 1.97 eV VGCH-class; PR #158) |
| ERR2 P1.b | Migrate CRYSTAL cluster (3 sites) from InvalidInput to InvalidCrystal (PR #159) |
| VGCH-2 Part A | Per-term energy trace on 7 heavy-atom systems — localizes residual to SCF dynamics not assembly; heavy-atom `|ΔE_1e|/|−ΔE_H|` ∈ [1.2, 2.9] linear-response δρ signature (PR #160) |
| GGAP Phase F-light | Wire Al/C/Cu/GaAs/NaCl/MgO PBE tests; Si PBE tolerance tightened; **full 16-cell VQEF matrix populated**. MgO PBE 6.5× better than LDA (largest functional-sensitive improvement) (PR #161) |
| TSEN | Add missing `-TS` smearing contribution to `total_energy` — closes 70-260 meV on every metal; **Al LDA + Al PBE flip GREEN**; Harris-Foulkes gets `-TS` symmetrically (PR #162) |
| ERR2 P1.c | Migrate PARAM cluster (9 sites); `InvalidInput` catch-all: 12 → 3 remaining (PR #163) |
| Si E_F diag | Diagnostic: 1.35 eV rigid offset, std 0.6 meV across 8 Si bands — proves pure V_loc(G=0) gauge mismatch; fix scoped for VGCH-SiEF B1 (PR #164) |
| BSUM | Band-sum identity gate `|Σ w_k f_ik ε_ik|` for all 16 VQEF cells; reveals heavy-atom ratio 1.5-3.3× partial-cancellation; Fe LDA exception at 0.99× (Hamiltonian-side outlier) (PR #165) |
| VGCH-SiEF-B1 | V_loc(G=0) gauge fix — stop zeroing Hamiltonian diagonal; deleted `with_g0_shift`; closes **Si E_F 1.35 eV → 6.4 meV GREEN** + new C E_F GREEN + Al E_F GREEN; E_total algebraic identity preserved (0.023 meV drift) (PR #166) |
| VGCH-2 Part B | Cu transplant experiment — **H3 mixer-basin REFUTED** (+16.34 eV gap at shared density on Cu); new leading suspect = Fermi-finder/smearing on dense d-manifold DOS (PR #167) |
| VGCH-MECH | Mechanism taxonomy — three-class split of remaining VQEF YELLOWs (Class A 8 / Class B 1 / Class C 2); Part B+C can run parallel to Class A (PR #168) |

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
