# Grooming log

Historical grooming-pass consolidation paragraphs rotated out of `INDEX.md`. `INDEX.md` keeps only the most recent pass; older ones live here for regression archaeology (e.g. "when did we archive DVSN?", "what was the VQEF scoreboard on 2026-04-19?").

## 2026-04-19 (GRM11 end-of-day consolidation — 59 PRs merged)

The day produced a full run of the GGA/PBE functional (Phases A+A.1+B+C+D+F-light all landed — Si PBE 12 meV GREEN, Al PBE 8 meV GREEN, Fe PBE retains M=2.16 μB at 1.97 eV VGCH-class), closed the VQEF matrix from `0 G / 8 Y / 8 R` → `4 G / 12 Y / 0 R` on E_total + `3 G` on E_F (Si/Al/C), and fractured the "VGCH heavy-atom residual" into three distinct mechanism classes (VGCH-MECH #168). Five root-cause hypotheses ruled out: β_q projectors (VGCH-1a/b), SAD initial density (VGCH-1c), V_loc(G=0) gauge (SiEF-B1 #166), smearing entropy reporting (TSEN #162), mixer basin (VGCH-2B #167). Archived: **MIXA**, **TSPL**, **LOGH**, **MOAD**, **VGCH** (parent) — all landed.

## 2026-04-18

Critical-Physics slate is empty. NCFX closed the 13.4 eV Si gap to 0.26 eV; PCFX closed the 1.2 eV per-component residual to 3.5e-11 eV by moving density symmetrization to G-space. The remaining ~23 meV Si residual vs QE is attributed to Monkhorst-Pack shifted-vs-Γ-centered grid convention (SYKP territory), below proposal-priority threshold.

## 2026-04-19 (GRUM grooming pass)

VNLM closed (PR #49, 2.9–4.5× V_NL speedup); MODR closed (all phases A–D landed as PRs #46/#50/#48/#47). ITEV moved to "Deferred — Blocked on upstream" pending a faer 0.24 `iterate_lanczos` reorthogonalization bug fix. WFRX elevated from Low to High under a new "High — Performance" subsection: technique 1 (subspace diag) is independent of ITEV and delivers 20–30 % SCF speedup on the current dense eigensolver.

## 2026-04-19 (Stack decision)

Engineering stack codified — observability stays on `log` + `env_logger` + `indicatif`; profiling adopts `samply` (PROF); benchmarks stay on `criterion`. The `tracing` ecosystem was considered and dropped; reopen triggers documented in PROF § "When to reconsider tracing". MIXL / LOGH / MOAD / DEAD / XCTH / PROF land independently — each is the simplest tool for its job rather than a piece of a unified observability framework.

## 2026-04-19 (GRM5 grooming pass — stale-scope review)

Archived **DVSN** (superseded by ITEV — faer's native partial solver), **SPRS** (sparsity doesn't fit plane-wave basis — Hamiltonian is structurally dense because V_eff is a G-space convolution), **HD5I** (3 days deferred with no user demand; re-open fresh when MD / geometry-optimization lands). Re-scoped **CFGN** post-CFGN1 (#114): priority medium→low, complexity large→medium, 10 knobs left. Nine PRs merged: MPSH #110, ERR2-AX #111, CLAU #113, CFGN1 #114, TPRF #115, PR #112 (observability), GRM4 #116, TPRFB #117, DCLN #118. New `Archived` section added between Deferred and Completed.

## 2026-04-19 (GRM4 grooming pass — post-MY_THOUGHTS.md review)

Added 5 new proposals drawn from user review: **ESPL** (split ElectronSettings + max_iter 100→50), **ECUT** (per-PP recommended ecutwfc from PseudoDojo table, drop hardcoded 204.09 eV), **DCLN** (strip 57 proposal-ID + 51 QE tokens from public rustdoc), **ELMN** (trim `atoms.rs` to a `pub use`), **TYPB** (revert premature i16 Miller narrowing, `fft_grid_size` → u32). DCLN blocks MOAD (which writes 14 new module headers — land DCLN first so those headers are clean). TPRF landed as PR #115 (test profile `opt-level = 3`; test wall 11 min → 95 s, 7×). PR #112 landed as squash #116 (PROF / DEAD / XCTH / LOGH / MIXL / MOAD + TRCE delete). MPSH empirically refuted the shift-convention prior for C / Al / Fe residuals — only Si E_total closed; C stall and Al 83 meV gap are NOT shift-related and need separate proposals (C mixer / ecut, Al ecut / Kerker). VGCH V_loc(G=0) eigenvalue zero-reference (1.35 eV Si E_F shift) is a distinct issue from VGCH's heavy-atom residual.

## 2026-04-19 (GRM2 grooming pass — 22 PRs #80–#101 merged)

WFRX Technique 1 landed (PR #99, opt-in `scf.subspace_diag`, 7 % at `n_pw = 725`; Technique 2 stays deferred on ITEV). Promoted WFRX to Completed. Added MLFX / QELK / UNTS / DOCX / CLNP to Completed as small / reactive landings (no proposal files). GGAP Phase A landed (PR #85). ALOC Finding F-5 landed (PR #100, alloc traffic 16.8 GB → 0 per SCF at production sizes); F-7 and F-12 remain. TRV2 F1 (PR #98) and F3 (PR #96) landed; F2 (CCMX extraction) deferred on WFRX / driver refactor; 10 Category 2–5 findings remain. ERR2 P0 landed (PR #86); ERR2-AX (`operations.rs` annotations) and P1 (`InvalidInput` split) remain. MAUD-AC still in flight. Machine-lock enforcement is now owner-scoped end-to-end (MLFX). Next strategic item: GGAP Phase B (PBE semilocal + gradient FFT helper) once MPSH drivers land; that unblocks 7 PBE validation cells in VQEF.
