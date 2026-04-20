# Proposal Index

Proposals use 4/5-letter IDs (e.g., `SIMP` or `UTRU4`) to avoid numbering conflicts when multiple agents work concurrently. Each proposal file is named `XXXXX-slug.md` with YAML frontmatter containing metadata.

## Active

### High — Validation (VQEF track)

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| PYQE | Route QE invocation/comparison through `pwdft_validation.qe` sub-package (phases A+B+C+D+E: runner + parser + lock + skill + CI gate) | large | medium | — | VQEF, VGCH-2, VGCH-MECH |
| RWHK | Reward-hacking audit (1 Critical GPU-vs-GPU tautology, 4 Major, 6 Minor). RWHK-FIX landed #171; remaining items open | small | low | — | VQEF |
| PZPW | Rewire LDA correlation from PZ-81 → Slater+PW92 to match QE's `SLA+PW` default. Predicted to close Fe LDA Class B outlier + compress Class A floor. PW92 helpers already exist in `src/potential/xc.rs` (used by PBE) | small | low | VGCH-2D | VQEF |
| CNLC | C diamond NLCC pin + NLCC-off ablation (primary Class A light-atom suspect per VGCH-2E). Template: Ga/As/O/Cl pins from VGCH-2F | small | low | VGCH-2E | VQEF |
| VNLM-CUD | Cu Γ-point per-m pin for the two-radial-d-projector sum in l=2 channel (H-C5 narrowed from VGCH-2F; unique vs Si in PP_DIJ structure) | small | low | VGCH-2F | VQEF |
| VGCH-2 | Heavy-atom residual hunt. Remaining open: H-C2 `n_bands` margin audit. (H1/H2/H3/H-C1/H-C3/H-C4 all CLEARED; H-C5 → VNLM-CUD) | medium | medium | — | VQEF |
| VGCH-MECH | Mechanism taxonomy — Class A (9 cells, tracked under VGCH-2 + PZPW + CNLC + VNLM-CUD), Class B (Fe LDA → PZPW), Class C emptied | medium | medium | VGCH-2, BSUM | VQEF |
| VQEF | Full LDA+PBE QE validation matrix (8 systems × 2 functionals). **Scoreboard: 4 GREEN / 12 YELLOW / 0 RED on E_total + 3 GREEN on E_F (Si/Al/C).** Caveat per RWHK C1: GPU-vs-CPU consistency signal structurally absent | medium | low | RWHK, VGCH-2, VGCH-MECH, PZPW, CNLC, VNLM-CUD, GGAP-E | — |

### Medium — Physics & Performance

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| PMTL | Split GPU backend into a new workspace crate `pwdft-metal`; define `GridAccelerator` trait in core; unblocks clean second-backend | medium | medium | — | CUCL |
| GGAP | GGA/PBE functional. Only Phase E remains (GPU PBE shader, deferred until CPU validates fully) | medium | medium | — | HYBR |
| HYBR | Hybrid functional (PBE0, HSE06) with ACE compression (phased 0–6; ~7–11 CE-weeks) | large | high | GGAP-E | — |
| ITEV | Iterative eigensolver — correctness defects closed (ITEV2 #140). Only step-4 remains: end-to-end SCF wall-time bench with WFRX active → decide default flip | small | low | — | — |
| GOPT | GPU audit — PR-A/B landed; remaining lever is F1+F2 chain fusion (blocked on `src/scf/driver.rs` quiescence) | medium | low-medium | — | — |
| ALOC | Per-iteration allocation audit — F-5 landed (16.8 GB → 0 per SCF); F-7 + F-12 remain | medium | low | — | — |
| ERR2 | Panic-free production — Phase 1 closed. Only P1.e `transplant.rs` 2-site mop-up remains (~20 min) | small | low | — | — |
| MADOC | Mathematical documentation push — A+B+C landed. Remaining phases D+E | medium | low | — | — |
| MAUD | Mathematical accuracy audit — MAUD-AC addressed top 2; 7 C-level items remain | small | low | MADOC | — |
| TRV2 | Fresh test-suite review — F1+F3 landed; 10 Categories 2–5 findings remain | medium | low | — | — |
| CFGN | Numerics knobs — CFGN1/DFLT/G2ZT landed; 10 knobs remain across Fermi-search / iterative-eigensolver / Ewald / floors | medium | low | — | — |
| ESPL | Split `ElectronSettings` — system physics vs convergence knobs; drop `scf.max_iter` default 100 → 50 | small | low | — | — |
| STYS | Settings type-sharpening audit — replace `String`-typed enum fields + `[usize;3]` flag-triples on `Settings` with their proper typed variants; surface serde errors for invalid values | medium | low | ESPL | — |

### Low / Deferred

| ID | Title | Complexity | Risk | Depends On | Blocks |
|----|-------|-----------|------|------------|--------|
| ROTI | Revert `SpaceGroupOp`/`SymmOp` rotation `[[i8;3];3]` → `[[i32;3];3]`; removes 4 `#[expect]` blocks + 27 `try_from` calls (~−80 LOC) | small | low | — | — |
| URES | Replace 82 `#[must_use]` annotations with `unused_results` lint | small | low | — | — |
| CUCL | CubeCL GPU kernels (deferred — explicit trigger conditions in proposal) | large | high | — | — |
| FLUP | Follow-up backlog — remaining unpromoted items: MXB2, EIGV/EIGW, FLP3, ITVF tracker | small | low | — | — |
| LTOB | Benchmark `lto = "fat"` vs `"thin"` on SCF hot path; adopt iff ≥2% production win and link <5 min | small | low | — | — |
| CISP | CI speedup follow-ups (mold, cargo-nextest, per-step job split). Pursue after TDBG lands real numbers | small | low | TDBG | — |

## Archived

Bodies retained in `proposals/completed/` with `status: archived` frontmatter.

| ID | Reason |
|----|--------|
| DVSN | Superseded by ITEV (faer native partial solver) |
| SPRS | Sparsity wrong for plane-wave DFT; V_eff convolution fills the Hamiltonian |
| HD5I | No user demand; re-open fresh when MD or geometry-opt lands |

## Reference Documents

- `psuedopotentials-sota.md` — survey of pseudopotential formalisms, libraries, and formats
