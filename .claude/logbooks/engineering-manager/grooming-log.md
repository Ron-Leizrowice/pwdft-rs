# EM grooming log

Rolling record of INDEX grooming passes. Paragraph-per-pass; keep the file lean (last ~5 passes is plenty — older history is reconstructable from `git log proposals/INDEX.md` and the `proposals/completed-log.md` shipping table).

## 2026-04-20 (GRM12 — EOD consolidation + INDEX slim)

7 PRs merged today + 1 direct-to-main micro-commit (72 cumulative). Three waves: **morning (physics)** — URES #172 + strict-stance revision #173 rejecting bare `let _ =` in production. **Midday (physics)** — **ERR2 Phase 1 CLOSED** (#175 P1.d; only P1.e `transplant.rs` 2-site mop-up remains); **DOCLEAN** (#176) stale `SI_REFERENCE_FERMI_EV` + two CLAUDE.md Conventions gotchas + `assert_band_sum_matches_qe` → `assert_one_electron_sum_matches_qe` rename (RWHK H3). Three VGCH Class B/C diagnostics: **VGCH-2D** (#177) designs the Fe LDA PZ-vs-PW92 discriminator and surfaces a global LDA functional mismatch (every QE LDA deck uses `SLA+PW` / PW92, pwdft-rs uses Slater+PZ-81 at `pwdft/pwdft-core/src/potential/xc.rs:143,394`); **VGCH-2E** (#178) reclassifies C diamond Class C → Class A light-atom (primary suspect NLCC ρ_core(G) on C); **VGCH-2F** (#179) refutes H-C4 (Ga/As/O/Cl ρ_core(G) all pinned < 1e-5 e/Å³), narrows H-C5 to a Cu two-radial-d-projector sum (filed as VNLM-CUD). **Afternoon (infra)** — **CICI** (#180) wires the clippy `-D warnings` GitHub Action on default + `--features gpu`; **DOCB2** (#181) bulk docs/proposal-file sync landing CISP/PMTL/PYQE/ROTI/TDBG/URES bodies + agent-def updates; **CINM** (`236c747` direct-to-main) renames the two CI workflow check names for branch-protection wiring. **EOD groom slim** — moved this grooming log out of `proposals/grooming-log.md` into the EM logbook; moved the Completed proposal table + Notes out of INDEX into `proposals/completed-log.md`; INDEX is now only Active + Archived tables + the Reference Documents pointer. **Next-day entry points** in priority order: (1) **PZPW** — rewire LDA to SLA+PW92 (Researcher spawned to draft proposal); (2) **CNLC** — C NLCC ρ_core pin + NLCC-off ablation (Researcher); (3) **VNLM-CUD** — Cu per-m two-radial-d diagnostic (Researcher); (4) **ERR2-P1.e** — `transplant.rs` 2-site mop-up.

## 2026-04-19 (GRM11 — end-of-day consolidation, 59 PRs merged)

Full run of the GGA/PBE functional (Phases A+A.1+B+C+D+F-light all landed — Si PBE 12 meV GREEN, Al PBE 8 meV GREEN, Fe PBE retains M=2.16 μB at 1.97 eV VGCH-class). Closed the VQEF matrix from `0 G / 8 Y / 8 R` → `4 G / 12 Y / 0 R` on E_total + `3 G` on E_F (Si/Al/C). Fractured the "VGCH heavy-atom residual" into three distinct mechanism classes (VGCH-MECH #168). Five root-cause hypotheses ruled out: β_q projectors (VGCH-1a/b), SAD initial density (VGCH-1c), V_loc(G=0) gauge (SiEF-B1 #166), smearing entropy reporting (TSEN #162), mixer basin (VGCH-2B #167). Archived: **MIXA**, **TSPL**, **LOGH**, **MOAD**, **VGCH** (parent) — all landed.

## 2026-04-19 (GRM5 — stale-scope review)

Archived **DVSN** (superseded by ITEV — faer's native partial solver), **SPRS** (sparsity doesn't fit plane-wave basis — Hamiltonian is structurally dense because V_eff is a G-space convolution), **HD5I** (3 days deferred with no user demand; re-open fresh when MD / geometry-optimization lands). Re-scoped **CFGN** post-CFGN1 (#114): priority medium→low, complexity large→medium, 10 knobs left. Nine PRs merged. New `Archived` section added between Deferred and Completed.

## 2026-04-19 (GRUM — post-VNLM/MODR)

VNLM closed (PR #49, 2.9–4.5× V_NL speedup); MODR closed (all phases A–D landed as PRs #46/#50/#48/#47). ITEV moved to "Deferred — Blocked on upstream" pending a faer 0.24 `iterate_lanczos` reorthogonalization bug fix (later resolved by GRM8 vendoring the fix). WFRX elevated from Low to High under a new "High — Performance" subsection: technique 1 (subspace diag) is independent of ITEV and delivers 20–30 % SCF speedup on the current dense eigensolver.

## 2026-04-19 (Stack decision)

Engineering stack codified — observability stays on `log` + `env_logger` + `indicatif`; profiling adopts `samply` (PROF); benchmarks stay on `criterion`. The `tracing` ecosystem was considered and dropped; reopen triggers documented in PROF § "When to reconsider tracing". MIXL / LOGH / MOAD / DEAD / XCTH / PROF land independently — each is the simplest tool for its job rather than a piece of a unified observability framework.
