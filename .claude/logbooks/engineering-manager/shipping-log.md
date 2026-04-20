# EM shipping log

Rolling record of shipped proposals and load-bearing notes that carry over beyond completion. Authoritative file bodies live in `proposals/completed/`; this is the skimmable roll-up. Grouped thematically, not chronologically — `git log proposals/` is the dated source of truth.

Cumulative landings through **2026-04-20**: **72 PRs**.

## Physics validation — CLOSED

- **PP→H assembly bit-correct vs QE** (VGCMP, VGC5, NCFX) — 13.4 eV Si gap → 0.26 eV
- **G-space density symmetrization** (PCFX) — Fd-3m τ=(1/4,1/4,1/4): per-component residual 1.2 eV → 3.5e-11 eV
- **Si/C/Al E_F GREEN** (VGCH-SiEF-B1) — V_loc(G=0) gauge; stopped zeroing H diagonal
- **TS smearing in total energy** (TSEN) — flipped Al LDA + Al PBE GREEN; 70–260 meV on every metal
- **MP Γ-centered default** (MPSH) — matches QE convention
- **Spin mixer** (CCMX) — `(ρ_total, m)` basis; Fe BCC FM limit cycle → 14 iters, \|HF-KS\|=1e-4 eV
- **PBE semilocal stack** (GGAP A/A.1/B/C/D/F-light) — Si 12 meV GREEN, Al 8 meV GREEN, Fe retains FM M=2.16 μB
- **Heavy-atom residual localization** (VGCH-2 A/B/C + VGCH-MECH) — 5 hypotheses ruled out (β_q, SAD, V_loc(G=0), mixer basin, smearing fn, Fermi-finder, ρ_core(G) for 5 metals); three-class mechanism taxonomy

## Perf wins

- **Eigensolver correctness** (ITEV2) — defects 1+2 closed; Si ecut=100 Dense↔Iterative \|ΔE\|=4.5e-12 eV
- **V_NL GEMM** (VNLM) — 27 ms → 5.2 ms at n_pw=725 (5.2×)
- **Subspace warm-start** (WFRX) — opt-in; 7% at n_pw=725
- **GPU upload pool** (GOPT-A) — 128³ hartree −34%, v_eff −40%
- **XC grid + spin parallelization** (XCPR) — 1.12–1.24×
- **FFT buffer reuse** (FFTB) — 23–33% on `fft/scf_iter_20x_*`
- **Fused multiply-add** (FMAD) — −3–4% on `lda_xc_grid_*`
- **Alloc traffic** (ALOC F-5) — 16.8 GB → 0 per SCF
- **Test profile** (TPRF, TPRFB, TSPL) — cargo test 11 min → 12 s Tier-1 (≈55× warm)

## Infra + quality

- **Error handling** (ERRH, ERR2 P0/AX/P1.a/b/c/d) — Phase 1 closed; only P1.e 2-site `transplant.rs` mop-up remains
- **Module refactor** (MODR) — god-modules split, 4 phases
- **Type audit** (TYPE, TYPB, CAST) — i8 rotations, i16 Miller, ~148-site cast triage
- **Machine lock** (MLFX, QELK) — owner-scoped acquire, QE-aware policy, 17-case shell test suite
- **Observability** (PROF, MIXL, LOGH, LOGH-2) — samply canonical; `log::info!` on SCF summary
- **Doc hygiene** (MADOC-A/B/C, MOAD, MOAD-2, DCLN, CLSS, DWGT, RDOC) — module headers, `# Errors`/`# Panics`, rustdoc `-D warnings` gate
- **CI** (CICI, CINM) — clippy `-D warnings` on default + gpu; branch-protection check names
- **Conventions** (UNTS, DOCX, DOCLEAN) — eV/Å internal units, `RUSTDOCFLAGS` recipe, Rust-vs-Python rounding + QE `wk` × `degspin` gotchas

## Audits + grooming

- **RWHK** — reward-hacking audit (1 Critical, 4 Major, 6 Minor); fixes via RWHK-FIX
- **TAUD, TACC, TRV2** — test-suite audits (≥13 findings landed)
- **MAUD-AC** — maths accuracy audit (top 2 findings done; 7 items remain)
- **GRM6–12** — grooming passes + INDEX hygiene

## Notes still live (beyond Active row summaries)

- **ITEV** — only Phase-5 step-4 remains: end-to-end SCF wall-time bench with WFRX active to decide the default flip.
- **WFRX Technique 2** — ~30 lines once ITEV production-flip unblocks; cache plumbing landed with Technique 1.
- **Fe BCC PBE** — retains M=2.16 μB vs QE 2.34; 1.97 eV residual tracked under VGCH Class A.
- **Fe BCC LDA** — Hamiltonian-side outlier (BSUM ratio 0.99×) tracked under VGCH Class B → PZPW.
- **C diamond** — Class C → Class A light-atom reclassification; NLCC ρ_core(G) primary suspect → CNLC.
- **HD5I** (archived) — referenced the deleted `src/input.rs`; if re-opened, rewrite against YAML `Settings`.
