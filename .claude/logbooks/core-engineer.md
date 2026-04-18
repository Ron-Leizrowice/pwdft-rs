# Core Engineer Logbook

Entries: date, proposal ID, what was done, what remains, anything surprising. Keep it brief.

## 2026-04-18 — MXBA submitted (PR #57)

Branch `MXBA/adaptive-beta`, rebased onto origin/main (post-MODR-B/C/D + VNLM + CAST + DEAD).

**Implementation:** Eyert 1996 §3.3 residual-norm monitor as a private `AdaptiveBeta` helper in `src/scf/mixing/mod.rs`. Hooked into `AndersonMixer::push_history` (top, after Kerker preconditioning) and `BroydenMixer::mix` (same spot). `PeriodicPulayMixer` inherits via inner Anderson. `Mixer::new` gains `adaptive_beta: bool`; `Mixer::current_beta()` exposes β to the driver. New `KerkerSetup<'a>` bundle struct was added to keep mixer constructors inside the `too_many_arguments` limit (the 3 Kerker params now travel together).

**Defaults:** growth_threshold=1.2, damp_factor=0.7, restore_threshold=0.5, restore_window=3, β_min=max(0.05·β_start, 0.01).

**`adaptive_beta` default = false.** Empirical finding: on Fe CCMX (the existing `test_ccmx_fe_free_magnetization_converges`), adaptive β on damps β to β_min≈0.017 during the initial Δρ≈0.34 plateau (first ~5 iters) before Anderson has built DIIS history. Starved DIIS cannot escape the plateau → ConvergenceFailure at iter 80. Fixed-β 0.3 converges in 14 iters (unchanged baseline). Documented failure pinned by `tests/mxba_adaptive_beta_fe.rs` (`#[ignore]`).

**Numbers:**
- Fe CCMX fixed β=0.3 (default): 14 iters, |HF-KS|=1.06e-4 eV, M=0 μB — bit-identical to CCMX logbook entry.
- Fe CCMX adaptive β on: ConvergenceFailure after 80 iters, Δρ=0.341.

**Tests status:** 273 CPU + 282 GPU pass, 0 failed. Clippy clean both feature sets. 9 new unit tests (`AdaptiveBeta` logic + per-mixer driven-residual trajectory + backward-compat bit-identical).

**Surprises:**
- Fe CCMX failure under adaptive β was not anticipated by the proposal text. Root cause: the proposal assumed DIIS is always present to "do most of the work"; in practice DIIS needs `max_history ≥ 2` BEFORE adaptive β sees a residual trajectory it can interpret. Early-iter plateaus confuse Eyert's monitor. This is exactly why VASP/ABINIT/QE don't adapt β inside `mix_rho` — they do it in user scripts with explicit restart logic.
- `scf::report::IterationReport` is a clean place to pipe β into the per-iteration log; updating both drivers was one edit to `report.rs` + one `.current_beta()` in each driver's `log_iteration` call.
- Rebase hit expected conflicts: MODR-B split `scf/mod.rs` into driver/driver_spin; my inline edits moved to the new files trivially. INDEX.md conflict was 2-way (CAST had landed; I needed to keep CAST and append MXBA).

**Flagged for follow-up:**
- Tune Eyert thresholds on a metallic case where adaptive β helps (blocked on C diamond @ 30 Ry plain mixing converging). A fat follow-up proposal, not trivial.
- DIIS warm-up window (suppress monitor for first `max_history` iters). May be a quick mitigation of the Fe CCMX failure — if a warm-up window fixes the Fe regression, adaptive could become the default.
- VNLM proposal file is still in `proposals/` but the code landed as PR #49 (INDEX shows it in Active + Completed-equivalents). Core Engineer hand-off to EM to reconcile.

## 2026-04-18 — CAST submitted (PR #56)

Branch `CAST/numeric-cast-audit` from `origin/main` (post DOCS #53, TACC-I #54, DEAD #55 landings). Enabled the three `cast_*` correctness lints in `Cargo.toml`; walked ~148 hits.

- 44 `#[allow(..., reason = "...")]` annotations across 18 files (most function-level covering a cluster).
- 4 assertion-guarded rewrites: `fft_grid_size` (n_max >= 0), `NonlocalPotential::new` (proj.l >= 0 loop), `real_sph_harmonics` (debug_assert lmax >= 0). These catch malformed UPF / negative Miller bounds that would otherwise silently blow up allocations. Not "wrong physics" bugs — "silent hang/abort instead of clean panic" defensiveness wins.
- ~7 stylistic rewrites `(x as f64)` → `f64::from(x)` in kpoints + symmetry/kpoints.

Rebase mid-session: VNLM/DEAD/TACC-I landed while I was working. Stashed, fast-forwarded local branch to origin/main, popped stash — one conflict (fe_debug.rs deleted by TACC-I, accepted deletion). After rebase, 9 new CAST hits appeared from VNLM's rewrite of `src/potential/nonlocal.rs` — all `l as usize` / `lmax+1 as usize`; resolved as the assertion-guarded rewrites above.

**Quality gate:** clippy default+gpu clean (was 174 + 205 warnings before); `cargo test --release` 258/0/9, `--features gpu` 267/0/9.

**Flagged for follow-up (add to PR body too):**
- `src/pseudopotential/upf/convert.rs` accepts negative `angular_momentum` verbatim. The new assert in `NonlocalPotential::new` is a safety net; parse-time validation error would be cleaner. Small proposal for Code Reviewer / Researcher.

## 2026-04-18 — TACC-I submitted (PR #54)

Branch `TACC-I/ignore-reasons-and-fe-debug` from origin/main. TACC findings #2 + #3 (finding #1 deferred until MXBA lands).

**Finding #2 — 6 `#[ignore]` reason strings in `tests/qe_validation.rs`:** all six previously cited archived VERF. Ran `--ignored --nocapture` once to harvest current pwdft-rs values, then rewrote:
- C (Z=6) → SYKP: stalls at Δρ ≈ 4.1e-6 / 80 iters.
- Al (Z=13) → SYKP: ~73 meV residual on 8×8×8 (E = -64.197 eV vs -64.269 eV; barely fails 50 meV tol).
- GaAs (Z=31+33) → VGCMP: ~33.6 eV.
- Cu (Z=29 semicore) → VGCMP: ~16.2 eV.
- NaCl (Cl Z=17) → VGCMP: ~7.7 eV.
- MgO (Mg semicore) → VGCMP: ~10.1 eV.
2× SYKP / 4× VGCMP. File-header VERF narrative refreshed too. No test newly passes — nothing un-ignored.

**Finding #3:** deleted `tests/fe_debug.rs` (223 LOC, 6 tests — 5 weak/dead, superseded by VGCMP Phases 1-4). Migrated `test_fe_ewald_energy` → `qe_validation.rs::test_fe_bcc_ewald_vs_qe` (new Ewald section, <0.01 eV tol vs -171.779_065_80 Ry). Verified green.

**Quality gate:** clippy default+gpu clean, `cargo test --release` + `--features gpu` green, migrated Ewald passes.

**TACC proposal stays active** (finding #1 open — deferred until MXBA lands). No production-code changes.

**Surprises:** Al at ≈73 meV is closest-to-passing of all ignored QE validation tests; if SYKP/MPSH lands it will almost certainly pass. Flagged for EM in PR body.

## 2026-04-18 — MODR-D landed (PR #47)

Branch `MODR-D/upf-folder` from `origin/main@9a5e9e9`. Pure-move: `src/pseudopotential/upf.rs` (467 LOC) → `src/pseudopotential/upf/{mod,xml,convert}.rs`. Git detected `upf.rs → upf/convert.rs` as 86% rename.

**Split:**
- `mod.rs` (25 LOC): `pub fn parse` facade → `convert::parse_body`.
- `xml.rs` (66 LOC): `pub(super)` text helpers. No tests (every existing test goes through `parse()`).
- `convert.rs` (395 LOC): `pub(super) parse_body` + every `#[test]` verbatim (NCFX/NLCC pins included).

**Visibility tightenings:** three XML helpers went from file-private `fn` to `pub(super)` (same effective scope). `parse_body` new, `pub(super)`. Dropped now-dead `use crate::consts::{RY_TO_EV, BOHR_TO_ANG}` from `pseudopotential/mod.rs` — convert.rs imports them directly.

**Tests:** 209 lib + all integration pass (CPU); 212 lib + all integration (GPU). Clippy clean both feature sets.

**Surprises:** None. Only one cross-module call site (`pseudopotential/mod.rs:73 upf::parse(&content)`) was unchanged.

**Flagged for follow-up (Technical Writer):** stale `src/pseudopotential/upf.rs` prose in `CLAUDE.md:64`, `docs/units.md:34`, `docs/nonlocal.md:88`, `LOGBOOK.md:92`, `src/scf/energy.rs:26`, `src/scf/mod.rs:52`.

MODR: A ✅ (PR #46), D ✅ (PR #47). B (split `scf/mod.rs`) and C (symmetry/density folder) still open — concurrent agents likely in flight.

## 2026-04-18 — MODR-C submitted (PR #48)

Branch `MODR-C/symmetry-density-folder` from origin/main (post-MODR-A).
Pure-move refactor: `src/symmetry/density.rs` (~770 LOC) → `src/symmetry/density/{mod,real_space,g_space}.rs`. Git detected `density.rs → g_space.rs` rename (61% similarity).

**Layout:**
- `mod.rs` (~85): facade with `pub use real_space::symmetrize_density;` + `pub use g_space::symmetrize_density_g;` + shared helpers `check_grid_compatibility` / `compatible_grid_dims`.
- `real_space.rs` (~240): legacy `#[deprecated]` symmetrizer + `frac_to_grid_idx` + 7 tests.
- `g_space.rs` (~435): G-space symmetrizer + Miller helpers + 6 tests (incl. PCFX 18³ idempotency).

**No visibility tightenings.** Considered narrowing to `pub(super)` on the sub-module symbols and re-exporting as `pub`, but the simpler "`pub fn` at source, `pub use` at facade" pattern matches Phase A and keeps `pub` as the single source of truth. `symmetry/mod.rs` untouched — call sites use `density::<name>` which resolves equivalently through a folder-with-mod.rs.

**Tests:** 14/14 density tests pass, correctly routed under `symmetry::density::{real_space,g_space}::tests`. Full suite 209 CPU + 212 GPU pass, clippy clean both feature sets. 9 pre-existing `#[ignore]`.

**Surprises:** none. `#[allow(deprecated)]` on the facade's `pub use` needed explicitly — the re-export alone without the attr fires the deprecation warning through the module boundary.

## 2026-04-18 — MODR-A landed (PR #46)

Branch `MODR-A/split-mixing` from `origin/main@eb54d69`. Pure-move refactor: `src/scf/mixing.rs` (971 LOC) → `src/scf/mixing/{mod,anderson,broyden,kerker,linalg}.rs`. `AndersonMixer` + `PeriodicPulayMixer` co-located in `anderson.rs` (periodic wrapper pokes Anderson private fields).

**Visibility tightenings vs pre-refactor:**
- `Mixer` enum, `AndersonMixer`, `BroydenMixer`, `PeriodicPulayMixer`: `pub` → `pub(crate)`. Wanted `pub(super)` on the structs but `private_interfaces` fires because `Mixer` is crate-public and names them as variants — `pub(crate)` is the tightest setting that compiles.
- `precondition_residual`, `auto_q_tf_squared`, `solve_linear_system`: module-private → `pub(super)` (consumed by sibling files).
- `MixingMode` stays `pub` (serde adapter + tests/*.rs).

**Tests:** 209 lib + all integration pass (CPU); 212 lib + all integration pass (GPU). Both clippy invocations clean. Every `#[test]` from old `mixing.rs` moved verbatim — no test dropped, no test merged.

**Surprises:** None. `src/scf/mod.rs` needed zero changes (folder-with-mod.rs vs file is transparent at use-sites). Git detected `mixing.rs → mixing/anderson.rs` as a rename (59% similarity).

Phase B (split `scf/mod.rs` into driver/driver_spin/report) and Phase C (symmetry/density folder) remain.

## 2026-04-18 — CCMX landed (PR #43)

Branch `CCMX/coupled-channel`, rebased onto origin/main (post-PRPL/MODR/NCFX).

**Implementation** (`src/scf/mod.rs::run_scf_spin`): replaced `mixer_up`/`mixer_down` with `mixer_total`/`mixer_mag`. Forward basis change `(ρ↑, ρ↓) → (ρ_total, m)` before mix call; inverse after. Matches QE `rhoz_or_updw` (scf_mod.f90:1360-1414). Kerker disabled on `mixer_mag`; match handles all 4 MixingMode variants (Plain, Kerker, Broyden, PeriodicPulay) by flipping kerker off while preserving period.

**Key numbers — Fe BCC 4×4×4 nspin=2 free-mag starting_mag=0.5, 15 Ry, Kerker:**
- Pre-CCMX: Δρ limit cycle at 0.254 for 200+ iters, |HF-KS| ≈ 13 eV.
- Post-CCMX: **14 iters**, Δρ=5.7e-4, |HF-KS|=1.06e-4 eV, M=0 μB.

8×8×8 trace: energy stable to 6 decimals by iter 11; |HF-KS|=5e-5 eV; Δρ then enters a numerical-noise floor (spikes between 1e-8 and 3e-5 as Anderson history becomes rank-deficient — the existing singular-pivot guard handles this gracefully).

**New regression test:** `tests/spin_polarization.rs::test_ccmx_fe_free_magnetization_converges`. Fails hard if mixer topology reverts (|HF-KS|<1e-3 eV, M<0.05, iters<80).

**Tests status:** 191 lib + all integration CPU pass, 194 lib + all GPU pass. 8 pre-existing ignored. Clippy clean both `--all-targets` feature sets.

**Unchanged:** Si nspin=2 regression still sub-μeV |HF-KS| (basis change is exact). `test_fe_ferromagnetic_fixed_moment` still correctly ConvergenceFailure (fixed-mag=2 not a fixed point for this PP, regardless of mixer).

**Flagged for follow-up:**
- `tests/qe_validation.rs::test_fe_bcc_fm_vs_qe` now converges cleanly but stays `#[ignore]` pending VGCMP (9.5 eV heavy-atom V_local gap).
- Δρ wobble at numerical-noise floor in tight-tolerance 8×8×8 Fe runs. Energy + |HF-KS| pinned, so benign; MXBA adaptive-beta could damp Anderson history sooner. Not a correctness issue.

## 2026-04-18 — PRPL nits (PR #39, rebased on NCFX)

Two APPROVE-WITH-NITS follow-ups applied on top of PR #39:
- `ScfParams::validate()` rejects `MixingMode::PeriodicPulay { period: 0, .. }` with `PwdftError::InvalidInput`; new unit test `validate_rejects_zero_pulay_period` in `src/scf/mod.rs`.
- `periodic_pulay_vs_plain_scf_convergence` in `src/scf/mixing.rs`: silent `(Err, Err)` arm replaced with explicit `panic!` (default expectation is both SCFs converge on Si Γ-only).

Rebased cleanly onto `origin/main` (post-NCFX). Clippy clean on default and `--features gpu`; `cargo test --release` 191 pass + all integration green; `--features gpu` 194 pass.

## 2026-04-18 — NCFX landed (critical-path fix)

Branch `NCFX/nlcc-core-density-fix` (rebased onto `origin/main` post-VGC5).
Both compounding bugs fixed as diagnosed in VGC5:

1. **`src/pseudopotential/upf.rs`** — PP_NLCC conversion: `/BOHR_TO_ANG` → `/BOHR3_TO_ANG3`. PP_NLCC stores bare ρ_core(r) in e/Bohr³ (not 4πr²·ρ in e/Bohr); QE `rhoc_mod.f90:107` explicitly multiplies by r² in its Bessel transform, confirming the convention.
2. **`src/scf/potentials.rs::compute_core_density`** — Integrand now has `r²` weight and `4π` prefactor, matching QE `init_tab_rhc`.

**Key numbers:**
- Si diamond (ecut=15 Ry, 4×4×4): E_total −218.18 → −231.87 eV (QE: −231.61). Gap: **13.43 → 0.26 eV** (52× reduction). E_xc residual dropped from +13.74 eV to −0.31 eV.
- Fe BCC (nspin=1, 4×4×4): E_total −3101.24 → −3051.89 eV. E_xc residual: −48.85 → +0.69 eV.
- Si partial core charge: 0.7399 e (pinned, expected 0.74 e for Si ONCVPSP).

**Residual 0.26 eV on Si** attributed to Monkhorst-Pack shifted-vs-Γ-centered grid (SYKP): QE uses `4 4 4 0 0 0` (Γ-centered), pwdft-rs hard-codes shifted MP-1976. Γ eigenvalues differ by ~1 eV consistent with different k-meshes. Test `test_si_diamond_vs_qe` remains `#[ignore]`; ignore message updated to point at SYKP/MPSH.

**Tests updated:**
- New unit test: `test_si_core_charge_integrates_to_partial_core` in `src/pseudopotential/upf.rs`.
- GPU pins in `tests/gpu_consistency.rs`: Si total −198.8926 → −213.0283 eV; Si E_F 6.969 → 6.709 eV (same −14.1 eV shift as CPU).
- VGC5 pins in `tests/vgc5_per_component_si.rs`: Si + Fe pins refreshed; pre-NCFX baselines retained inline as comments.
- Removed `PRE-NCFX` label from module header; test now a post-NCFX regression guard.

**Tests status:** 180 lib + all integration pass (CPU and GPU). Clippy clean on `--all-targets` and `--features gpu --all-targets`.

**Flagged for follow-up:**
- Si VGC5 self-check `Σ(components) − E_total = 1.18 eV` — was present pre-NCFX too (PCRS territory).
- Si one-electron residual vs QE: +2.00 eV (was +2.38 pre-NCFX). Small improvement; tracked by VGCMP/PCRS.
- Si E_hartree residual vs QE: −0.79 eV (was −1.51 pre-NCFX). Same — improvement but non-zero.
- Fe residual +8.27 eV on E_total: MP-shift + ecut convergence. Not NCFX domain.
- `tests/qe_validation.rs::test_fe_bcc_fm_vs_qe` (nspin=2 Kerker 8×8×8) does not converge in 80 iters post-NCFX. Pre-NCFX it passed at −3059.44 eV; this is likely CCMX-class behavior (independent ↑/↓ Anderson can't cope with the newly-exposed magnetic landscape). Reason on `#[ignore]` still fine; worth noting that Fe now needs CCMX to converge cleanly rather than benefiting from accidental XC cancellation.

## 2026-04-16 — Orientation + KBTF (PR #1)

- DDUP/SIMP/VERF/HRFK proposals verified accurate; ERRH needed line refresh; CFGN blocked on DDUP+SIMP.
- PR #1 `KBTF/kb-test-failures`: fixed 3 failing `tests/kb_projector_validation.rs` — test_09 rewritten (test bug: wrong D_ij assumption for UPF vs raw HGH h^l), test_07 `#[ignore]` (trapezoidal artifact, SIMP territory), test_vloc `#[ignore]` (bare Coulomb, VERF territory). Suite: 9 pass, 2 ignored.
- Pre-existing on main: Si energy 13.33 eV off, C diamond non-convergence (qe_validation.rs).

## 2026-04-16 — CLEN (PR #2)

5 cleanups: removed GPU staging double-copy; replaced `n_projectors`/`has_nlcc` fields with methods; dropped unused `_z_val` from `add_atomic_density_from_pp`; deduped `apply_rotation` in detect.rs (now uses `SymmOp::apply`). Net: 9 files, +31/−41, 159 tests pass.

## 2026-04-16 — HRFK (PR #6)

Harris-Foulkes energy diagnostic added. `harris_foulkes_energy()` in `energy.rs`; E_HF computed in both spin and non-spin SCF loops from INPUT-density quantities; new `ScfResult.harris_foulkes_energy` field; warning if `|HF-KS| > 0.01 eV` at convergence. Flagged subtle spin-path inconsistency (input `exc_r` mixed with output `rho_xc_total`) — later addressed by SPXC.

## 2026-04-16 — DDUP (commit e09513e)

5 dedup steps in SCF: `assemble_v_eff()` replaces inline V_eff build in spin loop (gains `par_iter`); FFT normalisation delegates to `density_r_to_g()`; `real_to_g_space` delegates likewise; `rho_core_half` hoisted out of spin loop; `compute_occupations()` helper replaces 5 identical blocks. Net: +49/−42 across 3 files.

## 2026-04-16 — SIMP (PR #7)

Simpson quadrature at all 4 radial integral sites via new `src/numerics.rs::simpson_integrate()` (matches QE `simpsn.f90` including even-mesh correction). Sites: V_local form factor, beta projectors, NLCC core density, SAD. Key outcome: **Fe BCC QE validation NOW PASSES** (was 45.4 eV off); Si still 13.43 eV off (remaining error not in quadrature). Un-ignored test_07 at threshold 0.12. VERF unblocked.

## 2026-04-16 — BROY (PR #8)

Modified Broyden second method (Johnson PRB 38 12807) added as `BroydenMixer` + `Mixer` enum dispatch. New `MixingMode::Broyden { kerker }` and YAML keys `broyden`/`broyden_kerker`. 8 new tests pass. Si energy matches Anderson within 0.01 eV. Adaptive-beta (Step 3) and periodic-Pulay (Step 4) deferred.

## 2026-04-17 — VERF/SPXC/QEVL agents aborted (isolation failure)

Three parallel agents escaped their worktrees and wrote to main checkout. Session reset to main. Key findings preserved in proposals:
- **VERF** erf subtraction does NOT close Si 13.4 eV gap (numerically equivalent to bare-Coulomb + Simpson). Critical flag lowered; follow-up `VGCMP` opened.
- **SPXC** fix drafted (~20-line change). Test convergence trouble at conv=1e-7; rerun at 1e-6.
- **QEVL** QE reference data generated for 8 Tier 1+2 systems (`/tmp/pwdft-rescue/qe_validation_data/`).

Isolation enforcement added: agents must `pwd`-check worktree at session start.

## 2026-04-17 — VERF finalized + VGCMP opened

Branch `VERF/vloc-erf-finalize`. Decision: **LAND** erf subtraction even though numerically a no-op on Si/Fe today.

- `src/pseudopotential/mod.rs` — replaced bare-Coulomb `G≠0` branch with QE's erf decomposition `[r·V(r) + Z·e²·erf(r)]·sin(Gr)/G` and analytic correction `−4π·Z·e²·exp(−G²/4)/(Ω·G²)`. G=0 unchanged.
- New `tests/vloc_erf_consistency.rs` — 2 tests. First 20 Si |G| shells: `max |Δ| = 5.91e-9 eV`, `max relative = 6.4e-6` — confirms erf ≡ bare-Coulomb on our log mesh post-SIMP. Clean regression guard.
- `proposals/VERF-vloc-erf-subtraction.md` → `proposals/completed/` with `outcome: landed-as-cosmetic`.
- New `proposals/VGCMP-vloc-g-cross-check.md` — Researcher-owned follow-up. Plan: Python reference for `V_local(G)` (Phase 1, 1 day), β_l(q) (Phase 2), D_ij (Phase 3), single Hamiltonian element (Phase 4). Phase 1 almost certainly isolates the 13.4 eV Si gap.
- `proposals/INDEX.md` — VERF moved to Completed, VGCMP added as critical (dependency of QEDX, QEVL).

**Rationale for landing VERF:** matches QE convention exactly (simplifies VGCMP Phase 1 by removing one axis of variation), bounded integrand at r=0 is an insurance policy for future high-Z PPs where bare-Coulomb hits billions of eV near the origin, and the regression test pins the mathematical equivalence so any future drift fires loudly.

**Numbers:**
- Si QE validation: still −218.18 eV (13.43 eV above QE) — unchanged. C diamond still doesn't converge. Both pre-existing; VERF was never expected to fix them.
- Fe BCC QE validation: still passes (−3059.44 eV).
- 10/11 kb_projector tests pass, 1 ignored (VGCMP territory).
- `vloc_erf_consistency`: 2/2 pass.
- Clippy clean.

**Isolation footgun:** this session initially wrote 4 files to the main checkout via the Write/Edit tool with `/Users/ronleizrowice/Documents/github/pwdft-rs/...` paths — both paths pointed at the main checkout, not the worktree. Caught it after `cargo test` reported "no test target" because the worktree didn't have the file. Cleaned up main checkout and re-applied to the worktree correctly. Lesson: **always verify Write/Edit paths resolve inside `.claude/worktrees/agent-*` before using them.**

## 2026-04-17 — SPXC implemented (attempt 2)

PR on `SPXC/spin-xc-consistency` branch. Implements the exact fix from the proposal: after `density_r_to_g`, recompute `(exc_r_out, vxc_up_r_out, vxc_down_r_out)` via `xc::lda_xc_spin_grid` from OUTPUT spin densities + half-core, use them for E_KS, keep INPUT-derived quantities for E_HF. 13-line logical change in `src/scf/mod.rs::run_scf_spin`.

**Empirical |E_HF - E_KS| for Fe BCC (4x4x4, 15 Ry, fixed mag=2, starting_mag=0.5, conv=1e-6, 300 max_iter):**
- Pre-fix: 22.24 eV  (E_KS=-3108.67 eV, 238 iters)
- Post-fix: 13.03 eV (E_KS=-3099.46 eV, 244 iters)
- Improvement: ~1.7x (not the 10x the task hoped for)

**Why not 10x:** nspin=2 convergence check uses only rho_total, not per-spin. Spin density (zeta) is never driven to self-consistency, so E_xc[rho_in, zeta_in] vs E_xc[rho_out, zeta_out] differ at ~13 eV level regardless of how long SCF runs. The SPXC fix removes the artificial extra ~10 eV from mixing zeta_in with rho_out, but the residual zeta-inconsistency gap remains. Worth a follow-up proposal: convergence criterion should include per-spin density difference.

**Regression test:** `test_fe_spin_xc_consistency_regression` in `tests/spin_polarization.rs`. Uses the Fe fixed-mag=2 setup with starting_mag=0.5 (required — without it Fe doesn't converge) and max_iter=300 (244 needed). Asserts `|HF-KS| < 18 eV`, which fails pre-fix (22.24 eV) and passes post-fix (13.03 eV). Takes ~3s in release mode.

**Side notes:**
- Existing `test_fe_ferromagnetic_fixed_moment` still does not converge (delta=8.46e-2 after 100 iters); its assertions are gated on `Ok`, so test passes vacuously. Pre-existing, unchanged by this fix.
- `qe_validation.rs` Si/C tests fail on main too; not caused by SPXC.
- Clippy clean, 173 lib + 44 integration (excl. qe_validation) tests pass.
- Tangential: nspin=2 convergence metric should include `density_diff(rho_up_in, rho_up_new)` and `density_diff(rho_down_in, rho_down_new)`, not just rho_total. Worth a proposal.

## 2026-04-17 — SPNC implemented

Branch `SPNC/per-spin-convergence`. Proposal + implementation + test update.

**Code change** (`src/scf/mod.rs::run_scf_spin`): replace the total-density delta at line 639 with `max(density_diff(rho_up_r, rho_up_sym), density_diff(rho_down_r, rho_down_sym))`. Added per-channel delta to info log (`Δρ=... (↑... ↓...)`) so future debugging sees channel behavior directly. Also added `env_logger::builder().is_test(true).try_init()` to the regression test so `RUST_LOG=info cargo test -- --nocapture` shows iterations.

**Empirical surprise — Fe fixed-mag=2 does NOT converge under per-spin criterion:**
- Pre-SPNC (total-only, SPXC in place): "converged" at iter 244 with Δρ_total=1e-6, |HF-KS|=13.03 eV.
- Post-SPNC (per-spin max): **ConvergenceFailure** — both channels pinned at `Δρ_up = Δρ_down = 0.254` from iter ~5 onwards (steady-state limit cycle). Total density stable (dE~3e-7), so the previous "convergence" was a +ε/−ε spin flip cancelling into the total.
- Root cause: fixed-mag=2 is not a stable SCF fixed point for this PP (LDA Fe ground state is non-magnetic; test comment already noted "this PP favours non-magnetic Fe"). Independent Anderson mixers on up/down cannot co-ordinate the inter-channel charge transfer that would close the gap. Fixing this needs a coupled-channel mixer or a different PP — out of SPNC scope.

**Regression test reframed:** repurposed `test_fe_spin_xc_consistency_regression` to **Si nspin=2** (starting_mag=0.2, free mag, 100 eV ecut). Si is non-magnetic, relaxes to M=0, and the spin XC machinery is fully exercised during SCF iterations. Demonstrates true O(Δρ²) quadratic convergence with SPXC+SPNC in place:
- **|HF-KS| = 7.19e-7 eV at conv=1e-6, 23 iters, M=0.0**.
- Assertion threshold: 1e-5 eV (~14× empirical, platform headroom).
- Test comment documents the Fe fixed-mag history and why it's now an unsuitable regression target (pathological mixer oscillation).

**Tests:** 220 pass, 9 ignored (pre-existing QE validation + kb_projector vloc), 0 fail. Clippy clean.

**Follow-up candidates (not in this PR, good fodder for future proposals):**
- Coupled-channel nspin=2 mixer (mix `(ρ_total, m)` instead of `(ρ_up, ρ_down)`) — QE does this via `mix_rho` which takes the full nspin-component vector as one residual. Would unblock Fe fixed-mag.
- Revisit `test_fe_ferromagnetic_fixed_moment` — passes vacuously (gated on `Ok`, SCF always returns Err). Either remove or swap to a ferromagnetic-stable PP.
- Adaptive spin mixing `mix_beta` separate from charge `mix_beta` — standard QE setting, damps spin oscillations independently of charge.

## 2026-04-17 — TAUD-A implemented

Branch `TAUD-A/silent-pass-fixes`. PR A of the test-quality audit: replaced 5 silent-pass `match`/`if let` patterns with hard `.expect()` + convergence guards.

**Files modified** (tests only, no `src/` changes, no tolerance changes):
- `tests/spin_polarization.rs` — finding 1.1 (`test_si_nspin2_matches_nspin1`): 4-arm match → `.expect()` + `n_iter < max_iter` guard for both nspin=1 and nspin=2.
- `tests/parallel_consistency.rs` — findings 1.3 + 1.4: `if let (Ok, Ok)` → `.expect()` + convergence guards for both `test_scf_serial_vs_parallel` and `test_scf_kerker_serial_vs_parallel`.
- `tests/gpu_consistency.rs` — finding 1.5: `(Err, Err)` arm now `panic!`s instead of `eprintln!`; `(Ok, Ok)` arm gained 5.2 convergence guards. Finding 1.7: `test_gpu_scf_kerker_converges` match → `.expect()` + convergence guard.

**Pre-existing compile fix incidental to 1.7:** line 440 of `gpu_consistency.rs` built an `ScfParams` literal that was missing the `energy_threshold`, `smearing_scheme`, `nspin`, `starting_magnetization`, `tot_magnetization` fields added since this test was written — the whole GPU test binary wouldn't compile. Added `..Default::default()` so the test compiles and runs. Documented in PR body.

**Test results post-fix (all pass strictly, no bugs unmasked):**
- `test_si_nspin2_matches_nspin1` — PASS
- `test_scf_serial_vs_parallel` — PASS
- `test_scf_kerker_serial_vs_parallel` — PASS
- `test_gpu_vs_cpu_scf_direct_comparison` — PASS (GPU path)
- `test_gpu_scf_kerker_converges` — PASS (GPU path)

Full suite: 220 CPU tests pass + 6 GPU consistency tests pass, 9 ignored (pre-existing VERF/VGCMP-gated). Clippy clean on `--all-targets`; `--features gpu` surfaces pre-existing warnings only.

**Rebase note:** branch rebased onto origin/main (ea1a300) before final commit to include XCPR merge.

**No follow-up proposals opened** — none of the 5 hardened tests revealed a regression. The silent-pass patterns were pure hygiene debt.

## 2026-04-17 — XCLN cleanup (docs + proposal paperwork)

Branch `XCLN/cleanup-ccmx-qedx-sykp`. Three bundled items, no behavioural code changes.

**1. CCMX proposal opened** (`proposals/CCMX-coupled-channel-mixer.md`). Captures the SPNC follow-up identified on 2026-04-17: independent `(ρ↑, ρ↓)` Anderson mixers can't drive Fe BCC fixed-mag=2 to a converged fixed point. Plan: mix `(ρ_total, m)` instead, mirroring QE's `rhoz_or_updw` basis change. Cited QE sources verified:
- `qe-7.5/PW/src/sum_band.f90:307` — `rhoz_or_updw(rho, 'r_and_g', '->rhoz')` called right after band-sum density accumulation.
- `qe-7.5/PW/src/scf_mod.f90:1360-1414` — the basis-change subroutine itself (`vi=1.0` forward, `vi=0.5` inverse).
- `qe-7.5/PW/src/v_of_rho.f90:320-360` — local back-conversion for XC evaluation.
`depends_on: [SPNC]`, priority high, complexity medium, risk medium. Mixer internals remain channel-agnostic — only `run_scf_spin` mixing step is touched.

**2. QEDX archived** (moved to `proposals/completed/QEDX-qe-energy-discrepancy.md`, status → `superseded`). Added completion note at top citing SIMP + VERF as done and VGCMP as the active investigation. INDEX.md updated: removed QEDX row from Critical (only VGCMP remains there), added to Completed table, updated the `## Notes` block.

**3. SYKP D1 done**. Updated 3 docstrings + 1 assertion comment:
- `src/kpoints.rs::monkhorst_pack` — now explicitly names the shifted MP-1976 convention, cites QE equivalence (`k1=k2=k3=1`, not default `0 0 0`), and points at the SYKP audit.
- `src/symmetry/kpoints.rs::reduce_kpoints` — convention section explaining the 10-vs-8 IBZ reduction story.
- `src/symmetry/kpoints.rs::mp_fractional` — same convention note, must-match callout to `kpoints::monkhorst_pack`.
- `src/symmetry/kpoints.rs` test comment for `test_si_4x4x4_reduces_to_8` — replaced the misleading "incomplete boundary handling" prose with the correct convention-mismatch explanation. Assertion range (`8..=10`) unchanged. Appended "2026-04-17 — D1 done" block to `proposals/SYKP-symmetry-ibz-audit.md` (the file still lives in active dir alongside INDEX.md's "Low/Deferred" row — task spec referenced `proposals/completed/` but SYKP hasn't been formally archived yet; I appended to its actual location).

**Quality gate:** `cargo test` all pass (spin_polarization 3/3 including `test_fe_spin_xc_consistency_regression` Si nspin=2, plus VGCMP/VLOC/VERF cross-checks). `cargo clippy -q --all-targets` clean. No behavioural code paths touched — purely docstring + test comment.

**Surprises:** none. The SPNC logbook entry already flagged CCMX as a follow-up candidate with the correct QE reference (`PW/src/mix_rho.f90`-via-`rhoz_or_updw`); I just verified the exact line numbers and sketched the scope. SYKP's `proposals/completed/` path was a minor spec mismatch (the file is still in active/`proposals/`) but doesn't affect correctness.

## 2026-04-17 — FFTB implemented

Branch `FFTB/fft-buffer-reuse`. `FFT3D` now owns two reusable `Array3<Complex64>` scratch buffers (`buf_a`, `buf_b`), allocated once in `new()`. `forward()`/`inverse()` no longer allocate per call — they `copy_from_slice` into `buf_a` and ping-pong through the three 1D transforms.

**Parallel-sharing risk:** investigated and cleared. Only two places run in parallel relative to FFTs:
- `density::compute_density` — already constructs a `FFT3D::new` per rayon worker inside `fold`, so each thread gets its own buffers. No restructuring needed.
- k-point eigensolve loops in `scf/mod.rs` — these don't touch `ctx.grid.fft` at all (pure linalg + `vnl_cache` read).

All `ctx.grid.fft` usages are sequential in the outer SCF loop. Signatures `&mut self` on `forward`/`inverse` were already present; only the implementations changed.

**Perf (criterion `fft/scf_iter_20x_NxNxN`, 20 forward+inverse on same FFT3D):**

| Grid | Before | After | Delta |
|---|---|---|---|
| 16³ | 1.60 ms | 1.07 ms | −33% |
| 20³ | 6.36 ms | 4.53 ms | −29% |
| 24³ | 5.22 ms | 3.76 ms | −28% |
| 32³ | 20.5 ms | 14.3 ms | −30% |
| 48³ | 90.7 ms | 69.8 ms | −23% |

Matches proposal's allocation-rate argument: per-call we save `data.to_vec()` (512 KB at 32³) + `Array3::zeros` (another 512 KB with zero-fill). Speedup consistent with avoided memcpy + page-fault traffic.

**Tests:** full `cargo test` green (173 unit + all integration passing). Clippy clean on `--all-targets`.

**Minor observation (not blocking):** during this session, in-progress edits to `src/fft.rs` and `benches/scf_benchmarks.rs` were silently reverted between bash invocations (the compiled bench binary lacked the newly-added `fft_scf_iteration` group despite the source having it minutes earlier). Re-applied the edits, verified hashes before each `cargo` invocation, and confirmed via `--list` that the bench binary contains the new symbols. Cost ~10 min. Worth investigating the worktree/hook setup separately.

## 2026-04-17 — TAUD B+C+D+E bundled (PR #26)

Branch `TAUD/prs-bcde-bundle`. Bundled remaining TAUD fix PRs (B/C/D/E) since all touch only `tests/` and are mechanical.

**Empirical margins read under lock (each test ran once to set data-driven tolerances):**
- 2.1 Si nspin=1 vs nspin=2 energy diff: ~0 (prints as "0.000000") → threshold 1e-5 eV (was 0.5)
- 2.2 Si M: ~0 (prints as "0.000000") → threshold |M| < 1e-4 μB (was M < 0.1, unsigned!)
- 2.4 GPU Hartree max rel err: 1.8e-7 → threshold 1e-5 (was 1e-3)
- 2.5 GPU Si E: -198.892595 eV → `|E - (-198.8926)| < 0.1` (was `∈ [-300, -100]`)
- 2.6 GPU Si E_F: **6.968947 eV** (not 5.97 as my first guess — all 4 bands occupied, E_F sits ~1 eV above HOMO) → `|E_F - 6.969| < 0.1` (was `∈ [-5, 10]`)
- 2.7 5σ Fermi tail: documented (5σ ≈ 0.67% FD occupation), no numeric change
- 2.8 Kerker GPU Si E: -198.892588 → same as 2.5

**PR B (1.2) — Fe test inverted.** Now asserts `matches!(result, Err(ConvergenceFailure { .. }))`. Empirical hit: 100 iters, delta=0.0846 (matches SPNC logbook line 148 exactly). `ScfResult` has no `Debug` derive, so used a 3-arm `match` with explicit eprintln/panic instead of `{result:?}` in assert format strings. Same idiom in PR E.

**PR D (3.1) — Real finding surfaced.** Un-ignored `test_vloc_comparison_with_qe` per plan; it **FAILED** post-VERF:
- V_local(G=0): ours +1.343 eV vs QE -1.003 eV (sign flip, diff 2.35 eV)
- |V_local(G=(1,0,0))|: ours 5.468 vs QE 6.968 (diff 1.50 eV)
- |V_local(G=(1,1,1))|: same 1.50 eV

This is the same class of convention mismatch VGCMP Phase 1 is targeting — re-ignored with specific numbers in the `#[ignore = "..."]` reason string and pointer to VGCMP Phase 1. **Did not open a new proposal** — would duplicate VGCMP. Listed in PR body "Flagged for follow-up" for EM visibility.

**PR E (5.3) — ConvergenceFailure variant.** Replaced `.is_err()` in `parallel_consistency.rs::test_scf_serial_vs_parallel` with explicit 3-arm match (ConvergenceFailure OK, Ok panics, other Err panics with `{other}`). Inline to avoid `clippy::items_after_statements` from a nested fn.

**Results:** 222 CPU pass + 231 GPU pass, 9 ignored (8 qe_validation + 1 kb_projector vloc), zero failures, clippy clean on both `--all-targets` and `--features gpu --all-targets`. Remaining warnings (benches/gpu_benchmarks.rs deprecated `black_box`, uninlined format args) are pre-existing, tracked under QLN2.

**Clippy gotchas encountered:** `const FOO: f64 = ...` inside a fn after let-bindings fires `clippy::items_after_statements`; switched to `let foo = 1.23_f64;`. Also `fn expect_...()` helper inside the test triggered same — inlined the match.

**Handoff notes for next session:**
- `test_vloc_comparison_with_qe` will need attention when VGCMP Phase 1 resolves the V_local convention. At that point, un-ignore and expect tolerance to need tightening from the current 0.5 eV (which was a wish at write-time, not empirical).
- `ScfResult` doesn't derive `Debug`. If future tests want to print Result values, either derive Debug on ScfResult (see `src/scf/mod.rs:131`) or use the 3-arm-match idiom from this PR (spin_polarization.rs Fe test, parallel_consistency.rs).
- Machine lock held 1007s total across baseline + final-test runs.

## 2026-04-17 — PRPL implemented

PR #39 created: `PRPL/periodic-pulay`. Periodic Pulay mixer (Banerjee et al., JCTC 12, 3053 (2016)) as a BROY follow-up.

**Implementation:**
- Refactored `AndersonMixer::mix` into `push_history` + `diis_step` (no behavior change; `mix` is now a thin wrapper). `push_history` applies Kerker preconditioning if enabled, appends to history, trims to `max_history`. `diis_step` requires `history_len >= 1`, does linear mixing for the first iteration then solves the DIIS system for iter ≥ 2.
- `PeriodicPulayMixer` wraps `AndersonMixer`, calls `push_history` unconditionally, and gates `diis_step` on `iteration.is_multiple_of(period) && history_len() >= 2`. On non-Pulay iterations it does `ρ + β·R` against the most recent (already-preconditioned) residual stored in history — this keeps Kerker preconditioning live on linear steps.
- `MixingMode::PeriodicPulay { period, kerker }` threaded through `Mixer::new/mix`. New `MixingModeType::{PeriodicPulay, PeriodicPulayKerker}` + `pulay_period: usize` (default 3) on `ElectronSettings`. New `MixingModeType::to_scf_mode(period)` helper; back-compat `From` uses default period 3.
- `src/scf/mod.rs` and `src/scf/energy.rs` untouched (VGC5 safe).

**Convergence numbers (Si Γ-only, ecut=100, 16³ grid, conv=1e-6):**
- Plain:          10 iters
- PeriodicPulay:   8 iters (period=3, ΔE≈1.15e-6 eV)

20% iteration reduction on a small insulator. Paper's strongest wins are metals/TMOs; this is still in the right direction.

**Tests added (10):** 8 in `src/scf/mixing.rs` (period=1 matches Anderson bit-for-bit, period=∞ matches plain linear bit-for-bit, history accumulation, Pulay fires on iter 3 of period=3, first-iter is linear, Kerker-finite, synthetic fixed-point convergence, Si SCF smoke test), 2 in `src/settings.rs` (YAML parse for `periodic_pulay` and `periodic_pulay_kerker`).

**Results:** 189 CPU + 192 GPU unit tests pass (+10 from PRPL), all integration tests pass, clippy clean on both feature sets. Rebase onto origin/main was clean (VGC5 had landed; no conflict since both PRs touch disjoint files).

## Flagged for follow-up
- `PeriodicPulayMixer::mix` non-Pulay path reaches into private fields `anderson.{beta, history_in, history_res}`. Fine within the module, but future refactors may want to move this into `AndersonMixer::apply_linear_step_from_last_history` to localize knowledge. Not worth a proposal on its own — notional cleanup.
- Paper recommends `period = 5–8` for metals; our default is 3 (insulator/semiconductor sweet spot). Consider auto-selecting based on system class (gap-detected via initial-density sloshing amplitude, or explicit `system.metallic: bool`). Likely a new proposal if anyone has a concrete metallic test case to tune against.
- No real metallic regression test for PRPL — the existing Γ-only Si test is an insulator. A follow-up might add a BCC Fe SCF test with PRPL vs Broyden, but BROY already covers Fe well enough that the marginal value is unclear.

## 2026-04-18 — PCFX landed (PR #44)

Branch `PCFX/g-space-symmetrization`. Moved density symmetrization from real-space (rounding-sensitive on non-symmorphic grids) to G-space (exact via phase factors).

**Convention — verified against QE line-by-line, worth pinning here:**

- `SpaceGroupOp::rotation` is the fractional-direct-space rotation `R`: atoms transform as `r' = R·r + τ`.
- Under the pullback `(S·ρ)(r) = ρ(S⁻¹ r)`, Miller indices rotate as `n → R^T · n` (NOT `R⁻¹`, NOT `R^{-T}`).
- Phase: `exp(-i·2π·n_dst·τ_S)` using the DESTINATION Miller, not source.
- `P² = P` proven analytically via `m·τ_{S₁·S₂} = m·τ_{S₁} + (R_{S₁}^T m)·τ_{S₂}` cancellation. Numerically verified on 12³.
- QE stores `s(:,:,ns)` as the *transpose* of our `R` (proof: `symm_base.f90:533` atoms rotate as `rau = s^T · xau`). So QE's `s(:,:,invs(ns))·g0 = R^{-T}·g0` matches our R^T under the `S → S⁻¹` relabel.

**Key numbers (Si, ecut=15, 4×4×4):**
- Per-component self-check: **1.204 eV → 3.5e-11 eV** (10 orders of magnitude; target was 1e-5 eV).
- E_total: −231.8653 → −231.8429 eV (23 meV; matches proposal's 17 meV estimate).
- Fe Im-3m (τ=0): unchanged as expected.

**Non-obvious landmine (hit this; avoid in future work):** the G-space formula ALSO requires the grid-compatibility condition `N·τ ∈ ℤ` UNLESS the input density is band-limited away from the Nyquist. Derivation: the projector proof uses `miller_to_flat(R^T m)` to reduce mod N. If R^T m wraps (i.e. `R^T m ∉ [−N/2, N/2]`), the coefficient read corresponds to Miller `k = R^T m − N·δ` for some integer `δ`. The inner application's phase `exp(-i·2π·k·τ)` differs from `exp(-i·2π·R^T m·τ)` by `exp(+i·2π·N·δ·τ)` which is 1 iff `N·τ ∈ ℤ`. For our SCF use case the density is band-limited (|G|² ≤ 4·ecutwfc; FftGrid sized with ≥ factor-2 margin), so rotations don't wrap. Unit tests had to use band-limited inputs; the naive `sin(i·0.37).abs()` generator broke idempotence by 5% on 18³ until I synthesized modes with |m| ≤ 2.

**Handoff notes:**
- `src/symmetry/density.rs` has an extensive docstring on the convention and band-limitation requirement. Read before modifying.
- GPU-resident symmetrization is out-of-scope (proposal explicit); density-grid FFT stays CPU-serial even with `--features gpu`. If someone wants to move it, they'd need to expose the FFT buffer to a GPU kernel and redo the phase sum in WGSL.
- The PCFX self-check now pins `|Σ − E_total| < 1e-5 eV` in `vgc5_per_component_si.rs`. Any regression in symmetrization (or the density reconstruction path) will fail this aggressively.

**Machine lock:** ~900s total (compile + tests + clippy + GPU tests).

## 2026-04-18 — MODR-B landed (PR #50)

Pure-move refactor of `src/scf/mod.rs`. Branch `MODR-B/split-scf-mod`.

**Split:**
- `run_scf` (non-spin hot loop) → `scf/driver.rs::run_scf_unpolarized` (pub(crate))
- `run_scf_spin` → `scf/driver_spin.rs::run_scf_spin` (pub(crate))
- logging helpers → `scf/report.rs` (IterationReport + log_iteration + log_convergence_summary + log_entropy + log_components, all pub(super))
- `EnergyComponents` → `scf/energy.rs`; re-exported via `pub use`
- Helper unit tests (real_to_g_space, assemble_v_eff, hartree_on_fft_grid, density_diff) moved to `scf::energy::tests` with their production code. `validate_rejects_zero_pulay_period` stays with `ScfParams` in mod.rs.

**mod.rs: 1292 → 271 LOC.** `mod driver`, `mod driver_spin`, `mod report` all fully private (tighter than Phase A's `pub mod mixing`).

**Landmine avoided:** first draft of report.rs collapsed both drivers into one `log_convergence_summary(..., ts, n_atoms)`. The original spin driver never emitted `Entropy (-TS):`; non-spin emits it when `|TS| > 1e-8`. That would have been a behavior change. Split into `log_convergence_summary` (both drivers) + `log_entropy` (non-spin only). Pinned by `log_entropy` docstring.

**Results:** 209 CPU + 212 GPU unit tests pass, all integration tests pass, clippy clean on both feature sets. Rebased onto origin/main after MODR-C (#48) and MODR-D (#47) landed mid-session — zero conflicts (disjoint file sets, as predicted in the phase plan).

## Flagged for follow-up
- The `+ ctx.v_local_g0 * ctx.n_electrons` G0-shift expression is now duplicated in `scf/driver.rs` (e_total + e_harris) and `scf/driver_spin.rs` (e_total + e_harris) — 4 call sites across 2 files. A `with_g0_shift()` helper in `scf/energy.rs` is the natural DRY. Proposal MODR flags as Core Engineer follow-up; did NOT do it here per explicit scope fence.
- `SpinIterationFields` in `scf::report` is currently pub(super) and could stay that way, but if a third driver variant ever lands (non-collinear spin? DFT+U?), the spin-extension pattern of "Option<SpinIterationFields>" on IterationReport will not generalize cleanly. Revisit if/when.
