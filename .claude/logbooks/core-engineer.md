# Core Engineer Logbook

Entries: date, proposal ID, what was done, what remains, anything surprising. Keep it brief.

## 2026-04-18 — HKIN: drop unused Option V_eff param (PR #123)

Deleted `build_hamiltonian` from `src/hamiltonian.rs` (29-line fn with dead `Option<&dyn Fn(usize, usize) -> Complex64>` branch). Kept `build_kinetic` as the only constructor. Dropped the matching `Option` from `compute_band_structure`. Touched 5 files, -26 net LOC.

**Mechanical refactor, zero behavior change** — the 12 free-electron band tests pass with byte-identical eigenvalues because the `Option` branch was inert at every one of the ~13 `None` call sites.

**Drive-by doc fix:** `src/scf/potentials.rs:151` docstring still pointed at the removed `build_hamiltonian`. Retargeted at `build_kinetic` (which is what the comment was really contrasting against anyway — "kept separate from" the kinetic-only public path).

**GitHub push footgun:** my first commit used `user.email=ron.leizrowice@protonmail.com` from the worktree's local git config, which GitHub's email-privacy protection rejects with GH007. Fixed with a one-shot `git -c user.email=174235040+...@users.noreply.github.com commit --amend --reset-author --no-edit` — did NOT touch global or repo config (per the "never update git config" rule). Future note: the main-branch recent commits use either `ron@pelanor.io` or the noreply form; if the Bash PreToolUse hook ever blocks the amend, fall back to the pelanor.io email.

Gate: clippy default 17 / gpu 23 (baseline unchanged), rustdoc clean, all tests pass.

## 2026-04-18 — MLFX: machine-lock hardening (PR #101)

Four bugs in `.claude/bin/{machine-lock,check-cargo-lock.sh}` flagged by WFRX #99. All fixed + 17-case shell test suite at `.claude/bin/tests/machine-lock.test.sh`. Shellcheck clean.

**Atomic primitive:** `mkdir` (not `flock`). BSD flock on macOS has different semantics from Linux's `flock(2)`; mkdir is POSIX-atomic on both.

**Lock format change:** flat file → directory `machine.lock.d/{agent,desc,ts,pid,worktree}`. Legacy flat-file locks detected by `_has_legacy_lock` and treated as stale (cleared on next acquire).

**PPID scoping landmine.** First draft used `$$` (the machine-lock script's own PID) as the "owner PID" — but the script exits immediately after `acquire`, so the PID looks dead to every subsequent caller, and the lock auto-stales. Fix: record `$PPID` (the parent shell that invoked machine-lock). The parent stays alive for the duration of the agent's session. Tests override via `$ML_OWNER_PID` so they can use a known-alive background sleep as the owner without relying on whatever PPID the test harness happens to have.

**Hook cwd extraction:** the PreToolUse input JSON has `tool_input.cwd` (falls back to top-level `cwd`, then `$PWD`). Needed to compare against the locked worktree root for the owner check.

**Shared lock path across worktrees:** `git rev-parse --git-common-dir` + parent resolves to the main repo regardless of which worktree invokes the script. Without this, each worktree would have its own isolated `.claude/locks/`, defeating the point of the machine lock.

**Test hang debug:** `sleep 600 &` for the live-PID fixture, then racers as `( ... ) &`. Using bare `wait` at the outer scope hangs forever because it waits for the sleep too. Fix: collect specific PIDs (`racer_pids+=($!)`) and `wait "$pid"` each one.

**Quality gate numbers:** cargo test+clippy+clippy-gpu+doc via the new `machine-lock run` took ~28 min wall (doc was quick; the WFRX subspace tests dominate at 118s). Two pre-existing `clippy::expect_used` warnings in `src/symmetry/operations.rs` landed with ALOC-F5 (PR #100), unrelated to MLFX — flagged to EM for separate cleanup.

**Flagged for follow-up:**
- `src/symmetry/operations.rs:120,151` — two `i8::try_from(v).expect(...)` sites trip the new ERR2 `clippy::expect_used` lint. Code Reviewer / Researcher to convert to fallible with a tighter input-domain assertion.

## 2026-04-18 — MXB3: AdaptiveBeta::update Fe-trajectory unit test (PR #71)

Direct unit test in `src/scf/mixing/mod.rs::adaptive_beta_tests` feeding a synthetic 80-iter residual sequence to reproduce the Fe CCMX β-floor failure mode documented by `tests/mxba_adaptive_beta_fe.rs`.

**Sequence:** 10-iter plateau at 0.34 (monitor silent — ratios ≈ 1.0 hit the hysteresis band), then 70 iters oscillating `[0.34·1.3, 0.34·0.9]`. The alternating ratios (1.444 damp / 0.692 band) chain damps every 2 iters while making the 3-iter < 0.5 restore streak architecturally unreachable.

**β schedule (hand-verified against the code):** 0.3 → 0.21 → 0.147 → 0.103 → 0.072 → 0.0504 → 0.0353 → 0.0247 → 0.0173 → clamped at 0.015 by iter 27; stays at 0.015 through iter 80.

**Defense-in-depth assertion:** `max_good_streak < restore_window` — if a future MXB2 tuning change relaxes `restore_threshold` from 0.5 to 0.8, this test will fail loudly instead of silently no-oping, reminding whoever's on the branch to rework the multipliers.

**Clippy landmine:** `assert_eq!(float, float)` trips `clippy::float_cmp` even on exact constants. Switched to `(a - b).abs() < 1e-15` pattern — project convention anyway per CLAUDE.md `approx::relative_eq!`, but for sanity-pin of compile-time constants the `< 1e-15` form is fine.

**Not touched:** `tests/mxba_adaptive_beta_fe.rs` (per task: integration test stays as-is, this is its cheap companion). `AdaptiveBeta` logic/API (per task: test-only).

## 2026-04-18 — VNMT: FLUP brief was architecturally wrong (PR #62)

Added `test_single_channel_l2_m_isolation` in `src/potential/nonlocal.rs::tests` with hand-computed reference.

**Surprise — FLUP shape was unreachable.** Brief said "single (l, m) channel via D_ij zap", but production `NonlocalPotential` ties all m values of a given (radial-projector-pair, same-l) slot to the same `D_{ij}` — physics of the KB form. True single-m isolation via D_ij alone needs a non-physical synthetic projector. Landed as **three-live-m pin** (m=0, +1, +2 nonzero at chosen G-vectors with φ=0; m=-1, -2 zero) with each Y_{2,m} hand-computed in a 90-line doc block. Still catches √2-on-one-m bugs — sanity-check by temporarily multiplying Y_{2,+2} by √2 produced residual 3.5e-7 ≫ 1e-10 tolerance.

## 2026-04-18 — MXBA: adaptive β *regressed* Fe CCMX (PR #57)

Eyert 1996 §3.3 residual-norm monitor hooked into Anderson + Broyden `push_history` (after Kerker). Default β_min = max(0.05·β_start, 0.01), growth=1.2, damp=0.7, restore=0.5.

**Surprise — adaptive_beta defaults to false.** On Fe CCMX (existing `test_ccmx_fe_free_magnetization_converges`), adaptive β damps to β_min≈0.017 during the initial Δρ≈0.34 plateau (first ~5 iters) before DIIS has built history. Starved DIIS cannot escape → ConvergenceFailure at iter 80. Fixed-β 0.3 converges in 14 iters. Root cause: proposal assumed DIIS always "does most of the work"; in practice DIIS needs `max_history ≥ 2` BEFORE monitor can interpret a trajectory. Early-iter plateaus confuse Eyert's monitor. This is why VASP/ABINIT/QE don't adapt β inside `mix_rho` — they do it in user scripts with explicit restart logic.

Documented failure pinned by `tests/mxba_adaptive_beta_fe.rs` (`#[ignore]`). RCA correction: proposal initially hypothesized a "flat residual damps β" mechanism; real mechanism is *oscillation* around the plateau confusing the monitor.

**New helper: `KerkerSetup<'a>`** bundle struct — the 3 Kerker params kept mixer constructors inside `too_many_arguments` limit. Useful pattern for future mixer variants.

**Flagged for follow-up:**
- DIIS warm-up window (suppress monitor for first `max_history` iters) may be a quick fix — if it lands, adaptive could become default.
- Tune Eyert thresholds on a metallic case where adaptive β *helps* (blocked on C diamond @ 30 Ry plain converging first).

## 2026-04-18 — CAST: parse-time UPF validation is missing (PR #56)

Enabled 3 `cast_*` correctness lints; walked ~148 hits. 44 `#[allow(reason=...)]` annotations, 4 assertion-guarded rewrites, ~7 stylistic.

**Surprise — `src/pseudopotential/upf/convert.rs` accepts negative `angular_momentum` verbatim.** The new assert in `NonlocalPotential::new` is a safety net, but a parse-time validation error would be cleaner. Small proposal for Code Reviewer / Researcher. **Warning:** any future changes to UPF parsing need to consider this; negative-l malformed PPs will hit `lmax+1 as usize` wrapping and cause silent huge allocations without the debug_assert.

**Rebase warning:** mid-session VNLM/DEAD/TACC-I landed; after rebase, 9 new CAST hits appeared from VNLM's rewrite of `nonlocal.rs`. Expected pattern during concurrent agent sessions — leave buffer time for follow-up sweeps after rebases.

## 2026-04-18 — TACC-I (PR #54): ignore-reasons + fe_debug deletion

Rewrote 6 `#[ignore]` strings in `tests/qe_validation.rs` after `--ignored` harvest: 2× SYKP / 4× VGCMP. Deleted `tests/fe_debug.rs` (223 LOC, 6 tests — 5 weak/dead, superseded by VGCMP Phases 1-4). Migrated `test_fe_ewald_energy` → `qe_validation.rs::test_fe_bcc_ewald_vs_qe` (<0.01 eV tol vs -171.77906580 Ry).

**Handoff to EM:** Al (Z=13) at ≈73 meV is closest-to-passing of all ignored QE validation tests — barely over 50 meV tol on 8×8×8. If SYKP/MPSH lands it almost certainly passes.

## 2026-04-18 — MODR A/B/C/D: pure-move SCF split

Four pure-move refactors, each ≤1 PR, no behavior changes:

- **MODR-A** (PR #46): `scf/mixing.rs` (971 LOC) → `scf/mixing/{mod,anderson,broyden,kerker,linalg}.rs`. `AndersonMixer` + `PeriodicPulayMixer` co-located (periodic wrapper pokes Anderson private fields). Visibility tightened `pub` → `pub(crate)` on mixer types.
- **MODR-B** (PR #50): `scf/mod.rs` 1292 → 271 LOC. Split into `driver.rs`/`driver_spin.rs`/`report.rs`/`energy.rs`; `mod driver` etc. fully private.
- **MODR-C** (PR #48): `symmetry/density.rs` → `density/{mod,real_space,g_space}.rs`.
- **MODR-D** (PR #47): `pseudopotential/upf.rs` → `upf/{mod,xml,convert}.rs`.

**Landmine avoided (MODR-B):** first draft of `report.rs` collapsed both drivers into one `log_convergence_summary(..., ts, n_atoms)`. Original spin driver never emitted `Entropy (-TS):`; non-spin emits it when `|TS| > 1e-8`. Would have been a behavior change — split into `log_convergence_summary` (both) + `log_entropy` (non-spin only). Pinned by docstring.

**Warning for future refactors:** `#[allow(deprecated)]` on a facade `pub use` is required explicitly — the re-export alone without the attr fires deprecation through the module boundary.

**Flagged for follow-up:** `ctx.v_local_g0 * ctx.n_electrons` G0-shift expression now duplicated across 4 call sites in driver.rs + driver_spin.rs. A `with_g0_shift()` helper in `scf/energy.rs` is the natural DRY (not done per explicit MODR scope fence).

**Flagged for TW:** stale `src/pseudopotential/upf.rs` prose in `CLAUDE.md:64`, `docs/units.md:34`, `docs/nonlocal.md:88`, `LOGBOOK.md:92`, `src/scf/energy.rs:26`, `src/scf/mod.rs:52`.

## 2026-04-18 — PCFX landed (PR #44) — G-space symmetrization

Moved density symmetrization from real-space (rounding-sensitive on non-symmorphic grids) to G-space (exact via phase factors). Fixes the PCRS 1.204 eV plateau.

**Convention (verified against QE line-by-line — worth pinning):**
- `SpaceGroupOp::rotation` is the fractional-direct-space rotation `R`: atoms transform as `r' = R·r + τ`.
- Under the pullback `(S·ρ)(r) = ρ(S⁻¹ r)`, Miller indices rotate as `n → R^T · n` (NOT `R⁻¹`, NOT `R^{-T}`).
- Phase: `exp(-i·2π·n_dst·τ_S)` using DESTINATION Miller.
- QE stores `s(:,:,ns)` as the *transpose* of our `R` (`symm_base.f90:533`: atoms rotate as `rau = s^T · xau`).

**Non-obvious landmine — band-limitation requirement.** G-space formula ALSO requires `N·τ ∈ ℤ` UNLESS the input density is band-limited away from the Nyquist. If `R^T m` wraps (i.e. `R^T m ∉ [−N/2, N/2]`), the coefficient read corresponds to `k = R^T m − N·δ`; the inner application's phase differs by `exp(+i·2π·N·δ·τ)` which is 1 iff `N·τ ∈ ℤ`. Unit tests had to use band-limited inputs — naive `sin(i·0.37).abs()` generator broke idempotence by 5% on 18³.

**Key numbers (Si ecut=15, 4×4×4):** per-component self-check 1.204 → 3.5e-11 eV; E_total −231.8653 → −231.8429 eV (23 meV shift).

## 2026-04-18 — CCMX landed (PR #43)

Replaced `mixer_up`/`mixer_down` with `mixer_total`/`mixer_mag` in `run_scf_spin`. Forward basis change `(ρ↑, ρ↓) → (ρ_total, m)` before mix; inverse after. Matches QE `rhoz_or_updw` (`scf_mod.f90:1360-1414`). Kerker disabled on `mixer_mag`.

**Fe BCC 4×4×4 nspin=2 free-mag, Kerker:**
- Pre-CCMX: Δρ limit cycle at 0.254 for 200+ iters, |HF-KS| ≈ 13 eV.
- Post-CCMX: **14 iters**, Δρ=5.7e-4, |HF-KS|=1.06e-4 eV, M=0 μB.

**8×8×8 Δρ wobble:** energy stable to 6 decimals by iter 11, |HF-KS|=5e-5 eV, but Δρ then spikes between 1e-8 and 3e-5 as Anderson history becomes rank-deficient. The existing singular-pivot guard handles gracefully — benign, MXBA adaptive-β could damp sooner.

**Flagged:** `test_fe_bcc_fm_vs_qe` now converges cleanly but stays `#[ignore]` pending VGCMP (9.5 eV V_local gap). `test_fe_ferromagnetic_fixed_moment` still correctly ConvergenceFailure.

## 2026-04-18 — NCFX landed — Si 13.43 → 0.26 eV

Two compounding NLCC bugs fixed (diagnosed in VGC5):
1. `pseudopotential/upf.rs` PP_NLCC conversion: `/BOHR_TO_ANG` → `/BOHR3_TO_ANG3`.
2. `scf/potentials.rs::compute_core_density` integrand now has r² weight + 4π prefactor (QE `init_tab_rhc`).

**Si residual 0.26 eV** attributed to MP-shifted-vs-Γ-centered grid (SYKP): QE uses `4 4 4 0 0 0` (Γ-centered), pwdft-rs hard-codes shifted MP-1976. Γ eigenvalues differ by ~1 eV consistent with different k-meshes. `test_si_diamond_vs_qe` stays `#[ignore]` pointing at SYKP/MPSH.

**Flagged (residual work after NCFX):**
- Si VGC5 self-check `Σ − E_total = 1.18 eV` pre-existing (→ PCRS → PCFX).
- `test_fe_bcc_fm_vs_qe` (nspin=2 Kerker 8×8×8) no longer converges in 80 iters post-NCFX — pre-NCFX it "passed" via accidental XC cancellation. Needs CCMX to converge cleanly post-NCFX.

## 2026-04-17 — SPNC: Fe fixed-mag=2 fails under per-spin convergence

Replaced total-density delta with `max(density_diff(up), density_diff(down))` in `run_scf_spin`.

**Surprise — pre-SPNC "convergence" was illusory.** Fe fixed-mag=2 (total-only, SPXC in place) "converged" at iter 244 with Δρ_total=1e-6, |HF-KS|=13 eV. Post-SPNC: **ConvergenceFailure**, both channels pinned at `Δρ_up = Δρ_down = 0.254` from iter ~5 (steady-state limit cycle). Total density stable (dE~3e-7), so previous "convergence" was a +ε/−ε spin flip cancelling into total. Fixed-mag=2 is not a stable SCF fixed point for this LDA Fe PP (ground state is non-magnetic) — independent Anderson mixers can't coordinate inter-channel charge transfer. Fixing needs coupled-channel mixer → opened CCMX.

**Regression test reframed to Si nspin=2** (non-magnetic, relaxes to M=0, exercises spin XC): |HF-KS|=7.19e-7 eV at conv=1e-6, 23 iters, M=0.0. Assertion threshold 1e-5 eV (14× empirical headroom).

## 2026-04-17 — SPXC attempt 2 implemented

Recompute `lda_xc_spin_grid` from OUTPUT spin densities for E_KS; keep INPUT-derived for E_HF.

**Fe fixed-mag=2 4×4×4:** |E_HF−E_KS| 22.24 → 13.03 eV (~1.7×, not 10× proposal hoped). Why not 10×: nspin=2 convergence check used only rho_total, not per-spin — spin density (zeta) never driven to self-consistency, so E_xc[ρ_in,ζ_in] vs E_xc[ρ_out,ζ_out] differ at ~13 eV level regardless of SCF length. SPXC removes the artificial extra ~10 eV; residual zeta-inconsistency gap needs per-spin convergence criterion → opened SPNC.

**Isolation footgun (from attempt 1):** three parallel agents (VERF/SPXC/QEVL) escaped their worktrees and wrote to main checkout via `/Users/.../pwdft-rs/...` absolute paths. Lesson: **always verify Write/Edit paths resolve inside `.claude/worktrees/agent-*` before using them.** Isolation hook was subsequently added; do not rely on memory.

**Surprising side-note from VERF:** erf subtraction does NOT close Si 13.4 eV gap — numerically equivalent to bare-Coulomb + Simpson on our log mesh (max |Δ|=5.91e-9 eV across 20 shells). Landed anyway as cosmetic match to QE convention (insurance for future high-Z PPs where bare-Coulomb hits billions of eV near origin).

## 2026-04-17 — TAUD A/B/C/D/E: silent-pass hardening

**PR A:** 5 silent-pass `match`/`if let` patterns → `.expect()` + convergence guards. No bugs unmasked — pure hygiene debt.

**PR D — real finding.** Un-ignored `test_vloc_comparison_with_qe`; it FAILED post-VERF:
- V_local(G=0): ours +1.343 eV vs QE −1.003 eV (sign flip, diff 2.35 eV)
- |V_local(G=(1,0,0))|: ours 5.468 vs QE 6.968 (diff 1.50 eV)

Re-ignored with specific numbers + pointer to VGCMP Phase 1. Did NOT open new proposal (duplicate of VGCMP).

**Empirical margins for PR B/C/D/E tolerances** (replaced wish-at-write-time numbers):
- Si nspin=1 vs nspin=2 E diff ~0 → 1e-5 eV (was 0.5)
- GPU Si E_F: **6.969 eV** (not 5.97 as first guessed — all 4 bands occupied, E_F sits ~1 eV above HOMO) → `|E_F−6.969|<0.1` (was `∈[-5,10]`)
- GPU Si E: −198.8926 → `|E−(−198.8926)|<0.1` (was `∈[-300,-100]`)

**Clippy gotchas:** `const FOO: f64 = ...` inside a fn after let-bindings fires `items_after_statements`; use `let foo = 1.23_f64;`. `fn expect_...()` helper inside a test triggers same — inline the match.

**Handoff:** `ScfResult` doesn't derive `Debug`. Future tests printing Result values either derive Debug (`src/scf/mod.rs:131`) or use the 3-arm-match idiom.

## 2026-04-17 — PRPL landed (PR #39)

Periodic Pulay mixer (Banerjee et al., JCTC 12, 3053 (2016)) as BROY follow-up. Refactored `AndersonMixer::mix` into `push_history` + `diis_step`; `PeriodicPulayMixer` wraps Anderson, gates `diis_step` on `iteration.is_multiple_of(period) && history_len() >= 2`.

**Si Γ-only, ecut=100, conv=1e-6:** plain 10 iters, PeriodicPulay 8 iters (20% reduction on an insulator). Paper's strongest wins are metals/TMOs.

**Approve-with-nits follow-up:** `ScfParams::validate()` now rejects `MixingMode::PeriodicPulay { period: 0, .. }` with `InvalidInput`.

**Flagged:** no real metallic regression test for PRPL (Γ-only Si is an insulator); paper recommends `period = 5-8` for metals, our default is 3 (insulator sweet spot). Consider auto-selecting based on system class if a concrete metallic test case emerges.

## 2026-04-17 — FFTB (PR shipped)

`FFT3D` now owns two reusable `Array3<Complex64>` scratch buffers. Forward/inverse `copy_from_slice` into buf_a and ping-pong through 3 1D transforms. Measured -23 to -33% across 16³-48³ grids on `fft/scf_iter_20x_NxNxN` — matches proposal's allocation-rate argument (save ~512 KB `to_vec` + 512 KB zero-fill per call at 32³).

**Parallel-sharing cleared:** only two parallel-FFT sites — `density::compute_density` already per-worker FFT3D, k-point eigensolve doesn't touch `ctx.grid.fft`.

**Worktree weirdness observed:** during this session, in-progress edits to `src/fft.rs` and `benches/scf_benchmarks.rs` were *silently reverted* between bash invocations (compiled bench binary lacked newly-added `fft_scf_iteration` group despite source having it). Cost ~10 min. Worth investigating worktree/hook setup separately.

## 2026-04-16 — Early proposals: KBTF, CLEN, HRFK, DDUP, SIMP, BROY

- **KBTF (#1):** fixed 3 kb_projector failures. test_09 rewrite (UPF D_ij vs raw HGH h^l misconception); test_07 `#[ignore]` (trapezoidal, → SIMP); test_vloc `#[ignore]` (bare-Coulomb, → VERF).
- **CLEN (#2):** 5 small cleanups. GPU staging double-copy removed; `n_projectors`/`has_nlcc` become methods; deduped `apply_rotation` via `SymmOp::apply`.
- **HRFK (#6):** Harris-Foulkes diagnostic. Warns if `|HF-KS| > 0.01 eV` at convergence. Flagged subtle spin-path inconsistency (input exc_r mixed with output rho_xc_total) → later addressed by SPXC.
- **DDUP (e09513e):** `assemble_v_eff()` replaces inline V_eff in spin loop (gains `par_iter`); FFT norm delegated to `density_r_to_g()`; `rho_core_half` hoisted; `compute_occupations()` helper.
- **SIMP (#7):** Simpson quadrature at all 4 radial sites, matches QE `simpsn.f90` even-mesh correction. **Fe BCC QE validation NOW PASSES** (was 45.4 eV off); Si still 13.43 eV off (→ VGCMP). Un-ignored test_07. VERF unblocked.
- **BROY (#8):** Modified Broyden second method (Johnson PRB 38 12807) as `BroydenMixer` + `Mixer` enum dispatch. Adaptive-β and periodic-Pulay deferred to MXBA/PRPL.
