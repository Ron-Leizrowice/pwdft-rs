# CLAUDE.md

Guidance for Claude Code (claude.ai/code) when working in this repository.

## What this is

A plane-wave density functional theory (PWDFT) solver in Rust, targeting macOS with Apple Metal GPU acceleration. Used for real research — correctness is paramount. Results are validated against Quantum ESPRESSO 7.5.

The repo is a Cargo workspace:

- `pwdft/pwdft-core/` — the Rust solver (library + binary, integration tests, benches)
- `pwdft/pwdft-validation/` — Python (uv) validation harness: QE invocation, reference-data parsing, diagnostic scripts
- `pwdft/faer/` — vendored `faer` v0.24.0 with one local patch (see § Vendored dependencies)
- `data/qe/`, `data/csv/` — reference data (QE outputs, CSV pins)
- `inputs/` — example YAML input decks
- `qe-7.5/` — Quantum ESPRESSO 7.5 source + build, symlinked during setup
- `pseudopotentials/` — NC / USPP / PAW libraries consumed by QE and pwdft-rs

## Build & run

```bash
cargo build --release                          # CPU only
cargo build --release --features gpu           # with Metal/Vulkan via wgpu

cargo run --release -- --input inputs/si_scf.yaml
cargo run --release -- --input inputs/si_free_electron.yaml -o bands.tsv
```

## Workflow skills

All common cargo / profiling / review operations are wrapped in skills under `.claude/skills/`. Prefer skills over raw bash — they handle the machine lock, arg parsing, and feature-flag invariants for you:

- `/quality-gate` — the full merge-criterion gate (fix clippy + default clippy + gpu clippy + rustdoc `-D warnings` + Tier-1 tests).
- `/test [args]` — tier-aware `cargo test`. Default Tier-1; `--tier2` runs `--ignored`; `--all` runs `--include-ignored`; anything else passes through.
- `/bench <bench-name>` — `cargo bench` under the lock.
- `/profile <cmd>` — `samply record <cmd>` under the lock.
- `/lint` — both clippy invocations (auto-fix + default + `--features gpu`).
- `/cargo <args>` — escape hatch: any cargo command, lock-wrapped.
- `/worktree-start PROP/slug` — fetch + branch from `origin/main` boilerplate.
- `/pr-submit` — quality gate + rebase + `gh pr create`.
- `/pr-review <N>` — EM review helper: metadata + diff + checklist.
- `/merge <N>` — EM merge trilogy.
- `/proposal {create,list,complete}` — proposal lifecycle.
- `/qe-runner` — Quantum ESPRESSO invocation recipes.

## Tests & tiers

Integration tests live in `pwdft/pwdft-core/tests/`: free-electron bands (Si, C diamond, BCC Fe), KB projector, non-local symmetry, parallel consistency, GPU vs CPU consistency, VGC5 per-component energies + MADOC band-sum identity, QE validation (8-system reference set), spin polarization, WFRX subspace consistency, ITEV eigensolver cross-checks.

The suite is split into two tiers via Rust's `#[ignore]` attribute with a `TSPL Tier-2: ...` reason string on each heavy case:

- **Tier 1** (default via `/test`) — unit tests + lightweight integration; every test is a pure unit check or a single-shot small-matrix operation. Target ≤ 2 min wall; currently ~12 s warm-cache on M3 Max.
- **Tier 2** (`/test --tier2`) — every case that runs a production-scale SCF loop (> 20 iters at `n_pw ≥ 100`, or two SCFs back-to-back). Target 5–10 min wall; currently ~58 s warm-cache.

**Tier-2 PR policy.** Any PR that touches `pwdft/pwdft-core/src/{scf,potential,symmetry,pseudopotential,eigensolver,gpu}/`, `basis.rs`, `fft.rs`, `ewald.rs`, `crystal.rs`, `kpoints.rs`, or bumps a numerics dep (faer / ndrustfft / nalgebra / ndarray) must run `/test --tier2` and report the outcome in the PR body. Doc-only, proposal-only, and lint-only PRs skip Tier 2. `/test --tier2` hits the pre-TSPL physics-blocker ignores (VGCH heavy-atom cells, Al ecut, C mixer, MXBA) that fail by design — the authoritative skip list is in the `#[ignore]` reason strings in `pwdft/pwdft-core/tests/qe_validation.rs` and `pwdft/pwdft-core/tests/mxba_adaptive_beta_fe.rs`.

**Profile split: local vs. CI.** Local `cargo test` uses the O3 `[profile.test]` from `Cargo.toml` (per TPRF) — warm-cache M3 Max is runtime-bound, and `opt-level=3` keeps the few hot-path unit tests fast (12 s warm Tier-1). CI is the other axis: cold-cache Ubuntu runners are compile-bound, so `.github/workflows/rust.yml` overrides Tier-1 to `cargo test --profile=dev -p pwdft-core` (per TDBG). Baseline on run `24663826105` was 343 s for the test step; the debug-profile override targets the 4-minute (240 s) budget enforced inline in the workflow. Tier-2 stays on the O3 profile everywhere — those SCF-loop tests would be glacial at `opt-level=0`.

## Observability stack

Three layers, three tools, no overlap:

- **Observability** ("what is the SCF doing right now?") → `log` + `env_logger` + `indicatif`. Always-on, human-readable.
- **Profiling** ("where does wall-time go?") → `samply` via `/profile`. Zero code overhead; emits a Firefox Profiler HTML artifact.
- **Benchmarks** ("did this change regress function X?") → `criterion` via `/bench`. Benches in `pwdft/pwdft-core/benches/`.

Do **not** reach for `cargo flamegraph`, `tracing-flame`, or hand-rolled `Instant::now()` timers — samply sees inside `faer` / `ndrustfft` / BLAS where annotation-based tools can't. Instruments.app is the fallback for Metal GPU timelines only. The `tracing` ecosystem stays out of the stack until (a) async code enters the codebase, (b) structured per-span post-mortem parsing becomes necessary, or (c) distributed tracing across a multi-process calculation is required.

## Architecture

**Entry point:** `pwdft/pwdft-core/src/main.rs` parses CLI args and YAML input (`settings.rs`), then either computes a free-electron band structure or runs SCF.

**SCF loop** (`scf::run_scf` in `pwdft/pwdft-core/src/scf/mod.rs` — thin dispatcher into `scf::driver::run_scf_unpolarized` for `nspin=1` or `scf::driver_spin::run_scf_spin` for `nspin=2`):

1. Build local pseudopotential V_local on the FFT grid (spherical Bessel transform). If any PP has NLCC (`core_correction="T"`), also build ρ_core(r) on the grid.
2. Initialize density via SAD (superposition of atomic densities).
3. Each iteration: Hartree (valence only) → LDA XC (on ρ_val + ρ_core if NLCC — Louie, Froyen, Cohen, PRB 26, 1738 (1982)) → assemble V_eff → Hamiltonian (kinetic + V_eff + KB non-local) → diagonalize (faer; `EigensolverKind::Dense` default) → Fermi-Dirac occupations → density → symmetrize (G-space phase factors, PCFX) → convergence check → density mixing. Mixers: Anderson/Pulay (DIIS), modified Broyden, Periodic Pulay; any can be Kerker-preconditioned. Spin driver uses coupled-channel (ρ_total, m) basis (CCMX).
4. Total energy (kinetic + local + non-local + Hartree + XC + Ewald), with NLCC double-counting subtraction applied to E_xc. Drivers return an `EnergyComponents` breakdown (VGC5 diagnostic).

All source paths below are relative to `pwdft/pwdft-core/src/`.

**Module groups:**

- **Crystal & basis:** `crystal.rs`, `basis.rs`, `kpoints.rs`, `atoms.rs`.
- **Pseudopotentials:** `pseudopotential/mod.rs` (`PseudopotentialData` + `v_local_of_g`); `pseudopotential/upf/` holds `mod.rs` (40-line `parse(&str)`), `xml.rs`, `convert.rs` (Ry→eV, Bohr→Å at the UPF boundary; assembles `PP_RHOATOM`, `PP_NLCC`).
- **Potentials:** `potential/xc.rs` (PZ LDA, spin-polarized variant, PBE GGA), `potential/local.rs`, `potential/nonlocal.rs` (Kleinman-Bylander, arbitrary l; V_NL via single GEMM per VNLM). Hartree is assembled inline in `scf::{energy,driver,driver_spin}` — no separate `hartree.rs`.
- **SCF internals:** `scf/{mod,driver,driver_spin,report,energy,density,initial_density,smearing,context,potentials,grid}.rs`.
- **Density mixing:** `scf/mixing/mod.rs` (`MixingMode` + `Mixer` dispatcher). Algorithms: `anderson.rs`, `broyden.rs` (Johnson PRB 38, 12807), `kerker.rs`, `linalg.rs`.
- **Numerics:** `fft.rs` (ndrustfft, zero unsafe), `eigensolver/{mod,dense,iterative}.rs`, `numerics.rs`, `ewald.rs`.
- **Symmetry:** `symmetry/{mod,operations,detect,kpoints}.rs`; `symmetry/density/` has `mod.rs` facade, `real_space.rs` (legacy `#[deprecated]`), `g_space.rs` (PCFX, what SCF uses).
- **GPU:** `gpu/mod.rs` (wgpu compute), `gpu/shaders/` (WGSL kernels for Hartree, LDA XC, V_eff assembly).

**GPU strategy:** Optional `gpu` feature flag. `GpuAccelerator` with `BufferPool` of pre-allocated f32 buffers. GPU kernels run in f32, CPU in f64, with conversion at boundaries. Falls back to CPU (rayon) when GPU unavailable.

**Eigensolver:** `EigensolverKind::Dense` (default) uses `faer::SelfAdjointEigen` — full O(n³). `EigensolverKind::Iterative` (opt-in via `scf.eigensolver: iterative`) uses faer's partial Arnoldi/Krylov-Schur. **Stay on `Dense` for now** — Iterative has open correctness issues tracked in ITEV: size-independent `n_request` padding drops 3-fold-degenerate valence clusters at `n_pw ≳ 725`, and the iterative dispatch bypasses WFRX warm-start. See `proposals/ITEV-*.md` § Status.

**Key types:** `Crystal`, `BasisSet`, `KPoint`, `PseudopotentialData`, `ScfParams`/`ScfResult`, `EnergyComponents`, `EigensolverKind`, `MixingMode`/`Mixer`, `NonlocalPotential`, `EigenResult`, `SymmetryInfo`/`SpaceGroupOp`.

## Vendored dependencies

`faer` v0.24.0 is vendored at `./pwdft/faer/` with a single local edit: `MAX_REORTH = 3` on `iterate_lanczos` reorthogonalization to prevent the upstream infinite-loop on near-null Krylov vectors. `Cargo.toml` wires it via `[patch.crates-io]` so the vendored copy takes over for both `faer` and `faer-traits`.

- **Editing faer is allowed and expected.** A proposal that needs to change the vendored solver edits files under `pwdft/faer/faer/src/...` directly, runs the pwdft-rs quality gate, and commits the faer-side + pwdft-rs-side changes together on its feature branch. No separate upstream PR is required for the local build.
- **Upstream-submit procedure** lives in `docs/FAER_ITERATE_LANCZOS_FIX.md`. Follow it when a vendored edit is ready to submit back to [`codeberg.org/sarah-quinones/faer`](https://codeberg.org/sarah-quinones/faer.git).
- **Build artifacts** (`pwdft/faer/target/`, `pwdft/faer/*.out`, stray `.git`) are in `.gitignore`. Running `cargo test` inside `pwdft/faer/` is safe.
- **Version bumps:** re-vendor by copying new source over `pwdft/faer/`, re-applying the patch (or dropping it if upstream absorbed the fix), and running the pwdft-rs gate. `pwdft/faer/Cargo.lock` is separate from pwdft-rs's `Cargo.lock`.

## Workflow

Branch-and-PR workflow with isolated worktrees and multiple concurrent agents. **Agents should read the shared protocols in `.claude/agents/shared/` at session start** — those files are the authoritative source for worktree, machine-lock, quality-gate, FLUP, and session-end rules. This section is a pointer, not a duplicate.

- **Branch from `origin/main`** (local `main` can lag). Branch name: `<PROPOSAL-ID>/<slug>`.
- **Rebase on `origin/main` before opening the PR.**
- **Never commit directly to main** — all work goes through PR.
- **One proposal per branch.** PR title: `<PROPOSAL-ID>: <description>`; commit messages the same. PR body must include Summary + Test Plan.
- **Quality gate before PR** — see § Code quality gate above; also in `.claude/agents/shared/quality-gate.md`.
- **Proposals drive work.** See `proposals/INDEX.md`. Each proposal has a 4/5-letter ID, frontmatter (priority/complexity/risk/dependencies), and an implementation plan. Propose first, implement after EM approval.
- **Logbooks** at `.claude/logbooks/<role>/` — one file per session, named `YYYY-MM-DD-<slug>.md`. Sub-agents write their entry inside their worktree and the file lands as part of their PR. See `.claude/logbooks/README.md` and `.claude/agents/shared/session-end.md`. Search (`rg <topic> .claude/logbooks/`) at session start before re-investigating a known area.

**Agent roles** (see `.claude/agents/*.md` for full prompts; ):

- **`engineering-manager`** — reviews PRs, manages proposals, merges to main. Never writes code. Runs on the main checkout. Specialists never spawn this agent; delegation flows one-way (see `.claude/agents/shared/flup.md`).
- **`core-engineer`** — implements proposals. Scientist-developer mindset.
- **`performance-engineer`** — benchmarks, profiles, optimizes. Measure first.
- **`researcher`** — owns physics / math correctness. Proposes features, validates against QE.
- **`code-reviewer`** — owns code quality. Hunts dead code, enforces idioms, improves tests.
- **`technical-writer`** — owns documentation. README, CLAUDE.md, docstrings, comments.

## Machine coordination

The machine lock serializes CPU-bound work to protect `cargo bench` wall-time measurements from contamination. See `.claude/agents/shared/machine-lock.md` for the protocol, and `.claude/bin/machine-lock` for the implementation. TL;DR:

- Acquire before any `cargo` subcommand or QE invocation (`pw.x`, `ph.x`, `pp.x`, `bands.x`, `projwfc.x`, `q2r.x`, `matdyn.x`, `dos.x`).
- Prefer the `run` one-liner: `.claude/bin/machine-lock run "<role>" "<desc>" -- <cmd>`.
- Lock is worktree-scoped (MLFX). Running cargo from a different worktree while the lock is held is denied by the Bash PreToolUse hook. Staleness is PID-based with a secondary 3-hour time cap.
- Lock state lives in `.claude/locks/machine.lock.d/` (gitignored). Never force-remove another agent's lock.
- Shell-test suite: `.claude/bin/tests/machine-lock.test.sh` — run when touching the lock scripts.

## Conventions

- `Cargo.toml` uses `>=` version specifiers (not `^` or exact).
- Rust edition 2024. Release profile: `opt-level = 3`, thin LTO.
- Floating-point comparisons in tests use `approx::relative_eq!`.
- **Internal units are eV / Å / e·Å⁻³.** Used by SCF arrays, `EnergyComponents`, Hamiltonian assembly, and every routine downstream of pseudopotential parsing. Conversion constants in `pwdft/pwdft-core/src/consts.rs` (`HA_TO_EV`, `RY_TO_EV`, `BOHR_TO_ANG`, `BOHR3_TO_ANG3`, `E2_COULOMB` in eV·Å, `HBAR2_OVER_2M` in eV·Å²). Ry/Bohr appear **only** at the UPF boundary in `pwdft/pwdft-core/src/pseudopotential/upf/convert.rs`.
- Input files are YAML in `inputs/`, parsed via `serde_yaml_ng` into `Settings`.
- Pure Rust stack: faer (eigensolver), ndrustfft (FFT), nalgebra (geometry), ndarray (grid ops).
- No system dependencies required for the default build. GPU requires the `gpu` feature flag.
- Validation harness in `pwdft/pwdft-validation/pwdft_validation/` (Python, uv environment): QE invocation, reference-data parsing, diagnostic scripts. Reference data in `data/qe/` and `data/csv/`.
- QE 7.5 source in `qe-7.5/` for reference during validation.
- **Rust `.round()` vs Python `round()`.** Rust rounds half-away-from-zero; Python 3 rounds half-to-even (banker's rounding). When porting a grid-index computation between Rust and a Python validation script, verify the rounding convention — mismatches appear as off-by-one errors on the ±0.5 boundary.
- **QE k-point weights pre-multiply `degspin`.** In QE output and in `charge-density.dat`, `wk` already includes the spin degeneracy factor (`degspin = 2` for `nspin=1`, `1` for `nspin=2`). Do not multiply by `degspin` again when cross-checking from Python — summing `wk * occ` directly yields the correct electron count.
