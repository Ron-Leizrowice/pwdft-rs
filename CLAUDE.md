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

## Tests & benchmarks

```bash
cargo test                                    # Tier 1 (default-fast)
cargo test --features gpu                     # with GPU tests
cargo test -- --ignored                       # Tier 2 heavy SCF suites
cargo test -- --include-ignored               # both tiers

cargo bench --bench scf_benchmarks
cargo bench --bench gpu_benchmarks --features gpu
```

Integration tests live in `pwdft/pwdft-core/tests/`: free-electron band validation (Si, C diamond, BCC Fe), KB projector validation, non-local symmetry, parallel consistency, GPU vs CPU consistency, VGC5 per-component energies + MADOC band-sum identity, QE validation (8-system reference set), spin polarization, WFRX subspace consistency, ITEV eigensolver cross-checks.

### Tier policy

The suite is split into two tiers via Rust's `#[ignore]` attribute with a `TSPL Tier-2: ...` reason string on each heavy case. Default `cargo test` runs Tier 1 only; opt into Tier 2 via `--ignored`.

- **Tier 1 — `cargo test`.** Unit tests + lightweight integration. Every test is either a pure unit check or a single-shot operation on a small matrix. Target ≤ 2 min wall on M3 Max; currently ~12 s warm-cache.
- **Tier 2 — `cargo test -- --ignored`.** Every case that runs a production-scale SCF loop (> 20 iters at `n_pw ≥ 100`, or two SCFs back-to-back). Target 5–10 min wall; currently ~58 s warm-cache.

**Tier-2 PR policy.** Any PR that touches `pwdft/pwdft-core/src/{scf,potential,symmetry,pseudopotential,eigensolver,gpu}/`, `basis.rs`, `fft.rs`, `ewald.rs`, `crystal.rs`, `kpoints.rs`, or bumps a numerics dep (faer / ndrustfft / nalgebra / ndarray) must run `cargo test -- --ignored` and report the outcome in the PR body. Doc-only, proposal-only, and lint-only PRs skip Tier 2.

`cargo test -- --ignored` also hits the pre-TSPL physics-blocker ignores (VGCH heavy-atom cells, Al ecut, C mixer, MXBA) which fail by design. The authoritative skip list is in the `#[ignore]` reason strings in `pwdft/pwdft-core/tests/qe_validation.rs` and `pwdft/pwdft-core/tests/mxba_adaptive_beta_fe.rs`.

## Observability, profiling, benchmarking

Three layers, three tools, no overlap:

- **Observability** — "what is the SCF doing right now?" → `log` crate (`info!` / `warn!` / `error!`) via `env_logger`, plus `indicatif` progress bars. Always-on, human-readable, low overhead.
- **Profiling** — "where does wall-time go?" → `samply` (canonical). Cross-platform, unprivileged, zero code overhead, emits a Firefox Profiler HTML artifact. Install once: `cargo install samply`. Always run under the machine lock — samply saturates CPU like `cargo bench`.
- **Benchmarks** — "did this change regress function X?" → `criterion` via `cargo bench`. Benches in `pwdft/pwdft-core/benches/`.

Do **not** use `cargo flamegraph`, `tracing-flame`, or hand-rolled `Instant::now()` timers for new profiling work. Samply sees inside `faer` / `ndrustfft` / BLAS where annotation-based tools cannot. Instruments.app is a fallback for Metal GPU timelines only (Xcode Metal debugger, GPU trace captures).

The `tracing` ecosystem is deliberately out of the stack. Reopen the question only if (a) async code enters the codebase, (b) structured per-span post-mortem analysis becomes necessary (regression watcher parsing SCF events as records), or (c) distributed tracing across a multi-process calculation becomes a requirement.

## Code quality gate

Run after finishing a batch of work — this is the merge criterion:

```bash
.claude/bin/machine-lock run "<role>" "quality gate" -- bash -c '
  cargo clippy -q --fix --allow-dirty --allow-staged --all-targets &&
  cargo clippy -q --all-targets &&
  cargo clippy -q --all-targets --features gpu &&
  RUSTDOCFLAGS="-D warnings" cargo doc --no-deps &&
  cargo test'
```

Both clippy invocations are required — without `--features gpu`, the `pwdft/pwdft-core/src/gpu/` tree and GPU-only test binaries are not linted. Rustdoc's `-D warnings` must travel through the `RUSTDOCFLAGS` env var; current cargo rejects it after `--`. Do not suppress clippy warnings like `too_many_arguments` — refactor instead. Do not `#[allow]` rustdoc warnings — fix the prose.

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

- **Always use a worktree.** Never modify the main checkout. Spawn sub-agents with `isolation: "worktree"`; hook-enforced via `.claude/bin/check-worktree.sh` PreToolUse.
- **Branch from `origin/main`** (local `main` can lag). Branch name: `<PROPOSAL-ID>/<slug>`.
- **Rebase on `origin/main` before opening the PR.** Prevents DIRTY/CONFLICTING states.
- **Never commit directly to main** — all work goes through PR.
- **One proposal per branch.** PR title: `<PROPOSAL-ID>: <description>`; commit messages the same. PR body must include Summary + Test Plan.
- **Quality gate before PR** — see § Code quality gate above; also in `.claude/agents/shared/quality-gate.md`.
- **Proposals drive work.** See `proposals/INDEX.md`. Each proposal has a 4-letter ID, frontmatter (priority/complexity/risk/dependencies), and an implementation plan. Propose first, implement after EM approval.
- **Logbooks** at `.claude/logbooks/<role>.md`. Read at session start, append at session end. Sub-agents in worktrees paste handoff text into the PR body (the hook blocks writes to the main checkout's logbook); the EM merges logbook entries on merge.

**Agent roles** (see `.claude/agents/*.md` for full prompts):

- **Engineering Manager** — reviews PRs, manages proposals, merges to main. Never writes code.
- **Core Engineer** — implements proposals. Scientist-developer mindset.
- **Performance Engineer** — benchmarks, profiles, optimizes. Measure first.
- **Researcher** — owns physics / math correctness. Proposes features, validates against QE.
- **Code Reviewer** — owns code quality. Hunts dead code, enforces idioms, improves tests.
- **Technical Writer** — owns documentation. README, CLAUDE.md, docstrings, comments.

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
