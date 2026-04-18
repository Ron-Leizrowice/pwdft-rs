# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

A plane-wave density functional theory (PWDFT) solver in Rust, targeting macOS with Apple Metal GPU acceleration. Used for real research — correctness is paramount. Results are validated against Quantum ESPRESSO 7.5.

## Build & Run

```bash
# Build (CPU only)
cargo build --release

# Build with GPU (Metal/Vulkan via wgpu)
cargo build --release --features gpu

# Run SCF calculation
cargo run --release -- --input examples/si_scf.yaml

# Run band structure
cargo run --release -- --input examples/si_free_electron.yaml -o bands.tsv
```

## Tests & Benchmarks

```bash
cargo test                                    # all tests (~265, ~24s)
cargo test --features gpu                     # with GPU tests (~268, ~28s)
cargo test test_name                          # single test by name
cargo test --test free_electron_bands         # single integration test file
cargo test -- --nocapture                     # with stdout

cargo bench --bench scf_benchmarks            # SCF benchmarks
cargo bench --bench gpu_benchmarks --features gpu  # GPU benchmarks
```

Integration tests in `tests/`: free-electron band validation (Si, C diamond, BCC Fe), KB projector validation, non-local symmetry, parallel consistency, GPU vs CPU consistency.

## Code Quality

After finishing a batch of work, always run:
```bash
cargo clippy -q --fix --allow-dirty --allow-staged --all-targets
cargo clippy -q --all-targets                  # default-feature warnings
cargo clippy -q --all-targets --features gpu   # GPU feature warnings
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps # rustdoc is error-clean
cargo test                                     # verify nothing broke
```
Both clippy invocations are required: without `--features gpu`, the `gpu/` source tree and the GPU-only test binaries are not linted, so warnings accumulate silently (QLN2 caught an `uninlined_format_args` violation in `tests/gpu_consistency.rs` that had slipped past the default-feature gate for weeks). Fix auto-fixable warnings, address remaining ones. Do not suppress codesmell warnings like `too_many_arguments` — refactor the code instead.

`RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` is now part of the gate (RDOC 2026-04-18). Any new docstring that breaks an intra-doc link, leaves a bracket unescaped, or links at a private item will fail the gate — either fix the prose or drop the link. Do not `#[allow]` rustdoc warnings. (Note: `cargo doc --no-deps -- -D warnings` is rejected by current cargo — the `-D warnings` flag must travel through the `RUSTDOCFLAGS` env var.)

## Architecture

**Entry point:** `src/main.rs` parses CLI args and YAML input (`src/settings.rs`), then either computes a free-electron band structure or runs SCF.

**SCF loop** (`scf::run_scf` in `src/scf/mod.rs` — a thin dispatcher that validates inputs and hands off to `scf::driver::run_scf_unpolarized` in `src/scf/driver.rs` for `nspin=1` or `scf::driver_spin::run_scf_spin` in `src/scf/driver_spin.rs` for `nspin=2`). The central computation pipeline:
1. Build local pseudopotential V_local on FFT grid (spherical Bessel transform). If any PP has NLCC (`core_correction="T"`), also build ρ_core(r) on the grid (same Bessel transform; see `scf::potentials::compute_core_density`).
2. Initialize density via SAD (superposition of atomic densities).
3. Each iteration: Hartree potential (valence density only) → LDA XC (on ρ_val + ρ_core if NLCC — Louie, Froyen, Cohen, PRB 26, 1738 (1982)) → assemble V_eff → build Hamiltonian (kinetic + V_eff + KB non-local) → diagonalize (faer, `EigensolverKind::Dense` by default or `Iterative` when opted in) → Fermi-Dirac occupations → reconstruct density → symmetrize (G-space phase factors, PCFX) → check convergence → density mixing. Available mixers: Anderson/Pulay (DIIS), modified Broyden (BROY), and Periodic Pulay (PRPL); any mixer can be combined with Kerker preconditioning. The spin driver uses the coupled-channel (ρ_total, m) basis (CCMX) rather than independent (ρ↑, ρ↓), so both channels share residual history.
4. Compute total energy (kinetic + local + non-local + Hartree + XC + Ewald), with the NLCC double-counting subtraction applied to E_xc when core correction is active. Each driver also returns an `EnergyComponents` breakdown (VGC5 diagnostic: per-term energies and the Harris-Foulkes stationary estimator).

**Module groups:**

- **Crystal & basis:** `crystal.rs` (lattice + atoms), `basis.rs` (G-vectors up to ecut), `kpoints.rs` (Monkhorst-Pack, band paths), `atoms.rs` (elements 1-92).
- **Pseudopotentials:** `pseudopotential/mod.rs` (`PseudopotentialData` + `v_local_of_g`); `pseudopotential/upf/` is a folder with `mod.rs` (40-line `parse(&str)` entry point), `xml.rs` (text helpers for UPF v2 XML), and `convert.rs` (Ry→eV, Bohr→Å unit conversion, `PP_RHOATOM`, `PP_NLCC` assembly).
- **Potentials:** `potential/xc.rs` (Perdew-Zunger LDA, spin-polarized variant), `potential/local.rs`, `potential/nonlocal.rs` (Kleinman-Bylander separable form, arbitrary l via recurrence; V_NL assembly via single GEMM, VNLM). The Hartree potential is assembled inline in `scf::energy` / `scf::driver` / `scf::driver_spin`; there is no standalone `potential/hartree.rs` module.
- **SCF internals:** `scf/mod.rs` (thin `run_scf` dispatcher + `ScfParams`/`ScfResult`), `scf/driver.rs` (non-spin hot loop), `scf/driver_spin.rs` (spin-polarized hot loop), `scf/report.rs` (per-iteration progress + final summary; `IterationReport`, `log_iteration`, `log_convergence_summary`), `scf/energy.rs` (total-energy assembly + `EnergyComponents`), `scf/density.rs`, `scf/initial_density.rs` (SAD), `scf/smearing.rs` (Fermi-Dirac, Gaussian, Methfessel-Paxton, cold), `scf/context.rs` (immutable per-calculation state), `scf/potentials.rs` (Hamiltonian assembly helpers), `scf/grid.rs` (FFT grid setup).
- **Density mixing:** `scf/mixing/mod.rs` exposes `MixingMode` (Plain / Kerker / Broyden / PeriodicPulay) and the `Mixer` dispatcher enum. Algorithms live in siblings: `anderson.rs` (Anderson/Pulay DIIS + PeriodicPulay wrapper), `broyden.rs` (modified Broyden, Johnson PRB 38, 12807), `kerker.rs` (preconditioner + Thomas-Fermi `q_TF` auto-estimate), `linalg.rs` (small Gauss-elimination solver for the DIIS/Broyden linear system).
- **Numerics:** `fft.rs` (3D FFT via ndrustfft, zero unsafe), `eigensolver/mod.rs` (`EigensolverKind::{Dense, Iterative}`), `eigensolver/dense.rs` (full faer Hermitian eigendecomposition), `eigensolver/iterative.rs` (ITEV — faer's partial Arnoldi/Krylov-Schur, shift-and-flip for lowest `n_bands`; experimental, opt-in, blocked on an upstream faer 0.24 `iterate_lanczos` reorthogonalization bug), `numerics.rs` (Simpson and radial quadrature), `ewald.rs` (ion-ion energy, erfc via puruspe).
- **Symmetry:** `symmetry/mod.rs` (`SymmetryInfo`), `symmetry/operations.rs` (`SpaceGroupOp` {R|τ}), `symmetry/detect.rs` (space group finder), `symmetry/kpoints.rs` (k-point reduction to IBZ), `symmetry/density/` (folder: `mod.rs` with the facade, `real_space.rs` — legacy `#[deprecated]` `nint`-rounding path, `g_space.rs` — PCFX phase-factor form; SCF always uses the G-space path).
- **GPU:** `gpu/mod.rs` (wgpu compute), `gpu/shaders/` (WGSL kernels for Hartree, LDA XC, V_eff assembly).

**GPU strategy:** Optional `gpu` feature flag. `GpuAccelerator` with `BufferPool` of pre-allocated f32 buffers. GPU kernels run in f32, CPU in f64, with conversion at boundaries. Falls back to CPU (rayon) when GPU unavailable. Three WGSL shaders handle the per-iteration grid operations.

**Eigensolver:** The SCF loop diagonalizes the Kohn-Sham Hamiltonian at each k-point once per iteration. `EigensolverKind::Dense` (the default) uses `faer::SelfAdjointEigen` — full O(n³) LAPACK-equivalent decomposition. `EigensolverKind::Iterative` (opt-in via `scf.eigensolver: iterative` in YAML) uses faer's partial Arnoldi/Krylov-Schur solver to compute only the lowest `n_bands` eigenpairs, with a shift-and-flip trick to map "algebraically lowest of H" to "largest-magnitude of σI − H". Iterative is typically 3-10× faster at `n_pw ≥ 200`, but it is currently gated behind a FIXME: faer 0.24's `iterate_lanczos` can spin on near-null Krylov vectors, so SCF loops on ill-conditioned inputs may hang. Until the upstream bug is fixed, stay on `Dense`.

**Key types:** `Crystal`, `BasisSet`, `KPoint`, `PseudopotentialData`, `ScfParams`/`ScfResult`, `EnergyComponents`, `EigensolverKind`, `MixingMode`/`Mixer`, `NonlocalPotential`, `EigenResult`, `SymmetryInfo`/`SpaceGroupOp`.

## Workflow

This project uses a branch-and-PR workflow. Multiple agents may work concurrently.

- **Always use a worktree.** Never modify files in the user's main checkout. Use `isolation: "worktree"` when spawning agents, or `EnterWorktree` for interactive work. This keeps the main checkout clean and allows concurrent agents.
- **Worktree isolation is hook-enforced.** `.claude/bin/check-worktree.sh` (PreToolUse on Edit/Write/MultiEdit) blocks writes to (a) the main checkout's `src/`, `tests/`, `benches/`, `build.rs` from any cwd, and (b) anything outside your worktree when your cwd is in a worktree. If a write is blocked, your `file_path` is wrong — fix the path, don't disable the hook. See `.claude/agents/*.md` for the full protocol.
- **Branch from `origin/main`, not local `main`.** Local `main` may be stale. Run `git fetch origin && git checkout -b <PROPOSAL-ID>/<slug> origin/main` from inside your worktree.
- **Rebase your branch onto `origin/main` before submitting the PR.** Avoids "DIRTY/CONFLICTING" PRs that the EM has to resolve manually.
- **Never commit directly to main.** All work happens on feature branches.
- **One proposal per branch.** Branch name: `<PROPOSAL-ID>/<slug>` (e.g., `SIMP/simpson-quadrature`).
- **PRs against main.** Title format: `<PROPOSAL-ID>: <description>`. PR body must reference the proposal and include a summary and test plan.
- **Quality gate before PR:** `cargo test` plus both `cargo clippy -q --all-targets` and `cargo clippy -q --all-targets --features gpu` must pass (see § Code Quality for why both clippy invocations are required).
- **Commit messages:** `<PROPOSAL-ID>: <imperative description>`.
- **Proposals drive work.** See `proposals/INDEX.md` for the backlog. Each proposal has a 4-letter ID (e.g., `SIMP`), frontmatter with priority/complexity/risk/dependencies, and an implementation plan.
- **Propose first, implement after approval.** All roles draft proposals and wait for EM approval before starting work.
- **Logbooks** at `.claude/logbooks/<role>.md` — every agent reads theirs at session start, appends findings at session end. Cross-read other roles' logbooks when relevant.
- **Agent roles** (see `.claude/agents/` for full definitions):
  - **Engineering Manager** — reviews PRs, manages proposals, merges to main. Never writes code.
  - **Core Engineer** — implements proposals (refactors, features, bug fixes). Scientist-developer mindset.
  - **Performance Engineer** — benchmarks and optimizes. Measure first, optimize second.
  - **Researcher** — owns physics/math correctness. Proposes features, validates against QE and literature.
  - **Code Reviewer** — owns code quality. Hunts dead code, enforces idioms, improves logging/tests.
  - **Technical Writer** — owns documentation. README, CLAUDE.md, docstrings, code comments.

## Machine Coordination

Multiple agents share this machine. The machine lock exists for **benchmark integrity**: it serializes CPU-bound work so that `cargo bench` wall-time measurements aren't polluted by background CPU load from another agent's test run, compile, or reference calculation. This is a benchmark-contamination concern, not a test-isolation concern — two heavy CPU jobs are usually correct when run concurrently; they just invalidate each other's timings.

**General principle:** any CPU-bound job that could contaminate a benchmark measurement must hold the lock while it runs. That includes all `cargo` subcommands *and* all Quantum ESPRESSO invocations.

**What requires the lock:**
- `cargo test`, `cargo bench`, `cargo build`, `cargo clippy` — any cargo command that compiles or runs code.
- **Any Quantum ESPRESSO run:** `pw.x`, `mpirun pw.x`, `ph.x`, `pp.x`, `bands.x`, `projwfc.x`, `q2r.x`, `matdyn.x`, `dos.x` — whether it's a one-off reference calculation, a validation against `qe_validation/`, or regeneration of reference data. QE is multi-threaded/multi-process and saturates the CPU; running it during a benchmark window corrupts the numbers.
- Any other long-running CPU-bound job (Python validation scripts under `scripts/` that spin up BLAS, etc.).

**What does NOT require the lock:** reading files, editing code in worktrees, writing proposals, git operations, `machine-lock status`.

```bash
# Check lock state:
.claude/bin/machine-lock status

# Acquire before cargo commands (from inside your worktree — scope is recorded):
.claude/bin/machine-lock acquire "Core Engineer" "cargo test"
cargo test
.claude/bin/machine-lock release

# Block until the lock is free (default 10 min, configurable):
.claude/bin/machine-lock acquire --wait --timeout=900 "Core Engineer" "long bench"

# One-liner (acquire --wait + run + release). `run` uses --wait internally,
# so two concurrent run calls serialize instead of racing:
.claude/bin/machine-lock run "Core Engineer" "cargo test" -- cargo test

# QE runs use the same lock — wrap the whole mpirun invocation.
# Auto-detect cores; do not hardcode rank counts:
NP=$(sysctl -n hw.ncpu)  # macOS; use $(nproc) on Linux
.claude/bin/machine-lock run "Researcher" "QE Si SCF validation" -- \
  gtimeout 600 mpirun -np "$NP" qe-7.5/build/bin/pw.x -in si.in
```

- **The lock is worktree-scoped (MLFX 2026-04-18).** `acquire` records your
  worktree root; the Bash PreToolUse hook compares every cargo command's cwd
  against that root. An agent running `cargo test` from a different worktree
  (or from the main checkout) while the lock is held is denied, with a
  message pointing at the owning worktree. Acquire from your worktree once
  at session start, run all cargo/QE commands from that same worktree.
- **Acquire is atomic.** `machine-lock` uses `mkdir` (POSIX-atomic) to claim
  the lock directory — two racing acquires cannot both win.
- **Staleness is PID-based first.** The acquirer's shell PID is recorded;
  a lock is only stale if that PID is dead (`kill -0` returns non-zero).
  A secondary 3-hour time cap protects against PID reuse on torn-down
  worktrees. This means long benches no longer get stolen at 30 min.
- **Always check/acquire before any cargo or QE command.** If blocked, wait
  (pass `--wait`) or retry — do not force-remove another agent's lock.
- **Release promptly.** Don't hold the lock while reading code or writing
  proposals. `trap 'machine-lock release' EXIT` at session start is a good
  pattern for multi-command sessions.
- Lock state: `.claude/locks/machine.lock.d/` (directory, gitignored). The
  lock records `agent`, `desc`, `ts`, `pid`, and `worktree`. Pre-MLFX flat
  files at `.claude/locks/machine.lock` are auto-cleared on next acquire.
- Shell-test suite: `.claude/bin/tests/machine-lock.test.sh` (17 cases).
  Run whenever touching the lock scripts.

## Conventions

- Cargo.toml uses `>=` version specifiers (not `^` or exact).
- Rust edition 2024. Release profile: opt-level 3, thin LTO.
- Floating-point comparisons use `approx::relative_eq!` in tests.
- Internal units are **eV for energies, Å for lengths, e/Å³ for densities** — used by SCF arrays, `EnergyComponents`, Hamiltonian assembly, and every routine downstream of pseudopotential parsing. Conversion constants in `consts.rs` (`HA_TO_EV`, `RY_TO_EV`, `BOHR_TO_ANG`, `BOHR3_TO_ANG3`, `E2_COULOMB` in eV·Å, `HBAR2_OVER_2M` in eV·Å²). Ry/Bohr appear **only** at the UPF unit boundary in `src/pseudopotential/upf/convert.rs`, which converts PP_LOCAL (Ry→eV), PP_DIJ (Ry→eV), PP_R/PP_RAB (Bohr→Å), PP_RHOATOM (e/Bohr→e/Å), and PP_NLCC (e/Bohr³→e/Å³) before the data enters the engine.
- Input files are YAML (see `examples/`), parsed via serde_yaml_ng into `Settings`.
- Pure Rust stack: faer (eigensolver), ndrustfft (FFT), nalgebra (geometry), ndarray (grid ops).
- No system dependencies required for default build. GPU requires wgpu feature flag.
- Validation scripts in `scripts/` use Python (uv environment).
- QE 7.5 source in `qe-7.5/` for reference during validation.
