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
cargo test                                    # all tests (~177, ~24s)
cargo test --features gpu                     # with GPU tests (~186, ~28s)
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
cargo test                                     # verify nothing broke
```
Both clippy invocations are required: without `--features gpu`, the `gpu/` source tree and the GPU-only test binaries are not linted, so warnings accumulate silently (QLN2 caught an `uninlined_format_args` violation in `tests/gpu_consistency.rs` that had slipped past the default-feature gate for weeks). Fix auto-fixable warnings, address remaining ones. Do not suppress codesmell warnings like `too_many_arguments` — refactor the code instead.

## Architecture

**Entry point:** `src/main.rs` parses CLI args and YAML input (`src/settings.rs`), then either computes a free-electron band structure or runs SCF.

**SCF loop** (`src/scf/mod.rs` — `run_scf()`): the central computation pipeline:
1. Build local pseudopotential V_local on FFT grid (spherical Bessel transform)
2. Initialize density via SAD (superposition of atomic densities)
3. Each iteration: Hartree potential → LDA XC → assemble V_eff → build Hamiltonian (kinetic + V_eff + KB non-local) → diagonalize (faer) → Fermi-Dirac occupations → reconstruct density → check convergence → Anderson/Pulay mixing (with optional Kerker preconditioning)
4. Compute total energy (kinetic + local + non-local + Hartree + XC + Ewald)

**Module groups:**

- **Crystal & basis:** `crystal.rs` (lattice + atoms), `basis.rs` (G-vectors up to ecut), `kpoints.rs` (Monkhorst-Pack, band paths), `atoms.rs` (elements 1-92)
- **Pseudopotentials:** `pseudopotential/upf.rs` (QE UPF v2). Parses into `PseudopotentialData` with local potential, beta projectors, D_ij matrix.
- **Potentials:** `potential/hartree.rs`, `potential/xc.rs` (Perdew-Zunger LDA), `potential/local.rs`, `potential/nonlocal.rs` (Kleinman-Bylander separable form, arbitrary l via recurrence)
- **SCF internals:** `scf/density.rs`, `scf/initial_density.rs` (SAD), `scf/mixing.rs` (Anderson/Pulay + Kerker preconditioning), `scf/smearing.rs` (Fermi-Dirac)
- **Numerics:** `fft.rs` (3D FFT via ndrustfft, zero unsafe), `eigensolver/dense.rs` (faer Hermitian eigendecomposition), `ewald.rs` (ion-ion energy, erfc via puruspe)
- **Symmetry:** `symmetry/detect.rs` (space group finder), `symmetry/kpoints.rs` (k-point reduction), `symmetry/density.rs` (density symmetrization)
- **GPU:** `gpu/mod.rs` (wgpu compute), `gpu/shaders/` (WGSL kernels for Hartree, LDA XC, V_eff assembly)

**GPU strategy:** Optional `gpu` feature flag. `GpuAccelerator` with `BufferPool` of pre-allocated f32 buffers. GPU kernels run in f32, CPU in f64, with conversion at boundaries. Falls back to CPU (rayon) when GPU unavailable. Three WGSL shaders handle the per-iteration grid operations.

**Key types:** `Crystal`, `BasisSet`, `KPoint`, `PseudopotentialData`, `ScfParams`/`ScfResult`, `NonlocalPotential`, `EigenResult`, `SymmetryInfo`/`SpaceGroupOp`.

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

Multiple agents share this machine. The machine lock serializes CPU-intensive commands to prevent test runs from contaminating benchmark results.

**What requires the lock:** `cargo test`, `cargo bench`, `cargo build`, `cargo clippy` — any command that compiles or runs code.

**What does NOT require the lock:** reading files, editing code in worktrees, writing proposals, git operations, `machine-lock status`.

```bash
# Check lock state:
.claude/bin/machine-lock status

# Acquire before cargo commands:
.claude/bin/machine-lock acquire "Core Engineer" "cargo test"
cargo test
.claude/bin/machine-lock release

# Or use the one-liner (acquire + run + release):
.claude/bin/machine-lock run "Core Engineer" "cargo test" -- cargo test
```

- **Always check/acquire before any cargo command.** If blocked, wait and retry — do not force-remove another agent's lock.
- **Release promptly.** Don't hold the lock while reading code or writing proposals.
- **Stale locks** (>30 min old) are auto-cleared on next acquire.
- Lock state: `.claude/locks/machine.lock` (gitignored).

## Conventions

- Cargo.toml uses `>=` version specifiers (not `^` or exact).
- Rust edition 2024. Release profile: opt-level 3, thin LTO.
- Floating-point comparisons use `approx::relative_eq!` in tests.
- Physical constants in `consts.rs` (Hartree atomic units: energies in Ry, lengths in Bohr).
- Input files are YAML (see `examples/`), parsed via serde_yaml_ng into `Settings`.
- Pure Rust stack: faer (eigensolver), ndrustfft (FFT), nalgebra (geometry), ndarray (grid ops).
- No system dependencies required for default build. GPU requires wgpu feature flag.
- Validation scripts in `scripts/` use Python (uv environment).
- QE 7.5 source in `qe-7.5/` for reference during validation.
