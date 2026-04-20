# pwdft-rs

A plane-wave density functional theory (DFT) solver written in pure Rust, with optional GPU acceleration via Apple Metal (wgpu). Results are validated against [Quantum ESPRESSO 7.5](https://www.quantum-espresso.org/).

Designed for real research on macOS (Apple Silicon). The code prioritizes correctness over speed — every energy component is cross-checked against QE reference values.

## Features

- **SCF solver** with four density mixers (Anderson/Pulay DIIS, modified Broyden, Periodic Pulay, plain linear) — any mixer composable with Kerker preconditioning.
- **Norm-conserving pseudopotentials** — reads QE UPF v2 directly; Kleinman-Bylander non-local projectors with arbitrary angular momentum.
- **Nonlinear core correction (NLCC)** — Louie-Froyen-Cohen prescription for pseudopotentials with `core_correction="T"`.
- **Exchange-correlation** — Perdew-Zunger LDA (Ceperley-Alder parametrization), spin-polarized LSDA variant.
- **Crystal symmetry** — automatic space-group detection, IBZ reduction, G-space phase-factor density symmetrization (exact for non-symmorphic groups), time reversal.
- **Collinear spin polarization** — two spin channels with coupled-channel `(ρ_total, m)` mixing so both channels share residual history.
- **Smearing schemes** — Fermi-Dirac, Gaussian, Methfessel-Paxton, Marzari-Vanderbilt cold, fixed occupations.
- **Dual eigensolvers** — full dense (faer `SelfAdjointEigen`, default) and experimental iterative (faer Arnoldi / Krylov-Schur, opt-in).
- **Free-electron band structures** — band paths along high-symmetry directions with TSV output.
- **GPU acceleration** (optional `gpu` feature) — Metal / Vulkan compute via wgpu for Hartree, LDA XC, and V_eff assembly; f32 on GPU, f64 on CPU, with automatic fallback.
- **Ewald summation** — ion-ion electrostatic energy via reciprocal-space Ewald.

## Repository layout

This is a Cargo workspace plus a Python validation package:

```text
pwdft/
  pwdft-core/              Rust solver: library + binary, integration tests, benches
  pwdft-validation/        Python (uv) validation harness: QE parsing + CSV pin generators
  faer/                    Vendored faer v0.24.0 with one local patch (see CLAUDE.md)
data/
  qe/                      Committed QE reference outputs (QE_REF_DIR)
  csv/                     CSV pin files for Rust unit tests (CSV_REF_DIR)
inputs/                    YAML input decks (si_scf.yaml, etc.)
pseudopotentials/          UPF libraries: nc/lda, nc/pbe, uspp/pbe, paw/pbe
qe-7.5/                    Quantum ESPRESSO source + build (symlink, created by setup.sh)
proposals/                 Engineering-proposal backlog (see proposals/INDEX.md)
docs/                      Topic notes on physics + numerics
.claude/                   Agent prompts, skills, and logbooks (agent-driven workflow)
```

## Quick start (Rust)

Requirements: Rust 2024 edition (1.85+). No system dependencies for the default CPU-only build; GPU build needs a Metal- or Vulkan-capable GPU.

```bash
cargo build --release                   # CPU only
cargo build --release --features gpu    # + Metal / Vulkan via wgpu

cargo run --release -- --input inputs/si_scf.yaml
cargo run --release -- --input inputs/si_free_electron.yaml -o bands.tsv
```

CLI:

```text
pwdft-rs --input <path>     Path to YAML input file (required)
         --output <path>    Output file for band structure TSV (default: stdout)
```

Logging:

```bash
RUST_LOG=info cargo run --release -- --input inputs/si_scf.yaml
```

## Quick start (Python validation harness)

`pwdft/pwdft-validation/` is a [uv](https://docs.astral.sh/uv/)-managed Python package that provides the `pwdft-validate` CLI for parsing QE outputs, generating CSV pin files, and cross-checking pwdft-rs against QE. Run `./setup.sh` once to bootstrap the uv environment and symlink QE 7.5.

```bash
uv sync                                      # install deps
uv run pwdft-validate --help                 # list sub-apps
uv run pwdft-validate energy --help          # per-term decomposition
uv run pwdft-validate fermi --help           # Fermi-level bisection reference
uv run pwdft-validate pbe --help             # PBE reference extraction
uv run pwdft-validate reference --help       # regenerate CSV pin files
uv run pwdft-validate density --help         # parse QE charge-density.dat
uv run pwdft-validate diag --help            # diagnostic scripts
```

All paths (`PROJECT_ROOT`, `DATA_DIR`, `QE_REF_DIR`, `CSV_REF_DIR`, `PSEUDO_DIR`, `INPUTS_DIR`) live in `pwdft/pwdft-validation/pwdft_validation/paths.py` — import them instead of deriving paths ad-hoc.

### Running QE itself

QE calculations are invoked directly (the Python harness consumes their output; it does not launch QE). Wrap every QE run in the machine lock so it doesn't contaminate `cargo bench` numbers:

```bash
NP=$(sysctl -n hw.ncpu)   # macOS; use $(nproc) on Linux
.claude/bin/machine-lock run "Researcher" "QE Si SCF reference" -- \
  gtimeout 600 mpirun -np "$NP" qe-7.5/build/bin/pw.x -in si.in > si.out 2>&1
```

See `.claude/skills/qe-runner/SKILL.md` for the full protocol (mandatory parameter limits, pre-flight checklist, output parsing).

## Input format

Input files are YAML. The `system` and `kpoints` sections are required; everything else has sensible defaults following QE conventions. See `inputs/si_scf.yaml` and `inputs/si_free_electron.yaml` for annotated examples covering SCF and band-structure calculations respectively.

## Tests and benchmarks

```bash
cargo test                               # Tier-1 (fast default)
cargo test --features gpu                # + GPU tests
cargo test -- --ignored                  # Tier-2 heavy SCF suites
cargo test -- --include-ignored          # both

cargo bench --bench scf_benchmarks
cargo bench --bench gpu_benchmarks --features gpu
```

Integration tests are in `pwdft/pwdft-core/tests/`. Tier policy: Tier 1 is fast (unit + lightweight integration); Tier 2 runs production-scale SCF loops. Any PR that touches the Tier-2 trigger list (see `CLAUDE.md`) must run `cargo test -- --ignored` and report the outcome.

## Architecture

**Entry point:** `pwdft/pwdft-core/src/main.rs` parses CLI args and YAML input, then computes a free-electron band structure or runs SCF.

**SCF pipeline** (`pwdft/pwdft-core/src/scf/`):

1. Build V_local on the FFT grid. If NLCC, also build ρ_core(r).
2. Initialize density via SAD (superposition of atomic densities).
3. Each iteration: Hartree → LDA XC (on ρ_val + ρ_core if NLCC) → assemble V_eff → Hamiltonian (kinetic + V_eff + KB non-local) → diagonalize → Fermi-Dirac occupations → reconstruct density → G-space symmetrize → convergence check → mix.
4. Total energy = E_kinetic + E_local + E_nonlocal + E_Hartree + E_xc + E_Ewald, with the NLCC double-counting subtraction when core correction is active.

See `CLAUDE.md` § Architecture for the full module tree and the `docs/` folder for topic-focused physics / numerics notes (`basis-and-fft`, `density`, `ewald`, `nonlocal`, `potentials`, `smearing`, `symmetry`, `total-energy`, `units`, etc.).

### Units

pwdft-rs uses **eV / Å / e·Å⁻³** internally. Ry / Bohr appear only at the UPF boundary in `pwdft/pwdft-core/src/pseudopotential/upf/convert.rs`. QE uses Ry / Bohr everywhere; cross-comparison requires conversion.

## Development workflow

Work on pwdft-rs is driven by specialist agents (Engineering Manager, Core Engineer, Performance Engineer, Researcher, Code Reviewer, Technical Writer) coordinating through a branch-and-PR workflow. Agents operate in isolated git worktrees; each PR ships with a per-session logbook entry.

Key entry points:

- `CLAUDE.md` — project conventions, architecture, workflow
- `proposals/INDEX.md` — backlog of engineering proposals
- `.claude/agents/*.md` — per-role prompts and shared protocols
- `.claude/skills/*/SKILL.md` — invocable workflow skills (`/quality-gate`, `/test`, `/bench`, `/profile`, `/lint`, `/pr-submit`, `/merge`, `/proposal`, `/qe-runner`, `/worktree-start`, `/pr-review`, `/cargo`)
- `.claude/logbooks/<role>/` — session handoff notes

## License

No license file is currently committed; contact the maintainers before redistribution.
