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
cargo clippy -q --all-targets  # check remaining warnings
cargo test                     # verify nothing broke
```
Fix auto-fixable warnings, address remaining ones. Do not suppress codesmell warnings like `too_many_arguments` — refactor the code instead.

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
