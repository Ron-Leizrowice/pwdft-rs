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
cargo run --release -- --input examples/si_scf.toml

# Run band structure
cargo run --release -- --input examples/si_free_electron.toml -o bands.tsv
```

LAPACK uses the Accelerate framework on macOS (`lapack-accelerate` feature in nalgebra-lapack).

## Tests & Benchmarks

```bash
cargo test                                    # all tests
cargo test test_name                          # single test by name
cargo test --test free_electron_bands         # single integration test file
cargo test -- --nocapture                     # with stdout

cargo bench --bench scf_benchmarks            # SCF benchmarks
cargo bench --bench gpu_benchmarks --features gpu  # GPU benchmarks
```

Integration tests in `tests/`: free-electron band validation, LAPACK smoke, GPU vs CPU consistency, KB projector validation, non-local symmetry, parallel consistency.

## Architecture

**Entry point:** `src/main.rs` parses CLI args and TOML input, then either computes a free-electron band structure or runs SCF.

**SCF loop** (`src/scf/mod.rs` — `run_scf()`): the central computation pipeline:
1. Build local pseudopotential V_local on FFT grid (spherical Bessel transform)
2. Initialize density via SAD (superposition of atomic densities)
3. Each iteration: Hartree potential → LDA XC → assemble V_eff → build Hamiltonian (kinetic + V_eff + KB non-local) → diagonalize (LAPACK zheev) → Fermi-Dirac occupations → reconstruct density → check convergence → Anderson/Pulay mixing
4. Compute total energy (kinetic + local + non-local + Hartree + XC + Ewald)

**Module groups:**

- **Crystal & basis:** `crystal.rs` (lattice + atoms), `basis.rs` (G-vectors up to ecut), `kpoints.rs` (Monkhorst-Pack, band paths), `atoms.rs` (elements 1-92)
- **Pseudopotentials:** `pseudopotential/upf.rs` (QE UPF v2), `pseudopotential/psp8.rs` (ABINIT/PseudoDojo). Both parse into `PseudopotentialData` with local potential, beta projectors, D_ij matrix.
- **Potentials:** `potential/hartree.rs`, `potential/xc.rs` (Perdew-Zunger LDA), `potential/local.rs`, `potential/nonlocal.rs` (Kleinman-Bylander separable form)
- **SCF internals:** `scf/density.rs`, `scf/initial_density.rs` (SAD), `scf/mixing.rs` (Anderson/Pulay), `scf/smearing.rs` (Fermi-Dirac)
- **Numerics:** `fft.rs` (3D FFT via rustfft, rayon-parallel batch 1D), `eigensolver/dense.rs` (LAPACK zheev wrapper), `ewald.rs` (ion-ion energy)
- **Symmetry:** `symmetry/detect.rs` (space group finder), `symmetry/kpoints.rs` (k-point reduction), `symmetry/density.rs` (density symmetrization)
- **GPU:** `gpu/mod.rs` (wgpu compute), `gpu/shaders/` (WGSL kernels for Hartree, LDA XC, V_eff assembly)

**GPU strategy:** Optional `gpu` feature flag. `GpuAccelerator` with `BufferPool` of pre-allocated f32 buffers. GPU kernels run in f32, CPU in f64, with conversion at boundaries. Falls back to CPU (rayon) when GPU unavailable. Three WGSL shaders handle the per-iteration grid operations.

**Key types:** `Crystal`, `BasisSet`, `KPoint`, `PseudopotentialData`, `ScfParams`/`ScfResult`, `NonlocalPotential`, `EigenResult`, `SymmetryInfo`/`SpaceGroupOp`.

## Conventions

- Cargo.toml uses `>=` version specifiers (not `^` or exact).
- Rust edition 2024. Release profile: opt-level 3, thin LTO.
- Floating-point comparisons use `approx::relative_eq!` in tests.
- Physical constants in `consts.rs` (Hartree atomic units: energies in Ry, lengths in Bohr).
- Input files are TOML (see `examples/`). Settings/config also supports YAML via serde_yaml_ng.
- Validation scripts in `scripts/` use Python (uv environment).
- QE 7.5 source in `qe-7.5/` for reference during validation.
