# pwdft-rs

A plane-wave density functional theory (DFT) solver written in pure Rust, with optional GPU acceleration via Apple Metal (wgpu). Results are validated against [Quantum ESPRESSO 7.5](https://www.quantum-espresso.org/).

Designed for real research on macOS (Apple Silicon). The code prioritizes correctness over speed -- every energy component is cross-checked against QE reference values.

## Features

- **Self-consistent field (SCF) solver** with Anderson/Pulay density mixing and optional Kerker preconditioning
- **Norm-conserving pseudopotentials** -- reads Quantum ESPRESSO UPF v2 format directly
- **Kleinman-Bylander non-local projectors** -- separable form with arbitrary angular momentum via spherical harmonic recurrence
- **Exchange-correlation** -- Perdew-Zunger LDA (Ceperley-Alder parametrization)
- **Crystal symmetry** -- automatic space group detection, irreducible Brillouin zone reduction, density symmetrization, time-reversal symmetry
- **Collinear spin polarization** -- two spin channels with configurable starting magnetization and optional fixed total magnetization
- **Smearing schemes** -- Fermi-Dirac, Gaussian, Methfessel-Paxton, Marzari-Vanderbilt cold smearing, and fixed occupations
- **Free-electron band structures** -- band paths along high-symmetry directions with TSV output
- **GPU acceleration** (optional) -- Metal/Vulkan compute via wgpu for Hartree potential, LDA XC, and V_eff assembly; f32 on GPU, f64 on CPU, with automatic fallback
- **Ewald summation** -- ion-ion electrostatic energy via reciprocal-space Ewald

## Quick Start

### Requirements

- Rust 2024 edition (1.85+)
- No system dependencies for the default (CPU-only) build
- GPU build requires a Metal- or Vulkan-capable GPU

### Build

```bash
# CPU only
cargo build --release

# With GPU acceleration (Metal on macOS, Vulkan elsewhere)
cargo build --release --features gpu
```

### Run an SCF calculation

```bash
cargo run --release -- --input examples/si_scf.yaml
```

Output (printed to stderr):

```
SCF converged in 18 iterations
Total energy: -215.583201 eV
Fermi energy: 6.422810 eV
```

### Compute a band structure

```bash
cargo run --release -- --input examples/si_free_electron.yaml -o bands.tsv
```

This writes a tab-separated file with columns `k_distance` and one column per band, suitable for plotting with gnuplot, matplotlib, or any TSV-aware tool.

### CLI reference

```
pwdft-rs --input <path>     Path to YAML input file (required)
         --output <path>    Output file for band structure TSV (default: stdout)
```

Logging is controlled via the `RUST_LOG` environment variable:

```bash
RUST_LOG=info cargo run --release -- --input examples/si_scf.yaml
```

## Input Format

Input files are YAML. The `system` and `kpoints` sections are required; everything else has sensible defaults following Quantum ESPRESSO conventions.

### SCF calculation

```yaml
# Self-consistent field calculation for bulk Si (FCC diamond, 2 atoms)

system:
  lattice:                              # Lattice vectors in Angstroms
    - [0.0, 2.7155, 2.7155]            #   a1
    - [2.7155, 0.0, 2.7155]            #   a2
    - [2.7155, 2.7155, 0.0]            #   a3
  atoms:                                # Fractional (crystal) coordinates
    - { symbol: Si, position: [0.0, 0.0, 0.0] }
    - { symbol: Si, position: [0.25, 0.25, 0.25] }

basis:
  ecutwfc: 204.09                       # Wavefunction cutoff in eV (= 15 Ry)
  ecutrho_ratio: 4                      # Charge density cutoff = 4 * ecutwfc
  fft_grid: [20, 20, 20]               # Optional: explicit FFT grid dimensions

kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]

scf:
  max_iter: 60                          # Maximum SCF iterations (default: 100)
  conv_threshold: 1.0e-7               # RMS density change in e/A^3 (default: 1e-6)
  energy_threshold: 1.0e-5             # Energy change in eV (default: 1e-5)
  n_bands: 8                           # Number of Kohn-Sham bands (default: auto)

electrons:
  mixing_beta: 0.3                      # Density mixing parameter (default: 0.3)
  mixing_ndim: 8                        # Anderson/Pulay history depth (default: 8)
  mixing_mode: plain                    # plain | kerker (default: plain)
  smearing: fermi_dirac                 # fermi_dirac | gaussian | methfessel_paxton
                                        #   | cold | fixed (default: fermi_dirac)
  smearing_width: 0.05                  # Smearing width in eV (default: 0.05)
  occupations: smearing                 # smearing | fixed (default: smearing)
  nspin: 1                              # 1 = unpolarized, 2 = collinear spin (default: 1)

xc:
  functional: pz                        # pz (LDA, default). pbe/pbe0/hse06 defined but not yet implemented.

symmetry:
  enabled: true                         # Detect and use crystal symmetry (default: true)
  time_reversal: true                   # k -> -k symmetry (default: true)
  tolerance: 1.0e-5                     # Symmetry detection tolerance (default: 1e-5)

pseudopotentials:
  Si: "path/to/Si.upf"                 # UPF v2 files, keyed by element symbol

output:
  verbosity: normal                     # low | normal | high (default: normal)
  write_density: false                  # Write converged density (default: false)
  write_bands: true                     # Write eigenvalues (default: true)
```

### Band structure

```yaml
system:
  lattice:
    - [0.0, 2.7155, 2.7155]
    - [2.7155, 0.0, 2.7155]
    - [2.7155, 2.7155, 0.0]
  atoms:
    - { symbol: Si, position: [0.0, 0.0, 0.0] }
    - { symbol: Si, position: [0.25, 0.25, 0.25] }

basis:
  ecutwfc: 200.0

scf:
  n_bands: 10

kpoints:
  type: band_path
  npoints: 50                           # Points per segment (default: 50)
  path:
    - { label: "L", frac: [0.5, 0.5, 0.5] }
    - { label: "G", frac: [0.0, 0.0, 0.0] }
    - { label: "X", frac: [0.5, 0.0, 0.5] }
    - { label: "K", frac: [0.375, 0.375, 0.75] }
    - { label: "G", frac: [0.0, 0.0, 0.0] }
```

## Architecture

### Source tree

```
src/
  main.rs               CLI entry point (clap). Dispatches to band structure or SCF.
  lib.rs                Public module exports.
  settings.rs           YAML input parsing (serde). All input configuration lives here.
  consts.rs             Physical constants and unit conversions (eV, Bohr, Hartree).

  crystal.rs            Lattice, atoms, and coordinate transforms.
  basis.rs              Plane-wave basis: G-vectors up to kinetic energy cutoff.
  kpoints.rs            Monkhorst-Pack grids and high-symmetry band paths.
  atoms.rs              Element data for Z = 1..92 (symbols, atomic numbers).

  pseudopotential/
    upf.rs              QE UPF v2 parser -> PseudopotentialData (local V, beta projectors, D_ij).

  potential/
    local.rs            Local pseudopotential on FFT grid (spherical Bessel transform).
    hartree.rs          Hartree potential from Poisson equation in reciprocal space.
    xc.rs               Perdew-Zunger LDA exchange-correlation.
    nonlocal.rs         Kleinman-Bylander separable non-local potential.

  scf/
    mod.rs              run_scf() -- the main SCF loop.
    density.rs          Charge density construction from wavefunctions.
    initial_density.rs  Superposition of atomic densities (SAD) for initial guess.
    mixing.rs           Anderson/Pulay mixing with optional Kerker preconditioning.
    smearing.rs         Fermi-Dirac, Gaussian, Methfessel-Paxton, cold smearing.
    energy.rs           Total energy components (kinetic, Hartree, XC, local, non-local, Ewald).
    context.rs          Per-iteration SCF state.
    grid.rs             FFT grid setup from cutoff or explicit dimensions.
    potentials.rs       Hamiltonian assembly (kinetic + V_eff + V_NL).

  eigensolver/
    dense.rs            Hermitian eigendecomposition via faer.

  ewald.rs              Ewald summation for ion-ion electrostatic energy.
  fft.rs                3D FFT wrapper (ndrustfft). Zero unsafe code.
  hamiltonian.rs        Hamiltonian matrix construction.
  bandstructure.rs      Free-electron band structure computation.

  symmetry/
    detect.rs           Space group detection from crystal structure.
    kpoints.rs          k-point reduction to irreducible Brillouin zone.
    density.rs          Charge density symmetrization.

  gpu/                  (behind "gpu" feature flag)
    mod.rs              wgpu compute pipeline, buffer pool, GPU/CPU dispatch.
    shaders/
      hartree.wgsl      Hartree potential kernel.
      lda_xc.wgsl       LDA exchange-correlation kernel.
      v_eff_add.wgsl    Effective potential assembly kernel.
```

### SCF loop

The self-consistent field loop in `src/scf/mod.rs` (`run_scf`) implements the standard Kohn-Sham DFT algorithm:

```
1. Build V_local on FFT grid (spherical Bessel transform of pseudopotential)
2. Initialize electron density via SAD (superposition of atomic densities)
3. For each iteration:
   a. Solve Poisson equation -> V_Hartree
   b. Evaluate XC functional -> V_xc (Perdew-Zunger LDA)
   c. Assemble V_eff = V_local + V_Hartree + V_xc
   d. Build Hamiltonian H = T_kinetic + V_eff + V_NL (Kleinman-Bylander)
   e. Diagonalize H at each k-point (faer, parallelized over k with rayon)
   f. Determine Fermi energy and occupations (Fermi-Dirac smearing)
   g. Reconstruct density from occupied wavefunctions
   h. Check convergence (both density RMS and energy change)
   i. Mix input/output densities (Anderson or Pulay, with optional Kerker)
4. Compute total energy = E_kinetic + E_local + E_nonlocal + E_Hartree + E_xc + E_Ewald
```

Convergence requires both the density change (RMS in e/A^3) and the energy change (eV) to fall below their respective thresholds.

### GPU strategy

The `gpu` feature flag enables wgpu-based GPU acceleration for per-grid-point operations that dominate wall time in larger calculations:

- **Hartree potential** -- reciprocal-space Poisson solve
- **LDA XC** -- Perdew-Zunger exchange-correlation on the real-space grid
- **V_eff assembly** -- summing local + Hartree + XC potentials

GPU kernels run in f32 (Metal/Vulkan compute shaders in WGSL). CPU code uses f64 throughout. Precision conversion happens at the GPU/CPU boundary. A `BufferPool` pre-allocates GPU buffers to avoid per-iteration allocation overhead.

When no GPU is available or the `gpu` feature is disabled, all operations fall back to CPU with rayon parallelism.

### Units

The code uses eV for energies and Angstroms for lengths in user-facing quantities (input, output, logging). Internally, the SCF pipeline works in these same units. Physical constants and conversion factors (Hartree, Rydberg, Bohr) are defined in `src/consts.rs`.

## Dependencies

| Crate | Purpose |
|-------|---------|
| `faer` | Hermitian eigendecomposition (Kohn-Sham diagonalization) |
| `ndrustfft` | 3D FFT (charge density, potentials) |
| `ndarray` | N-dimensional arrays for grid operations |
| `nalgebra` | Linear algebra, lattice geometry |
| `num-complex` | Complex arithmetic for wavefunctions |
| `rayon` | Parallel iteration over k-points |
| `puruspe` | Complementary error function (Ewald summation) |
| `mendeleev` | Periodic table data (elements 1--92) |
| `clap` | CLI argument parsing |
| `serde` + `serde_yaml_ng` | YAML input file parsing |
| `thiserror` | Error type definitions |
| `indicatif` | SCF progress bar |
| `log` + `env_logger` | Structured logging |
| `approx` | Floating-point comparisons in tests |
| `wgpu` | GPU compute (optional, behind `gpu` feature) |
| `pollster` | Async runtime for wgpu (optional) |
| `bytemuck` | Safe transmutes for GPU buffers (optional) |
| `criterion` | Benchmarking (dev dependency) |

All dependencies are pure Rust. No system libraries, BLAS, LAPACK, or FFTW required.

## Tests and Benchmarks

### Running tests

```bash
cargo test                                    # All tests (~177, ~24s)
cargo test --features gpu                     # Include GPU tests (~186, ~28s)
cargo test test_name                          # Single test by name
cargo test --test free_electron_bands         # Single integration test file
cargo test -- --nocapture                     # Show stdout/stderr
```

### Integration test coverage

| Test file | What it validates |
|-----------|-------------------|
| `free_electron_bands.rs` | Free-electron band energies for Si, C diamond, BCC Fe against analytic results |
| `kb_projector_validation.rs` | Kleinman-Bylander non-local projector construction and matrix elements |
| `nonlocal_symmetry.rs` | Non-local potential respects crystal symmetry |
| `parallel_consistency.rs` | Rayon parallelism produces identical results to serial |
| `gpu_consistency.rs` | GPU kernels match CPU results within f32 tolerance |
| `qe_validation.rs` | SCF results compared against Quantum ESPRESSO 7.5 reference values |
| `spin_polarization.rs` | Collinear spin-polarized SCF correctness |
| `fe_debug.rs` | Iron-specific debugging and edge cases |
| `lapack_smoke.rs` | Eigensolver sanity checks |

### Benchmarks

```bash
cargo bench --bench scf_benchmarks                        # SCF iteration benchmarks
cargo bench --bench gpu_benchmarks --features gpu          # GPU kernel benchmarks
```

Benchmarks use [criterion](https://github.com/bheisler/criterion.rs) and generate HTML reports in `target/criterion/`.

## License

*License not yet specified.*
