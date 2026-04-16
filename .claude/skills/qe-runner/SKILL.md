---
name: qe-runner
description: |
  Run Quantum ESPRESSO 7.5 calculations and extract results for validating a Rust planewave DFT implementation. Use this skill whenever: creating QE input files for pw.x / ph.x / pp.x / bands.x / projwfc.x, running any QE calculation, validating DFT results (energies, eigenvalues, forces, stress, phonons, band structure, DOS) against QE, comparing planewave DFT outputs, or parsing QE output. Trigger on any mention of "Quantum ESPRESSO", "QE", "pw.x", "validate against QE", "reference calculation", "run SCF", "band structure", "phonon calculation", "density of states", or requests to check DFT results against a known-good implementation.
---

# Quantum ESPRESSO 7.5 Runner

## Paths

Relative to the pwdft-rs project root:

- **Executables**: `qe-7.5/build/bin/` — `pw.x`, `ph.x`, `pp.x`, `bands.x`, `projwfc.x`, `q2r.x`, `matdyn.x`, `dos.x`
- **Installed pw.x**: `$HOME/qe/bin/pw.x` (optimized build with libxc + veclibfort/Accelerate)
- **Pseudopotentials**: `pseudopotentials/` (primary — NC, USPP, PAW libraries)
- **Working dirs**: create under `qe-7.5/runs/<material>/<calc_type>/`
- **Benchmark data**: `$HOME/qe-bench/results/` (baseline.json, optimized.json, comparison.md)

### Runtime configuration

For **small systems** (< 16 atoms), pure MPI is fastest:
```bash
export OMP_NUM_THREADS=1 LC_ALL=C LANG=C
export OMPI_MCA_btl=self,vader OMPI_MCA_pml=ob1
ulimit -s unlimited
mpirun -np 12 qe-7.5/build/bin/pw.x -in input.in > output.out 2>&1
```

For **larger systems** (16+ atoms), the installed optimized build may benefit from hybrid MPI+OMP:
```bash
export OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close
export VECLIB_MAXIMUM_THREADS=1 LC_ALL=C LANG=C
mpirun -np 6 $HOME/qe/bin/pw.x -in input.in > output.out 2>&1
```

Always set `outdir = './tmp'` in input files.

## Input file reference

Read the relevant reference file before constructing an input:

| Program | Reference | Purpose |
|---------|-----------|---------|
| pw.x | `references/pw-input.md` | SCF, NSCF, bands, relax, vc-relax. All namelists (`&CONTROL`, `&SYSTEM`, `&ELECTRONS`, `&IONS`, `&CELL`), all cards (`ATOMIC_SPECIES`, `ATOMIC_POSITIONS`, `K_POINTS`, `CELL_PARAMETERS`), `ibrav` table, XC functionals. |
| ph.x | `references/ph-input.md` | Phonons via DFPT at single q-points or q-grids. Dispersion workflow with `q2r.x` + `matdyn.x`. |
| pp.x | `references/pp-input.md` | Post-processing: charge density, potentials, wavefunctions. `plot_num` values, output formats, Cube file format. |
| bands.x | `references/bands-input.md` | Extract band energies from pw.x bands calculation. Output file formats. High-symmetry k-paths for FCC/BCC/HEX. |
| projwfc.x | `references/projwfc-input.md` | Projected DOS onto atomic orbitals. k-resolved PDOS for fat bands. |

## Pseudopotentials

All pseudopotentials live in `pseudopotentials/` at the project root. Set `pseudo_dir` in QE input files to the appropriate subdirectory.

### Directory layout

```
pseudopotentials/
  nc/pbe/     — Norm-conserving, PBE (72 elements: H–Zr)
  nc/lda/     — Norm-conserving, LDA (74 elements, includes Fe_dalcorso, Ga_oncv, N_oncv, Si_hgh variants)
  uspp/pbe/   — Ultrasoft, PBE (41 files, SSSP selection)
  paw/pbe/    — PAW, PBE (12 files)
  SSSP_1.3.0_PBE_efficiency.tar.gz   — full SSSP efficiency archive
  SSSP_1.3.0_PBE_precision.tar.gz    — full SSSP precision archive
```

### Which to use

- **For pwdft-rs validation**: use `nc/pbe/` or `nc/lda/` (pwdft-rs implements norm-conserving only). Files are named `<Element>.upf`.
- **For QE-only reference calculations**: `uspp/pbe/` or `paw/pbe/` are fine. File naming varies (SSSP conventions).
- **Legacy PPs** in `qe-7.5/pseudo/` (Si_r.upf, Au_ONCV_PBE_FR_.upf, etc.) still work but prefer the organized `pseudopotentials/` tree.

### Example pseudo_dir in input files

```
pseudo_dir = '<project_root>/pseudopotentials/nc/pbe'
```

## Multi-step workflows

| Goal | Pipeline |
|------|----------|
| Band structure | pw.x (`scf`) → pw.x (`bands`) → bands.x |
| Density of states | pw.x (`scf`) → pw.x (`nscf`, dense k-grid) → projwfc.x |
| Phonons | pw.x (`scf`) → ph.x |
| Phonon dispersion | pw.x (`scf`) → ph.x (`ldisp`) → q2r.x → matdyn.x |
| Charge density | pw.x (`scf`) → pp.x (`plot_num=0`) |

Each step reads from the previous step's `outdir`. Same `prefix` and `outdir` throughout.

## Units

QE uses **Rydberg atomic units**:

| Quantity | QE unit | Conversion |
|----------|---------|------------|
| Energy | Ry | 1 Ry = 0.5 Ha = 13.6057 eV |
| Length | Bohr | 1 Bohr = 0.529177 Å |
| Force | Ry/Bohr | |
| Stress | Ry/Bohr³ | also printed in kbar |
| Eigenvalues | eV (in text output) | Ha in XML |
| Phonon freq | THz and cm⁻¹ | |

## Output parsing

### Key grep patterns (pw.x)

```bash
grep '!    total energy'       output.out   # converged total energy (Ry)
grep 'one-electron'            output.out   # kinetic + local + nonlocal pseudo
grep 'hartree contribution'    output.out   # Hartree energy
grep 'xc contribution'         output.out   # XC energy
grep 'ewald contribution'      output.out   # Ewald ion-ion energy
grep 'convergence has been'    output.out   # verify SCF converged
grep 'force ='                 output.out   # forces (Ry/Bohr)
grep 'total   stress'          output.out   # stress tensor
```

Eigenvalues appear per-k-point when `verbosity = 'high'`:
```
          k = 0.0000 0.0000 0.0000 (   749 PWs)   bands (ev):
    -5.7032   6.2555   6.2555   6.2555
```

### XML (full precision)

`tmp/<prefix>.save/data-file-schema.xml` contains all energies, eigenvalues, k-points, occupations, and cell data at full machine precision. Prefer this for automated validation over text scraping.

### Phonon output (ph.x)

```
     freq (    1) =      15.298080 [THz] =     510.374390 [cm-1]
```

## Build instructions

Two build configurations are available. Use the **fast build** for small validation benchmarks and the **optimized build** when libxc or OpenMP threading is needed.

### Fast build (OpenBLAS, no OpenMP — best for small systems)

```bash
cd qe-7.5 && rm -rf build && mkdir build && cd build
cmake .. -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="/opt/homebrew/opt/openblas;/opt/homebrew/opt/fftw;/opt/homebrew/opt/scalapack" \
  -DBLA_VENDOR=OpenBLAS \
  -DLAPACK_LIBRARIES="/opt/homebrew/opt/openblas/lib/libopenblas.dylib" \
  -DBLAS_LIBRARIES="/opt/homebrew/opt/openblas/lib/libopenblas.dylib" \
  -DQE_ENABLE_MPI=ON -DQE_ENABLE_OPENMP=OFF -DQE_ENABLE_SCALAPACK=ON
make -j$(sysctl -n hw.ncpu) pw pp ph
```

### Optimized build (Accelerate/AMX, OpenMP, libxc — best for larger systems)

Requires: `brew install veclibfort libxc`

**Note**: QE 7.5 CMakeLists.txt requests libxc >= 5.1.2 but uses major-version matching. If libxc 7.x is installed, temporarily patch lines 578+581 in the top-level CMakeLists.txt: change `5.1.2` to `7.0.0` in both `find_package(Libxc ...)` calls. Revert after configure.

**Note**: LTO (`-flto`) fails on macOS Tahoe due to LLVM version mismatch between Apple Clang (C files, LLVM 22.x) and gfortran's linker (LLVM 17.x). Do not use.

```bash
cd qe-7.5 && rm -rf build && mkdir build && cd build
cmake .. -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpifort \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_FLAGS="-O3 -mcpu=apple-m3" \
  -DCMAKE_Fortran_FLAGS="-O3 -mcpu=apple-m3 -mtune=apple-m3 -ffree-line-length-none -funroll-loops" \
  -DQE_ENABLE_OPENMP=ON -DQE_ENABLE_MPI=ON -DQE_ENABLE_SCALAPACK=OFF \
  -DQE_ENABLE_LIBXC=ON -DQE_ENABLE_HDF5=OFF \
  -DBLAS_LIBRARIES=/opt/homebrew/lib/libveclibfort.dylib \
  -DLAPACK_LIBRARIES=/opt/homebrew/lib/libveclibfort.dylib \
  -DCMAKE_PREFIX_PATH="/opt/homebrew/opt/libxc" \
  -DCMAKE_INSTALL_PREFIX=$HOME/qe
make -j$(sysctl -n hw.ncpu) pw pp ph
cp bin/pw.x $HOME/qe/bin/pw.x   # install
```

### Performance notes (M3 Max, Si 2-atom, ecutwfc=50 Ry, 8x8x8 k-grid)

| Config | Build | Wall time |
|--------|-------|-----------|
| MPI=12, OMP=1 | Fast (OpenBLAS) | **1.19s** |
| MPI=12, OMP=1 | Optimized (Accelerate) | 1.77s |
| Serial | Fast (OpenBLAS) | 6.51s |

For small systems, OpenBLAS outperforms Accelerate/AMX because AMX coprocessor startup latency exceeds the tiny BLAS call durations. The crossover favoring Accelerate is around 30+ atoms.
