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

### Machine lock — MANDATORY for every QE run

**Every QE invocation must be wrapped in the machine lock.** QE runs (`pw.x`, `mpirun pw.x`, `ph.x`, `pp.x`, `bands.x`, `projwfc.x`, `q2r.x`, `matdyn.x`, `dos.x`) are CPU- and memory-intensive. The machine lock serializes them against `cargo bench` windows so that benchmark wall-time measurements aren't polluted by background CPU load.

**The concern is benchmark integrity, not test isolation.** Two QE calculations can run concurrently against each other without correctness issues — the lock exists so that the Performance Engineer's bench numbers stay clean. Skipping the lock corrupts someone else's measurements silently.

Use the `machine-lock run` one-liner, which handles acquire / release / cleanup on failure:

```bash
# Wrap the full mpirun pw.x invocation — lock covers compile-free runtime only,
# so acquire is cheap and the lock window is as short as the run itself.
.claude/bin/machine-lock run "Researcher" "QE Si SCF reference" -- \
  gtimeout 600 mpirun -np 12 qe-7.5/build/bin/pw.x -in si.in > si.out 2>&1
```

For multi-step pipelines (e.g. scf → bands → bands.x), either wrap the whole shell pipeline in one `machine-lock run` (preferred — holds the lock across the full workflow), or acquire / release explicitly around each step if there's user-facing waiting time between them:

```bash
.claude/bin/machine-lock acquire "Researcher" "QE Si band structure"
gtimeout 600 mpirun -np 12 qe-7.5/build/bin/pw.x -in si.scf.in  > si.scf.out
gtimeout 600 mpirun -np 12 qe-7.5/build/bin/pw.x -in si.bands.in > si.bands.out
gtimeout 600 qe-7.5/build/bin/bands.x        -in si.bandsx.in   > si.bandsx.out
.claude/bin/machine-lock release
```

If the lock is held by another agent, wait and retry — never force-remove another agent's lock. See the "Machine Coordination" section in the project root `CLAUDE.md` for the full policy.

### Runtime configuration

**Every QE invocation must be bounded by a timeout — 10 minutes max.** Use both QE's internal `max_seconds` (graceful stop, writes checkpoint) and an external wall-clock timeout (hard kill) as belt-and-suspenders. If a validation run doesn't finish in 10 min, the parameters are wrong — cut them down rather than extend the timeout.

On macOS, `timeout` is not built in. Install once: `brew install coreutils` (provides `gtimeout`). The examples below assume `gtimeout` is available; fall back to `perl -e 'alarm shift; exec @ARGV' 600 ...` if not.

For **small systems** (< 16 atoms), pure MPI is fastest:
```bash
export OMP_NUM_THREADS=1 LC_ALL=C LANG=C
export OMPI_MCA_btl=self,vader OMPI_MCA_pml=ob1
ulimit -s unlimited
.claude/bin/machine-lock run "Researcher" "QE small-system pw.x" -- \
  gtimeout 600 mpirun -np 12 qe-7.5/build/bin/pw.x -in input.in > output.out 2>&1
```

For **larger systems** (16+ atoms), the installed optimized build may benefit from hybrid MPI+OMP:
```bash
export OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close
export VECLIB_MAXIMUM_THREADS=1 LC_ALL=C LANG=C
.claude/bin/machine-lock run "Researcher" "QE large-system pw.x" -- \
  gtimeout 600 mpirun -np 6 $HOME/qe/bin/pw.x -in input.in > output.out 2>&1
```

Always set `outdir = './tmp'` in input files.

## Runtime safety — mandatory parameter limits

Validation runs must be fast. The goal is a correct reference value, not a production-quality calculation. Use the smallest parameters that still give a physically meaningful result. Excessively long runs block other agents (machine lock), waste wall time, and rarely improve validation precision.

### Always use the PP's suggested cutoff

**Do not guess `ecutwfc`.** Every UPF file carries a recommended cutoff from the author. Use it — don't pad it. pw.x cost scales roughly as `ecutwfc^{1.5}`, so inflating from 40 → 80 Ry nearly triples wall time with no meaningful precision gain for validation.

Extract it from the UPF header (substitute the element you're working with):

```bash
grep -iE 'suggested|cutoff|wfc_cutoff|rho_cutoff' pseudopotentials/nc/pbe/<Element>.upf | head -20
# or for UPF v2 XML:
grep -oE 'wfc_cutoff="[0-9.]+"|rho_cutoff="[0-9.]+"' pseudopotentials/nc/pbe/<Element>.upf
```

Typical values (Ry) once located:

| PP type | ecutwfc | ecutrho | Notes |
|---------|---------|---------|-------|
| ONCV NC (pseudo-dojo) | 30–50 | 4× ecutwfc | Most common for pwdft-rs validation |
| NC HGH / older NC | 40–60 | 4× ecutwfc | |
| USPP (SSSP efficiency) | 25–40 | 8–10× ecutwfc | High ecutrho is mandatory |
| PAW | 30–50 | 8–12× ecutwfc | |

If the UPF doesn't list a suggested cutoff, run a quick convergence sweep (e.g. 20, 30, 40, 50 Ry on a 2-atom cell) and pick the smallest value where total energy changes by <1 meV/atom between steps. Cache the result in a comment at the top of the input file.

### Hard ceilings (exceed only with explicit user approval)

| Parameter | Ceiling | Rationale |
|-----------|---------|-----------|
| External wall timeout | **10 min** (`gtimeout 600`) | Hard kill; prevents runaway jobs |
| `max_seconds` in `&CONTROL` | **540** (9 min) | QE stops gracefully ~60s before hard kill, writing checkpoint |
| `electron_maxstep` in `&ELECTRONS` | **60** | Well-posed SCF converges in 15–30 iters; >60 signals bad mixing/smearing |
| `nstep` in `&CONTROL` (relax) | **30** | Ionic steps; bail early if forces not decreasing |
| Atoms in unit cell | **8** (routine), 16 (absolute max) | pw.x scales ~O(N³); 16 atoms already risks the 10-min cap |
| `ecutwfc` (Ry) | **PP's suggested cutoff, +5 Ry at most** | See section above; never a blanket 80 Ry |
| `ecutrho` (Ry) | **4× ecutwfc** NC, **8–12× ecutwfc** USPP/PAW | Match PP type; too low gives ringing, too high wastes FFT grid |
| k-grid (Monkhorst-Pack) | **≤ 6×6×6** metals, **≤ 4×4×4** insulators | Use symmetry; irreducible k-point count is what matters |
| NSCF dense grid (DOS) | **≤ 12×12×12** | Check runtime first; reduce if >5 min |
| Phonon q-grid (`ldisp`) | **≤ 2×2×2** | Each q-point ≈ one SCF; 4³ routinely busts the cap |

### Required `&CONTROL` settings for every run

```fortran
&CONTROL
  calculation     = 'scf'           ! or bands/nscf/relax/vc-relax
  prefix          = '<material>'
  outdir          = './tmp'
  pseudo_dir      = '<project_root>/pseudopotentials/nc/pbe'
  max_seconds     = 540             ! MANDATORY — QE graceful stop before 10 min hard kill
  disk_io         = 'low'           ! skip unnecessary wavefunction writes
  verbosity       = 'default'       ! use 'high' only when eigenvalues per k are needed
  tstress         = .false.         ! enable only if stress is the validation target
  tprnfor         = .false.         ! enable only if forces are the validation target
/
```

### Required `&ELECTRONS` settings

```fortran
&ELECTRONS
  electron_maxstep = 60
  conv_thr         = 1.0d-8         ! do not tighten below 1d-10 for validation
  mixing_beta      = 0.4            ! 0.7 insulators, 0.3 metals; Kerker helps metals
  mixing_mode      = 'plain'        ! 'local-TF' for charged/metallic
  diagonalization  = 'david'
/
```

### Pre-flight checklist (run before every QE invocation)

1. **Atom count** — ≤8 for routine validation; 16 absolute max.
2. **Pseudopotential** — sourced from `pseudopotentials/` (never `qe-7.5/pseudo/` or elsewhere). XC functional matches `&SYSTEM`.
3. **ecutwfc** — pulled from the UPF's suggested cutoff (see section above), not guessed. `ecutrho` matches the PP type.
3. **k-grid** — smallest grid that still resolves the physics (4³ insulators, 6³ metals for validation).
4. **max_seconds = 540** present in `&CONTROL`? External `gtimeout 600` in the launch command?
5. **disk_io = 'low'** unless a later step needs the wavefunctions?
6. **Estimated cost** — ballpark on 12 MPI ranks: 2-atom sp-bonded insulator, ecutwfc=30, 4³ k ≈ 1–3s; 8-atom cell, ecutwfc=40, 6³ k ≈ 20–60s. Heavier elements, d/f electrons, spin-polarized, magnetic, or metallic systems cost 3–10× more. If projected runtime exceeds 3 min, cut parameters further.

If any item fails, shrink the calculation before launching. When a run is killed by timeout, **do not re-run with a longer timeout** — diagnose why it was slow (too many k-points, SCF not converging, wrong smearing for a metal) and fix the input.

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

**Always pull pseudopotentials from the `pseudopotentials/` tree at the project root.** Do not reach into `qe-7.5/pseudo/`, `$HOME/qe/pseudo/`, or any other location — those are legacy/untracked and may differ from what pwdft-rs validates against. If the element you need isn't in `pseudopotentials/`, stop and tell the user rather than substituting an unvetted PP.

- **For pwdft-rs validation**: use `nc/pbe/` or `nc/lda/` (pwdft-rs implements norm-conserving only). Files are named `<Element>.upf`.
- **For QE-only reference calculations** (no pwdft-rs comparison): `uspp/pbe/` or `paw/pbe/` are fine. File naming follows SSSP conventions.
- **Match the XC functional between PP and `&SYSTEM`** — mixing a PBE PP with LDA input or vice versa gives silently wrong energies.

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
