---
name: qe-runner
description: |
  Run Quantum ESPRESSO 7.5 calculations and extract results for validating pwdft-rs. Use this skill whenever: creating QE input files for pw.x / ph.x / pp.x / bands.x / projwfc.x, running any QE calculation, validating DFT results (energies, eigenvalues, forces, stress, phonons, band structure, DOS) against QE, or parsing QE output. Trigger on "Quantum ESPRESSO", "QE", "pw.x", "validate against QE", "reference calculation", "run SCF", "band structure", "phonon calculation", "density of states", or requests to check pwdft-rs against QE.
---

# Quantum ESPRESSO 7.5 runner

## Paths (relative to project root)

- **Executables:** `qe-7.5/build/bin/` — `pw.x`, `ph.x`, `pp.x`, `bands.x`, `projwfc.x`, `q2r.x`, `matdyn.x`, `dos.x`. `qe-7.5/` is a symlink created by `./setup.sh`; if missing, fix the symlink rather than hard-coding absolute paths.
- **Installed optimized `pw.x`:** `$HOME/qe/bin/pw.x` (libxc + veclibfort / Accelerate; optional, used for larger systems).
- **Pseudopotentials:** `pseudopotentials/` at the project root — see § Pseudopotentials.
- **QE reference outputs:** `data/qe/` — committed `.out` logs that pwdft-rs validates against. `pwdft/pwdft-validation/pwdft_validation/paths.py` exposes this as `QE_REF_DIR`.
- **Input decks:** `inputs/` for pwdft-rs YAML; QE `.in` files for reference calculations live next to their outputs under `data/qe/`.
- **Working dirs:** create under `qe-7.5/runs/<material>/<calc_type>/`. Set `outdir = './tmp'` inside the input.
- **Build instructions:** `docs/qe-build.md` (not this skill — you build QE once, not per validation).

## Python validation harness — `pwdft-validate`

The canonical way to extract and compare QE results is the `pwdft-validate` CLI (`pwdft/pwdft-validation/pwdft_validation/`), installed as a uv entry point. Prefer it over hand-rolled grep where a sub-command exists:

```bash
uv run pwdft-validate --help          # list sub-apps
uv run pwdft-validate energy --help   # per-term decomposition (VGC5 / VGCH-2)
uv run pwdft-validate fermi --help    # Fermi-level bisection reference
uv run pwdft-validate pbe --help      # PBE reference extraction from QE logs
uv run pwdft-validate reference --help  # CSV pin file generators (vloc, beta, D_ij, hamiltonian, NLCC, SAD)
uv run pwdft-validate density --help  # parse QE charge-density.dat → binary
uv run pwdft-validate diag --help     # cross-check diagnostics
```

`paths.py` holds `PROJECT_ROOT`, `DATA_DIR`, `QE_REF_DIR`, `CSV_REF_DIR`, `PSEUDO_DIR`, `INPUTS_DIR` — scripts import these rather than deriving paths, so reorganizations don't break them.

## Machine lock — mandatory for every QE run

Wrap every invocation in the machine lock — QE saturates CPU and silently corrupts anyone else's `cargo bench` numbers otherwise. Use the `run` one-liner:

```bash
NP=${NP:-$(sysctl -n hw.ncpu)}  # macOS; use $(nproc) on Linux
.claude/bin/machine-lock run "Researcher" "QE Si SCF reference" -- \
  gtimeout 600 mpirun -np "$NP" qe-7.5/build/bin/pw.x -in si.in > si.out 2>&1
```

For multi-step pipelines (scf → bands → bands.x), wrap the whole pipeline in one `machine-lock run` so the lock covers the full workflow.

If the lock is held by another agent, wait (`acquire --wait`) or retry. Never force-remove another agent's lock. See `.claude/agents/shared/machine-lock.md` for the full policy.

## Runtime configuration

**10-minute hard cap per run.** Use both QE's `max_seconds = 540` in `&CONTROL` (graceful stop, writes checkpoint at ~9 min) and an external `gtimeout 600` (hard kill) as belt-and-suspenders. If a validation run doesn't finish in 10 min, the parameters are wrong — shrink them, don't extend the timeout. On macOS, `timeout` is not built in; `brew install coreutils` provides `gtimeout`.

**Detect core count once per shell session.** Never hard-code MPI ranks or OMP threads. On macOS: `NP=$(sysctl -n hw.ncpu)`; on Linux: `NP=$(nproc)`.

For **small systems** (< 16 atoms), pure MPI is fastest:

```bash
NP=${NP:-$(sysctl -n hw.ncpu)}
export OMP_NUM_THREADS=1 LC_ALL=C LANG=C
export OMPI_MCA_btl=self,vader OMPI_MCA_pml=ob1
ulimit -s unlimited
.claude/bin/machine-lock run "Researcher" "QE small-system pw.x" -- \
  gtimeout 600 mpirun -np "$NP" qe-7.5/build/bin/pw.x -in input.in > output.out 2>&1
```

For **larger systems** (16+ atoms), hybrid MPI+OMP with the optimized build may help. Split `$NP` into `RANKS × OMP`:

```bash
NP=${NP:-$(sysctl -n hw.ncpu)}
OMP=${OMP:-2}
RANKS=$(( NP / OMP ))
export OMP_NUM_THREADS="$OMP" OMP_PLACES=cores OMP_PROC_BIND=close
export VECLIB_MAXIMUM_THREADS=1 LC_ALL=C LANG=C
.claude/bin/machine-lock run "Researcher" "QE large-system pw.x" -- \
  gtimeout 600 mpirun -np "$RANKS" $HOME/qe/bin/pw.x -in input.in > output.out 2>&1
```

## Mandatory parameter limits

Validation runs are for correct reference values, not production-quality calculations. Use the smallest parameters that still give a meaningful result.

**Cutoffs — always use the PP's suggested cutoff.** pw.x cost scales ≈ `ecutwfc^1.5`; inflating 40 → 80 Ry nearly triples wall time for no precision gain.

```bash
grep -oE 'wfc_cutoff="[0-9.]+"|rho_cutoff="[0-9.]+"' pseudopotentials/nc/pbe/<Element>.upf
```

| PP type | ecutwfc (Ry) | ecutrho | Notes |
|---------|-------------|---------|-------|
| ONCV NC (pseudo-dojo) | 30–50 | 4 × ecutwfc | Most common for pwdft-rs |
| NC HGH / older NC | 40–60 | 4 × ecutwfc | |
| USPP (SSSP efficiency) | 25–40 | 8–10 × ecutwfc | High ecutrho mandatory |
| PAW | 30–50 | 8–12 × ecutwfc | |

**Hard ceilings** (exceed only with explicit user approval):

| Parameter | Ceiling | Rationale |
|-----------|---------|-----------|
| External wall timeout | **10 min** (`gtimeout 600`) | Hard kill |
| `max_seconds` in `&CONTROL` | **540** (9 min) | QE graceful stop before hard kill |
| `electron_maxstep` in `&ELECTRONS` | **60** | Good SCFs converge in 15–30 iters |
| `nstep` in `&CONTROL` (relax) | **30** | Bail early if forces not decreasing |
| Atoms in unit cell | **8** routine, 16 absolute max | pw.x scales ~O(N³) |
| `ecutwfc` | **PP's suggested cutoff, +5 Ry at most** | |
| `ecutrho` | **4× ecutwfc NC, 8–12× USPP/PAW** | Match PP type |
| k-grid (Monkhorst-Pack) | **≤ 6³ metals, ≤ 4³ insulators** | |
| NSCF dense grid (DOS) | **≤ 12³** | Check runtime first |
| Phonon q-grid (`ldisp`) | **≤ 2³** | Each q-point ≈ one SCF |

### Required `&CONTROL` settings

```fortran
&CONTROL
  calculation     = 'scf'           ! or bands/nscf/relax/vc-relax
  prefix          = '<material>'
  outdir          = './tmp'
  pseudo_dir      = '<project_root>/pseudopotentials/nc/pbe'
  max_seconds     = 540             ! MANDATORY graceful-stop before hard kill
  disk_io         = 'low'           ! skip unnecessary wavefunction writes
  verbosity       = 'default'       ! 'high' only when per-k eigenvalues needed
  tstress         = .false.         ! enable only if stress is the validation target
  tprnfor         = .false.         ! enable only if forces are the validation target
/
```

### Required `&ELECTRONS` settings

```fortran
&ELECTRONS
  electron_maxstep = 60
  conv_thr         = 1.0d-8         ! do not tighten below 1d-10 for validation
  mixing_beta      = 0.4            ! 0.7 insulators, 0.3 metals
  mixing_mode      = 'plain'        ! 'local-TF' for charged / metallic
  diagonalization  = 'david'
/
```

### Pre-flight checklist

1. **Atom count** ≤ 8 for routine validation; 16 absolute max.
2. **Pseudopotential** from `pseudopotentials/` only. XC matches `&SYSTEM`.
3. **ecutwfc** from the UPF's suggested cutoff. `ecutrho` matches PP type.
4. **k-grid** smallest that still resolves the physics.
5. **`max_seconds = 540`** in `&CONTROL`? External `gtimeout 600` in the launch command?
6. **`disk_io = 'low'`** unless a later step needs the wavefunctions?
7. **Estimated cost** — on 12 MPI ranks: 2-atom sp-bonded insulator ecut 30, 4³ k ≈ 1–3 s; 8-atom cell ecut 40, 6³ k ≈ 20–60 s. Heavy elements, d/f electrons, spin-polarized, magnetic, and metallic systems cost 3–10× more. Projected runtime > 3 min ⇒ shrink further.

If any item fails, shrink before launching. A run killed by timeout means the inputs were wrong — diagnose (too many k-points, SCF not converging, wrong smearing for a metal), don't extend the timeout.

## Pseudopotentials

All PPs live in `pseudopotentials/` at the project root. Set `pseudo_dir` in QE input files to the appropriate subdirectory.

```text
pseudopotentials/
  nc/pbe/     — Norm-conserving, PBE (72 elements: H–Zr)
  nc/lda/     — Norm-conserving, LDA (74 elements, includes Fe_dalcorso, Ga_oncv, N_oncv, Si_hgh variants)
  uspp/pbe/   — Ultrasoft, PBE (41 files, SSSP selection)
  paw/pbe/    — PAW, PBE (12 files)
  SSSP_1.3.0_PBE_{efficiency,precision}.tar.gz — full archives
```

- **For pwdft-rs validation:** use `nc/pbe/` or `nc/lda/` — pwdft-rs implements norm-conserving only. Files are named `<Element>.upf`.
- **For QE-only reference calculations:** `uspp/pbe/` or `paw/pbe/` are fine.
- **Match the XC functional between PP and `&SYSTEM`** — mixing PBE PP with LDA input (or vice versa) gives silently wrong energies.
- **Never** reach into `qe-7.5/pseudo/`, `$HOME/qe/pseudo/`, or any other location — those may differ from what pwdft-rs validates against. If an element isn't in `pseudopotentials/`, stop and tell the user.

## Input-file reference

Read the relevant reference file before constructing an input:

| Program | Reference | Purpose |
|---------|-----------|---------|
| pw.x | `references/pw-input.md` | SCF, NSCF, bands, relax, vc-relax. All namelists, cards, `ibrav`, XC functionals. |
| ph.x | `references/ph-input.md` | Phonons via DFPT. Dispersion via `q2r.x` + `matdyn.x`. |
| pp.x | `references/pp-input.md` | Charge density, potentials, wavefunctions. `plot_num` values, Cube format. |
| bands.x | `references/bands-input.md` | Band energies from pw.x bands run. High-symmetry k-paths for FCC/BCC/HEX. |
| projwfc.x | `references/projwfc-input.md` | Projected DOS onto atomic orbitals. |

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

QE uses Rydberg atomic units:

| Quantity | QE unit | Conversion |
|----------|---------|------------|
| Energy | Ry | 1 Ry = 0.5 Ha = 13.6057 eV |
| Length | Bohr | 1 Bohr = 0.529177 Å |
| Force | Ry/Bohr | |
| Stress | Ry/Bohr³ | also printed in kbar |
| Eigenvalues | eV (text output) | Ha (XML) |
| Phonon freq | THz and cm⁻¹ | |

pwdft-rs's internal units are eV / Å. All conversion happens at the UPF boundary in `pwdft/pwdft-core/src/pseudopotential/upf/convert.rs`.

## Output parsing

Prefer `pwdft-validate` sub-commands (§ Python validation harness). Where no sub-command exists, use `data-file-schema.xml` for full precision — never hand-scrape the text output for production validation.

### Text output (pw.x) — quick grep patterns

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

Eigenvalues with `verbosity = 'high'`:

```text
          k = 0.0000 0.0000 0.0000 (   749 PWs)   bands (ev):
    -5.7032   6.2555   6.2555   6.2555
```

### XML (full precision)

`tmp/<prefix>.save/data-file-schema.xml` — all energies, eigenvalues, k-points, occupations, cell data at machine precision. The `pwdft-validate` parsers read the XML where available; text scraping is fallback only.

### Phonon output (ph.x)

```text
     freq (    1) =      15.298080 [THz] =     510.374390 [cm-1]
```
