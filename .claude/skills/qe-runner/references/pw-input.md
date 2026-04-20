# pw.x Input Reference (QE 7.5)

pw.x is the main program for self-consistent calculations, structural optimization, and molecular dynamics.

## Input file structure

```text
&CONTROL
  ... control parameters ...
/
&SYSTEM
  ... system parameters ...
/
&ELECTRONS
  ... electronic iteration parameters ...
/
&IONS          (required for relax, vc-relax, md)
  ... ionic parameters ...
/
&CELL          (required for vc-relax)
  ... cell parameters ...
/
ATOMIC_SPECIES
ATOMIC_POSITIONS
K_POINTS
CELL_PARAMETERS   (required if ibrav = 0)
```

---

## &CONTROL namelist

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `calculation` | string | `'scf'` | `'scf'`, `'nscf'`, `'bands'`, `'relax'`, `'vc-relax'`, `'md'` |
| `prefix` | string | `'pwscf'` | Prepended to output/temporary filenames |
| `pseudo_dir` | string | `'$ESPRESSO_PSEUDO'` | Directory containing pseudopotential files |
| `outdir` | string | `'$ESPRESSO_TMPDIR'` | Directory for large temporary/output files |
| `restart_mode` | string | `'from_scratch'` | `'from_scratch'` or `'restart'` |
| `verbosity` | string | `'low'` | `'low'` or `'high'`. High prints eigenvalues, symmetry ops, etc. |
| `tprnfor` | logical | `.false.` | Print forces on atoms |
| `tstress` | logical | `.false.` | Print stress tensor |
| `forc_conv_thr` | real | `1.0d-3` | Force convergence threshold (Ry/Bohr) for relax |
| `etot_conv_thr` | real | `1.0d-4` | Energy convergence threshold (Ry) for relax |
| `disk_io` | string | `'medium'` | `'none'`, `'low'`, `'medium'`, `'high'` — controls what is written to disk |
| `max_seconds` | real | `1.0d+7` | Max wall time in seconds before checkpointing |
| `nstep` | integer | 50/1 | Max ionic/MD steps (50 for relax, 1 for scf) |

---

## &SYSTEM namelist

### Crystal structure

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ibrav` | integer | — | **Required.** Bravais lattice index (0–14). See table below. |
| `celldm(1)` | real | — | Lattice parameter `a` in **Bohr** (required for ibrav ≠ 0) |
| `celldm(2)` | real | — | `b/a` ratio |
| `celldm(3)` | real | — | `c/a` ratio |
| `celldm(4)` | real | — | cos(α) for ibrav=14, cos(γ) for ibrav=12,13 |
| `celldm(5)` | real | — | cos(β) |
| `celldm(6)` | real | — | cos(γ) for ibrav=14 |
| `A` | real | — | Lattice parameter `a` in **Angstrom** (alternative to celldm) |
| `B, C` | real | — | Lattice parameters in Angstrom |
| `cosAB, cosAC, cosBC` | real | — | Cosines of angles |
| `nat` | integer | — | **Required.** Number of atoms in the unit cell |
| `ntyp` | integer | — | **Required.** Number of atomic types |

### Plane-wave basis

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ecutwfc` | real | — | **Required.** KE cutoff for wavefunctions (Ry) |
| `ecutrho` | real | `4*ecutwfc` | KE cutoff for charge density (Ry). Use 4× for NC PPs, 8–12× for US/PAW |

### Electronic structure

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nbnd` | integer | auto | Number of bands. Default: enough for all electrons + some empty. Set explicitly for NSCF/bands. |
| `occupations` | string | `'fixed'` | `'fixed'` (insulators), `'smearing'` (metals), `'tetrahedra'`, `'tetrahedra_lin'`, `'tetrahedra_opt'` |
| `smearing` | string | `'gaussian'` | `'gaussian'`, `'methfessel-paxton'`, `'marzari-vanderbilt'` (cold), `'fermi-dirac'` |
| `degauss` | real | `0.0` | Smearing width (Ry). Typical: 0.01–0.02 for metals |
| `nspin` | integer | `1` | 1 = unpolarized, 2 = LSDA (collinear spin), 4 = noncollinear |
| `noncolin` | logical | `.false.` | Noncollinear magnetism |
| `lspinorb` | logical | `.false.` | Spin-orbit coupling (requires fully-relativistic PP) |
| `starting_magnetization(i)` | real | `0.0` | Initial magnetization for species `i` (-1 to 1) |
| `tot_charge` | real | `0.0` | Total charge of the system |
| `tot_magnetization` | real | `-1` | Total magnetization (only nspin=2). -1 = unconstrained. |
| `input_dft` | string | — | Override XC functional. Examples: `'PBE'`, `'PZ'`, `'BLYP'`, `'HSE'` |

### Symmetry

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nosym` | logical | `.false.` | Disable spatial symmetry |
| `noinv` | logical | `.false.` | Disable time-reversal symmetry |
| `no_t_rev` | logical | `.false.` | Disable time reversal for noncollinear |

### Other

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `assume_isolated` | string | `'none'` | `'2D'`, `'esm'`, `'makov-payne'` for isolated/2D systems |
| `force_symmorphic` | logical | `.false.` | Force symmorphic symmetry operations only |
| `use_all_frac` | logical | `.false.` | Use all fractional translations |

---

## &ELECTRONS namelist

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `electron_maxstep` | integer | `100` | Max SCF iterations |
| `conv_thr` | real | `1.0d-6` | SCF convergence threshold on total energy (Ry). Use 1e-10 or tighter for validation. |
| `mixing_beta` | real | `0.7` | Mixing factor for SCF. Reduce to 0.1–0.3 for difficult convergence. |
| `mixing_mode` | string | `'plain'` | `'plain'`, `'TF'` (Thomas-Fermi screening), `'local-TF'` |
| `mixing_ndim` | integer | `8` | Number of iterations used in mixing (Broyden) |
| `diagonalization` | string | `'david'` | `'david'` (Davidson), `'cg'` (conjugate gradient), `'ppcg'`, `'paro'` |
| `startingpot` | string | `'atomic'` | `'atomic'` or `'file'` (read from previous calc) |
| `startingwfc` | string | `'atomic+random'` | `'atomic'`, `'atomic+random'`, `'random'`, `'file'` |

---

## &IONS namelist (relax, vc-relax, md only)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ion_dynamics` | string | `'bfgs'` | `'bfgs'` (relax), `'damp'`, `'verlet'` (MD), `'langevin'` |
| `upscale` | real | `100.0` | Factor to reduce conv_thr during optimization |
| `trust_radius_max` | real | `0.8` | Max BFGS displacement (Bohr) |

---

## &CELL namelist (vc-relax only)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cell_dynamics` | string | `'bfgs'` | `'bfgs'`, `'damp-pr'`, `'damp-w'` |
| `press` | real | `0.0` | Target pressure (kbar) |
| `press_conv_thr` | real | `0.5` | Pressure convergence threshold (kbar) |
| `cell_dofree` | string | `'all'` | Which cell parameters to optimize: `'all'`, `'x'`, `'y'`, `'z'`, `'xy'`, `'shape'`, `'volume'`, `'ibrav'` |

---

## Cards

### ATOMIC_SPECIES

```text
ATOMIC_SPECIES
  Symbol  Mass  PseudopotentialFile
```

Example:

```text
ATOMIC_SPECIES
  Si  28.086  Si_r.upf
  C   12.011  C.UPF
```

### ATOMIC_POSITIONS

```text
ATOMIC_POSITIONS {units}
  Symbol  x  y  z  [if_pos(1) if_pos(2) if_pos(3)]
```

Units: `{alat}` (in units of celldm(1)), `{bohr}`, `{angstrom}`, `{crystal}` (fractional coordinates).

Optional `if_pos` flags (0 or 1) control which directions are relaxed.

### K_POINTS

```text
K_POINTS {type}
```

Types:

- **`{automatic}`**: `nk1 nk2 nk3 sk1 sk2 sk3` — Monkhorst-Pack grid with shift
- **`{gamma}`**: Gamma point only (optimized)
- **`{crystal}`**: Explicit list in crystal coordinates. First line = number of k-points. Each line: `kx ky kz weight`
- **`{crystal_b}`**: Band path in crystal coordinates. First line = number of high-symmetry points. Each line: `kx ky kz n_intermediate_points`. The last point has n=1.
- **`{tpiba}`**: Explicit list in units of 2π/a
- **`{tpiba_b}`**: Band path in units of 2π/a

### CELL_PARAMETERS

Required when `ibrav = 0`:

```text
CELL_PARAMETERS {units}
  v1x  v1y  v1z
  v2x  v2y  v2z
  v3x  v3y  v3z
```

Units: `{alat}` (default, scaled by celldm(1)), `{bohr}`, `{angstrom}`.

---

## ibrav reference

| ibrav | Lattice | Required celldm | Lattice vectors (in units of a) |
|-------|---------|------------------|---------------------------------|
| 0 | Free | — | Specified in CELL_PARAMETERS |
| 1 | Cubic P (sc) | (1) | a₁=(1,0,0), a₂=(0,1,0), a₃=(0,0,1) |
| 2 | Cubic F (fcc) | (1) | a₁=(-½,0,½), a₂=(0,½,½), a₃=(-½,½,0) |
| 3 | Cubic I (bcc) | (1) | a₁=(½,½,½), a₂=(-½,½,½), a₃=(-½,-½,½) |
| -3 | Cubic I (bcc) alt | (1) | a₁=(-½,½,½), a₂=(½,-½,½), a₃=(½,½,-½) |
| 4 | Hexagonal | (1),(3) | a₁=(1,0,0), a₂=(-½,√3/2,0), a₃=(0,0,c/a) |
| 5 | Trigonal R (3-fold c) | (1),(4) | cos(α)=celldm(4) |
| -5 | Trigonal R (3-fold ⟨111⟩) | (1),(4) | Alternative orientation |
| 6 | Tetragonal P | (1),(3) | a₁=(1,0,0), a₂=(0,1,0), a₃=(0,0,c/a) |
| 7 | Tetragonal I (bct) | (1),(3) | Body-centered tetragonal |
| 8 | Orthorhombic P | (1),(2),(3) | a₁=(a,0,0), a₂=(0,b,0), a₃=(0,0,c) |
| 9 | Orthorhombic C (bco) | (1),(2),(3) | Base-centered |
| -9 | Orthorhombic C alt | (1),(2),(3) | Alternative base-centered |
| 10 | Orthorhombic F | (1),(2),(3) | Face-centered |
| 11 | Orthorhombic I | (1),(2),(3) | Body-centered |
| 12 | Monoclinic P (c unique) | (1),(2),(3),(4) | celldm(4)=cos(γ) |
| -12 | Monoclinic P (b unique) | (1),(2),(3),(5) | celldm(5)=cos(β) |
| 13 | Monoclinic C (c unique) | (1),(2),(3),(4) | Base-centered |
| -13 | Monoclinic C (b unique) | (1),(2),(3),(5) | Base-centered |
| 14 | Triclinic | (1)–(6) | celldm(4)=cos(α), (5)=cos(β), (6)=cos(γ) |

---

## Common XC functionals (input_dft values)

| Value | Functional | Notes |
|-------|-----------|-------|
| `'PZ'` | LDA Perdew-Zunger | Local density approximation |
| `'PW'` | LDA Perdew-Wang | |
| `'PBE'` | PBE GGA | Most common GGA |
| `'PBESOL'` | PBEsol | Optimized for solids |
| `'BLYP'` | BLYP GGA | Becke exchange + LYP correlation |
| `'BP'` | BP86 | Becke exchange + P86 correlation |
| `'revPBE'` | revised PBE | Zhang-Yang |
| `'HSE'` | HSE06 | Hybrid (expensive) |
| `'PBE0'` | PBE0 | Hybrid (expensive) |

If `input_dft` is not set, QE uses the XC functional specified in the pseudopotential file. For validation, check the output line `Exchange-correlation=` to confirm which functional is actually used.

---

## Complete examples

### Metallic system (Al fcc)

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'al'
  pseudo_dir  = '../../../pseudo'
  outdir      = './tmp'
  verbosity   = 'high'
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 7.65
  nat       = 1
  ntyp      = 1
  ecutwfc   = 30.0
  occupations = 'smearing'
  smearing    = 'marzari-vanderbilt'
  degauss     = 0.02
/
&ELECTRONS
  conv_thr = 1.0d-10
/
ATOMIC_SPECIES
  Al  26.982  Al.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS {alat}
  Al  0.00  0.00  0.00
K_POINTS {automatic}
  12 12 12  0 0 0
```

### Structural relaxation

```fortran
&CONTROL
  calculation   = 'relax'
  prefix        = 'sic'
  pseudo_dir    = '../../../pseudo'
  outdir        = './tmp'
  tprnfor       = .true.
  forc_conv_thr = 1.0d-5
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 8.24
  nat       = 2
  ntyp      = 2
  ecutwfc   = 40.0
/
&ELECTRONS
  conv_thr = 1.0d-10
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
  Si  28.086  Si_r.upf
  C   12.011  C.UPF
ATOMIC_POSITIONS {alat}
  Si  0.00  0.00  0.00
  C   0.25  0.25  0.25
K_POINTS {automatic}
  6 6 6  1 1 1
```

### No symmetry (for validation without symmetry reduction)

Add to &SYSTEM:

```fortran
  nosym  = .true.
  noinv  = .true.
```

And specify the full k-point grid. With `K_POINTS {automatic}` and nosym, QE will use all nk1×nk2×nk3 k-points without reduction.
