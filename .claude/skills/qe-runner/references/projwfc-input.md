# projwfc.x Input Reference (QE 7.5)

projwfc.x projects wavefunctions onto orthogonalized atomic orbitals to compute the projected density of states (PDOS) and atom/orbital-resolved band character. It requires a prior pw.x NSCF calculation (or SCF, but NSCF with a denser k-grid gives smoother DOS).

## Prerequisites

1. Run pw.x with `calculation = 'scf'`
2. Run pw.x with `calculation = 'nscf'` using a denser k-grid and `nbnd` set to include enough empty states. Use the same prefix/outdir.

For k-resolved PDOS (fat bands), use the same k-path as the bands calculation instead of a dense grid.

## Input file structure

```text
&PROJWFC
  ... parameters ...
/
```

---

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `prefix` | string | `'pwscf'` | Must match pw.x prefix |
| `outdir` | string | `'./'` | Must match pw.x outdir |
| `filpdos` | string | `prefix` | Root name for PDOS output files |
| `ngauss` | integer | `0` | Gaussian broadening type: 0=ordinary, 1=Methfessel-Paxton, -1=Fermi-Dirac, -99=tetrahedra |
| `degauss` | real | `0.0` | Broadening width (Ry). Typical: 0.01–0.02 |
| `Emin` | real | auto | Minimum energy for DOS (eV) |
| `Emax` | real | auto | Maximum energy for DOS (eV) |
| `DeltaE` | real | `0.01` | Energy grid spacing (eV) |
| `kresolveddos` | logical | `.false.` | Compute k-resolved PDOS (for fat bands) |
| `filproj` | string | `' '` | File for projection coefficients |
| `lsym` | logical | `.true.` | Symmetrize projections |
| `tdosinboxes` | logical | `.false.` | Compute local DOS in specified real-space boxes |
| `plotboxes` | logical | `.false.` | Output box geometries for visualization |

---

## Example: Total and projected DOS

**Step 1 — NSCF with dense k-grid:**

```fortran
&CONTROL
  calculation = 'nscf'
  prefix      = 'si'
  pseudo_dir  = '../../../pseudo'
  outdir      = './tmp'
  verbosity   = 'high'
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 10.20
  nat       = 2
  ntyp      = 1
  ecutwfc   = 40.0
  nbnd      = 12
  occupations = 'tetrahedra'
/
&ELECTRONS
  conv_thr = 1.0d-10
/
ATOMIC_SPECIES
  Si  28.086  Si_r.upf
ATOMIC_POSITIONS {alat}
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS {automatic}
  12 12 12  0 0 0
```

**Step 2 — projwfc.x:**

```fortran
&PROJWFC
  prefix   = 'si'
  outdir   = './tmp'
  filpdos  = 'si'
  ngauss   = -99
  DeltaE   = 0.01
/
```

Using `ngauss = -99` with `occupations = 'tetrahedra'` in the NSCF step gives the most accurate DOS integration (no broadening artifacts). For Gaussian broadening instead, use `ngauss = 0` and set `degauss`.

---

## Example: k-resolved PDOS (fat bands)

Use the same k-path as a bands calculation:

```fortran
&PROJWFC
  prefix         = 'si'
  outdir         = './tmp'
  filpdos        = 'si.k'
  ngauss         = 0
  degauss        = 0.01
  DeltaE         = 0.01
  kresolveddos   = .true.
  filproj        = 'si.proj.dat'
/
```

---

## Output files

### PDOS files

projwfc.x produces several output files:

| File | Contents |
|------|----------|
| `filpdos.pdos_atm#N(element)_wfc#M(orbital)` | PDOS for atom N, orbital M |
| `filpdos.pdos_tot` | Total DOS |

Example filenames for Si with 2 atoms:

```text
si.pdos_atm#1(Si)_wfc#1(s)
si.pdos_atm#1(Si)_wfc#2(p)
si.pdos_atm#2(Si)_wfc#1(s)
si.pdos_atm#2(Si)_wfc#2(p)
si.pdos_tot
```

### File format

Each PDOS file has columns:

```text
# E (eV)   ldos(E)   pdos(E)
-12.000    0.0000    0.0000
-11.990    0.0001    0.0001
...
```

For spin-polarized calculations, there are separate columns for up and down.

### Projection file (filproj)

If `filproj` is set, projwfc.x writes the projection coefficients |⟨φ_i|ψ_nk⟩|² for every band `n` and k-point `k` onto each atomic orbital `φ_i`. This is useful for orbital-resolved band character analysis.

---

## Parsing for validation

For comparing DOS against pwdft-rs:

- The `.pdos_tot` file gives the total DOS on a uniform energy grid
- Individual PDOS files give atom- and orbital-resolved contributions
- Column format is simple: energy (eV) followed by DOS values (states/eV)

For band character comparison:

- The `filproj` file contains |⟨φ_i|ψ_nk⟩|² — these should sum to ~1 for each (n,k) if the projection basis is complete
