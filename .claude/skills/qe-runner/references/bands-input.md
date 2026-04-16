# bands.x Input Reference (QE 7.5)

bands.x is a post-processing tool that reads band energies from a prior pw.x `calculation = 'bands'` run and writes them to a file suitable for plotting or comparison.

## Prerequisites

1. Run pw.x with `calculation = 'scf'` (self-consistent charge density)
2. Run pw.x with `calculation = 'bands'` using the same prefix/outdir, with `K_POINTS {crystal_b}` or `{tpiba_b}` specifying the k-path

## Input file structure

```
&BANDS
  ... parameters ...
/
```

---

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `prefix` | string | `'pwscf'` | Must match pw.x prefix |
| `outdir` | string | `'./'` | Must match pw.x outdir |
| `filband` | string | `'bands.out'` | Output file for band data |
| `lsym` | logical | `.true.` | Classify bands by symmetry at high-symmetry k-points |
| `lsigma(1)` | logical | `.false.` | Write spin components of expectation value of σ_x |
| `lsigma(2)` | logical | `.false.` | σ_y component |
| `lsigma(3)` | logical | `.false.` | σ_z component |
| `spin_component` | integer | `0` | For LSDA: 1=up, 2=down |
| `no_overlap` | logical | `.true.` | Don't compute overlap matrix for band ordering |
| `plot_2d` | logical | `.false.` | 2D band structure plot |

---

## Example

```
&BANDS
  prefix  = 'si'
  outdir  = './tmp'
  filband = 'si_bands.dat'
  lsym    = .true.
/
```

---

## Output files

### filband

The main output contains eigenvalues in eV for each k-point along the path. Format:

```
 &plot nbnd=   8, nks=  81 /
          -0.500000  0.500000  0.500000
  -2.0025   5.7234   5.7234   5.7234   8.5671   8.5671   8.5671  13.4892
          -0.475000  0.475000  0.475000
  -1.9418   4.4321   5.8901   5.8901   8.4192   8.5038   8.5038  13.3067
  ...
```

First line is a namelist with `nbnd` (number of bands) and `nks` (number of k-points).
Then alternating: k-point coordinates (Cartesian, 2π/a) followed by eigenvalues (eV) for that k-point.

### filband.gnu

A gnuplot-ready file with two columns:
```
k_position  eigenvalue_eV
```

Each band is a separate block separated by blank lines. Plot with:
```gnuplot
plot 'si_bands.dat.gnu' with lines
```

### filband.rap

If `lsym = .true.`, this file contains the symmetry representations of each band at high-symmetry points.

---

## Parsing for validation

To extract eigenvalues programmatically from the `.gnu` file:
- Column 1: cumulative k-path distance (in 2π/a units)
- Column 2: eigenvalue (eV)
- Blank line separates bands

For direct numerical comparison, the XML file at `tmp/<prefix>.save/data-file-schema.xml` contains eigenvalues at all k-points with full precision (Hartree units in the XML).

---

## High-symmetry k-points for common lattices

### FCC (Si, GaAs, Al, Cu, etc.)

```
K_POINTS {crystal_b}
5
  0.500  0.500  0.500  20  ! L
  0.000  0.000  0.000  20  ! Gamma
  0.500  0.000  0.500  10  ! X
  0.625  0.250  0.625  20  ! U|K
  0.000  0.000  0.000   1  ! Gamma
```

### BCC (Fe, W, Na, etc.)

```
K_POINTS {crystal_b}
4
  0.000  0.000  0.000  20  ! Gamma
  0.500 -0.500  0.500  20  ! H
  0.250  0.250  0.250  20  ! P
  0.000  0.000  0.000   1  ! Gamma
```

### Hexagonal (graphene, BN, etc.)

```
K_POINTS {crystal_b}
4
  0.000  0.000  0.000  20  ! Gamma
  0.500  0.000  0.000  20  ! M
  0.333  0.333  0.000  20  ! K
  0.000  0.000  0.000   1  ! Gamma
```
