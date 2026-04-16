# pp.x Input Reference (QE 7.5)

pp.x extracts and processes data from a prior pw.x calculation — charge density, potentials, wavefunctions, etc. It reads from `outdir/prefix.save/` and writes to plot files.

## Input file structure

```
&INPUTPP
  ... what to extract ...
/
&PLOT
  ... how to format the output ...
/
```

---

## &INPUTPP parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `prefix` | string | `'pwscf'` | Must match pw.x prefix |
| `outdir` | string | `'./'` | Must match pw.x outdir |
| `plot_num` | integer | — | **Required.** What quantity to extract (see table below) |
| `filplot` | string | `'tmp.pp'` | Intermediate output file |
| `spin_component` | integer | `0` | 0 = total, 1 = up, 2 = down (for nspin=2) |
| `kpoint(1)` | integer | `1` | k-point index (for plot_num=7) |
| `kband(1)` | integer | `1` | Band index (for plot_num=7) |

### plot_num values

| Value | Quantity | Notes |
|-------|----------|-------|
| 0 | Charge density ρ(r) | Most common for validation |
| 1 | Total potential V_tot(r) | V_local + V_Hartree + V_xc |
| 2 | Local ionic potential V_local(r) | From pseudopotential |
| 3 | Local density of states at E_F | |
| 4 | Local density of electronic entropy | |
| 5 | Spin polarization ρ↑-ρ↓ | nspin=2 only |
| 6 | Spin polarization (absolute) | nspin=2 only |
| 7 | |ψ(r)|² for specific k,band | Specify kpoint, kband |
| 8 | Electron localization function (ELF) | |
| 9 | Charge density minus superposition of atomic densities | Δρ |
| 10 | ILDOS (integrated local DOS) in energy range | Specify emin, emax |
| 11 | Electrostatic potential V_bare + V_Hartree | No XC |
| 12 | Sawtooth potential (for E-field) | |
| 13 | Noncollinear magnetization (all) | |
| 17 | All-electron charge density (PAW) | PAW only |
| 19 | Reduced density gradient | |
| 20 | Electrostatic potential on 1D grid | |
| 21 | XC potential V_xc(r) | |
| 22 | All-electron valence charge (PAW) | PAW only |

---

## &PLOT parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `iflag` | integer | — | **Required.** Dimensionality: 0=1D, 1=1D planar avg, 2=2D, 3=3D, 4=2D polar |
| `output_format` | integer | — | **Required.** Output format (see table below) |
| `fileout` | string | — | Output file name |
| `interpolation` | string | `'fourier'` | `'fourier'` or `'bspline'` |

### output_format values

| Value | Format | iflag | Notes |
|-------|--------|-------|-------|
| 0 | gnuplot 1D | 0,1 | x, f(x) columns |
| 1 | gnuplot contour | 2 | x, y, f(x,y) |
| 2 | gnuplot 2D | 2 | |
| 3 | XCrySDen XSF (3D) | 2,3 | Viewable in VESTA, XCrySDen |
| 4 | obsolete | — | |
| 5 | gnuplot 3D (cartesian) | 3 | Plain text grid |
| 6 | Gaussian Cube | 3 | **Best for validation** — simple format |
| 7 | gnuplot 1D (spherical avg) | 0 | |

For 3D data (iflag=3), also specify grid dimensions:
```
&PLOT
  iflag         = 3
  output_format = 6
  fileout       = 'charge.cube'
  nx = 50, ny = 50, nz = 50
/
```

If `nx, ny, nz` are omitted, pp.x uses the FFT grid dimensions.

---

## Examples

### Total charge density (Gaussian Cube)

```
&INPUTPP
  prefix   = 'si'
  outdir   = './tmp'
  plot_num = 0
  filplot  = 'si_charge.pp'
/
&PLOT
  iflag         = 3
  output_format = 6
  fileout       = 'si_charge.cube'
/
```

### Total potential

```
&INPUTPP
  prefix   = 'si'
  outdir   = './tmp'
  plot_num = 1
  filplot  = 'si_vtot.pp'
/
&PLOT
  iflag         = 3
  output_format = 6
  fileout       = 'si_vtot.cube'
/
```

### Wavefunction |ψ|² at k=1, band=4

```
&INPUTPP
  prefix    = 'si'
  outdir    = './tmp'
  plot_num  = 7
  kpoint(1) = 1
  kband(1)  = 4
  filplot   = 'si_psi2.pp'
/
&PLOT
  iflag         = 3
  output_format = 6
  fileout       = 'si_psi2.cube'
/
```

### XC potential

```
&INPUTPP
  prefix   = 'si'
  outdir   = './tmp'
  plot_num = 21
  filplot  = 'si_vxc.pp'
/
&PLOT
  iflag         = 3
  output_format = 6
  fileout       = 'si_vxc.cube'
/
```

---

## Gaussian Cube file format

The Cube format is the easiest to parse programmatically for validation:

```
Comment line 1
Comment line 2
N_atoms  origin_x  origin_y  origin_z          (Bohr)
N1  voxel1_x  voxel1_y  voxel1_z              (axis 1)
N2  voxel2_x  voxel2_y  voxel2_z              (axis 2)
N3  voxel3_x  voxel3_y  voxel3_z              (axis 3)
Z_atom  charge  x  y  z                        (per atom)
data...                                        (N1 × N2 × N3 values)
```

Values are on the real-space FFT grid. Each line has up to 6 values.
