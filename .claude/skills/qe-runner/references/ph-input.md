# ph.x Input Reference (QE 7.5)

ph.x computes phonon frequencies and eigenvectors using density-functional perturbation theory (DFPT). It requires a prior pw.x SCF calculation with the same `prefix` and `outdir`.

## Input file structure

```text
Title line (free text, ignored)
&INPUTPH
  ... parameters ...
/
qx qy qz          (q-point in Cartesian coords, units 2π/a — unless ldisp=.true.)
```

For `ldisp = .true.` (grid of q-points), the q-point line is omitted.

---

## &INPUTPH parameters

### Required / essential

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `prefix` | string | `'pwscf'` | Must match the pw.x prefix |
| `outdir` | string | `'./'` | Must match the pw.x outdir |
| `tr2_ph` | real | `1.0d-12` | Convergence threshold for DFPT self-consistency. Use 1e-14 for validation. |
| `fildyn` | string | `'matdyn'` | Output file for dynamical matrices |
| `amass(i)` | real | from PP | Atomic mass for species `i` in amu. Override to match your implementation. |

### Dispersion (q-point grid)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ldisp` | logical | `.false.` | Compute phonons on a grid of q-points |
| `nq1, nq2, nq3` | integer | `0` | q-point grid dimensions (with ldisp) |

### Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `epsil` | logical | `.false.` | Compute dielectric constant and Born effective charges (only at q=Gamma) |
| `zeu` | logical | `epsil` | Compute Born effective charges |
| `trans` | logical | `.true.` | Compute phonons. Set `.false.` if you only want dielectric properties. |
| `electron_phonon` | string | `' '` | `'interpolated'`, `'simple'`, `'ahc'` for el-ph coupling |
| `recover` | logical | `.false.` | Restart from interrupted calculation |
| `alpha_mix(i)` | real | `0.7` | Mixing factor for iteration `i` of DFPT |
| `niter_ph` | integer | `100` | Max iterations for DFPT self-consistency |
| `search_sym` | logical | `.true.` | Use symmetry to reduce representations |

---

## Examples

### Gamma-point phonons

```text
Phonons at Gamma
&INPUTPH
  prefix   = 'si'
  outdir   = './tmp'
  tr2_ph   = 1.0d-14
  fildyn   = 'si_gamma.dyn'
  amass(1) = 28.086
  epsil    = .true.
/
0.0 0.0 0.0
```

### Single q-point (zone boundary X)

```text
Phonons at X
&INPUTPH
  prefix   = 'si'
  outdir   = './tmp'
  tr2_ph   = 1.0d-14
  fildyn   = 'si_X.dyn'
  amass(1) = 28.086
/
1.0 0.0 0.0
```

The q-point is in Cartesian coordinates, units of 2π/a.

### Full phonon dispersion (q-grid)

```text
Phonon dispersion on 4x4x4 grid
&INPUTPH
  prefix   = 'si'
  outdir   = './tmp'
  tr2_ph   = 1.0d-14
  fildyn   = 'si.dyn'
  amass(1) = 28.086
  ldisp    = .true.
  nq1 = 4, nq2 = 4, nq3 = 4
/
```

No q-point line when `ldisp = .true.`.

---

## Key output

### Phonon frequencies

```text
     freq (    1) =      -0.136809 [THz] =      -4.563360 [cm-1]
     freq (    2) =      -0.136809 [THz] =      -4.563360 [cm-1]
     freq (    3) =      -0.136809 [THz] =      -4.563360 [cm-1]
     freq (    4) =      15.298080 [THz] =     510.374390 [cm-1]
     freq (    5) =      15.298080 [THz] =     510.374390 [cm-1]
     freq (    6) =      15.298080 [THz] =     510.374390 [cm-1]
```

At Gamma, the first 3 modes are acoustic (near zero). For Si diamond, the optical modes at Gamma should be ~15.5 THz (~520 cm⁻¹) with PBE.

### Dynamical matrix file

The `fildyn` file contains the dynamical matrix in Cartesian coordinates. For dispersion, one file per irreducible q-point is produced (fildyn1, fildyn2, ...).

### Grep patterns

```bash
grep 'freq' output.out                    # phonon frequencies
grep 'Dielectric constant' output.out     # dielectric tensor (if epsil=.true.)
grep 'Effective charges' output.out       # Born charges
```

---

## Workflow for phonon dispersion

1. **SCF** with pw.x — use a dense k-grid (e.g., 8×8×8)
2. **ph.x** with `ldisp = .true.` — compute dynamical matrices on q-grid
3. **q2r.x** — Fourier-transform dynamical matrices to real-space force constants
4. **matdyn.x** — Interpolate to get frequencies along arbitrary q-paths

q2r.x input:

```text
&INPUT
  fildyn = 'si.dyn'
  zasr   = 'simple'
  flfrc  = 'si.fc'
/
```

matdyn.x input:

```text
&INPUT
  asr    = 'simple'
  flfrc  = 'si.fc'
  flfrq  = 'si.freq'
  flvec  = 'si.modes'
  q_in_band_form = .true.
/
5
0.500 0.500 0.500  20   ! L
0.000 0.000 0.000  20   ! Gamma
0.500 0.000 0.500  20   ! X
0.750 0.250 0.750  20   ! W
0.500 0.500 0.500   1   ! L
```
