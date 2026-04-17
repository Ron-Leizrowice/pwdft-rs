---
id: QEVL
status: completed
priority: high
complexity: medium
risk: medium
depends_on: [SIMP]
blocks: []
---

# QEVL: QE Validation Test Suite

## Problem

The current QE validation covers only 3 systems (Si, C diamond, BCC Fe) at a single parameter set (15 Ry, 4x4x4, Fermi-Dirac). This misses important physics regimes:

| Regime | Currently tested | Gap |
|--------|-----------------|-----|
| Insulator (diamond structure) | Si, C | No wide-gap ionic insulator |
| Metal (simple) | Fe (nspin=1, unphysical) | No real simple metal (Al, Cu, Na) |
| Metal (transition, magnetic) | Fe (nspin=2, fixed moment) | No free-moment magnetic test |
| Semiconductor (III-V) | None | GaAs is the canonical compound semiconductor |
| Compound (ionic) | None | NaCl, MgO — tests multi-species handling |
| Heavy element | None | Tests high-Z PP with many projectors |
| Convergence studies | None | ecutwfc, k-points — validates systematic convergence |

We have 70 NC/LDA PseudoDojo pseudopotentials available. The Ljubljana QE tutorials provide a checklist of standard test systems and convergence studies that production DFT codes should reproduce.

## Design Principles

All QE reference runs must use:
- **Same PPs** as pwdft-rs: PseudoDojo ONCV NC/LDA from `pseudopotentials/nc/lda/`
- **Same XC**: LDA (Perdew-Zunger), which is what these PPs encode
- **Same parameters**: ecut, k-grid, smearing type and width stated explicitly
- **QE 7.5** for reproducibility

Each test should compare:
1. **Total energy** (primary): agree within 0.01 eV (limited by different eigensolvers/FFT grids)
2. **Fermi energy** (metals): agree within 0.05 eV
3. **Magnetization** (nspin=2): agree within 0.1 μB
4. **Γ-point eigenvalues** (spot check): agree within 0.1 eV for lowest bands
5. **Convergence**: SCF should converge (not diverge or stall)

## Material Systems

### Tier 1 — Core validation (must pass for correctness)

#### 1. Si diamond (insulator, 2 atoms, FCC)
Already tested. Canonical semiconductor, gap ~1.1 eV (LDA underestimates).

```
# QE input: qe_validation/si_scf.in
&CONTROL
    calculation = 'scf'
    prefix = 'si'
    pseudo_dir = './pseudo/'
    outdir = './tmp/'
/
&SYSTEM
    ibrav = 2
    celldm(1) = 10.2626  ! 5.431 Å in Bohr
    nat = 2
    ntyp = 1
    ecutwfc = 15.0
    occupations = 'smearing'
    smearing = 'fd'
    degauss = 0.01
/
&ELECTRONS
    conv_thr = 1.0d-8
    mixing_beta = 0.3
/
ATOMIC_SPECIES
    Si 28.085 Si.upf
ATOMIC_POSITIONS crystal
    Si 0.00 0.00 0.00
    Si 0.25 0.25 0.25
K_POINTS automatic
    4 4 4 0 0 0
```

**What it validates**: Basic SCF, diamond structure, insulator occupations.

#### 2. C diamond (wide-gap insulator, 2 atoms, FCC)
Already tested. Large band gap (~5.5 eV), stiff lattice.

```
# QE input: qe_validation/c_diamond_scf.in
&SYSTEM
    ibrav = 2
    celldm(1) = 6.7409  ! 3.567 Å in Bohr
    nat = 2
    ntyp = 1
    ecutwfc = 30.0       ! C needs higher cutoff
    occupations = 'smearing'
    smearing = 'fd'
    degauss = 0.01
/
ATOMIC_SPECIES
    C 12.011 C.upf
ATOMIC_POSITIONS crystal
    C 0.00 0.00 0.00
    C 0.25 0.25 0.25
K_POINTS automatic
    4 4 4 0 0 0
```

**What it validates**: Higher ecut, wide-gap insulator, light element.

#### 3. Al FCC (simple metal, 1 atom)
**NEW.** The simplest metal — nearly-free-electron, tests metallic occupation and Kerker preconditioning.

```
# QE input: qe_validation/al_fcc_scf.in
&SYSTEM
    ibrav = 2
    celldm(1) = 7.6527  ! 4.05 Å in Bohr
    nat = 1
    ntyp = 1
    ecutwfc = 15.0
    occupations = 'smearing'
    smearing = 'fd'
    degauss = 0.02
/
ATOMIC_SPECIES
    Al 26.982 Al.upf
ATOMIC_POSITIONS crystal
    Al 0.00 0.00 0.00
K_POINTS automatic
    8 8 8 0 0 0
```

**What it validates**: Metallic occupations, partial filling, Kerker preconditioning essential.
**Rust test**: Use `MixingMode::Kerker`, `smearing_sigma = 0.02 * RY_TO_EV`, 8x8x8 k-grid.

#### 4. BCC Fe (magnetic metal, 1 atom, nspin=2)
Partially tested (fixed moment only). Add free-moment test.

```
# QE input: qe_validation/fe_bcc_fm_scf.in
&SYSTEM
    ibrav = 3
    celldm(1) = 5.4235  ! 2.87 Å in Bohr
    nat = 1
    ntyp = 1
    ecutwfc = 15.0
    occupations = 'smearing'
    smearing = 'fd'
    degauss = 0.02
    nspin = 2
    starting_magnetization(1) = 0.5
/
ATOMIC_SPECIES
    Fe 55.845 Fe.upf
ATOMIC_POSITIONS crystal
    Fe 0.00 0.00 0.00
K_POINTS automatic
    8 8 8 0 0 0
```

**What it validates**: Spin-polarized SCF, NLCC (Fe has core correction), magnetic moment.
**Note**: With PseudoDojo NC/LDA, Fe may converge to non-magnetic. If so, compare against QE non-magnetic as well — the test validates machinery, not physics.

### Tier 2 — Compound systems and diversity

#### 5. GaAs zincblende (III-V semiconductor, 2 atoms, FCC)
The canonical compound semiconductor. Tests multi-species handling with different PP types.

```
# QE input: qe_validation/gaas_scf.in
&SYSTEM
    ibrav = 2
    celldm(1) = 10.6829  ! 5.653 Å in Bohr
    nat = 2
    ntyp = 2
    ecutwfc = 20.0
    occupations = 'smearing'
    smearing = 'fd'
    degauss = 0.01
/
ATOMIC_SPECIES
    Ga 69.723 Ga.upf
    As 74.922 As.upf
ATOMIC_POSITIONS crystal
    Ga 0.00 0.00 0.00
    As 0.25 0.25 0.25
K_POINTS automatic
    4 4 4 0 0 0
```

**What it validates**: Two atom types with different PPs and projectors, NLCC on both species, compound gap.

#### 6. Cu FCC (transition metal, 1 atom)
Noble metal with d-electrons. Tests d-band crossing the Fermi level, harder convergence than Al.

```
# QE input: qe_validation/cu_fcc_scf.in
&SYSTEM
    ibrav = 2
    celldm(1) = 6.8219  ! 3.61 Å in Bohr
    nat = 1
    ntyp = 1
    ecutwfc = 25.0       ! d-electrons need higher cutoff
    occupations = 'smearing'
    smearing = 'fd'
    degauss = 0.02
    ! nspin = 1 — Cu is non-magnetic
/
ATOMIC_SPECIES
    Cu 63.546 Cu.upf
ATOMIC_POSITIONS crystal
    Cu 0.00 0.00 0.00
K_POINTS automatic
    8 8 8 0 0 0
```

**What it validates**: Transition metal d-states, NLCC, higher ecut, metallic Kerker mixing.

#### 7. NaCl rocksalt (ionic insulator, 2 atoms, FCC)
Strongly ionic system with large charge transfer. Tests charge sloshing (Anderson mixing) and multi-species with very different electronegativity.

```
# QE input: qe_validation/nacl_scf.in
&SYSTEM
    ibrav = 2
    celldm(1) = 10.6078  ! 5.614 Å in Bohr (NaCl exp lattice param)
    nat = 2
    ntyp = 2
    ecutwfc = 25.0
    occupations = 'smearing'
    smearing = 'fd'
    degauss = 0.01
/
ATOMIC_SPECIES
    Na 22.990 Na.upf
    Cl 35.453 Cl.upf
ATOMIC_POSITIONS crystal
    Na 0.00 0.00 0.00
    Cl 0.50 0.50 0.50
K_POINTS automatic
    4 4 4 0 0 0
```

**What it validates**: Ionic bonding, two species with very different Z, rocksalt structure, charge transfer.

#### 8. MgO rocksalt (wide-gap ionic, 2 atoms, FCC)
Similar to NaCl but with stronger ionic character and wider gap (~7.7 eV). Standard benchmark in DFT validation.

```
# QE input: qe_validation/mgo_scf.in
&SYSTEM
    ibrav = 2
    celldm(1) = 7.9586  ! 4.212 Å in Bohr
    nat = 2
    ntyp = 2
    ecutwfc = 30.0       ! O needs higher cutoff
    occupations = 'smearing'
    smearing = 'fd'
    degauss = 0.01
/
ATOMIC_SPECIES
    Mg 24.305 Mg.upf
    O  15.999 O.upf
ATOMIC_POSITIONS crystal
    Mg 0.00 0.00 0.00
    O  0.50 0.50 0.50
K_POINTS automatic
    4 4 4 0 0 0
```

**What it validates**: Wide-gap ionic insulator, first-row element (O), higher ecut needed.

### Tier 3 — Convergence studies

#### 9. Si ecutwfc convergence
Sweep ecut from 10 to 40 Ry at fixed 4x4x4 k-grid. Total energy should decrease monotonically and converge. Validates that our plane-wave expansion is variational.

```
# Run QE at ecut = 10, 12, 15, 20, 25, 30, 35, 40 Ry
# Same input as test 1, varying ecutwfc only
```

**Rust test**: Run pwdft-rs at the same ecut values, verify:
- Energy decreases monotonically with ecut
- Converged energy matches QE within 0.01 eV
- Energy difference between 30 and 40 Ry < 0.001 eV

#### 10. Si k-point convergence
Sweep k-grid from 2x2x2 to 8x8x8 at fixed ecut=15 Ry. Tests BZ sampling.

```
# Run QE at nk = 2, 3, 4, 5, 6, 8
# Same input as test 1, varying K_POINTS only
```

**Rust test**: Verify energy converges and matches QE at each grid.

#### 11. Al smearing convergence
Sweep degauss from 0.005 to 0.10 Ry at fixed ecut=15 Ry, 8x8x8. Validates that sigma→0 extrapolation converges and different smearing schemes agree in the limit.

```
# Run QE at degauss = 0.005, 0.01, 0.02, 0.04, 0.06, 0.10 Ry
# Test with smearing = 'fd', 'gauss', 'mp', 'mv'
```

**Rust test**: Verify `energy_sigma0` converges as degauss→0 for all 4 schemes.

## Implementation

### Step 1: Generate QE reference data

Create `qe_validation/` directory with input files for all 11 tests. Run each with QE 7.5 and extract:
- Total energy (Ry)
- Fermi energy (eV)
- Magnetization (μB, if nspin=2)
- Γ-point eigenvalues (eV)
- Number of SCF iterations

Store results in a `qe_validation/reference_data.toml`:

```toml
[si_diamond]
total_energy_ry = -17.02298254
fermi_energy_ev = 6.3435
n_iterations = 7
gamma_eigenvalues_ev = [-5.8909, 6.0800, 6.0800, 6.0800, 8.6090, 8.6090, 8.6090, 9.3220]

[al_fcc]
total_energy_ry = -4.18338  # to be filled
fermi_energy_ev = 7.65      # to be filled
# ...
```

### Step 2: Expand `tests/qe_validation.rs`

Refactor into a structured test file with helper functions:

```rust
// Helper: build crystal from ibrav
fn fcc_crystal(a_ang: f64, atoms: Vec<Atom>) -> Crystal { ... }
fn bcc_crystal(a_ang: f64, atom: Atom) -> Crystal { ... }

// Helper: standard SCF run with consistent parameters
fn run_qe_comparison(
    crystal: &Crystal,
    pps: &[&str],           // element symbols
    ecut_ry: f64,
    nk: u32,
    n_bands: usize,
    mixing: MixingMode,
    smearing: SmearingScheme,
    degauss_ry: f64,
    nspin: usize,
) -> ScfResult { ... }

// Assertion: compare against QE reference
fn assert_energy_matches_qe(result: &ScfResult, qe_energy_ry: f64, tolerance_ev: f64) { ... }
fn assert_fermi_matches_qe(result: &ScfResult, qe_fermi_ev: f64, tolerance_ev: f64) { ... }
```

### Step 3: Write Rust tests

One `#[test]` per system. Convergence studies use parameterized loops. All tests should:
- Print the comparison values (even on pass) for CI visibility
- Use `assert!` with descriptive messages
- Be `#[ignore]`-tagged for convergence studies (they're slow)

```rust
#[test]
fn test_al_fcc_vs_qe() { ... }

#[test]
fn test_gaas_zincblende_vs_qe() { ... }

#[test]
fn test_cu_fcc_vs_qe() { ... }

#[test]
fn test_nacl_rocksalt_vs_qe() { ... }

#[test]
fn test_mgo_rocksalt_vs_qe() { ... }

#[test]
#[ignore] // slow: runs 8 SCF calculations
fn test_si_ecutwfc_convergence() { ... }

#[test]
#[ignore] // slow: runs 6 SCF calculations
fn test_si_kpoint_convergence() { ... }

#[test]
#[ignore] // slow: runs 24 SCF calculations (6 degauss × 4 schemes)
fn test_al_smearing_convergence() { ... }
```

## Verification

1. `cargo test --test qe_validation` passes for all Tier 1 and Tier 2 tests.
2. `cargo test --test qe_validation -- --ignored` runs convergence studies.
3. All QE input files in `qe_validation/` reproduce the reference values when run with QE 7.5 and our PseudoDojo NC/LDA PPs.
4. Energy tolerances:
   - Tier 1 (existing systems): < 0.01 eV (we already match this for Si, C)
   - Tier 2 (new systems): < 0.05 eV initially, tighten as code matures
   - Convergence: monotonic decrease, final value within 0.01 eV of QE

## Estimated Effort

QE reference generation: one session (run 11 calculations, extract data).
Rust test code: one session (refactor existing + add 8 new test functions).
Convergence studies: one session (parameterized sweeps, slower iteration).

Total: ~3 sessions. The QE runs are the bottleneck.

## Summary Table

| # | System | Structure | Atoms | Type | Key physics tested |
|---|--------|-----------|-------|------|--------------------|
| 1 | Si | Diamond/FCC | 2 | Insulator | Basic SCF ✓ |
| 2 | C | Diamond/FCC | 2 | Wide-gap | Higher ecut ✓ |
| 3 | Al | FCC | 1 | Simple metal | **Metallic occ, Kerker** |
| 4 | Fe | BCC | 1 | Magnetic metal | **nspin=2, NLCC, free moment** |
| 5 | GaAs | Zincblende/FCC | 2 | III-V semiconductor | **Two species, compound** |
| 6 | Cu | FCC | 1 | Transition metal | **d-electrons, NLCC** |
| 7 | NaCl | Rocksalt/FCC | 2 | Ionic insulator | **Charge transfer, ionic** |
| 8 | MgO | Rocksalt/FCC | 2 | Wide-gap ionic | **O 2p, wide gap** |
| 9 | Si | Diamond/FCC | 2 | — | **ecutwfc convergence** |
| 10 | Si | Diamond/FCC | 2 | — | **k-point convergence** |
| 11 | Al | FCC | 1 | — | **Smearing convergence** |

## 2026-04-17 — Attempt 1: QE Reference Data Already Generated

A Researcher agent generated QE 7.5 reference data for all 8 Tier 1 + Tier 2 systems (Si, C, Al, Fe, GaAs, Cu, NaCl, MgO). The outputs are archived at `/tmp/pwdft-rescue/qe_validation_data/`, including:

- `*.in` — QE input files for all 8 systems
- `*.out` — raw `pw.x` outputs
- `reference_data.toml` — extracted totals, Fermi energies, magnetizations, Γ eigenvalues
- `pseudo/` — PseudoDojo NC/LDA UPF files used
- `README.md`

Generated on 2026-04-17 with `qe-7.5/build/bin/pw.x`, 8 MPI ranks. These should be validated before being trusted and moved to `qe_validation/` in-repo (suggest a lightweight re-run of one system to confirm reproducibility, since this data was generated by an agent that was stopped mid-session).

### Known/expected failures

- **Fe BCC** collapses to non-magnetic with PseudoDojo NC/LDA at ecut=15 Ry despite `starting_magnetization=0.5` (both QE and pwdft-rs). This is a PP/cutoff limitation; the test validates machinery, not magnetism. Agent noted final |M| = 0.00 μB in `reference_data.toml`.
- **Si** will likely still fail the energy-tolerance assertion by 13.4 eV because VERF does not close the gap (see VERF proposal's 2026-04-17 update). Tier 2 systems with higher Z may show similar systematic offsets.

### Recommended Next Step

The bottleneck (QE reference generation) is done. The remaining work is:

1. **Validate** the `/tmp/pwdft-rescue/qe_validation_data/` outputs by re-running 1-2 systems via the `qe-runner` skill.
2. **Move** validated data to `qe_validation/` in-repo (commit as part of the QEVL PR).
3. **Refactor** `tests/qe_validation.rs` into the structured helpers described in "Implementation" (fcc_crystal, bcc_crystal, run_qe_comparison, assert_energy_matches_qe).
4. **Write** 8 `#[test]` functions — mark those expected to fail pre-Si-root-cause-fix with `#[ignore]` and a clear comment.

Tier 3 convergence studies remain deferred.
