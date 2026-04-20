# QE Validation Reference Data

Reference results from Quantum ESPRESSO 7.5 used to validate pwdft-rs SCF
results in `tests/qe_validation.rs`. All 8 systems in this directory use the
same PseudoDojo ONCV NC/LDA pseudopotentials as pwdft-rs (`pseudopotentials/nc/lda/`),
exposed here as `./pseudo/` via symlink so QE inputs can reference
`pseudo_dir = './pseudo'` without duplicating content.

## Provenance and trust

Reference values were generated on 2026-04-17 (M3 Max, `qe-7.5/build/bin/pw.x`,
8 MPI ranks) and spot-checked by re-running Si and Al on the same build:
both systems reproduced the archived `total energy` bit-for-bit
(Si: `-17.02299344 Ry`, Al: `-4.72371790 Ry`). The `pseudo/` symlink
points at the canonical `pseudopotentials/nc/lda/` so both pwdft-rs and
QE read identical pseudopotential files.

## Layout

| File | Purpose |
|------|---------|
| `reference_data.toml` | Machine-readable reference numbers (energies, Fermi, Γ-eigenvalues). |
| `<system>_scf.in` | QE pw.x input files. |
| `<system>_scf.out` | Raw pw.x output (committed for auditability). |
| `pseudo/` | Symlink to `../pseudopotentials/nc/lda/` (canonical PP source). |

## Re-running a reference calculation

From the project root:

```bash
cd qe_validation
NP=$(sysctl -n hw.ncpu)  # macOS; use $(nproc) on Linux. Do not hardcode the rank count.
export OMP_NUM_THREADS=1 LC_ALL=C LANG=C
export OMPI_MCA_btl=self,vader OMPI_MCA_pml=ob1
ulimit -s unlimited
perl -e 'alarm shift; exec @ARGV' 600 \
    mpirun -np "$NP" ../qe-7.5/build/bin/pw.x \
    -in si_scf.in > si_scf.out 2>&1
```

The same invocation works for any `<system>_scf.in`. If `gtimeout` is
installed (`brew install coreutils`) it may replace the `perl -e 'alarm'`
wrapper.

## Systems covered (Tier 1 + Tier 2 per proposal QEVL)

| # | System | Structure | Atoms | Physics tested |
|---|--------|-----------|-------|----------------|
| 1 | Si | Diamond/FCC | 2 | Basic SCF, insulator |
| 2 | C | Diamond/FCC | 2 | Higher ecut, wide-gap |
| 3 | Al | FCC | 1 | Simple metal, Kerker mixing |
| 4 | Fe | BCC | 1 | nspin=2 (collapses to NM — see note) |
| 5 | GaAs | Zincblende | 2 | Two species, III-V |
| 6 | Cu | FCC | 1 | Transition metal d-states |
| 7 | NaCl | Rocksalt | 2 | Ionic insulator, charge transfer |
| 8 | MgO | Rocksalt | 2 | Wide-gap ionic |

Tier 3 convergence studies (tests 9-11 in the proposal) are deferred.

## Fe magnetism note

With PseudoDojo NC/LDA pseudopotentials and ecutwfc = 15 Ry, BCC Fe
(a = 2.87 Å) collapses to a non-magnetic state despite
`starting_magnetization = 0.5`. This is a known limitation of the
pseudopotential/cutoff combination (real Fe is ferromagnetic). The test
validates nspin=2 machinery, not the magnetic physics. QE and pwdft-rs
should agree on the non-magnetic energy.
