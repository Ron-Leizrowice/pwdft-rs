# Units and Constants

## Convention

This code works in **eV and Angstroms** throughout (not Hartree/Bohr or Rydberg).

| Quantity | Unit |
|----------|------|
| Energy | eV |
| Length | Angstrom |
| Density | e/Angstrom^3 |
| Potential | eV |
| Wavevector | 1/Angstrom |
| hbar^2/2m | eV*Angstrom^2 (= 3.81 eV*Angstrom^2) |
| e^2 | 14.4 eV*Angstrom (Coulomb constant) |

## Physical Constants (CODATA 2018)

| Constant | Symbol | Value | Code |
|----------|--------|-------|------|
| Hartree to eV | HA_TO_EV | 27.211386245988 | `consts.rs` |
| Rydberg to eV | RY_TO_EV | 13.605693122994 | `consts.rs` |
| Bohr to Angstrom | BOHR_TO_ANG | 0.529177210903 | `consts.rs` |
| Bohr^3 to Angstrom^3 | BOHR3_TO_ANG3 | BOHR_TO_ANG^3 = 0.14818... | `consts.rs` |
| Coulomb constant | E2_COULOMB | 14.399645351950548 | `consts.rs` |
| hbar^2/2m | HBAR2_OVER_2M | ~3.81 | `consts.rs` |
| G=0 threshold | G2_ZERO_THRESHOLD | 1e-12 | `consts.rs` |

All values verified against CODATA 2018 recommended values.

## UPF Unit Conversion Chain

UPF files (Quantum ESPRESSO format) use Rydberg atomic units internally.
Conversion happens at parse time in `src/pseudopotential/upf/convert.rs`:

| Quantity | UPF unit | Conversion | Internal unit |
|----------|----------|------------|---------------|
| Radial grid r | Bohr | x BOHR_TO_ANG | Angstrom |
| Grid weights rab | Bohr | x BOHR_TO_ANG | Angstrom |
| V_local(r) | Ry | x RY_TO_EV | eV |
| Beta projectors chi(r) | Bohr^{-1/2} | / sqrt(BOHR_TO_ANG) | Angstrom^{-1/2} |
| D_ij matrix | Ry | x RY_TO_EV | eV |
| Atomic density 4pi r^2 rho | e/Bohr | / BOHR_TO_ANG | e/Angstrom |
| Core charge 4pi r^2 rho_core | e/Bohr | / BOHR_TO_ANG | e/Angstrom |

### Beta projector conversion detail

UPF stores `chi(r) = r * beta(r)` where beta is in Bohr^{-3/2}, so chi
is in Bohr^{-1/2}. Energy dimension enters only through D_ij.

`Bohr^{-1/2} -> Angstrom^{-1/2}`: divide by `sqrt(BOHR_TO_ANG)`.

## XC Unit Conversion

The XC functional is evaluated in Hartree atomic units internally:

1. Input: rho in e/Angstrom^3
2. Convert: `rho_Bohr = rho * BOHR3_TO_ANG3` (multiply by ~0.148)
3. Compute: epsilon_xc, V_xc in Hartree
4. Convert: multiply by HA_TO_EV to get eV

**Code:** `src/potential/xc.rs:80-96` (exchange), `xc.rs:105-147` (correlation)

## Known Issue: Duplicate E2

The value `14.399645351950548` is defined in 4 places. See Proposal 33 for
consolidation plan.

## Audit Status

| Item | Status |
|------|--------|
| All CODATA 2018 values | CORRECT |
| HBAR2_OVER_2M derivation | CORRECT (within 0.5%) |
| UPF r, rab conversion | CORRECT |
| UPF V_local conversion | CORRECT |
| UPF beta conversion | CORRECT |
| UPF D_ij conversion | CORRECT |
| UPF rho_atom conversion | CORRECT |
| UPF core_charge conversion | CORRECT |
| XC unit chain (Angstrom -> Bohr -> Ha -> eV) | CORRECT |
