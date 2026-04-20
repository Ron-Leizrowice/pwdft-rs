# Proposal 25: Fe BCC Energy Discrepancy — Missing NLCC

**Status: ROOT CAUSE IDENTIFIED**

## Problem

BCC Fe total energy is ~210 eV off QE. Diagnostic tests show every eigenvalue is shifted by a constant ~15.2 eV.

## Root Cause

The Fe pseudopotential (`Fe.pz-n-nc.UPF`) has **nonlinear core correction (NLCC)** enabled:

```text
nlcc=.true.
core_correction="T"
<PP_NLCC> ... 1191 radial values ... </PP_NLCC>
```

Our code completely ignores the `PP_NLCC` data. NLCC adds the frozen core charge density to the valence density when evaluating exchange-correlation. Without it, the XC potential is evaluated on ρ_valence alone instead of ρ_valence + ρ_core, producing a systematic shift in all eigenvalues.

Si has `core_correction="F"` — no NLCC — which is why it works fine.

## Evidence

Diagnostic test (`tests/fe_debug.rs`) shows:

- Basis size matches QE (79 PWs) ✓
- Kinetic eigenvalues correct ✓
- V_NL Hermitian, D_ij correct ✓
- Every eigenvalue shifted by ~15.2 eV (constant offset = V_local or XC error)
- V_local(G=0) = 21.16 eV (matches manual calculation from PP data)

The constant offset rules out basis set, kinetic, or V_NL bugs. It's the XC potential that's wrong because it doesn't include the core charge.

## Fix

### Step 1: Parse PP_NLCC in UPF parser

Add `core_charge: Vec<f64>` to `PseudopotentialData`. Parse `PP_NLCC` block (same radial grid as V_local). Convert from e/Bohr³ to e/ų.

### Step 2: Add core density to XC evaluation

In the SCF loop, before computing XC:

```rust
// Add core density for NLCC
let rho_for_xc: Vec<f64> = if has_nlcc {
    rho_r.iter().zip(rho_core_r.iter()).map(|(v, c)| v + c).collect()
} else {
    rho_r.to_vec()
};
let (exc_r, vxc_r) = xc::lda_xc_grid(&rho_for_xc);
```

The core density needs to be computed on the FFT grid via Bessel transform, similar to the SAD initial density.

### Step 3: Core energy correction

The XC energy has a correction term:

```text
E_xc[ρ_val + ρ_core] - E_xc[ρ_val]
```

This is handled automatically by passing the augmented density to `lda_xc_grid`.

## Impact

Fixes any pseudopotential with NLCC — common for transition metals (Fe, Ni, Co, Cu), some main-group elements, and most production PPs from PseudoDojo/SSSP.
