# Proposal 18: Fix PSP8 D_ij Matrix (Remaining Work)

**Status:** Partially complete. UPF rho_atom units were fixed (commit d8f9687). PSP8 D_ij remains all zeros.

## Problem

`src/pseudopotential/psp8.rs`, line 144:

```rust
// TODO: read ekb values properly for off-diagonal terms
let dij = vec![0.0; n_projectors * n_projectors];
```

The D_ij matrix is always zero for PSP8 pseudopotentials, meaning the non-local energy contribution `V_NL = Σ |β_i⟩ D_{ij} ⟨β_j|` is identically zero. Any calculation using PSP8 format pseudopotentials has no non-local contribution and gives wrong results.

PSP8 encodes KB energies (`ekb`) in each projector block header. For norm-conserving PPs, D_ij is diagonal: `D_ii = ekb_i`.

## Fix

Parse `ekb` from each projector block header and populate the diagonal of D_ij. See the original proposal in `completed/18-psp8-dij-fix.md` for detailed implementation steps and unit conversion notes.

## Acceptance Criteria

1. PSP8 D_ij has non-zero diagonal entries matching `ekb` values.
2. PSP8-based Si SCF matches UPF-based Si SCF within 0.01 eV.
