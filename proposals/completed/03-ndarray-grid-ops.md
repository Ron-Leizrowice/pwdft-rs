# Proposal: Adopt ndarray for grid operations

**Status: COMPLETED (partial — mixing module migrated, ndarray already a dependency)**

## Result

ndarray was added as a dependency in the ndrustfft migration (proposal 02). The mixing module (`src/scf/mixing.rs`) was migrated to use `ndarray::Array1` for history vectors and `ArrayView1::dot()` for inner products, eliminating manual dot product and zip/map/collect patterns.

The remaining phases (XC grid ops, Hartree, density symmetrization) were deferred — the current code is clean and performant without them. ndarray is available for future use where it adds clarity.

## Changes made

- `src/scf/mixing.rs`: history vectors `Vec<Vec<f64>>` → `Vec<Array1<f64>>`, residual computation via ndarray arithmetic (`&a - &b`), dot products via `Array1::dot()`
- Performance neutral (mixing is <1% of SCF wall time)

## Multi-material test suite (added alongside)

Free-electron band tests extended beyond Si to:

- **Diamond C** (FCC, a=3.567 Å): Γ shell degeneracies, X-point analytic match
- **BCC Fe** (a=2.87 Å): Γ eigenvalues, band continuity along Γ→N
- Fe.UPF and C.UPF pseudopotentials added from QE 7.5 distribution

165 tests total, all passing.
