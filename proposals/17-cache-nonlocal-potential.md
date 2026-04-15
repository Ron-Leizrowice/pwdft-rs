# Proposal 17: Cache NonlocalPotential Across SCF Iterations

## Problem

`NonlocalPotential::new()` is called inside the per-k-point parallel loop on every SCF iteration (`src/scf/mod.rs`, line 277):

```rust
let kpoint_results: Vec<_> = kpoints
    .par_iter()
    .map(|kp| {
        let mut h = build_hamiltonian_with_v_eff(basis, &kp.k, &v_eff_fft, grid.dims);
        let vnl = NonlocalPotential::new(crystal, basis, &kp.k, pseudopotentials);
        vnl.add_to_hamiltonian(&mut h, crystal, basis, &kp.k);
        dense::diagonalize_lowest(&h, params.n_bands)
    })
    .collect();
```

`NonlocalPotential::new()` computes Bessel-transform form factors `β_i(|k+G|)` for every projector and every G-vector (`src/potential/nonlocal.rs`, lines 53-128). This involves numerical integration over the radial grid per (projector, G) pair — O(n_proj × n_pw × n_radial) floating-point operations.

The form factors depend only on `|k+G|`, which is determined by the k-point and basis set — neither changes between SCF iterations. The non-local potential is 15-25% of SCF wall time, and this work is repeated 10-20× unnecessarily.

## Implementation

### Step 1: Precompute before SCF loop

In `src/scf/mod.rs`, before the SCF loop:

```rust
// Precompute non-local projectors (k-dependent, iteration-independent)
let vnl_per_kpoint: Vec<NonlocalPotential> = kpoints
    .par_iter()
    .map(|kp| NonlocalPotential::new(crystal, basis, &kp.k, pseudopotentials))
    .collect();
```

### Step 2: Use cached projectors inside the loop

Replace the per-iteration construction with a reference:

```rust
let kpoint_results: Vec<_> = kpoints
    .par_iter()
    .zip(vnl_per_kpoint.par_iter())
    .map(|(kp, vnl)| {
        let mut h = build_hamiltonian_with_v_eff(basis, &kp.k, &v_eff_fft, grid.dims);
        vnl.add_to_hamiltonian(&mut h, crystal, basis, &kp.k);
        dense::diagonalize_lowest(&h, params.n_bands)
    })
    .collect();
```

### Step 3: Pre-add non-local to Hamiltonian (optional further optimization)

Since `add_to_hamiltonian` itself has an O(n_pw²) loop over structure factors and angular terms, and the structure factors also don't change between iterations, the entire non-local Hamiltonian contribution `V_NL(k)` could be precomputed as a matrix:

```rust
let vnl_matrices: Vec<faer::Mat<Complex64>> = kpoints
    .par_iter()
    .map(|kp| {
        let vnl = NonlocalPotential::new(crystal, basis, &kp.k, pseudopotentials);
        let mut h_nl = faer::Mat::<Complex64>::zeros(n_pw, n_pw);
        vnl.add_to_hamiltonian(&mut h_nl, crystal, basis, &kp.k);
        h_nl
    })
    .collect();
```

Then inside the SCF loop, Hamiltonian construction becomes a simple matrix addition:

```rust
let mut h = build_hamiltonian_with_v_eff(basis, &kp.k, &v_eff_fft, grid.dims);
h += &vnl_matrices[ik]; // O(n_pw²) add instead of O(n_proj × n_pw² × n_radial)
```

This trades O(n_kpoints × n_pw²) memory for eliminating the entire non-local computation from the SCF loop.

## Performance Impact

For Si with n_pw=200, n_bands=8, 10 SCF iterations, 1 k-point:
- Current: `NonlocalPotential::new()` called 10× = 10 × O(n_proj × n_pw × n_radial)
- Step 1: Called 1× (10× speedup on non-local setup)
- Step 3: Non-local Hamiltonian add goes from O(n_proj × n_pw² × n_atoms) to O(n_pw²) per iteration

Net SCF speedup: ~15-25% (the full non-local fraction of wall time, minus one-time setup).

## Acceptance Criteria

1. **Identical results:** SCF converged energy matches current code to within 1e-12 eV.
2. **Measurable speedup:** Wall time for 10+ iteration SCF is reduced by >10%.
3. **Memory bounded:** Additional memory is O(n_kpoints × n_pw²), documented and acceptable for target system sizes.
4. **Parallel correctness:** Cached `NonlocalPotential` is `Sync` and can be shared across rayon threads.
