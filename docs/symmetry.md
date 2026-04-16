# Symmetry

## Space Group Detection

Finds all crystallographic symmetry operations `{R|τ}` where R is a
rotation/reflection and τ is a fractional translation.

### Algorithm

1. Compute metric tensor `M = L^T L` from lattice vectors
2. Enumerate all 3×3 integer matrices R satisfying:
   - `det(R) = ±1` (unimodular — preserves volume)
   - `R^T M R = M` (metric-preserving — preserves distances)
3. For each rotation R, find compatible translations τ by trying all atoms
   as mapping targets, then verifying ALL atoms map correctly

### Group verification

`verify_group_closure()` confirms the found operations satisfy:
- Closure: composition of any two operations is in the set
- Inverse: every operation has an inverse in the set
- Tolerance: 1e-5 for fractional coordinate matching

**Code:** `src/symmetry/detect.rs`

### Test coverage

- FCC Si: 48 operations (Fd-3m, O_h point group) — 24 proper, 24 improper
- BCC Fe: 48 operations (Im-3m)
- Triclinic P1: 1 operation (identity only)
- Triclinic P-1: 2 operations (identity + inversion)

## K-Point Reduction (IBZ)

Reduces the full Monkhorst-Pack grid to the irreducible Brillouin zone.

### Monkhorst-Pack formula

```
k_j = (2i_j - N_j + 1) / (2N_j)    for i_j ∈ [0, N_j)
```

- Odd N: includes Gamma point
- Even N: does not include Gamma
- Weight per point: `w = 1 / (N₁ × N₂ × N₃)`

**Code:** `src/kpoints.rs:30-32`

### Reduction algorithm

1. Generate full grid in fractional reciprocal coordinates
2. For each k-point, apply all symmetry operations using `(R^{-1})^T` transformation
3. Apply time-reversal symmetry: k → -k
4. Mark visited points, count orbit size
5. Weight = |orbit| / N_total (preserves normalization: Σ w = 1.0)

**Code:** `src/symmetry/kpoints.rs`

### Known limitation: BZ boundary handling

At Brillouin zone boundaries (k near ±0.5), rounding tolerance (1e-6) may
miss some equivalences. The implementation produces a valid superset of the
minimal IBZ (integration is still correct, weights sum to 1.0, but may use
more k-points than necessary).

Example: Si 4×4×4 gives 10 irreducible points vs QE's 8. Both produce
correct energies.

## Density Symmetrization

Enforces crystal symmetry on the real-space density:

```
ρ_sym(r) = (1/N_ops) Σ_S ρ(S^{-1} r)
```

### Algorithm

For each symmetry operation S = {R|τ}:
1. Compute inverse: `S^{-1} = {R^{-1} | -R^{-1}τ}`
2. Map each FFT grid point: `f' = R^{-1}·f + τ_inv` (mod 1)
3. Accumulate: `ρ_sym(r) += ρ(mapped_point)`
4. Divide by N_ops

### Grid compatibility

Before symmetrizing, verifies `R_{ij} × n_j ≡ 0 (mod n_i)` for all
operations, ensuring rotation maps grid points to grid points exactly.

**Code:** `src/symmetry/density.rs`

### Properties

- Preserves integral: `Σ ρ_sym = Σ ρ`
- Idempotent: `S(S(ρ)) = S(ρ)`
- Uniform density unchanged

## Audit Status

| Item | Status |
|------|--------|
| Metric preservation criterion | CORRECT |
| Group closure verification | CORRECT |
| Monkhorst-Pack formula | CORRECT |
| K-point (R^{-1})^T transformation | CORRECT |
| Time-reversal symmetry | CORRECT |
| Weight conservation (Σw = 1) | CORRECT |
| Density symmetrization formula | CORRECT |
| Grid compatibility check | CORRECT |
| FFT grid wrapping (periodic BC) | CORRECT |
| **BZ boundary k-point folding** | **MINOR ISSUE** — valid but non-minimal |
