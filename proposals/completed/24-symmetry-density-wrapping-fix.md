# Proposal 24: Fix Density Symmetrization Grid Wrapping

## Problem

`src/symmetry/density.rs` has an inconsistency between how fractional coordinates are mapped to grid indices in two different code paths:

### `symmetrize_density()` uses `round()` (line 65)

```rust
let jx = ((fp[0] * nx as f64).round() as i64 % nx as i64 + nx as i64) as usize % nx;
```

### `check_grid_compatibility()` / `frac_to_grid_index()` uses `floor()` (line 95)

```rust
f -= (f + 0.5).floor();  // wrap to [-0.5, 0.5)
```

`round()` and `floor()` disagree at half-integer boundaries. For a 10×10×10 grid, fractional coordinate `f = 0.05` maps to grid index:

- `round(0.05 × 10) = round(0.5) = 1` (rounds up)
- `floor(0.05 × 10) = floor(0.5) = 0` (rounds down)

This means the symmetrized density can pick up values from the wrong grid point at boundary cases. The error is small (one grid point) but breaks the exact symmetry that symmetrization is supposed to enforce.

### Additional issue: no `rho_orig` copy

`symmetrize_density()` accumulates into `rho` while reading from `rho_orig` (a clone made at line 33). This is correct. But the identity operation (always present) copies `rho[dst] += rho_orig[src]` where `src == dst`, meaning the identity contribution is correct. No bug here, but the averaging by `n_ops` at line 81-83 means the original density contributes with weight `1/n_ops` from the identity, plus `(n_ops-1)/n_ops` from the other operations. This is correct symmetrization.

## References

- QE source: `symmetrize_rho.f90` — uses `nint()` (nearest integer, equivalent to `round()`)
- VASP: uses nearest-grid-point mapping
- Convention: most codes use `round()` / `nint()` for grid mapping, but the key requirement is internal consistency

## Implementation

### Step 1: Use consistent `round()` everywhere

In `symmetrize_density()`, the existing `round()` is the standard convention. Update `frac_to_grid_index()` and `check_grid_compatibility()` to use the same convention:

```rust
/// Map fractional coordinate to grid index using nearest-grid-point (round).
fn frac_to_grid_idx(frac: f64, n: usize) -> usize {
    let scaled = frac * n as f64;
    ((scaled.round() as i64 % n as i64) + n as i64) as usize % n
}
```

Replace all three inline wrapping expressions in `symmetrize_density()` (lines 65-70) with calls to this helper:

```rust
let jx = frac_to_grid_idx(fp[0], nx);
let jy = frac_to_grid_idx(fp[1], ny);
let jz = frac_to_grid_idx(fp[2], nz);
```

### Step 2: Extract grid index helper

The wrapping logic `((x.round() as i64 % n as i64) + n as i64) as usize % n` appears 3 times in `symmetrize_density` and is error-prone. Extract to a well-tested helper function:

```rust
/// Map a fractional coordinate to the nearest grid index in [0, n).
///
/// Uses nearest-integer mapping: index = round(frac * n) mod n.
/// Handles negative fractional coordinates correctly via double-modulo.
fn frac_to_grid_idx(frac: f64, n: usize) -> usize {
    let ni = n as i64;
    let idx = (frac * n as f64).round() as i64;
    ((idx % ni) + ni) as usize % n
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_frac_to_grid_idx() {
        assert_eq!(frac_to_grid_idx(0.0, 10), 0);
        assert_eq!(frac_to_grid_idx(0.5, 10), 5);
        assert_eq!(frac_to_grid_idx(0.99, 10), 10 % 10); // rounds to 10 → wraps to 0
        assert_eq!(frac_to_grid_idx(-0.1, 10), 9); // -1 + 10 = 9
        assert_eq!(frac_to_grid_idx(1.0, 10), 0); // wraps
        assert_eq!(frac_to_grid_idx(0.05, 10), 1); // round(0.5) = 1
    }
}
```

### Step 3: Update `check_grid_compatibility()`

The compatibility check should verify that symmetry operations map grid points to grid points. With the `round()` convention, this is guaranteed when `R_{ij} × n_j` is divisible by `n_i` (so rotated fractional coordinates land exactly on grid points, not between them):

```rust
pub fn check_grid_compatibility(dims: [usize; 3], ops: &[SpaceGroupOp]) -> bool {
    for op in ops {
        for i in 0..3 {
            for j in 0..3 {
                let product = op.rotation[i][j] as i64 * dims[j] as i64;
                if product % dims[i] as i64 != 0 {
                    return false;
                }
            }
        }
    }
    true
}
```

This is already the existing logic. The fix is ensuring that `symmetrize_density` and `check_grid_compatibility` agree on which grid points are "exact" by using the same rounding convention.

## Acceptance Criteria

1. **Internal consistency:** `frac_to_grid_idx()` is used in both `symmetrize_density()` and `check_grid_compatibility()`.
2. **Symmetry preserved:** After symmetrization, the density has the full space group symmetry: `ρ(Rr + τ) = ρ(r)` for all operations, verified by re-applying symmetrization and checking `max|ρ_sym - ρ| < 1e-12`.
3. **Electron count conserved:** `Σ ρ(r) × dvol` is unchanged by symmetrization to within 1e-12.
4. **Edge cases tested:** Fractional coordinates at 0.0, 0.5, 1.0, -0.5, and near half-integer boundaries produce correct grid indices.
5. **Helper function tested:** `frac_to_grid_idx` has unit tests covering positive, negative, and boundary cases.
