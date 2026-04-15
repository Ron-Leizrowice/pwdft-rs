# Proposal 25: Ewald Structure Factor Cleanup and Numerical Hardening

## Problem

The Ewald summation (`src/ewald.rs`) has several code quality issues that, while not producing wrong results for typical inputs, reduce clarity and robustness:

### 1. Manual sin/cos structure factor accumulation (lines 62-69)

```rust
let mut s_re = 0.0;
let mut s_im = 0.0;
for (i, pos) in positions.iter().enumerate() {
    let phase = g.dot(pos);
    s_re += charges[i] * phase.cos();
    s_im += charges[i] * phase.sin();
}
let s_sq = s_re * s_re + s_im * s_im;
```

This should use `Complex64::cis()` and `.norm_sqr()` for clarity and to avoid redundant trigonometric decomposition. For systems with many atoms, this also accumulates floating-point error through naive summation.

### 2. Self-interaction threshold (line 95)

```rust
if r_norm < 1e-10 { continue; }
```

This hardcoded threshold skips the self-interaction term (same atom, same cell). For the `i == j, T == 0` case, `r_norm` should be exactly zero (or very close). The threshold works, but the intent is clearer with a logical check:

```rust
if atom_i == atom_j && t_is_zero { continue; }
```

### 3. Missing edge case tests

The test suite verifies NaCl Madelung energy and basic properties but lacks:
- Zero charges (should give zero energy)
- Charged supercell (non-zero background correction)
- Highly anisotropic cells (one very short axis)
- Single-atom cell (self-energy only)

## Implementation

### Step 1: Refactor structure factor with Complex64

```rust
// Before (lines 62-69):
let mut s_re = 0.0;
let mut s_im = 0.0;
for (i, pos) in positions.iter().enumerate() {
    let phase = g.dot(pos);
    s_re += charges[i] * phase.cos();
    s_im += charges[i] * phase.sin();
}
let s_sq = s_re * s_re + s_im * s_im;

// After:
let s: Complex64 = positions.iter().enumerate()
    .map(|(i, pos)| charges[i] * Complex64::cis(g.dot(pos)))
    .sum();
let s_sq = s.norm_sqr();
```

This is clearer, uses a single `cis()` call per atom, and `sum()` on `Complex64` is a standard operation.

### Step 2: Clarify self-interaction skip

```rust
// Before (lines 91-95):
let r = &positions[j] - &positions[i] + &t_vec;
let r_norm = r.norm();
if r_norm < 1e-10 { continue; }

// After:
let r = &positions[j] - &positions[i] + &t_vec;
let r_norm = r.norm();
// Skip self-interaction (same atom in same cell)
if r_norm < 1e-10 { continue; }
```

The threshold approach is actually fine here — the alternative (tracking `i == j` and `t == [0,0,0]`) is more complex and the threshold correctly handles numerical noise. Just add a comment explaining the intent.

### Step 3: Add near-zero G² guard in reciprocal sum

```rust
// Line 71, after skipping G=0:
// Guard against numerically near-zero G² from floating-point noise
if g2 < 1e-12 { continue; }
```

Currently the code skips exact `G == [0,0,0]` via integer comparison (line 55-56), which is correct. But for pathological reciprocal lattices, a G-vector with very small but non-zero `|G|²` could produce a huge `1/G²` term. This is unlikely in practice but costs nothing to guard against.

### Step 4: Add missing test cases

```rust
#[test]
fn test_ewald_zero_charges() {
    let lattice = Lattice::new(/* simple cubic */);
    let positions = vec![Vector3::new(0.0, 0.0, 0.0)];
    let charges = vec![0.0];
    let e = ewald_energy(&lattice, &positions, &charges);
    assert!(e.abs() < 1e-12, "Zero charges should give zero energy: {e}");
}

#[test]
fn test_ewald_single_atom() {
    // Single atom in a box: only self-energy and background
    let a = 5.0;
    let lattice = Lattice::new(
        Vector3::new(a, 0.0, 0.0),
        Vector3::new(0.0, a, 0.0),
        Vector3::new(0.0, 0.0, a),
    );
    let positions = vec![Vector3::new(0.0, 0.0, 0.0)];
    let charges = vec![4.0]; // Si-like
    let e = ewald_energy(&lattice, &positions, &charges);
    // Self-energy should be negative (attractive self-interaction correction)
    assert!(e < 0.0, "Single atom Ewald should be negative: {e}");
}

#[test]
fn test_ewald_anisotropic_cell() {
    // Slab geometry: a=b=10Å, c=2Å
    let lattice = Lattice::new(
        Vector3::new(10.0, 0.0, 0.0),
        Vector3::new(0.0, 10.0, 0.0),
        Vector3::new(0.0, 0.0, 2.0),
    );
    let positions = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(5.0, 5.0, 1.0),
    ];
    let charges = vec![1.0, 1.0];
    let e = ewald_energy(&lattice, &positions, &charges);
    // Should converge without overflow or NaN
    assert!(e.is_finite(), "Anisotropic cell energy should be finite: {e}");
}
```

## Acceptance Criteria

1. **Structure factor uses `Complex64`:** No manual `s_re`/`s_im` accumulation.
2. **Self-interaction skip documented:** Comment explains the `r_norm < 1e-10` threshold.
3. **G² guard present:** Near-zero `G²` values skipped in reciprocal sum.
4. **New tests pass:** Zero charges, single atom, and anisotropic cell all produce correct/finite results.
5. **NaCl test unchanged:** Existing Madelung energy test passes with identical result.
6. **No performance regression:** Refactored structure factor is not slower (benchmark if concerned).
