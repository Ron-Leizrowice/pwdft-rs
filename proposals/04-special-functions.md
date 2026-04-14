# Proposal: Replace hand-rolled special functions with validated libraries

## Motivation

The codebase contains three hand-implemented mathematical functions that are correct for the current use cases but limited in scope:

1. **Spherical Bessel functions** (`src/potential/nonlocal.rs`, lines 234-248): `spherical_bessel_j(l, x)` implemented for l = 0, 1, 2, 3 only. Panics for l > 3.

2. **Legendre polynomials** (`src/potential/nonlocal.rs`, lines 252-260): `legendre_p(l, x)` implemented for l = 0, 1, 2, 3 only. Panics for l > 3.

3. **Complementary error function** (`src/ewald.rs`, lines 120-129): `erfc(x)` using the Abramowitz & Stegun polynomial approximation (equation 7.1.26), accurate to ~1e-6.

These limits block:
- **d-electron pseudopotentials** (l=2 projectors are supported, but l=3 for f-electrons is the hard cutoff)
- **PAW/ultrasoft PPs** which may have higher angular momentum channels
- **High-accuracy Ewald sums** where 1e-6 erfc precision limits energy accuracy

## Options

### Option A: `puruspe` (pure Rust special functions)

```toml
puruspe = ">=0.3"
```

Provides `erfc`, `erf`, `gamma`, `ln_gamma`, `beta`, `regularized_incomplete_beta`. Pure Rust, no system dependencies. Covers the Ewald needs with full f64 precision.

Does **not** provide spherical Bessel or Legendre functions.

### Option B: `sphrs` (spherical harmonics)

```toml
sphrs = ">=0.2"
```

Provides spherical harmonics Y_lm for arbitrary l and m. Useful if you need full angular decomposition beyond the P_l(cos theta) sum-over-m shortcut currently used.

Does not directly provide spherical Bessel functions or Legendre polynomials as standalone functions.

### Option C: `special` crate

```toml
special = ">=0.10"
```

Provides Bessel functions (regular and modified), error functions, gamma functions, beta functions, hypergeometric functions. Covers `erfc` and Bessel functions but not spherical Bessel directly (would need `j_l(x) = sqrt(pi/2x) * J_{l+1/2}(x)` conversion).

### Option D: Keep hand-rolled, extend to arbitrary l (recommended)

The current implementations are straightforward and performant. The main limitation is the l <= 3 hardcoded match. A recurrence relation handles arbitrary l:

**Spherical Bessel:**
```rust
fn spherical_bessel_j(l: i32, x: f64) -> f64 {
    if x.abs() < 1e-10 {
        return if l == 0 { 1.0 } else { 0.0 };
    }
    if l == 0 { return x.sin() / x; }
    if l == 1 { return x.sin() / (x * x) - x.cos() / x; }
    // Upward recurrence: j_{l+1}(x) = (2l+1)/x * j_l(x) - j_{l-1}(x)
    let mut jlm1 = x.sin() / x;
    let mut jl = x.sin() / (x * x) - x.cos() / x;
    for n in 1..l {
        let jlp1 = (2 * n + 1) as f64 / x * jl - jlm1;
        jlm1 = jl;
        jl = jlp1;
    }
    jl
}
```

**Legendre:**
```rust
fn legendre_p(l: i32, x: f64) -> f64 {
    if l == 0 { return 1.0; }
    if l == 1 { return x; }
    // Bonnet's recurrence: (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)
    let mut plm1 = 1.0;
    let mut pl = x;
    for n in 1..l {
        let plp1 = ((2 * n + 1) as f64 * x * pl - n as f64 * plm1) / (n + 1) as f64;
        plm1 = pl;
        pl = plp1;
    }
    pl
}
```

**erfc:** Replace the Abramowitz & Stegun 5-term polynomial with a higher-precision rational approximation, or use `puruspe::erfc` for validated full-precision.

## Recommended Approach

**Combine Option D (recurrence) with Option A (`puruspe` for erfc):**

1. Extend `spherical_bessel_j` and `legendre_p` to arbitrary l via recurrence (lines 234-260 in `nonlocal.rs`). These are hot-path functions called O(n_pw^2) times, and the recurrence is faster than a library call with function pointer overhead.

2. Replace `erfc` in `ewald.rs` with `puruspe::erfc` for full f64 precision. Ewald is computed once per SCF, so performance is irrelevant — correctness matters.

## Scope of Changes

### File: `src/potential/nonlocal.rs`

Lines 234-260: Replace the match-based implementations with recurrence versions shown above. Remove the `panic!` branches. Add tests for l = 4, 5, 6 against known values.

### File: `src/ewald.rs`

Lines 120-129: Replace `fn erfc(x: f64) -> f64` with:
```rust
fn erfc(x: f64) -> f64 {
    puruspe::erfc(x)
}
```

Or inline the call at the two call sites (lines 98 and 108).

### Tests

- Add `test_spherical_bessel_j_l4` through `l6` using known values
- Add `test_legendre_p_l4` through `l6` using known polynomial values
- Add `test_erfc_high_precision` comparing to known values at more decimal places
- Validate that SCF total energy is unchanged after the switch

## Risks

- Upward recurrence for spherical Bessel functions is numerically unstable for large l and small x. For DFT pseudopotentials, l <= 3 is typical and l <= 6 covers all practical cases. The instability threshold (l >> x) is not reached.
- Adding `puruspe` for a single function may feel like dependency overhead. Alternatively, use a higher-precision polynomial approximation for erfc (7-term rational, ~1e-15 accuracy) and keep it self-contained.

## Expected Impact

- **Capability:** Enables f-electron pseudopotentials (l=3 projectors) and future PAW support.
- **Accuracy:** Full-precision erfc improves Ewald energy by ~1e-6 eV (negligible for most purposes, but removes a known approximation).
- **Code quality:** Removes hardcoded match arms and `panic!` branches.
