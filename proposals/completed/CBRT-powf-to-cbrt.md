---
id: CBRT
status: active
priority: medium
complexity: trivial
risk: low
depends_on: []
blocks: []
---

# CBRT: Replace powf(1/3) with cbrt()

## Problem

`f64::cbrt()` is a dedicated hardware-optimized function (3-5x faster than `powf(1.0/3.0)`, which routes through `exp(ln(x)/3)`). The XC functional evaluation — called once per grid point per SCF iteration — uses `powf` for every cube root. For a 32^3 grid with 20 iterations, that's ~1.3M unnecessary `powf` calls.

12 occurrences across 3 files:

| File | Line | Expression |
|------|------|------------|
| `src/potential/xc.rs` | 90 | `(3.0 * rho_bohr / PI).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 111 | `(3.0 / (4.0 * PI * rho_bohr)).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 237 | `(6.0 * rho_up_bohr / PI).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 242 | `(6.0 * rho_down_bohr / PI).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 274 | `(3.0 / (4.0 * PI * rho_bohr)).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 283 | `2.0_f64.powf(4.0 / 3.0)` |
| `src/potential/xc.rs` | 286 | `(1.0 + zeta).max(0.0).powf(4.0 / 3.0)` |
| `src/potential/xc.rs` | 287 | `(1.0 - zeta).max(0.0).powf(4.0 / 3.0)` |
| `src/potential/xc.rs` | 292 | `(1.0 + zeta).max(0.0).powf(1.0 / 3.0)` |
| `src/potential/xc.rs` | 293 | `(1.0 - zeta).max(0.0).powf(1.0 / 3.0)` |
| `src/ewald.rs` | 52 | `(n_atoms as f64 * PI / omega).powf(1.0 / 3.0)` |
| `src/scf/mixing.rs` | 195 | `(3.0 * PI² * rho_bohr).powf(1.0 / 3.0)` |

## Implementation

Mechanical replacements:

- `x.powf(1.0 / 3.0)` → `x.cbrt()`
- `x.powf(4.0 / 3.0)` → `{ let c = x.cbrt(); c * x }` (or inline where readable)
- `2.0_f64.powf(4.0 / 3.0)` → `2.0_f64.cbrt() * 2.0`

For the `(1+ζ)^{4/3}` and `(1+ζ)^{1/3}` pairs in `pz_correlation_spin`, compute the cbrt once and reuse:

```rust
let cbrt_op = (1.0 + zeta).max(0.0).cbrt();
let cbrt_om = (1.0 - zeta).max(0.0).cbrt();
let op_zeta = cbrt_op * cbrt_op * cbrt_op * cbrt_op;  // or: cbrt_op * (1.0 + zeta).max(0.0)
let om_zeta = cbrt_om * cbrt_om * cbrt_om * cbrt_om;
// ...
let dfz = (4.0 / 3.0) * (cbrt_op - cbrt_om) / f_denom;
```

## Verification

```bash
cargo test                         # correctness unchanged
cargo test --test qe_validation    # physics unchanged
cargo bench --bench scf_benchmarks # measure speedup
```

The existing `test_lsda_unpolarized_limit` test validates that the unpolarized and spin-polarized paths agree to 1e-10, so any numerical drift from the change will be caught.

## Estimated Effort

Under 30 minutes. Purely mechanical.
