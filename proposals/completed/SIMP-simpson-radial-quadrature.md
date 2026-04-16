---
id: SIMP
status: completed
priority: critical
complexity: medium
risk: medium
depends_on: []
blocks: [VERF, QEVL]
---

# SIMP: Simpson's Rule for Radial Integrals

## Problem

All radial integrals in the codebase use a plain sum `integral += f(r_i) * rab[i]`, which is O(h^2) trapezoidal-like accuracy. Quantum ESPRESSO uses Simpson's rule throughout (`qe-7.5/upflib/simpsn.f90`), which is O(h^4) — roughly 100-10000x more accurate on the same grid for smooth integrands.

This is the most likely root cause of the 13-45 eV energy discrepancy documented in Proposal 30. The per-electron error scales with Z_valence (Si Z=4: 1.66 eV/el, Fe Z=16: 2.84 eV/el), consistent with quadrature error in integrals whose integrands grow with Z.

### Affected integrals

| Location | Integral | Impact |
|----------|----------|--------|
| `src/pseudopotential/mod.rs:104-132` | V_local(G) form factor | **Critical** — largest potential contribution |
| `src/potential/nonlocal.rs:224-238` | Beta projector F(q) | **High** — enters every eigenvalue |
| `src/scf/potentials.rs:83-97` | NLCC core density Bessel transform | Medium — affects XC for NLCC elements |
| `src/scf/initial_density.rs:162-173` | Atomic density Bessel transform | Low — only affects initial guess |

All four use the identical pattern:

```rust
for i in 0..n {
    integral += f(r_i) * rab[i];
}
```

### QE reference

QE's Simpson implementation (`qe-7.5/upflib/simpsn.f90:9-54`) computes:

```
result = (1/3) × [c_1 f_1 rab_1 + c_2 f_2 rab_2 + ... + c_n f_n rab_n]
```

where c_i alternates 1, 4, 2, 4, 2, ..., 4, 1 (standard composite Simpson). For even mesh counts, a correction formula is applied at the boundary.

### Why this causes eigenvalue degeneracy breaking

The sin(Gr)/(Gr) oscillation in the V_local integrand samples different parts of the radial function for different |G| values. Lower-accuracy quadrature produces G-dependent errors in V_local(G). Since the Hamiltonian uses V_eff(G-G'), errors at different |G| shells distort the potential asymmetrically, numerically breaking the crystal symmetry that requires exact eigenvalue degeneracies.

## Research

Verified by reading QE source code line-by-line:

- `qe-7.5/upflib/vloc_mod.f90:147` — `CALL simpson(msh, aux, rab, tab_vloc)` for V_local
- `qe-7.5/upflib/beta_mod.f90:115` — `CALL simpson(kkbeta, aux, rab, vqint)` for beta projectors
- `qe-7.5/upflib/simpsn.f90:9-54` — The Simpson implementation itself (54 lines)

## Implementation

### Step 1: Add Simpson's rule utility

Port QE's `simpsn.f90` to Rust. ~20 lines of logic:

```rust
/// Simpson's rule integration: ∫ f(r) dr ≈ Σ c_i f_i rab_i
///
/// Matches QE's `simpsn.f90`. Weights alternate 1/3, 4/3, 2/3
/// with endpoint corrections. Handles both odd and even mesh sizes.
pub fn simpson_integrate(func: &[f64], rab: &[f64]) -> f64 {
    let n = func.len();
    assert_eq!(n, rab.len());
    if n < 3 {
        // Fallback to trapezoidal for tiny grids
        return func.iter().zip(rab).map(|(&f, &dr)| f * dr).sum();
    }

    let mut sum = 0.0;
    for i in 1..n - 1 {
        // Weight: 4 for even index (1-based), 2 for odd index (1-based)
        // In 0-based: i=1 is even in 1-based → weight 4
        let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * func[i] * rab[i];
    }

    if n % 2 == 1 {
        // Odd mesh: standard Simpson
        (sum + func[0] * rab[0] + func[n - 1] * rab[n - 1]) / 3.0
    } else {
        // Even mesh: boundary correction (matches QE/DFTK formula)
        (sum + func[0] * rab[0]
            - func[n - 3] * rab[n - 3] * 0.25
            + func[n - 2] * rab[n - 2]
            + func[n - 1] * rab[n - 1] * 1.25)
            / 3.0
    }
}
```

Place in `src/pseudopotential/mod.rs` (private) or a new `src/numerics.rs` if reuse across modules is cleaner.

### Step 2: Update V_local form factor

In `src/pseudopotential/mod.rs:98-133`, replace the inline sum with Simpson:

```rust
pub fn v_local_of_g(&self, g_norm: f64, omega: f64) -> f64 {
    const E2: f64 = 14.399_645_351_950_548;
    let n = self.r_grid.len();

    // Build integrand array
    let mut integrand = vec![0.0; n];

    if g_norm < 1e-12 {
        for i in 0..n {
            let r = self.r_grid[i];
            let v_short = self.v_local[i] + self.z_valence * E2 / r.max(1e-20);
            integrand[i] = r * r * v_short;
        }
        let integral = simpson_integrate(&integrand, &self.rab);
        4.0 * std::f64::consts::PI / omega * integral
    } else {
        for i in 0..n {
            let r = self.r_grid[i];
            let gr = g_norm * r;
            let v_short = self.v_local[i] + self.z_valence * E2 / r.max(1e-20);
            let sinc = if gr < 1e-10 { 1.0 - gr * gr / 6.0 } else { gr.sin() / gr };
            integrand[i] = r * r * v_short * sinc;
        }
        let integral = simpson_integrate(&integrand, &self.rab);
        4.0 * std::f64::consts::PI / omega * integral
            - 4.0 * std::f64::consts::PI * self.z_valence * E2 / (omega * g_norm * g_norm)
    }
}
```

### Step 3: Update beta projector form factors

In `src/potential/nonlocal.rs:217-238`:

```rust
fn bessel_transform_projector(
    r_grid: &[f64], rab: &[f64], r_beta: &[f64], l: i32, q: f64,
) -> f64 {
    let n = r_grid.len();
    let mut integrand = vec![0.0; n];

    for i in 0..n {
        let r = r_grid[i];
        let qr = q * r;
        let jl = spherical_bessel_j(l, qr);
        integrand[i] = r_beta[i] * jl * r;
    }

    4.0 * PI * simpson_integrate(&integrand, rab)
}
```

### Step 4: Update NLCC core density

In `src/scf/potentials.rs:83-97`, replace the inline sum with a Simpson call over the integrand `rho_c * j0`.

### Step 5: Update atomic density (SAD)

In `src/scf/initial_density.rs:162-173`, same pattern for the `rho_atom * j0` integral.

## Verification

```bash
cargo test                                                # no regressions
cargo test --release --test qe_validation -- --nocapture   # compare energies vs QE
cargo clippy -q --all-targets                              # no new warnings
```

**Success criteria:**
- Si total energy within 1.0 eV of QE (currently 13.3 eV off). If within 0.1 eV, Proposal 39 is unnecessary.
- Eigenvalue degeneracies at Gamma improve (bands 2-4 splitting < 0.1 eV, currently ~5 eV).
- All existing tests pass unchanged (the integrals converge to the same values on fine grids, just more accurately).

## Estimated Effort

1-2 hours. The Simpson utility is ~20 lines. The 4 integration sites are mechanical refactors (extract integrand array, call utility). No mathematical changes.
