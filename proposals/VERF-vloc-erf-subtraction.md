---
id: VERF
status: active
priority: critical
complexity: small
risk: medium
depends_on: [SIMP]
blocks: [QEVL]
---

# VERF: V_local erf Coulomb Subtraction

> **Prerequisite:** Implement Proposal 38 (Simpson's rule) first. This proposal is only needed if Simpson alone doesn't bring Si within 0.1 eV of QE.

## Problem

Our V_local(G) form factor computation subtracts the full Coulomb tail `Ze^2/r` from V_local(r), creating a divergent integrand at r=0. QE uses a smoother decomposition with `erf(r)/r`, producing a bounded integrand everywhere.

### Current approach (`src/pseudopotential/mod.rs:115-132`)

For G != 0:
```
integrand = r^2 × [V_local(r) + Z·e2/r] × sin(Gr)/(Gr)
analytical_correction = -4π Z e2 / (Ω G^2)
```

The term `V_local(r) + Z·e2/r` diverges as `Z·e2/r` at r→0 because pseudopotentials are smooth at the origin (V_local(0) is finite, not -Z/r). The r^2 factor makes the integral convergent, but the integrand can be very large at the first few grid points (billions of eV for Fe at r ~ 10^-6 Å).

### QE approach (`qe-7.5/upflib/vloc_mod.f90:136-148`)

For G != 0 (q > 0):
```fortran
aux(ir) = (r*vloc(ir) + Zp*e2*erf(r)) * sin(q*r)/q
```

Since `erf(r)/r → 2/√π` as r→0, the integrand is bounded:
```
integrand = [r·V_local(r) + Z·e2·erf(r)] × sin(qr)/q
```

The analytical Coulomb correction becomes:
```fortran
vloc(igl) = vloc(igl) - fpi*Zp*e2*exp(-gl*tpiba2*0.25d0)/gl * (1/omega)
```

i.e., `-4π Z e2 exp(-G^2/4) / (Ω G^2)`, which is the FT of `-Z·e2·erfc(r)/r`.

For G=0, both approaches use the same formula (full Coulomb subtraction). See QE `vloc_mod.f90:158-163`.

### Mathematical equivalence

Both decompositions give the same V_local(G):

```
V_local(r) = [V_local(r) + Z·e2/r] - Z·e2/r        (ours)
V_local(r) = [V_local(r) + Z·e2·erf/r] - Z·e2·erf/r (QE)
```

The FT of the long-range parts differ:
- Ours: FT[-Ze2/r] = -4πZe2/G^2
- QE: FT[-Ze2·erf(r)/r] = -4πZe2[1-exp(-G^2/4)]/G^2

But the short-range integrals compensate, yielding identical totals. The difference is purely in numerical quality.

## Implementation

### Step 1: Modify G != 0 branch in `v_local_of_g`

In `src/pseudopotential/mod.rs`, replace lines 115-132:

```rust
} else {
    // G ≠ 0: use erf subtraction (QE convention) for smoother integrand.
    // Decomposition: V(r) = [V(r) + Ze2·erf(r)/r] - Ze2·erf(r)/r
    // The bracketed term is smooth at r=0 (erf(r)/r → 2/√π).
    // FT of erf(r)/r → 4πZe2[1-exp(-G²/4)]/G², so we subtract
    // the erfc part: -4πZe2·exp(-G²/4)/(ΩG²)
    let mut integrand = vec![0.0; n];
    for i in 0..n {
        let r = self.r_grid[i];
        let gr = g_norm * r;
        // erf(r)/r is finite at r=0; use series for small r
        let erf_over_r = if r < 1e-20 {
            2.0 / std::f64::consts::PI.sqrt()  // limit of erf(r)/r as r→0
        } else {
            puruspe::erf(r) / r
        };
        let v_short = self.v_local[i] + self.z_valence * E2 * erf_over_r;
        let sinc = if gr < 1e-10 {
            1.0 - gr * gr / 6.0
        } else {
            gr.sin() / gr
        };
        // Integrand: r × v_short × sin(Gr)/G = r^2 × v_short × sinc(Gr)
        integrand[i] = r * r * v_short * sinc;
    }
    let integral = simpson_integrate(&integrand, &self.rab);
    let g2 = g_norm * g_norm;
    4.0 * PI / omega * integral
        - 4.0 * PI * self.z_valence * E2 * (-g2 / 4.0).exp() / (omega * g2)
}
```

### Step 2: Keep G=0 unchanged

The G=0 branch already uses the full Coulomb subtraction, matching QE's convention at `vloc_mod.f90:158-163`:
```fortran
aux(ir) = r * (r*vloc(ir) + Zp*e2)
```

No change needed.

### Step 3: Check erf availability

`puruspe` (already a dependency) provides `puruspe::erf()`. If not available, use `libm::erf()` or a simple polynomial approximation.

## Verification

```bash
cargo test                                                # no regressions
cargo test --release --test qe_validation -- --nocapture   # compare energies vs QE
cargo clippy -q --all-targets
```

**Key test:** Compare V_local(G) for the first 20 G-shells for Si against QE reference values. The two codes should now agree to ~1e-6 eV.

**Success criteria:**
- Si total energy within 0.1 eV of QE (-231.61 eV)
- Fe total energy within 0.1 eV of QE (-3059.46 eV)
- All eigenvalue degeneracies at Gamma exact to 1e-4 eV

## Estimated Effort

1-2 hours. Only modifies the G != 0 branch of a single function. The erf function is already available via `puruspe`.
