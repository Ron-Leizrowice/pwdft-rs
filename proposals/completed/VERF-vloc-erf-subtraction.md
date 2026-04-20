---
id: VERF
status: completed
priority: critical
complexity: small
risk: medium
depends_on: [SIMP]
blocks: []
completed: 2026-04-17
completed_by: core-engineer
outcome: landed-as-cosmetic
---

# VERF: V_local erf Coulomb Subtraction

> **Completion note (2026-04-17):** Landed the erf-subtraction convention even though it is numerically equivalent to the previous bare-Coulomb form on the current log mesh. Rationale (engineer's call): (1) matches QE's convention exactly, which simplifies the upcoming VGCMP V_local(G) cross-check; (2) produces a bounded integrand near r=0, which will matter for future high-Z pseudopotentials where the bare-Coulomb integrand can reach very large values at the first radial grid points. A regression test at `tests/vloc_erf_consistency.rs` pins the mathematical equivalence for Si's first 20 |G| shells (< 1e-6 eV absolute) so any future deviation fires loudly. **This change does NOT close the 13.4 eV Si gap — see the new VGCMP proposal for that investigation.**
>
> **Prerequisite:** Implement Proposal 38 (Simpson's rule) first. This proposal is only needed if Simpson alone doesn't bring Si within 0.1 eV of QE.

## Problem

Our V_local(G) form factor computation subtracts the full Coulomb tail `Ze^2/r` from V_local(r), creating a divergent integrand at r=0. QE uses a smoother decomposition with `erf(r)/r`, producing a bounded integrand everywhere.

### Current approach (`src/pseudopotential/mod.rs:115-132`)

For G != 0:

```text
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

```text
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

```text
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

## 2026-04-17 — Attempt 1: Negative Result

A Core Engineer agent implemented the exact change proposed above (full patch archived at `/tmp/pwdft-rescue/VERF-attempt.diff`, ~70 lines in `src/pseudopotential/mod.rs`). Full test suite ran with the change applied:

| System | QE (eV) | Ours before | Ours after VERF | Δ vs QE |
|--------|---------|-------------|-----------------|---------|
| Si diamond | -231.61 | -218.18 | -218.18 | **13.43 eV (unchanged)** |
| Fe BCC | -3059.46 | -3059.44 | -3059.44 | 0.02 eV (no regression) |
| C diamond | — | does not converge | does not converge | — |

All other 200+ unit tests pass; clippy clean. The change is algebraically correct but numerically identical to the current bare-Coulomb approach.

### Why VERF alone doesn't help

The two decompositions are mathematically equivalent by construction. The numerical advantage of erf subtraction requires that the radial quadrature be stressed by the `1/r` integrand near the origin. Post-SIMP (Simpson's rule over log-mesh) the current bare-Coulomb integrand is already handled to ~1e-6 precision — the `r²` factor is enough to tame the `1/r` growth at the first few grid points on the log mesh. The erf form gives the same integral to machine precision.

### Implication for the Si 13.4 eV gap

**VERF is NOT the root cause.** The Si/QE discrepancy must come from elsewhere. Candidates to investigate in order of likelihood:

1. **Kleinman–Bylander non-local projectors** — `D_ij` handling for UPF files where QE diagonalized the raw `h^l_ij` block and absorbed the rotation. Tests 05/06/08/10 pass, but cross-k-point or cross-l mixing may differ.
2. **NLCC (core density)** — Fe passes because NLCC was already cross-checked for that PP; Si has no NLCC, so this is unlikely to be the cause.
3. **Ewald parameters** — alpha, G-cutoff, real-space shell — but Fe passes, which has far more Ewald contribution per-atom, so unlikely.
4. **Kinetic G-set truncation** — cutoff-sphere vs cutoff-FFT-grid mismatch. Worth checking that |G|² ≤ 2·ecut for wavefunctions is applied identically to QE.
5. **Local pseudopotential tail** — comparing `V_local(G)` values for the first 20 G-shells against QE numerically is the smoking-gun test.

### Recommended Next Step

Before keeping the VERF code change, add a test that directly compares `V_local(G)` values for Si's first 20 G-shells against QE's `vloc.dat` output. If they agree to ~1e-6 Ry (expected), VERF is confirmed cosmetic and the proposal should be **archived as not-needed**. If they disagree, then there is a deeper bug in either quadrature or the PP data ingestion, and VERF belongs on the fix path.

Either way, **the Si 13.4 eV gap requires a different proposal**. Suggest opening `VGCMP — V_local(G) cross-validation vs QE` to pinpoint the source.
