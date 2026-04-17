---
id: PCFX
status: active
priority: high
complexity: medium
risk: medium
depends_on: [VGC5]
blocks: []
---

# PCFX: Density Symmetrization in G-Space (Fix for Non-Symmorphic τ)

## Origin

`proposals/PCRS-per-component-residual.md` traced a **1.204 eV**
per-component energy identity residual in Si to a real-space density
symmetrization that cannot exactly represent the Fd-3m fractional
translation τ = (1/4, 1/4, 1/4) on the current 18³ FFT grid.

The same symmetrization also injects a **~17 meV** error into the
converged total energy `E_KS` itself (the `e_local`, `e_hartree`,
`e_xc` integrals all use ρ_sym while the diagonalization sees a
slightly different effective potential built from the
mixer-input density).

For symmorphic groups (Fe Im-3m, Al Fm-3m, etc.) the plateau is
~45–75 meV — smaller but still above our target tolerance for QE
validation.

## Mathematical background

The exact reciprocal-space symmetrization formula is

```
ρ_sym(G) = (1 / N_ops) Σ_S  exp(i G · τ_S) · ρ(R_S^{-1} G)
```

where `S = {R_S | τ_S}` ranges over all space-group operations.
The phase factor `exp(i G · τ_S)` is analytic in τ_S, so the
symmetrization is exact for any fractional translation —
**no FFT-grid discretization error**. QE uses this form
(`PW/src/symme.f90::sym_rho`, `symme_module` module).

The real-space form is equivalent only when
`τ_{S,i} · n_i ∈ ℤ` for every operation S and grid axis i, i.e.
the grid captures every τ on an exact grid point. For Si Fd-3m
on an 18-grid this fails (18·¼ = 4.5 is not integer), producing
an O(1 eV) bias.

## Plan

### Option A — G-space symmetrization (preferred)

Add `src/symmetry/density_g.rs` with:

```rust
/// Symmetrize ρ on the FFT grid via its G-space representation.
///
/// ρ_sym(G) = (1/N_ops) Σ_S exp(i G · τ_S) · ρ(R_S⁻¹ G)
pub fn symmetrize_density_g(
    rho_r: &mut [f64],
    dims: [usize; 3],
    recip: &Matrix3<f64>,
    fft: &mut FFT3D,
    symmetry: &SymmetryInfo,
) {
    // 1. FFT ρ(r) -> ρ(G).
    // 2. For each G, sum exp(i G·τ_S) · ρ(R_S⁻¹ G) over all S.
    // 3. Divide by N_ops.
    // 4. IFFT back to real space.
}
```

Replace `symmetrize_density` call sites in `src/scf/mod.rs:395` and
`:743-744` with the new G-space version.

Performance: one extra forward + inverse FFT per iteration, plus an
O(N_grid · N_ops) scalar loop. On the Si 18³ × 48 ops × 10 iters budget
that's ~16 Mflops — negligible.

Memory: none beyond a reused `rho_g` buffer (already present via
`density_r_to_g`).

### Option B — enforce compatible grid only (stop-gap)

Extend `check_grid_compatibility` to also check translations:

```rust
pub fn check_grid_compatibility(dims, symmetry) -> bool {
    // … existing rotation check …
    for op in &symmetry.operations {
        for i in 0..3 {
            let shift = op.translation[i] * dims[i] as f64;
            if (shift - shift.round()).abs() > 1e-10 {
                return false;
            }
        }
    }
    true
}
```

And in `ScfContext::new`, round the FFT grid up to the nearest
compatible size before calling the FFT factory. For Si this would
force 18 → 20 or 24 (both divisible by 4).

Downside: actually increases grid size, costing FFT + memory. And
empirically at 20³ and 24³ the Anderson mixer loses conditioning
(see PCRS notes); so this option needs a mixer-robustness fix too.

### Recommendation

Do Option A. It's the correct physics, matches QE, has no grid-size
downside, and the extra FFT is cheap. Option B should be kept as a
`debug_assert!` fallback: when Option A is active, skipping it when
the grid is compatible gives bit-identical results to the current
real-space path.

## Verification

1. **Regression:** `tests/vgc5_per_component_si.rs` identity check
   `Σ(components) − E_total` drops from 1.204 eV to ≤ 1e-6 eV at
   conv_threshold = 1e-8.
2. **Convergence to QE:** Si total energy at ecut=15 Ry, 4×4×4 MP
   should shift by ~17 meV (removing the symmetrization-induced
   bias in `e_xc - e_vxc`), closing part of the 13.43 eV gap to QE.
3. **Free-electron invariance:** band structure on an empty lattice
   with symmetry ON must stay bit-identical to symmetry OFF
   (already tested in `tests/free_electron_bands.rs`).
4. **Round-trip:** symmetrize a purely-symmetric test density
   (`ρ_test(r) = Σ_i c_i · exp(i G_i · r)` with G_i in an already
   symmetric star) and assert
   `max |ρ_sym − ρ| < 1e-12 · ||ρ||_∞`.

Add a new `tests/symmetrize_rho_g.rs` with unit coverage for
the four point-group types represented in our validation suite:
Fd-3m (Si, diamond), Fm-3m (Al, NaCl), Im-3m (Fe, Cu), P6_3/mmc
(no current test case but cheap to synthesize).

## Scope

- `src/symmetry/density.rs`: add `symmetrize_density_g` and its
  helpers; leave the real-space `symmetrize_density` in place as
  a lower-level helper for the free-electron test.
- `src/symmetry/mod.rs`: export the new entry point.
- `src/scf/mod.rs`: two one-line call-site swaps (non-spin line 395,
  spin lines 743-744).
- `src/scf/context.rs`: no change.
- `check_grid_compatibility`: document as "real-space-only
  sanity check; do not use for Option-A path".
- `tests/symmetrize_rho_g.rs`: new integration test (see
  §Verification).

## Risk

- `check_grid_compatibility` is exposed as `pub`; other crates in
  the workspace (not currently) could depend on its behavior. Check
  before touching.
- `symmetrize_density` (real-space) is used by a couple of unit
  tests; preserve it as a direct-access helper.
- GPU feature flag: the new FFT pair is CPU-only in the first
  iteration; a GPU-resident symmetrization is out of scope here.

## References

- QE 7.5: `qe-7.5/PW/src/symme.f90` — see `sym_rho_init_shells`
  and `sym_rho` for the G-space phase loop.
- Martin, "Electronic Structure," §C.2 (space group action on plane
  waves).
- pwdft-rs `src/symmetry/operations.rs` — existing `R^{-1}`, τ data.
- PCRS proposal, tables & residual-scan script.
