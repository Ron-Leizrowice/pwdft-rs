---
id: PCRS
status: completed
priority: medium
complexity: small
risk: low
depends_on: [VGC5]
blocks: []
---

# PCRS: Per-Component Energy Residual Investigation

## Status

**COMPLETED — 2026-04-17.** Root cause identified; follow-up fix
proposal opened as **PCFX** (density symmetrization at FFT grids
incompatible with space-group fractional translations).

## Origin

VGC5 Phase 5 (PR #34, 2026-04-17) instrumented per-component SCF
energies and observed that for Si diamond, the identity

```
E_total = E_kinetic + E_local + E_local_G0_shift + E_nonlocal
        + E_Hartree + E_xc + E_ewald
```

closes to only **1.20 eV** at `conv_threshold = 1e-8`. Fe BCC shows
0.045 eV at the same conditions.

The first-pass hypothesis (from the VGC5 session) was that a
convergence criterion gap or a rescale step was responsible. This
proposal ran a conv_threshold sweep to separate SCF noise from a
structural bookkeeping bug.

## Findings (2026-04-17)

### 1. QE identity closes to bit-exact

From `qe_validation/si_scf.out`:

```
one_electron + hartree + xc + ewald = -17.02213518 Ry
internal E (= F - (-TS))            = -17.02213518 Ry
Δ = 0.000e+00 Ry
```

QE's decomposition is exact to the printed precision (10 digits
Ry ≈ 1e-9 eV). **Any residual in pwdft-rs is ours, not an identity
artifact.**

### 2. pwdft-rs conv_threshold sweep (Si, ecut=15 Ry, 4×4×4 MP)

Running `scripts/validate/pcrs_residual_scan.py`:

| conv_thr | iters | final Δρ | E_KS (eV) | E_sum (eV) | \|Δ\| (eV) |
|----------|-------|----------|-----------|------------|-----------|
| sym ON  1e-6  | 10   | 1.18e-09 | -218.180513 | -216.976374 | **1.204** |
| sym ON  1e-8  | 11   | 3.15e-11 | -218.180513 | -216.976374 | **1.204** |
| sym ON  1e-10 | 159  | 9.59e-11 | -218.180513 | -216.976374 | **1.204** |
| sym OFF 1e-6  | 11   | 1.95e-09 | -218.163250 | -218.163250 | 1.07e-07  |
| sym OFF 1e-8  | 12   | 1.10e-09 | -218.163250 | -218.163250 | 1.03e-08  |
| sym OFF 1e-10 | 128  | 9.73e-11 | -218.163250 | -218.163250 | 1.39e-09  |

- **With symmetry ON:** residual is a hard plateau at 1.204 eV. Tightening
  `conv_threshold` by 4 orders of magnitude does not move it one digit.
  This proves the residual is not SCF noise.
- **With symmetry OFF:** residual scales with Δρ at the O(Δρ) rate
  expected for self-consistency noise — falls from 1e-7 eV at conv=1e-6
  to 1.4e-9 eV at conv=1e-10.
- **Ratio sym-ON / sym-OFF ≈ 3 × 10⁷.**

E_KS itself also differs: **−218.180 eV (sym ON) vs −218.163 eV (sym OFF)**.
The 17 meV gap is the same symmetrization artifact leaking into the total
energy via the double-counting correction terms.

### 3. Fe BCC spot check

| Config            | Residual  |
|-------------------|-----------|
| Fe, sym ON  (CLI, IBZ=10 kpts) | 7.30e-02 eV |
| Fe, sym OFF (CLI, 64 kpts)     | 1.42e-07 eV |

Same pattern as Si, but smaller magnitude. Fe's Im-3m is symmorphic
(all τ=0) so the fractional-translation artifact is absent; the
residual comes from a different, smaller source (TBD, possibly pure
rotation aliasing or k-point weight rounding).

### 4. Root cause — non-symmorphic symmetrization on an incompatible FFT grid

The SCF density accumulator produces ρ_ψ(r) from the wavefunction sum.
VGC5 builds `EnergyComponents` by computing:

- `e_kinetic`, `e_nonlocal` directly from ψ coefficients (i.e. consistent
  with ρ_ψ — the **un-symmetrized** density).
- `e_local`, `e_hartree`, `e_xc` from grid integrals over `rho_r_new`,
  which is **symmetrized** (`src/scf/mod.rs:395`).

If ρ_ψ == ρ_sym to machine precision, the two sets match and the
identity closes. In practice:

- `src/symmetry/density.rs::symmetrize_density` applies each space-group
  operation `S = {R | τ}` by mapping grid indices via a `round()`-based
  `nint` at `src/symmetry/density.rs:12-16` / `:75-77`.
- For Si Fd-3m (space group 227), non-symmorphic generators have
  fractional translation **τ = (1/4, 1/4, 1/4)**.
- Our FFT grid for Si at ecut=15 Ry is **18 × 18 × 18**, and
  **18 is not divisible by 4**. So `nint(n · τ_i)` rounds a half-integer
  (`4.5 → 5` in Rust), and the inverse operation rounds back to a
  **different** grid point. The symmetry average therefore smears
  rather than copies the density.

- `src/symmetry/density.rs::check_grid_compatibility` only validates
  **rotation** compatibility (`R_{ij} · n_j ≡ 0 mod n_i`), not
  translations (`τ_i · n_i ∈ ℤ`). It returns `true` for the Si 18-grid
  even though the grid cannot represent the (1/4, 1/4, 1/4) shift
  exactly. And the SCF path does not consult the function at all —
  the grid is sized purely from `ecutrho_ratio` / explicit `fft_grid`.

- For Fe BCC (Im-3m, space group 229) all operations are symmorphic
  (τ = 0) so the fractional-translation artifact vanishes and the
  residual drops by ~2 orders of magnitude.

QE sidesteps this by **symmetrizing ρ(G) in reciprocal space** using
`ρ_sym(G) = (1/N_ops) Σ_S ρ(R^{-1} G) · exp(i G · τ)`. Phase factors are
exact for any τ; there is no grid-discretization error. That is why QE's
per-term identity closes bit-exactly.

### 5. Not a bug in `total_energy` itself

Subtracting the sum of components from `total_energy` and back-solving
for the implied `e_vxc` reproduces the correct self-consistent identity:

```
e_vxc_from_total_energy_bookkeeping    = -85.3128 eV  (∫ρ_sym · V_xc dr)
e_vxc_implied_by_band_identity          = -86.5169 eV  (∫ρ_ψ · V_xc dr)
difference                               = +1.2041 eV  ≡  the residual
```

The residual is precisely `∫(ρ_ψ - ρ_sym) · V_eff dr`. That integrates
~1 eV for Si at ecut=15 Ry, 18-grid, Fd-3m. All three `total_energy`,
`EnergyComponents`, and `Harris-Foulkes` are internally consistent; the
flaw is that `ρ_sym` ≠ `ρ_ψ` at the current grid+symmetry combination.

## Verdict

**Real bug identified; fix requires a source-code change.** The
investigation is complete and should be marked completed.

## Follow-up

Open new proposal **PCFX — fix density symmetrization for
non-symmorphic groups on incompatible FFT grids** (see
`proposals/PCFX-symmetrize-rho-g-space.md`). Options:

1. **Preferred:** symmetrize ρ in G-space using the QE approach
   (`ρ_sym(G) = (1/N_ops) Σ_S exp(i G · τ_S) · ρ(R_S^{-1} G)`). Exact
   for any τ; no discretization error.
2. **Stop-gap:** extend `check_grid_compatibility` to also require
   `τ_i · n_i ∈ ℤ` for every operation and round the FFT grid up to
   satisfy it. Keeps real-space symmetrization but with compatible
   grids.
3. **Fallback:** subsitute a trilinear interpolation for sub-grid shifts
   in the existing real-space loop. Cheap but still inexact.

Until PCFX lands, anyone validating per-component energies against QE
must either disable symmetry or inflate the identity tolerance to
the observed plateau (~1.2 eV for Si-like systems, ~50 meV for
symmorphic ones).

## Artifacts

- `scripts/validate/pcrs_residual_scan.py` — conv_threshold sweep, with
  sym-on / sym-off comparison. `uv run` compatible.

## References

- `proposals/VGC5-per-component-energy-accounting.md` — parent.
- `proposals/PCFX-symmetrize-rho-g-space.md` — follow-up fix.
- `src/scf/mod.rs:395` — `symmetrize_density` call in non-spin SCF.
- `src/symmetry/density.rs:26-91` — real-space symmetrization.
- `src/symmetry/density.rs:100-113` — `check_grid_compatibility`
  (rotations only).
- QE `PW/src/symme.f90::sym_rho` — reciprocal-space symmetrization.
- `proposals/completed/24-symmetry-density-wrapping-fix.md` — prior
  partial fix that addressed `round()` vs `floor()` consistency but
  didn't tackle non-symmorphic τ.
