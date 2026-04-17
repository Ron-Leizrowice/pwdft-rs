---
id: NCFX
status: active
priority: critical
complexity: small
risk: low
depends_on: [VGC5]
blocks: []
owner: core-engineer
---

# NCFX: NLCC Core-Density Unit and Radial-Weight Fix

## Origin

VGC5 (VGCMP Phase 5) per-component energy accounting (this branch) localized
the Si 13.4 eV vs QE gap almost entirely to the **XC term**:

| system | ΔE_xc (ours − QE) | ΔE_total | dominant? |
|--------|-------------------|----------|-----------|
| Si diamond    | +13.74 eV | +13.43 eV | yes |
| Fe BCC (nspin=1) | −48.85 eV | −41.08 eV | yes |

Other terms (Ewald, kinetic, local, non-local, Hartree) match QE to
≤ 2.4 eV total, with Ewald matching to ≤ 0.012 eV. The XC component
drifts far more than any other term, and the drift **scales with NLCC
core-density magnitude** (Si: moderate NLCC; Fe: large NLCC).

## Root cause

Two compounding bugs in how PP_NLCC is read and Fourier-transformed.

### Bug 1: Wrong unit conversion on parse

`src/pseudopotential/upf.rs:94-105` divides the raw PP_NLCC block by
`BOHR_TO_ANG` (= 0.529177 Å/Bohr). The comment claims PP_NLCC stores
`4πr²·ρ_core(r) in e/Bohr`, so the conversion to e/Å mirrors the
PP_RHOATOM treatment.

**PP_NLCC actually stores the bare `ρ_core(r)` in e/Bohr³.** QE's
`upflib/rhoc_mod.f90:107-108` explicitly multiplies `upf%rho_atc` by
`rgrid%r2` in the Bessel transform, confirming it is NOT pre-weighted.
Numeric check: the first PP_NLCC value in `pseudopotentials/nc/lda/Si.upf`
is `0.229` at `r ≈ 0`, inconsistent with `4π·r²·ρ` (which vanishes at
r = 0) but consistent with a bare density.

The conversion factor should therefore be `1/BOHR_TO_ANG³ ≈ 6.748`, not
`1/BOHR_TO_ANG ≈ 1.890`. **The stored values are currently low by a
factor of `BOHR_TO_ANG² ≈ 0.280`** (about 3.6×).

### Bug 2: Missing radial weight and 4π in the Bessel transform

`src/scf/potentials.rs:92-107` computes the G-space core charge as

```rust
integrand = rho_c · j₀(Gr)
integral  = simpson_integrate(integrand, rab)
rho_core(G) = struct_factor · integral / Ω
```

i.e. `ρ_core(G) = (1/Ω) ∫ ρ_core(r) · j₀(Gr) dr`.

The correct radial FT of a spherically symmetric density is

```
ρ_core(G) = (4π / Ω) ∫ ρ_core(r) · j₀(Gr) · r² dr
```

(see QE's `init_tab_rhc` at `upflib/rhoc_mod.f90:107-115`: aux = rho·r²·j₀,
final multiply by `fpi/omega`). **The `r²` weight and the `4π` prefactor
are both missing from pwdft-rs.**

### Combined error signature

Together, the two bugs give a G=0 core charge roughly

```
ρ_core(G=0)[pwdft-rs] ≈ (1/BOHR_TO_ANG) · ∫ρ_core(r) dr / Ω
ρ_core(G=0)[QE]      = (4π/BOHR_TO_ANG³) · ∫ρ_core(r)·r² dr / Ω
```

The integrals differ dimensionally (Å vs Å³) and by a factor of 4π. The
net discrepancy scales with the magnitude and spatial extent of the core
charge, matching the observed Si/Fe asymmetry.

## Fix

Both bugs sit in two small code blocks:

1. **`src/pseudopotential/upf.rs`** (≈ lines 94-109):
   - Update the comment to state that PP_NLCC is bare `ρ_core(r)` in
     `e/Bohr³`.
   - Change the parse-time conversion from `/BOHR_TO_ANG` to
     `/BOHR_TO_ANG.powi(3)` so `PseudopotentialData.core_charge` is in
     `e/Å³`.

2. **`src/scf/potentials.rs::compute_core_density`** (≈ lines 92-107):
   - Change the integrand from `rho_c · j₀(Gr)` to
     `rho_c · r² · j₀(Gr)` (all in Å, consistent units).
   - Multiply the integral by `4π` before dividing by `Ω`.

No API changes; only two files edited. Both changes are local and
mechanical.

## Verification

1. Rebuild and rerun `cargo test --test vgc5_per_component_si --
   --nocapture`. Expect `ΔE_xc` for Si to drop from +13.74 eV to below
   ~0.1 eV. Expect Fe `ΔE_xc` to drop correspondingly.

2. Remove `#[ignore]` from `test_si_diamond_vs_qe` (and any Tier 2 tests
   that unblock) in `tests/qe_validation.rs`. Si should match QE total
   energy to within 0.05 eV.

3. Update the VGC5 integration-test pins (`tests/vgc5_per_component_si.rs`)
   after the fix lands; bump Si XC to QE-match numbers.

4. Add a new Python-referenced unit test analogous to VGCMP Phase 1/2/3:
   compute ρ_core(G) with `scipy.integrate.simpson` on the UPF mesh for a
   few G shells, compare to the Rust implementation to ≤ 1e-4 e/Å³.

## Reproduction data (VGC5 baseline, pre-fix)

From `tests/vgc5_per_component_si.rs` pins (a161221 + VGC5 patch):

Si diamond, ecut=15 Ry, 4×4×4 MP, FD smearing σ=0.01 Ry:
- E_kinetic       =   82.866 eV      (QE part of "one-electron")
- E_local(G≠0)    =  -58.468 eV
- E_local(G=0)·N  =   10.745 eV
- E_nonlocal      =   33.465 eV
- E_hartree       =   13.593 eV      (QE: 15.104,  Δ = −1.51)
- E_xc            =  -70.658 eV      (QE: −84.396, Δ = **+13.74**)
- E_ewald         = -228.519 eV      (QE: −228.530, Δ = +0.011)
- E_total         = -218.181 eV      (QE: −231.610, Δ = +13.43)

Fe BCC (nspin=1, 4×4×4 MP, σ=0.02 Ry) shows the same XC-dominated pattern
with Δ_xc = −48.85 eV.

## References

- QE source: `qe-7.5/upflib/rhoc_mod.f90:101-120` (`init_tab_rhc`)
- QE source: `qe-7.5/PW/src/set_rhoc.f90:29-125` (`set_rhoc`)
- Louie, Froyen, Cohen, PRB **26**, 1738 (1982) — NLCC theory
- VGC5 measurement: `tests/vgc5_per_component_si.rs`,
  `scripts/validate/vgc5_per_component.py`,
  `scripts/validate/vgc5_qe_{si,fe}_components.csv`
- Pre-VGC5 VGCMP chain: see `proposals/VGCMP-vloc-g-cross-check.md`
  and `proposals/VGC5-per-component-energy-accounting.md`
