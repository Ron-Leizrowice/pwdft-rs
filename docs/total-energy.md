# Total Energy

## Kohn-Sham Total Energy Functional

```
E[ρ] = T_s[ρ] + E_H[ρ] + E_xc[ρ] + E_ext[ρ] + E_ion-ion
```

In practice, using the band energy identity (post-TSEN, matches QE's
`! total energy` which is the Mermin free energy `F = E − TS`):

```
E_total = E_band - E_H + E_xc - E_vxc + E_ewald + V_local(G=0) × N_el − TS
```

The `−TS` term is the Mermin smearing-entropy contribution (zero for
insulators and for the `SmearingScheme::Fixed` path). The pre-TSEN
internal energy `E_internal = E_total + TS` is still recoverable from
`ScfResult::total_energy + entropy_ts`.

**Code:** `src/scf/energy.rs` (assembly in `total_energy` and per-component `EnergyComponents`), and the final pass in `src/scf/driver.rs` / `src/scf/driver_spin.rs` which adds `V_local(G=0) * N_el`.

### Component definitions

| Term | Formula | Code |
|------|---------|------|
| E_band | `Σ_{n,k} f_{n,k} w_k ε_{n,k}` | `energy.rs:15-31` |
| E_H | `(Ω/2) Σ_{G≠0} \|ρ(G)\|² 4πe²/\|G\|²` | `energy.rs:34-49` |
| E_xc | `∫ ε_xc(r) ρ(r) dr` | `xc.rs:64-72` |
| E_vxc | `∫ V_xc(r) ρ_val(r) dr` | `energy.rs:66-70` |
| E_ewald | Ewald ion-ion energy | `ewald.rs` |
| V_local(G=0) × N_el | Constant shift from pseudopotential | `context.rs:81-83` |
| −TS | Smearing-entropy Mermin correction | `smearing.rs::entropy_ts` |

### Why this works

E_band counts E_H once (through V_H in V_eff) and E_vxc once (through V_xc).
We subtract E_H (which over-counts by 1/2) and replace E_vxc with E_xc (the
true XC energy, not the potential integral).

## NLCC Modification

With nonlinear core correction (Louie, Froyen, Cohen, PRB 26, 1738, 1982):

```
E_xc = ∫ ε_xc[ρ_val + ρ_core] × (ρ_val + ρ_core) dr
E_vxc = ∫ V_xc[ρ_val + ρ_core] × ρ_val dr
```

E_xc uses total density in both functional and integration measure.
E_vxc uses only valence density because that is what the eigenvalues contain.

**Code:** `src/scf/driver.rs` / `src/scf/driver_spin.rs` (core density added via `add_core_density` before each XC evaluation); `src/scf/energy.rs` `xc_energy_corrected` (NLCC double-counting).

### Audit status

| Item | Status | Verified against |
|------|--------|-----------------|
| Total energy formula | CORRECT | QE `electrons.f90` |
| Band energy | CORRECT | Standard definition |
| Hartree energy | CORRECT | Positive-definite, G=0 excluded |
| XC double-counting | CORRECT | QE `v_of_rho.f90` |
| NLCC prescription | CORRECT | Louie et al. 1982 |
| V_local(G=0) handling | CORRECT | QE `setlocal.f90` |
