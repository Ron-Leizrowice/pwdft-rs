# Total Energy

## Kohn-Sham Total Energy Functional

```
E[ρ] = T_s[ρ] + E_H[ρ] + E_xc[ρ] + E_ext[ρ] + E_ion-ion
```

In practice, using the band energy identity:

```
E_total = E_band - E_H + E_xc - E_vxc + E_ewald + V_local(G=0) × N_el
```

**Code:** `src/scf/energy.rs:76-83` (assembly), `src/scf/mod.rs:278-283` (V_local(G=0) addition)

### Component definitions

| Term | Formula | Code |
|------|---------|------|
| E_band | `Σ_{n,k} f_{n,k} w_k ε_{n,k}` | `energy.rs:15-31` |
| E_H | `(Ω/2) Σ_{G≠0} \|ρ(G)\|² 4πe²/\|G\|²` | `energy.rs:34-49` |
| E_xc | `∫ ε_xc(r) ρ(r) dr` | `xc.rs:64-72` |
| E_vxc | `∫ V_xc(r) ρ_val(r) dr` | `energy.rs:66-70` |
| E_ewald | Ewald ion-ion energy | `ewald.rs` |
| V_local(G=0) × N_el | Constant shift from pseudopotential | `context.rs:81-83` |

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

**Code:** `src/scf/mod.rs:275-276` (core density addition), `energy.rs:55-73` (corrected XC)

### Audit status

| Item | Status | Verified against |
|------|--------|-----------------|
| Total energy formula | CORRECT | QE `electrons.f90` |
| Band energy | CORRECT | Standard definition |
| Hartree energy | CORRECT | Positive-definite, G=0 excluded |
| XC double-counting | CORRECT | QE `v_of_rho.f90` |
| NLCC prescription | CORRECT | Louie et al. 1982 |
| V_local(G=0) handling | CORRECT | QE `setlocal.f90` |
