# Pitfalls (Lessons Learned)

Known issues encountered during development and their fixes.

## 1. NLCC Omission (Fe Bug)

**Symptom:** Constant eigenvalue shift (~15 eV for Fe) relative to QE.

**Cause:** `PP_NLCC` data ignored; XC evaluated on `rho_val` instead of
`rho_val + rho_core`. Si works because it has no NLCC.

**Fix:** Parse `PP_NLCC`, Bessel-transform to FFT grid, add to density before
XC evaluation. See Proposal 25.

**Code:** `src/scf/driver.rs` / `src/scf/driver_spin.rs` (core density added via `add_core_density` before each XC call).

## 2. V_local(G=0) Convention

**Symptom:** All eigenvalues shifted by a large constant.

**Cause:** `V_local(G=0)` is pseudopotential-dependent and arbitrary.

**Fix:** Exclude from Hamiltonian, add `V_local(G=0) * N_electrons` to total
energy. Matches QE convention.

**Code:** `src/scf/context.rs:80-83`

## 3. FFT Normalization

**Symptom:** Energies off by factors of N_grid or Omega.

**Cause:** Inconsistent 1/N vs 1/Omega normalization between density, potential,
and energy formulas.

**Convention:** Forward FFT is unnormalized; callers divide by N_grid to get
Fourier coefficients. `rho(G=0)` is the spatial average density. Hartree
energy includes explicit `* Omega` factor.

**Code:** `src/fft.rs` (convention), `src/scf/energy.rs:106-109` (normalization)

## 4. Spin Factor in Occupations

**Symptom:** Wrong electron count or doubled energies in spin-polarized.

**Cause:** Occupation function returns `f in [0, spin_factor]` where
`spin_factor = 2/nspin`. Must be consistent in density, energy, and Fermi
energy search.

**Convention:** `occupation()` includes spin_factor. Fermi search uses the same
function. Density accumulates `f * w_k * |psi|^2`.

## 5. E_xc vs E_vxc with NLCC

**Symptom:** Small but systematic energy error with NLCC pseudopotentials.

**Cause:** Using `rho_total` in both `E_xc` and `E_vxc`.

**Fix:** `E_xc = integral(epsilon_xc * rho_total)` but
`E_vxc = integral(V_xc * rho_val)`. The double-counting correction uses
valence-only because that is what the eigenvalues contain.

## 6. Radial Quadrature Accuracy

**Symptom:** 13-45 eV total energy discrepancy vs QE, broken eigenvalue
degeneracies at high-symmetry points. Error scales with Z_valence
(Si: 1.66 eV/el, Fe: 2.84 eV/el).

**Cause:** Two compounding issues:

1. Plain sum `integral += f * dr` (O(h^2)) vs QE's Simpson's rule (O(h^4))
2. V_local Coulomb subtraction adds `Ze^2/r` which diverges at r=0, while
   QE uses `erf(r)/r` which stays finite

The G-dependent quadrature error breaks crystal symmetry numerically, splitting
eigenvalues that should be degenerate.

**Fix:** Proposals 38 (Simpson's rule) and 39 (erf subtraction).

**QE reference:** `upflib/simpsn.f90`, `upflib/vloc_mod.f90:138`

See [Radial Integration](radial-integration.md) for details.
