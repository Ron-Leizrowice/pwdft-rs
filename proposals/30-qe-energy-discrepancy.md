# Proposal 30: Systematic Energy Discrepancy vs Quantum ESPRESSO

## Problem

pwdft-rs SCF total energies differ from QE 7.5 by 13–45 eV depending on the material, using the same pseudopotentials (PseudoDojo ONCV LDA), same ecut (15 Ry), same k-grid (4×4×4), and same smearing (Fermi-Dirac σ=0.01 Ry). Eigenvalue degeneracies that should be exact by symmetry are broken.

### Reproduction

```bash
# Run QE references (results already in qe-7.5/runs/{si,c,fe}_pseudodojo/)
cargo test --release --test qe_validation -- --nocapture
```

### Measured discrepancies

| System | Atoms | Z_val | QE energy (eV) | Our energy (eV) | ΔE (eV) | ΔE/el (eV) | Status |
|--------|-------|-------|-----------------|------------------|---------|-------------|--------|
| Si diamond | 2 | 4 | -231.61 | -218.28 | 13.3 | 1.66 | Converges (11 iter) |
| C diamond | 2 | 4 | -312.51 | — | — | — | Does not converge |
| Fe BCC | 1 | 16 | -3059.46 | -3104.83 | 45.4 | 2.84 | Converges (16 iter) |

### Eigenvalue degeneracy breaking

In Si (FCC) at Gamma, bands 2-4 should be triply degenerate. QE gives 6.080, 6.080, 6.080 eV. We give -1.04, 3.75, 3.76 eV — completely wrong pattern and magnitude.

In Fe (BCC) at Gamma, bands 2-4 (3p semicore) should be triply degenerate. QE gives -46.50, -46.50, -46.50 eV. We give -44.86, -44.86, -38.56 eV — one state is split by 6 eV.

## Research: What's validated, what's suspect

### Validated components

| Component | Test | Status |
|-----------|------|--------|
| **Ewald energy** | Fe: ours = -2337.167 eV, QE = -2337.173 eV | **Match (0.006 eV)** |
| **Basis size** | Fe Gamma: 79 PWs, matches QE exactly | **Match** |
| **V_NL Hermiticity** | Max violation = 2.2e-16 | **Correct** |
| **Free-electron eigenvalues** | Band 0 = 0, band 1 = 36.52 eV (BCC Fe) | **Correct** |
| **PZ correlation constants** | Verified against Table I of PZ 1981 | **Correct** |
| **Smearing/occupations** | Extensive unit tests, nspin=1/2 consistency | **Correct** |

### Suspect components (ordered by likelihood)

**1. V_local Fourier transform — HIGH SUSPICION**

`v_local_of_g()` in `src/pseudopotential/mod.rs:97-132` computes the spherical Bessel transform of V_local(r) analytically using Simpson-like integration on the PP radial grid. This is the single largest potential contribution and a subtle numerical integral.

Evidence: Si V_local(G=0) = 5.17 eV, Fe V_local(G=0) = 5.17 eV. These are suspiciously similar for very different elements. Fe with Z=16 should have a much larger (more negative?) short-range potential integral.

To verify: extract V_local(G) from QE's XML output and compare point-by-point. QE stores this in `tmp/<prefix>.save/charge-density.hdf5` or can be printed with pp.x.

**2. XC double-counting with NLCC — MEDIUM SUSPICION**

The energy expression is:
```
E_total = E_band - E_H + (E_xc[ρ_val+ρ_core] - E_vxc) + E_ewald + V_local(G=0)·N_el
```

In `xc_energy_corrected()` (`src/scf/energy.rs:55-73`), E_vxc uses `rho_val` (without core), but `vxc_r` was computed from `rho_val + rho_core`. This is correct per the NLCC prescription: the potential is from total density, but the double-counting integral uses only valence density. However, the `exc_r` used in E_xc is also from `rho_val + rho_core`, meaning:

```
E_xc = ∫ (ρ_val + ρ_core) · ε_xc(ρ_val + ρ_core) dr
E_vxc = ∫ ρ_val · V_xc(ρ_val + ρ_core) dr
```

QE does the same (see `v_of_rho.f90`). So this is likely correct but should be verified by comparing individual energy components.

**3. Projector unit conversion — MEDIUM SUSPICION**

UPF stores β projectors as `χ(r) = r·β(r)` in Bohr^{-1/2}. We convert to Å^{-1/2} by dividing by `√(BOHR_TO_ANG)` (`src/pseudopotential/upf.rs:68-71`). D_ij converts from Ry to eV (`dij_ry * RY_TO_EV`). The dimensional chain:

```
V_NL matrix element ~ F_i × D_ij × F_j / Ω
F_i = 4π ∫ χ(r) j_l(qr) r dr  [Å^{-1/2} · Å · Å = Å^{3/2}]
D_ij [eV]
Ω [ų]
V_NL ~ [ų] × [eV] / [ų] = [eV] ✓
```

But the actual numerical prefactor matters. QE's `init_us_2.f90` applies `(4π/Ω) × tpiba` factors that may differ from our convention. A 2π or √(2π) factor error would shift all eigenvalues uniformly.

**4. FFT normalization or grid mismatch — LOW-MEDIUM SUSPICION**

Our FFT convention: forward is unnormalized, we divide by N manually. QE uses the same convention. But if there's a mismatch in how V_eff(G-G') is indexed or normalized, it could shift eigenvalues.

The broken degeneracies suggest the potential doesn't respect the full crystal symmetry. Possible causes:
- FFT grid dimensions not respecting the point group (e.g., different dims along different axes for a cubic cell)
- Numerical noise in the radial integrals for V_local or projectors breaking the angular symmetry

**5. V_local(G=0) convention — LOW SUSPICION**

Our code excludes V_local(G=0) from the Hamiltonian and adds `V_local(G=0) × N_el` to the total energy, following QE convention. This was verified in previous debugging. The value itself (5.17 eV for both Si and Fe) needs cross-checking.

## Implementation: Diagnostic Steps

### Step 1: Add per-component energy printout at convergence

In `src/scf/mod.rs`, after computing `e_total`, also log the individual components:

```
info!("  E_band  = {:.6} eV", e_band);
info!("  E_H     = {:.6} eV", e_hartree);
info!("  E_xc    = {:.6} eV", e_xc);
info!("  E_vxc   = {:.6} eV", e_vxc);
info!("  E_ewald = {:.6} eV", e_ewald);
info!("  V_G0·N  = {:.6} eV", v_local_g0 * n_el);
```

Compare each against QE's output (which prints exactly these components in Ry).

Files: `src/scf/mod.rs`, `src/scf/energy.rs`

### Step 2: Compare V_local(G) point-by-point against QE

Write a test that computes V_local(G) for the first ~20 G-vectors and compares against QE's values (extractable from `pp.x` with `plot_num=1`).

Files: `tests/fe_debug.rs`, new QE pp.x run

### Step 3: Compare projector form factors F_i(q) against QE

Compute `F_i(|k+G|)` for the first few G-vectors and compare against QE's `init_us_2` output (available in XML with `verbosity='debug'` or by adding print statements to QE source).

Files: `tests/fe_debug.rs`

### Step 4: Verify XC on a known density

Construct a uniform electron gas at a known density (e.g., r_s = 2 Bohr), evaluate our LDA XC, and compare against published values from PZ 1981 Table I.

Files: `src/potential/xc.rs` (new test)

### Step 5: Run QE with `verbosity='debug'` to extract V_eff

This dumps the effective potential at each SCF step, enabling direct comparison of V_H, V_xc, V_local in G-space.

### Step 6: Binary search — V_local only, then V_local+V_NL

Run our code with V_NL disabled (D_ij = 0). If eigenvalues match QE's "kinetic + V_local" component, the bug is in V_NL. If they don't, the bug is in V_local.

This is the most information-efficient diagnostic: one bit tells us which half of the code to focus on.

## Verification

1. Si total energy within 0.1 eV of QE (-231.61 eV)
2. C diamond converges and matches QE within 0.1 eV
3. Fe total energy within 0.1 eV of QE (-3059.46 eV)
4. All eigenvalue degeneracies exact to 1e-6 eV at Gamma
5. `cargo test --release --test qe_validation` passes with all assertions < 0.5 eV

## Estimated Effort

Multi-session investigation. Step 1 (energy printout) and Step 6 (binary search) are the highest-value diagnostics — each takes ~30 minutes and together will localize the bug to V_local or V_NL. Steps 2-5 are follow-up validation once the root cause is identified.
