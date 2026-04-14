# Proposal 11: Entropy and Free Energy (Mermin Functional)

## Problem

For any system with finite smearing (sigma > 0), the DFT total energy depends on sigma with an O(sigma^2) error:

```
E(sigma) = E(0) + gamma * sigma^2 + O(sigma^3)
```

The physically meaningful quantities are:

- **F = E - T*S** — the Helmholtz free energy (Mermin functional). Forces and stresses are derivatives of F, not E.
- **E_0 = (E + F) / 2** — the sigma->0 extrapolated energy. Eliminates the leading O(sigma^2) error for Fermi-Dirac and Gaussian smearing.

The current code (`src/scf/mod.rs`, lines 465-518) computes only `E_total = E_band - E_H + E_xc - E_vxc + E_ewald` with no entropy contribution. For metals with sigma = 0.1-0.2 eV, the entropy term T*S can be 10-100 meV/atom — larger than typical energy differences of interest.

VASP reports three energies (TOTEN=F, energy without entropy=E, sigma->0 energy=E_0). QE includes the `-T*S` contribution in its total energy. Without this, energy differences between metallic phases are unreliable.

## References

- Mermin, N.D., Phys. Rev. 137, A1441 (1965) — finite-temperature DFT
- Gillan, M.J., J. Phys.: Condens. Matter 1, 689 (1989) — sigma->0 extrapolation
- De Vita, A. & Gillan, M.J., J. Phys.: Condens. Matter 3, 6225 (1991)
- VASP wiki: [ISMEAR](https://www.vasp.at/wiki/index.php/ISMEAR) — three energies explained
- VASP wiki: [smearing technique](https://www.vasp.at/wiki/index.php/Smearing_technique)

## Implementation

### Step 1: Compute entropy for Fermi-Dirac smearing

Add to `src/scf/smearing.rs`:

```rust
/// Compute the electronic entropy for Fermi-Dirac smearing.
///
/// S = -k_B * sum_{n,k} w_k * [f * ln(f) + (1-f) * ln(1-f)]
///
/// where f includes the spin factor of 2 (so f/2 is the per-spin occupation).
/// Returns entropy in eV (i.e., k_B * S in natural units where k_B = 1
/// and sigma is in eV).
pub fn fermi_dirac_entropy(
    eigenvalues: &[Vec<f64>],
    kpoints: &[KPoint],
    fermi_energy: f64,
    sigma: f64,
) -> f64 {
    if sigma < 1e-15 {
        return 0.0;
    }

    let mut entropy = 0.0;
    for (evs, kp) in eigenvalues.iter().zip(kpoints.iter()) {
        for &e in evs {
            let f = fermi_dirac(e, fermi_energy, sigma);
            // f is [0, 2] with spin factor; per-spin occupation is f/2
            let f1 = (f / 2.0).clamp(1e-30, 1.0 - 1e-30);
            entropy -= kp.weight * 2.0 * (f1 * f1.ln() + (1.0 - f1) * (1.0 - f1).ln());
        }
    }

    // entropy is dimensionless; multiply by sigma to get T*S in eV
    // (since x = (e - E_F) / sigma, the natural "temperature" is sigma)
    entropy * sigma
}
```

Note: The factor of 2 outside accounts for spin degeneracy. The `clamp` prevents log(0).

### Step 2: Extend ScfResult

In `src/scf/mod.rs`:

```rust
pub struct ScfResult {
    pub total_energy: f64,       // E (Kohn-Sham energy, no entropy)
    pub free_energy: f64,        // F = E - T*S (Mermin free energy)
    pub energy_sigma0: f64,      // E_0 = (E + F) / 2 (sigma->0 extrapolated)
    pub entropy_ts: f64,         // T*S (entropy contribution in eV)
    pub eigenvalues: Vec<Vec<f64>>,
    pub fermi_energy: f64,
    pub n_iterations: usize,
    pub rho_g: Vec<Complex64>,
}
```

### Step 3: Compute and report

After convergence in `run_scf`, compute entropy and the three energies:

```rust
let entropy_ts = smearing::fermi_dirac_entropy(
    &eigenvalues_all, kpoints, fermi_energy, params.smearing_sigma,
);
let free_energy = total_energy - entropy_ts;
let energy_sigma0 = (total_energy + free_energy) / 2.0;

info!("Energy without entropy (E):  {total_energy:.6} eV");
info!("Free energy (F = E - TS):    {free_energy:.6} eV");
info!("Energy sigma->0 (E_0):       {energy_sigma0:.6} eV");
info!("Entropy contribution (-TS):  {:.6} eV ({:.3} meV/atom)",
    -entropy_ts, -entropy_ts * 1000.0 / crystal.atoms.len() as f64);

return Ok(ScfResult {
    total_energy,
    free_energy,
    energy_sigma0,
    entropy_ts,
    eigenvalues: eigenvalues_all,
    fermi_energy,
    n_iterations: iter + 1,
    rho_g: rho_g_basis,
});
```

### Step 4: Use free energy for convergence monitoring

When Proposal 10 (energy convergence) is implemented, use the free energy F for the energy convergence check rather than E, since F is the variational quantity (it is stationary at self-consistency):

```rust
let energy_for_convergence = free_energy;  // not total_energy
```

## Acceptance Criteria

1. **Entropy is computed** via `fermi_dirac_entropy` and included in `ScfResult`.
2. **Three energies are reported:** E (KS), F (free), E_0 (sigma->0 extrapolated).
3. **Entropy magnitude test:** For an insulator (Si) with sigma = 0.01 eV, T*S should be < 1e-6 eV. For a metal (Al) with sigma = 0.1 eV, T*S should be on the order of 10-100 meV.
4. **Sigma extrapolation test:** Run the same system at sigma = 0.05, 0.1, 0.2 eV. E varies quadratically with sigma, but E_0 should be approximately constant (spread < O(sigma^3)).
5. **Zero-temperature limit:** When sigma < 1e-15, entropy is identically 0 and all three energies are equal.
6. **QE comparison:** For a metallic system, compare F and T*S against QE's reported values.
