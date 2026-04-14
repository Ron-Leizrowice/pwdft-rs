# Proposal 12: Gaussian and Methfessel-Paxton Smearing

## Problem

The code (`src/scf/smearing.rs`) implements only Fermi-Dirac smearing. While physically motivated, Fermi-Dirac has the worst sigma-dependence of the common schemes — the total energy error is O(sigma^2) with a large prefactor, requiring small sigma (and hence dense k-meshes) for accurate energies.

Production codes offer multiple schemes because different systems need different trade-offs:

| Scheme | Energy error | Occupations | Best for |
|--------|-------------|-------------|----------|
| Fermi-Dirac | O(sigma^2) | Always [0,2], monotonic | Physical temperature effects |
| Gaussian | O(sigma^2) | Always [0,2], monotonic | Safe default, semiconductors |
| Methfessel-Paxton 1 | O(sigma^4) | Can be negative | Metals (forces, relaxation) |
| Marzari-Vanderbilt | O(sigma^2)* | Always positive | Metals (safe M-P alternative) |

*M-V has the same O(sigma^2) scaling as Gaussian but with a ~10x smaller prefactor due to error cancellation.

M-P order 1 is the default in VASP (`ISMEAR=1`) and the recommended scheme for metallic relaxations. It allows using sigma = 0.1-0.2 eV (instead of 0.01 eV for F-D) while maintaining the same energy accuracy, reducing the required k-mesh density significantly.

## References

- Methfessel, M. & Paxton, A.T., Phys. Rev. B 40, 3616 (1989)
- Marzari, N. et al., Phys. Rev. Lett. 82, 3296 (1999) — cold smearing
- VASP wiki: [ISMEAR](https://www.vasp.at/wiki/index.php/ISMEAR)
- QE docs: `smearing` parameter in [pw.x input](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
- [arXiv:2212.07988](https://arxiv.org/abs/2212.07988) — Fermi energy determination for advanced smearing

## Implementation

### Step 1: Smearing type enum

In `src/scf/smearing.rs`:

```rust
/// Available smearing schemes for occupation numbers.
#[derive(Debug, Clone, Copy, Default, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SmearingType {
    /// Fermi-Dirac: f(x) = 1 / (1 + exp(x)). Physical temperature.
    #[default]
    FermiDirac,
    /// Gaussian: f(x) = (1/2) erfc(x). Safe general-purpose default.
    Gaussian,
    /// Methfessel-Paxton order 1: O(sigma^4) energy error. Best for metals.
    /// Warning: occupations can be negative.
    MethfesselPaxton,
    /// Marzari-Vanderbilt cold smearing: O(sigma^2) with small prefactor.
    /// Always-positive occupations. Good M-P alternative.
    ColdSmearing,
}
```

### Step 2: Occupation functions

```rust
/// Gaussian smearing: f(x) = (1/2) erfc(x)
/// Complementary error function gives smooth step.
pub fn gaussian_occupation(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 { return zero_temp(energy, fermi_energy); }
    let x = (energy - fermi_energy) / sigma;
    // erfc(x) = 1 - erf(x); using the approximation from ewald.rs or puruspe
    2.0 * 0.5 * erfc(x)  // spin factor of 2, erfc goes from 2 to 0
}

/// Methfessel-Paxton order 1.
/// delta_1(x) = (1/sqrt(pi)) * (2 - sqrt(2)*x) * ... wait, that's M-V.
/// M-P order 1: delta_1(x) = delta_0(x) + A_1 * H_2(x) * exp(-x^2)
///   where delta_0(x) = (1/sqrt(pi)) * exp(-x^2)
///   A_1 = -1/4sqrt(pi), H_2(x) = 4x^2 - 2
/// Step function: f_1(x) = f_0(x) + A_1 * H_1(x) * exp(-x^2)
///   where f_0(x) = (1/2)erfc(x), H_1(x) = 2x
pub fn methfessel_paxton_occupation(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 { return zero_temp(energy, fermi_energy); }
    let x = (energy - fermi_energy) / sigma;
    let gauss = (-x * x).exp() / PI.sqrt();
    let f0 = erfc(x) / 2.0;
    // Order-1 correction: A_1 * H_1(x) * exp(-x^2)
    // A_1 = -1/(4*sqrt(pi)), H_1(x) = 2x
    let f1 = f0 - 0.5 * x * gauss;
    2.0 * f1  // spin factor
}

/// Marzari-Vanderbilt cold smearing.
/// delta(x) = (1/sqrt(pi)) * (2 - sqrt(2)*x) * exp(-(x - 1/sqrt(2))^2)
/// f(x) = (1/2) [sqrt(2)*erfc(x_shifted) + ...] — integrated form.
/// Simpler: f(x) = (1/2) + erf(x+1/sqrt(2))/2 + (1/sqrt(2*pi)) * exp(-(x+1/sqrt(2))^2)
pub fn cold_smearing_occupation(energy: f64, fermi_energy: f64, sigma: f64) -> f64 {
    if sigma < 1e-15 { return zero_temp(energy, fermi_energy); }
    let x = (energy - fermi_energy) / sigma;
    let sq2_inv = 1.0 / 2.0_f64.sqrt();
    let arg = x + sq2_inv;
    let f = 0.5 * (1.0 + erf(arg)) + (1.0 / (2.0 * PI).sqrt()) * (-arg * arg).exp();
    // This goes from ~1 (x << 0) to ~0 (x >> 0), but it's (1-f) relative to standard convention
    2.0 * (1.0 - f)  // spin factor, flip convention
}

fn zero_temp(energy: f64, fermi_energy: f64) -> f64 {
    if energy < fermi_energy { 2.0 }
    else if (energy - fermi_energy).abs() < 1e-12 { 1.0 }
    else { 0.0 }
}
```

### Step 3: Unified dispatch

```rust
/// Compute occupation for a given smearing scheme.
pub fn occupation(
    smearing: SmearingType,
    energy: f64,
    fermi_energy: f64,
    sigma: f64,
) -> f64 {
    match smearing {
        SmearingType::FermiDirac => fermi_dirac(energy, fermi_energy, sigma),
        SmearingType::Gaussian => gaussian_occupation(energy, fermi_energy, sigma),
        SmearingType::MethfesselPaxton => methfessel_paxton_occupation(energy, fermi_energy, sigma),
        SmearingType::ColdSmearing => cold_smearing_occupation(energy, fermi_energy, sigma),
    }
}
```

### Step 4: Entropy for each scheme

Each smearing scheme has a different entropy formula. Add `entropy_weight` functions:

```rust
/// Entropy contribution per state for Gaussian smearing.
/// S_gauss = (1/sqrt(pi)) * exp(-x^2) * sigma
pub fn gaussian_entropy_per_state(x: f64) -> f64 {
    (-x * x).exp() / PI.sqrt()
}

/// Entropy for M-P order 1 (from derivative of free energy).
/// S_mp1 = (1/sqrt(pi)) * (1/2 - x^2) * exp(-x^2) * sigma
pub fn mp1_entropy_per_state(x: f64) -> f64 {
    (0.5 - x * x) * (-x * x).exp() / PI.sqrt()
}
```

Total entropy: `T*S = sigma * sum_{n,k} w_k * s(x_{n,k})` where `x = (e - E_F) / sigma`.

### Step 5: Input configuration

In `src/input.rs`:

```rust
/// Smearing scheme: "fermi_dirac" | "gaussian" | "methfessel_paxton" | "cold_smearing"
#[serde(default)]
pub smearing_type: SmearingType,
```

```toml
[scf]
smearing_type = "methfessel_paxton"
smearing_sigma = 0.1
```

### Step 6: Update `find_fermi_energy` and callers

`find_fermi_energy` (line 36) and the occupation computation in `run_scf` (lines 265-272) need to accept `SmearingType` and dispatch accordingly. The bisection algorithm in `find_fermi_energy` is generic — only the occupation function changes.

**M-P Fermi energy caveat:** Methfessel-Paxton occupations are non-monotonic, which means the electron count N(E_F) is not monotonic either. Bisection may find the wrong root. For M-P, the safest approach is to start from the Gaussian Fermi energy as initial guess, then refine with a few Newton-Raphson steps using the M-P delta function as the derivative. Alternatively, scan all roots and choose the one closest to the Gaussian E_F.

## Acceptance Criteria

1. **All four smearing schemes produce correct occupations:** Test with known eigenvalues that fully occupied states get f~2, empty states get f~0, and states near E_F get the correct intermediate values for each scheme.
2. **M-P negative occupations:** Verify that M-P can produce f < 0 for states well above E_F (this is expected, not a bug). Add a warning log when this occurs.
3. **Electron count conservation:** `find_fermi_energy` correctly integrates to N_electrons for all four schemes, including the non-monotonic M-P case.
4. **Entropy consistency:** For each scheme, verify that T*S -> 0 as sigma -> 0.
5. **Sigma extrapolation:** Run the same metallic system at sigma = 0.05, 0.1, 0.2 eV with each scheme. M-P E_0 should vary by < O(sigma^4); Gaussian/F-D E_0 should vary by < O(sigma^3).
6. **Input parsing:** TOML `smearing_type` field correctly selects the scheme. Default is `fermi_dirac` for backward compatibility.
