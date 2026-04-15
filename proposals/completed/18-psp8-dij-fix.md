# Proposal 18: Fix PSP8 D_ij Matrix and Audit UPF rho_atom Units

## Problem

### PSP8 D_ij is always zero

`src/pseudopotential/psp8.rs`, line 144:

```rust
// TODO: read ekb values properly for off-diagonal terms
let dij = vec![0.0; n_projectors * n_projectors];
```

The D_ij matrix couples projectors in the Kleinman-Bylander non-local potential: `V_NL = Σ |β_i⟩ D_{ij} ⟨β_j|`. With `D_ij = 0`, the non-local contribution is identically zero. Any calculation using PSP8 format pseudopotentials produces wrong results — the non-local energy is missing entirely.

PSP8 format encodes the KB energies (`ekb`) in the header line of each projector block. For norm-conserving pseudopotentials, D_ij is diagonal: `D_ii = ekb_i`. The `ekb` values are already parsed at line 85 but discarded.

### UPF rho_atom unit conversion is questionable

`src/pseudopotential/upf.rs`, lines 81-94:

```rust
// UPF stores 4π r² ρ(r) in e/Bohr units on the radial grid.
// ...
rho_raw.iter().zip(r_grid_bohr.iter())
    .map(|(&rho, &_r)| rho / BOHR_TO_ANG.powi(3))
    .collect()
```

The comment says UPF stores `4πr²ρ(r)` in `e/Bohr`, and the code divides by `Bohr³` to convert to `e/ų`. But `4πr²ρ(r)` has units of `e/Bohr` (charge per length), not `e/Bohr³` (charge per volume). Dividing by `Bohr³` gives `e/(Bohr × Å³)`, which is dimensionally inconsistent.

The downstream consumer (`src/scf/initial_density.rs:158-178`) uses the Bessel transform:

```
ρ_atom(G) = (1/Ω) × ∫ [4πr²ρ(r)] × j₀(|G|r) × dr
```

For this integral, `4πr²ρ(r)` must be in `e/length` and `dr` in the same length unit as `r` and `1/|G|`. If `r_grid` was converted to Å in UPF parsing, then `4πr²ρ(r)` should be in `e/Å`, not divided by `Bohr³`.

## References

- PSP8 format specification: ABINIT documentation, `dev_psp8_format.txt`
- UPF format specification: QE `upftools/` documentation
- QE source: `init_us_1.f90` (projector initialization), `atomic_rho.f90` (SAD density)

## Implementation

### Step 1: Fix PSP8 D_ij

In `src/pseudopotential/psp8.rs`, store `ekb` values during projector parsing and build a diagonal D_ij:

```rust
// During projector block parsing (around line 85):
let mut ekb_values: Vec<f64> = Vec::new();

for i_proj in 0..nproj_l {
    // The header line for each projector contains: l, i, ekb
    let header: Vec<f64> = lines[line_idx].split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    let ekb = header[2]; // KB energy in Hartree
    ekb_values.push(ekb * RY_TO_EV / 2.0); // Convert Hartree → eV (or Ry→eV if already Ry)
    // ... parse radial data ...
}

// Build diagonal D_ij from ekb values
let mut dij = vec![0.0; n_projectors * n_projectors];
for (i, &ekb) in ekb_values.iter().enumerate() {
    dij[i * n_projectors + i] = ekb;
}
```

The unit conversion depends on whether `ekb` is stored in Hartree or Ry in PSP8 format — verify against ABINIT's `m_psp8.F90`.

### Step 2: Audit and fix UPF rho_atom

Trace the full unit chain from UPF file through to the Bessel transform:

1. UPF stores `4πr²ρ(r)` on a radial grid `r` in Bohr. Units: `e/Bohr`.
2. `r_grid` is converted from Bohr to Å in UPF parsing (line 42): `r * BOHR_TO_ANG`.
3. `rab` (dr) is converted from Bohr to Å (line 46): `rab * BOHR_TO_ANG`.
4. For the Bessel transform `∫ f(r) j₀(Gr) dr` with `r` in Å and `G` in 1/Å, `f(r)` must be in `e/Å`.

So the conversion should be:
```rust
// 4πr²ρ(r) [e/Bohr] → [e/Å]
rho_raw.iter().map(|&rho| rho / BOHR_TO_ANG).collect()
```

Not `/ BOHR_TO_ANG.powi(3)`. Verify by checking that the integrated density `∫ 4πr²ρ(r) dr` over the full grid gives `z_valence`.

### Step 3: Add unit test for density integral

```rust
#[test]
fn test_rho_atom_integrates_to_z_valence() {
    let pp = load_test_pp("Si.upf");
    let integral: f64 = pp.rho_atom.iter()
        .zip(pp.rab.iter())
        .map(|(&rho, &dr)| rho * dr)
        .sum();
    // Should integrate to z_valence (4 for Si with 4 valence electrons)
    assert!((integral - pp.z_valence).abs() < 0.01,
        "rho_atom integral = {integral}, expected {}", pp.z_valence);
}
```

### Step 4: Add test for PSP8 non-local energy

```rust
#[test]
fn test_psp8_dij_nonzero() {
    let pp = load_test_pp("Si.psp8");
    let dij_max: f64 = pp.dij.iter().map(|d| d.abs()).fold(0.0, f64::max);
    assert!(dij_max > 0.01, "PSP8 D_ij should not be all zeros: max = {dij_max}");
}
```

## Acceptance Criteria

1. **PSP8 D_ij is populated:** Loading a PSP8 file produces non-zero diagonal D_ij entries matching the `ekb` values in the file header.
2. **UPF rho_atom integrates correctly:** `∫ rho_atom × dr` gives `z_valence ± 0.01` for Si, C, and Fe test pseudopotentials.
3. **SCF with PSP8 converges:** A PSP8-based Si SCF gives the same total energy as UPF-based Si SCF within 0.01 eV.
4. **SAD initial guess quality:** With corrected rho_atom units, the initial density is closer to self-consistency (fewer SCF iterations to converge compared to uniform initial density).
5. **No regression on UPF:** Existing UPF test cases produce unchanged results (if the rho_atom fix changes SAD but not final converged energy, document the intermediate behavior change).
