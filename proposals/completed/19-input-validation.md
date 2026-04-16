# Proposal 19: Input Validation and Defensive Checks

> **Note:** Line numbers reference the pre-ScfContext codebase (src/scf/mod.rs was ~1127 lines, now ~709). Verify locations before implementing.

## Problem

The codebase accepts invalid physics parameters without error, leading to silent wrong results or panics deep in computation. There is no validation at the boundary between user input and internal computation.

Specific gaps:

| Parameter | Current behavior | Failure mode |
|-----------|-----------------|--------------|
| `ecut = 0` or negative | Empty basis set | Silent garbage results |
| `conv_threshold = 0` | Never converges | Infinite loop until `max_iter` |
| `mixing_beta > 1` or `≤ 0` | Divergent SCF or NaN | Silent wrong results or panic |
| `n_bands = 0` | Empty eigenvalue arrays | Panic in occupations |
| `smearing_sigma = 0` | Division by zero in F-D | Panic or NaN |
| Lattice with zero volume | Division by zero in reciprocal | Panic |
| Empty atoms list | Empty density | Silent wrong results |
| K-point grid dims = 0 | No k-points generated | Panic |
| K-point weights ≠ 1.0 | Wrong occupations/energy | Silent wrong results |
| Unknown element symbol | `panic!()` at `input.rs:125` | Crash |
| `ecutrho_ratio < 1` | FFT grid smaller than basis | Aliasing errors |
| `Lattice::volume()` negative | Negative `omega` in normalization | Wrong sign in energies |

Additionally, the Anderson mixer silently returns uniform coefficients when its matrix is singular (`mixing.rs:217-218`), masking convergence problems without any diagnostic output.

## Implementation

### Step 1: Validate `ScfParams` at construction

Add a validation method called at the start of `run_scf()`:

```rust
// src/scf/mod.rs
impl ScfParams {
    pub fn validate(&self) -> crate::error::Result<()> {
        if self.n_bands == 0 {
            return Err(PwdftError::InvalidInput("n_bands must be > 0".into()));
        }
        if self.conv_threshold <= 0.0 {
            return Err(PwdftError::InvalidInput("conv_threshold must be positive".into()));
        }
        if self.mixing_beta <= 0.0 || self.mixing_beta > 1.0 {
            return Err(PwdftError::InvalidInput(
                format!("mixing_beta must be in (0, 1], got {}", self.mixing_beta)
            ));
        }
        if self.smearing_sigma < 0.0 {
            return Err(PwdftError::InvalidInput("smearing_sigma must be non-negative".into()));
        }
        if self.ecutrho_ratio < 1.0 {
            return Err(PwdftError::InvalidInput(
                format!("ecutrho_ratio must be >= 1.0, got {}", self.ecutrho_ratio)
            ));
        }
        Ok(())
    }
}
```

### Step 2: Validate crystal and basis at `run_scf()` entry

```rust
pub fn run_scf(...) -> Result<ScfResult> {
    params.validate()?;

    let omega = crystal.lattice.volume().abs();
    if omega < 1e-10 {
        return Err(PwdftError::InvalidInput("lattice has zero or near-zero volume".into()));
    }
    if crystal.atoms.is_empty() {
        return Err(PwdftError::InvalidInput("at least one atom is required".into()));
    }
    if kpoints.is_empty() {
        return Err(PwdftError::InvalidInput("at least one k-point is required".into()));
    }

    let weight_sum: f64 = kpoints.iter().map(|kp| kp.weight).sum();
    if (weight_sum - 1.0).abs() > 1e-6 {
        return Err(PwdftError::InvalidInput(
            format!("k-point weights sum to {weight_sum:.6}, expected 1.0")
        ));
    }

    // ... existing code ...
}
```

### Step 3: Fix `Lattice::volume()` to return absolute value

In `src/crystal.rs`, line 42:

```rust
pub fn volume(&self) -> f64 {
    self.a.cross(&self.b).dot(&self.c).abs()
}
```

This ensures positive `omega` regardless of lattice vector handedness. The `reciprocal()` method (line 46) computes `2π / volume()`, which requires a signed triple product to get the correct reciprocal vector directions. So either:
- Keep `volume()` signed and add a separate `volume_abs()`, or
- Use the signed triple product directly inside `reciprocal()`:

```rust
pub fn volume(&self) -> f64 {
    self.a.cross(&self.b).dot(&self.c).abs()
}

pub fn reciprocal(&self) -> Self {
    let triple = self.a.cross(&self.b).dot(&self.c); // signed
    let factor = 2.0 * PI / triple;
    Self {
        a: self.b.cross(&self.c) * factor,
        b: self.c.cross(&self.a) * factor,
        c: self.a.cross(&self.b) * factor,
    }
}
```

### Step 4: Replace `panic!` in input parsing with errors

In `src/input.rs`, line 124-125:

```rust
// Before:
let elem = crate::atoms::Element::from_symbol(&ai.symbol)
    .unwrap_or_else(|| panic!("unknown element: {}", ai.symbol));

// After:
let elem = crate::atoms::Element::from_symbol(&ai.symbol)
    .ok_or_else(|| PwdftError::InvalidInput(format!("unknown element: {}", ai.symbol)))?;
```

Same pattern in `src/settings.rs:343`.

### Step 5: Add warning on Anderson mixer singular fallback

In `src/scf/mixing.rs`, line 217:

```rust
if pivot.abs() < 1e-15 {
    log::warn!(
        "Anderson mixer: singular overlap matrix (pivot={:.2e}), \
         falling back to uniform coefficients. This may indicate \
         linearly dependent mixing history.",
        pivot
    );
    return vec![1.0 / (n + 1) as f64; n];
}
```

### Step 6: Add warning on near-zero density integral

In `src/scf/density.rs` and `src/scf/initial_density.rs`, where density normalization is skipped:

```rust
if integral.abs() > 1e-15 {
    let scale = n_electrons / integral;
    for v in &mut rho_r { *v *= scale; }
} else {
    log::warn!(
        "Density integral near zero ({:.2e}), skipping normalization. \
         This indicates a problem with the wavefunction or initial density.",
        integral
    );
}
```

## Acceptance Criteria

1. **Invalid parameters rejected:** Each invalid parameter in the table above produces a descriptive `PwdftError::InvalidInput` error, not a panic.
2. **Existing valid inputs unchanged:** All current test cases and examples pass without modification.
3. **Diagnostics present:** Anderson mixer singular fallback and density normalization skip produce `warn!` log messages.
4. **Left-handed lattice:** A crystal with left-handed lattice vectors produces the same energy as right-handed (volume is always positive).
5. **Error messages actionable:** Each error message names the parameter, its value, and the valid range.
