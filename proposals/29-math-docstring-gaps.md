# Proposal 29: Math Documentation and Constant Hygiene

## Problem

A full math and unit audit found all formulas correct and units consistent (eV + Å throughout). However, four functions lack docstrings explaining the math they implement, one comment is misleading, and the Coulomb constant `E2` is duplicated in four places.

### Documentation gaps

| Location | Function | Missing |
|----------|----------|---------|
| `src/crystal.rs:51` | `reciprocal()` | No docstring at all. Should state b_i = 2π(a_j × a_k)/Ω |
| `src/scf/mixing.rs:147-148` | Anderson coefficient solve | Doesn't explain the constraint-embedding technique: solves for m-1 coefficients, derives the last as α_last = 1 - Σα_prev to enforce Σc_i = 1 |
| `src/scf/mixing.rs:187-189` | `auto_q_tf_squared()` | Formula is documented but not cited. Should reference Ashcroft & Mermin or note it derives from the free-electron Thomas-Fermi screening length |
| `src/scf/energy.rs:75-83` | `total_energy()` | Docstring says `E_band - E_H + (E_xc - E_vxc) + E_ewald` but omits the `V_local(G=0) × N_el` correction added at both call sites (mod.rs:283, mod.rs:549) |

### Duplicate Coulomb constant

`E2 = 14.399645351950548` (eV·Å) is defined independently in four places:

| Location | Declaration |
|----------|-------------|
| `src/consts.rs:13` | `pub const E2_COULOMB: f64 = 14.399645351950548;` |
| `src/potential/hartree.rs:11` | `pub const E2: f64 = 14.399645351950548;` |
| `src/pseudopotential/mod.rs` | `const E2: f64 = 14.399645351950548;` |
| `src/pseudopotential/upf.rs` | `let e2 = 14.399645351950548;` |

All four have the same value, so there is no bug today. But if one is ever updated without the others, it would silently break energy calculations.

### Misleading XC conversion comment

`src/potential/xc.rs:85` — the comment describing the density unit conversion reads backwards. The code is correct (multiplies by `BOHR3_TO_ANG3 = 0.529177³ ≈ 0.149` to convert e/ų → e/Bohr³), but the comment phrasing is confusing.

## Implementation

### 1. `crystal.rs:51` — Add docstring to `reciprocal()`

```rust
/// Reciprocal lattice vectors via b_i = 2π(a_j × a_k) / Ω.
///
/// Uses the convention where plane waves are exp(iG·r) with G in
/// units of Å⁻¹, so the 2π factor is included here (not in the
/// plane wave definition).
pub fn reciprocal(&self) -> Self {
```

### 2. `mixing.rs:147-148` — Add comment explaining constraint embedding

```rust
// Solve the reduced DIIS system for m-1 coefficients. The constraint
// Σ c_i = 1 is embedded by eliminating the last coefficient:
//   α_last = 1 - Σ α_prev
// The matrix A[i,j] = ΔR_i · ΔR_j (where ΔR_i = R_i - R_last) and
// b[i] = -ΔR_i · R_last. This is equivalent to minimizing |Σ c_i R_i|²
// subject to Σ c_i = 1.
let alpha_prev = solve_linear_system(&a_mat, &b_vec, mm);
let alpha_last = 1.0 - alpha_prev.iter().sum::<f64>();
```

### 3. `mixing.rs:187-189` — Add citation to `auto_q_tf_squared()`

```rust
/// Auto-estimate Thomas-Fermi screening wavevector squared from average density.
///
/// q_TF² = 4(3π²ρ)^{1/3} / π  (in a.u., then convert from Bohr⁻² to Å⁻²)
///
/// This is the free-electron Thomas-Fermi screening length from
/// Ashcroft & Mermin, Solid State Physics (1976), Ch. 17.
fn auto_q_tf_squared(n_electrons: f64, omega: f64) -> f64 {
```

### 4. `energy.rs:75` — Document the full energy formula

```rust
/// Total Kohn-Sham energy (without V_local(G=0) correction).
///
///   E = E_band - E_H + (E_xc - E_vxc) + E_ewald
///
/// Callers must add `V_local(G=0) × N_el` to account for the G=0
/// local pseudopotential term that is zeroed in the Hamiltonian
/// (see `run_scf` and `run_scf_spin` in mod.rs).
pub(crate) fn total_energy(
```

### 5. Centralize `E2` constant

Remove the local definitions in `hartree.rs`, `pseudopotential/mod.rs`, and `upf.rs`. Import from `consts.rs` instead:

```rust
// In each file, replace local E2 with:
use crate::consts::E2_COULOMB;
```

Rename uses from `E2` to `E2_COULOMB`, or add a local alias if brevity matters:

```rust
use crate::consts::E2_COULOMB as E2;
```

### 6. `xc.rs:85` — Fix misleading conversion comment

```rust
// Before:
// rho [e/ų] → rho [e/Bohr³] = rho × Bohr_to_Å³ = rho × 0.529177³

// After:
// rho [e/ų] → rho [e/Bohr³]: multiply by (Bohr/Å)³ = BOHR_TO_ANG³ ≈ 0.149
// (1 ų holds more Bohr³, so density in Bohr⁻³ is smaller)
```

## Verification

```bash
cargo test            # constant centralization doesn't change values
cargo clippy -q       # no new warnings from imports
cargo doc --no-deps   # docstrings render correctly
```

## Estimated Effort

Docstrings: trivial. E2 centralization: straightforward find-and-replace across 3 files, verify with `cargo test`.
