---
id: HKIN
status: active
priority: low
complexity: trivial
risk: low
depends_on: []
blocks: []
---

# HKIN: Drop `Option<&dyn Fn>` from `build_hamiltonian` and `compute_band_structure`

## Problem

`src/hamiltonian.rs:35` exposes `build_hamiltonian` with an `Option<&dyn Fn(usize, usize) -> Complex64>` parameter that is **never** `Some` at any call site:

| Call site | V_eff arg |
|---|---|
| `src/bandstructure.rs:37` (inside `compute_band_structure`) | forwards caller's `Option` |
| `tests/free_electron_bands.rs:114` | `None` |
| `tests/free_electron_bands.rs:310` | `None` |
| `tests/free_electron_bands.rs:345` | `None` |
| `tests/free_electron_bands.rs:379` | `None` |
| `tests/free_electron_bands.rs:423` | `None` |
| `tests/free_electron_bands.rs:466` | `None` |
| `tests/free_electron_bands.rs:508` | `None` |
| `tests/free_electron_bands.rs:541` | `None` |

The one caller that forwards the `Option` — `compute_band_structure` — is itself called with `None` at every one of its five call sites:

| `compute_band_structure` call site | V_eff arg |
|---|---|
| `src/main.rs:58` (production, `KPointSettings::BandPath`) | `None` |
| `src/bandstructure.rs:118` (unit test) | `None` |
| `src/bandstructure.rs:151` (unit test) | `None` |
| `tests/free_electron_bands.rs:273` | `None` |

So: in the entire repo, **zero** call sites exercise the `Some` branch of either function. The SCF path does not use `build_hamiltonian` at all — it uses the separate `build_hamiltonian_with_v_eff` at `src/scf/potentials.rs:140`, which takes `v_eff` as a non-optional `&[Complex64]` on the FFT grid.

The `Option` is pure dead flexibility. CLAUDE.md is explicit:

> Don't add error handling, fallbacks, or validation for scenarios that can't happen. Trust internal code and framework guarantees.
> Don't design for hypothetical future requirements.

Every production and test caller is running a **free-electron (kinetic-only) band structure** — a validation/diagnostic mode, not a physics calculation. The analytic nearly-free-electron spectrum is the reference: `tests/free_electron_bands.rs` asserts eigenvalues of `build_hamiltonian(&basis, &k, None)` match `analytic_eigenvalues(&basis, &k, N_BANDS)` for Si, C diamond, and BCC Fe. The `None` branch is the point of those tests; the `Some` branch is vestigial.

If and when we later add a post-SCF band-structure workflow (read converged V_eff, compute bands along a path — QE's `bands.x` model), that feature will need a non-optional `v_eff` parameter *and* a different V_eff representation (FFT grid, not closure), so the current `Option<&dyn Fn>` is not a reusable hook for it anyway.

## Implementation

Three mechanical changes. No physics, no algorithms.

### 1. Delete `build_hamiltonian`; keep `build_kinetic`

`build_kinetic` in `src/hamiltonian.rs:12` is already `pub` and does the kinetic-only work. Drop `build_hamiltonian` entirely.

```rust
// src/hamiltonian.rs — delete lines 25-52 (the `build_hamiltonian` fn + its docstring).
// Keep `build_kinetic` as the only public constructor in this module.
```

### 2. Drop `v_eff: Option<...>` from `compute_band_structure`

```rust
// src/bandstructure.rs:26-51 — before:
pub fn compute_band_structure(
    basis: &BasisSet,
    kpoints: &[KPoint],
    distances: &[f64],
    n_bands: usize,
    v_eff: Option<&dyn Fn(usize, usize) -> num_complex::Complex64>,
) -> Result<BandStructure> {
    ...
    let h = hamiltonian::build_hamiltonian(basis, &kp.k, v_eff);
    ...
}

// after:
pub fn compute_band_structure(
    basis: &BasisSet,
    kpoints: &[KPoint],
    distances: &[f64],
    n_bands: usize,
) -> Result<BandStructure> {
    ...
    let h = hamiltonian::build_kinetic(basis, &kp.k);
    ...
}
```

Update the docstring to state plainly that this routine computes the **free-electron** band structure (kinetic only). Anything beyond free-electron goes through the SCF driver.

### 3. Fix the 4 call-site types

| File | Change |
|---|---|
| `src/main.rs:58` | drop the trailing `, None` |
| `src/bandstructure.rs:118` | drop the trailing `, None` |
| `src/bandstructure.rs:151` | drop the trailing `, None` |
| `tests/free_electron_bands.rs:273` | drop the trailing `, None` |

And for the 8 direct `build_hamiltonian(..., None)` sites in `tests/free_electron_bands.rs` (lines 114, 310, 345, 379, 423, 466, 508, 541) replace with `build_kinetic(&basis, &k)`.

### Optional — rename `compute_band_structure`

Consider `compute_free_electron_bands` to make the mode explicit at call sites. Not strictly required for this proposal; flag for the Engineering Manager during review. If renamed, update `src/main.rs:45`'s `KPointSettings::BandPath` arm to call the new name.

## Verification

Mechanical refactor — no behavior change. Gate:

1. `cargo test` — all free-electron band tests still pass with identical eigenvalues (the `build_hamiltonian(..., None)` path is definitionally equivalent to `build_kinetic` because the `Option` branch is inert when `None`).
2. `cargo clippy -q --all-targets` and `cargo clippy -q --all-targets --features gpu` — clean.
3. `cargo doc --no-deps -- -D warnings` — clean; the docstring on `build_hamiltonian` referenced `None`/`Some` which no longer exist, so the updated `compute_band_structure` docstring must not leave dangling intra-doc links.
4. Run the free-electron band example end-to-end: `cargo run --release -- --input examples/si_free_electron.yaml -o /tmp/bands.tsv` — output should be byte-identical to a pre-HKIN run.

No QE validation needed — this touches zero physics.
