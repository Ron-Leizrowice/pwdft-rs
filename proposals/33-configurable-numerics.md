# Proposal 33: Expose Hardcoded Numerics as Configurable Settings

## Problem

An audit of the codebase found ~60 hardcoded numeric values controlling algorithm behavior — convergence tolerances, cutoff radii, density floors, iteration limits, buffer sizes — scattered across 13 source files. These are all reasonable defaults, but an advanced researcher has no way to override them without editing source code. This blocks:

- **Sensitivity studies:** Can't sweep rho_floor or Ewald cutoff to assess numerical stability
- **Debugging:** Can't relax/tighten thresholds to isolate convergence issues
- **System-specific tuning:** Metals vs insulators, large vs small cells, extreme smearing regimes all benefit from different numerics
- **Hardware tuning:** GPU workgroup size is hardware-dependent (M2 vs A100)

Several values are also **inconsistent** — `potential/xc.rs` uses `1e-30` as a density floor in 4 independent locations while `consts.rs` defines `RHO_FLOOR = 1e-20`.

This proposal adds the settings structs and threads them through call sites. It depends on proposal 32 (YAML migration) having added the schema fields to `Settings` in `settings.rs`, but the plumbing work is independent and large enough to warrant its own scope.

## Inventory

### A. Ewald summation — `src/ewald.rs`

| Value | Line | Default | What it controls |
|---|---|---|---|
| eta | 52 | `(N·π/Ω)^{1/3}` | Real/reciprocal space partition |
| cutoff multiplier | 57, 89 | `10.0` | `g_max = 10η`, `r_max = 10/η`; QE uses 4–5x |
| self-interaction threshold | 106 | `1e-10` | Distance below which atom pair is treated as self-interaction |

### B. Fermi energy search — `src/scf/smearing.rs`

| Value | Line | Default | What it controls |
|---|---|---|---|
| bisection bounds factor | 66–67 | `10.0 * σ.max(0.1)` | Eigenvalue range expansion for Fermi search |
| bisection max iterations | 69 | `200` | Iteration cap |
| bisection convergence | 85 | `1e-14` | Fermi energy tolerance (eV) |
| zero-sigma threshold | 114+ | `1e-15` | σ below this → step-function occupations (all 4 schemes) |
| degenerate-state tolerance | 116+ | `1e-12` | \|E - E_F\| below this → occupation = 0.5 (all 4 schemes) |
| entropy overflow cutoff | 215 | `30.0` | \|x\| > this saturates FD entropy |
| entropy occupation floor | 219 | `1e-30` | Minimum occupation for -f·ln(f) entropy |

### C. Density and occupation — `src/consts.rs`, `src/potential/xc.rs`, `src/scf/density.rs`

| Value | File:Line | Default | What it controls |
|---|---|---|---|
| `RHO_FLOOR` | `consts.rs:19` | `1e-20` | Global density floor for XC (e/ų) |
| XC density floor | `potential/xc.rs:27,177,237,270` | `1e-30` | Density below which XC is skipped — **4 independent literals, inconsistent with RHO_FLOOR** |
| `G2_ZERO_THRESHOLD` | `consts.rs:17` | `1e-12` | \|G\|² treated as zero in Hartree/local potential |
| occupation skip | `scf/density.rs:61` | `1e-15` | Occupation below which band excluded from density |

### D. Mixing — `src/scf/mixing.rs`

| Value | Line | Default | What it controls |
|---|---|---|---|
| Kerker q_tf | 69–72 | auto from Thomas-Fermi | Screening vector; `MixingMode::Kerker { q_tf: Option }` exists but isn't in YAML |
| DIIS pivot tolerance | 231 | `1e-15` | Gauss elimination pivot in Anderson/Pulay |

### E. Initial density — `src/scf/initial_density.rs`

| Value | Line | Default | What it controls |
|---|---|---|---|
| Gaussian sigma | 28 | `1.0` Å | Width of Gaussian model charges for SAD |
| normalization tolerance | 103 | `1e-15` | Density integral normalization check |
| small-gr Bessel threshold | 166–167 | `1e-10` | j₀ Taylor expansion cutoff |

### F. Non-local potential — `src/potential/nonlocal.rs`

| Value | Line | Default | What it controls |
|---|---|---|---|
| q-norm threshold | 175–178 | `1e-12` | \|q\| below which cos(θ) defaults to 1.0 |
| projector sum cutoff | 187 | `1e-20` | Projector magnitude below which non-local term is skipped |
| Bessel small-x threshold | 250 | `1e-10` | Spherical Bessel j_l Taylor expansion cutoff |

### G. Symmetry — `src/symmetry/detect.rs`, `src/symmetry/kpoints.rs`

| Value | File:Line | Default | What it controls |
|---|---|---|---|
| search radius factor | `detect.rs:100` | `1.5` | Candidate lattice vector search multiplier |
| k-point rounding tolerance | `kpoints.rs:100` | `1e-6` | k-point grid snapping after symmetry rotation |

### H. GPU — `src/gpu/mod.rs`

| Value | Line | Default | What it controls |
|---|---|---|---|
| workgroup size | 65 | `256` | Compute shader threads per workgroup (64–512) |
| complex buffer pool | 130 | `5` | Pre-allocated complex f32 buffers |
| real buffer pool | 142 | `3` | Pre-allocated real f32 buffers |

### Values NOT exposed (with justification)

| Value | Location | Why hardcoded |
|---|---|---|
| PZ LDA coefficients (γ, β₁, β₂, a, b, c, d) | `potential/xc.rs` | Published physical parameters — changing them produces a different functional |
| Smearing math (M-P coefficient 0.5, exchange 3/4) | `potential/xc.rs`, `smearing.rs` | Algorithm-defined constants from the original papers |
| Zeta clamp bounds (-1, 1) | `potential/xc.rs:275` | Mathematical constraint: spin polarization ∈ [-1, 1] |
| FFT-friendly factors [2, 3, 5] | `fft.rs:120-121` | **Performance:** ndrustfft optimized for these primes |
| Grid Nyquist formula `2*n_max + 1` | `fft.rs:105` | Sampling theorem, not tunable |

## Implementation

The approach is to carry settings structs through the SCF context rather than reading global constants. Each module group below is an independent unit of work.

### Step 1: Ewald

Add `&EwaldSettings` parameter to `ewald_energy()`:

```rust
pub fn ewald_energy(
    crystal: &Crystal,
    pseudopotentials: &[&PseudopotentialData],
    settings: &EwaldSettings,
) -> f64 {
    let eta = settings.eta.unwrap_or_else(|| {
        (crystal.atoms.len() as f64 * PI / omega).powf(1.0 / 3.0)
    });
    let g_max = settings.cutoff_multiplier * eta;
    let r_max = settings.cutoff_multiplier / eta;
    // ... replace 1e-10 with settings.self_interaction_threshold
}
```

Update caller in `scf/energy.rs` (`total_energy`) to pass settings.

**Files:** `src/ewald.rs`, `src/scf/energy.rs`, `src/scf/mod.rs`

### Step 2: Numerics → potentials and energy

Thread `NumericsSettings` fields into XC, Hartree, and non-local code.

**XC (`potential/xc.rs`):** Replace all 4 independent `1e-30` density floor literals with a `rho_floor` parameter. The function signatures gain a `rho_floor: f64` argument:

```rust
pub fn lda_xc(rho: &[f64], rho_floor: f64) -> (Vec<f64>, Vec<f64>) {
    // replace: if rho_i < 1e-30 { ... }
    // with:    if rho_i < rho_floor { ... }
}
```

**Energy (`scf/energy.rs`):** Replace `G2_ZERO_THRESHOLD` constant reads with a `g2_zero_threshold` parameter in `hartree_energy()`, `hartree_on_fft_grid()`, etc.

**Non-local (`potential/nonlocal.rs`):** Pass `projector_sum_cutoff`, `bessel_small_x_threshold`, `q_norm_threshold` into `compute_nonlocal_potential()` or via a small `NonlocalThresholds` struct.

**`consts.rs` cleanup:** Keep `G2_ZERO_THRESHOLD` and `RHO_FLOOR` as constants referenced by `NumericsSettings::default()`, but call sites read from the settings struct.

**Files:** `src/potential/xc.rs`, `src/potential/nonlocal.rs`, `src/scf/energy.rs`, `src/scf/mixing.rs`, `src/consts.rs`

### Step 3: Electrons → smearing and density

Thread Fermi search and occupation parameters through the smearing module.

**Smearing (`scf/smearing.rs`):** `find_fermi_energy()` gains parameters (or a `FermiSearchParams` struct):

```rust
pub struct FermiSearchParams {
    pub bounds_factor: f64,       // default: 10.0
    pub max_iter: usize,          // default: 200
    pub tol: f64,                 // default: 1e-14
    pub zero_sigma: f64,          // default: 1e-15
    pub degenerate_tol: f64,      // default: 1e-12
    pub entropy_cutoff: f64,      // default: 30.0
    pub entropy_floor: f64,       // default: 1e-30
}
```

Each of the 4 occupation functions (`fermi_dirac_occupation`, `gaussian_occupation`, etc.) takes `zero_sigma` and `degenerate_tol` instead of hardcoded literals.

**Density (`scf/density.rs`):** Pass `occupation_threshold` into `construct_density()`.

**Mixing (`scf/mixing.rs`):** Pass `diis_pivot_tolerance` into the Anderson solver.

**Files:** `src/scf/smearing.rs`, `src/scf/density.rs`, `src/scf/mixing.rs`

### Step 4: Initial density

Pass `InitialDensitySettings` into `initial_density()`:

```rust
pub fn initial_density(
    crystal: &Crystal,
    basis: &BasisSet,
    grid: &FftGrid,
    settings: &InitialDensitySettings,
) -> Vec<f64> {
    // use settings.gaussian_sigma instead of DEFAULT_GAUSSIAN_SIGMA
    // use settings.normalization_tol instead of 1e-15
}
```

**Files:** `src/scf/initial_density.rs`

### Step 5: Symmetry

Pass `search_radius_factor` into `SymmetryInfo::from_crystal()` and `kpoint_rounding_tol` into `reduce_kpoints()`.

**Files:** `src/symmetry/detect.rs`, `src/symmetry/kpoints.rs`

### Step 6: GPU

Pass `GpuSettings` into `GpuAccelerator::new()` for workgroup size and buffer pool sizing. The WGSL shaders use a `WORKGROUP_SIZE` override constant — this must be set dynamically via shader preprocessing or `naga` constant override.

**Files:** `src/gpu/mod.rs`

### Step 7: Wire through `run_scf()`

The SCF entry point `run_scf()` needs access to the full `Settings` (or at least the sub-structs). Options:

**Option A — pass `&Settings` directly:** Simplest, but couples `run_scf` to the config format.

**Option B — expand `ScfParams`:** Add `NumericsSettings`, `EwaldSettings`, `InitialDensitySettings`, `GpuSettings` as fields on `ScfParams`. Keeps the existing API shape.

**Option C — new `ScfContext` struct:** Bundle `ScfParams` + numerics + ewald + initial_density + gpu into a single context. Clean but requires touching every `run_scf` call site.

Recommended: **Option B** — it's the least disruptive. `ScfParams` already carries mixing/smearing config; adding the remaining settings structs is consistent.

**Files:** `src/scf/mod.rs`, `src/main.rs`

## Verification

1. **Existing tests pass unchanged:** All defaults match current hardcoded values, so behavior is identical
2. **Non-default override test:** Parse a YAML with `rho_floor: 1e-14`, run XC evaluation, verify it uses the override (not the old `1e-30` literal)
3. **Ewald override test:** Explicit `ewald.eta` produces a different Ewald energy vs auto; verify both are physically reasonable
4. **Smearing edge case test:** Set `fermi_bisection_bounds: 50.0` and `fermi_bisection_tol: 1e-16`, verify Fermi search still converges
5. **Consistency fix verified:** After unifying XC density floors, `grep -n '1e-30' src/potential/xc.rs` returns 0 matches
6. **Full SCF unchanged:** Si SCF energy with all-default settings matches the value before this change to within machine epsilon
7. **Clippy clean:** `cargo clippy -q --all-targets`

## Estimated Effort

Two sessions.

- **Session 1:** Steps 1–4 (Ewald, numerics/potentials, smearing/density, initial density). These are the highest-value changes and touch the deepest call chains. ~200 lines across 10 files.
- **Session 2:** Steps 5–7 (symmetry, GPU, wiring through `run_scf`). Lighter touch, plus verification tests. ~100 lines across 5 files.
