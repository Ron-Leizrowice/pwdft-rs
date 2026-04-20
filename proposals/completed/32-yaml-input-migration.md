# Proposal 32: Migrate to Pure YAML Input

## Problem

The codebase has two parallel input parsing paths that describe the same physical calculation:

1. **`src/input.rs`** — TOML-based `InputFile`, used by `main.rs` and all 4 example files (`examples/*.toml`)
2. **`src/settings.rs`** — YAML-based `Settings`, with 1 example file (`examples/si_scf_settings.yaml`) but **never used by `main.rs`**

The TOML path is the only one wired to the CLI. The YAML path is strictly better:

| Concern | TOML (`InputFile`) | YAML (`Settings`) |
|---|---|---|
| Sections | 3 (`system`, `kpoints`, `scf`) | 9 (`system`, `basis`, `kpoints`, `scf`, `electrons`, `xc`, `symmetry`, `pseudopotentials`, `output`) |
| Defaults | Partial — `ecut` lives in `[system]`, mixing in `[scf]` | Full — every optional section has `#[serde(default)]` with QE-convention defaults |
| Smearing types | None — only `smearing_sigma: f64` | Enum: `fermi_dirac`, `gaussian`, `methfessel_paxton`, `cold`, `fixed` |
| XC functional | Not configurable | Enum: `pz`, `pbe`, `pbe0`, `hse06` |
| Symmetry control | Not exposed | `enabled`, `time_reversal`, `tolerance` |
| Output control | Not exposed | `verbosity`, `write_density`, `write_bands` |
| Conversion helpers | `to_crystal()`, `to_high_sym_path()` | All of those plus `to_scf_params()`, `to_symmetry_info()`, `pseudopotential_path()`, `ecutwfc()`, `mp_grid()` |
| Test coverage | 1 test | 18 tests (roundtrip, defaults, partial overrides, error cases) |
| `n_bands` location | `system.n_bands` (physics leak) | `scf.n_bands` (correct section) |
| `ecut` location | `system.ecut` (mixed concerns) | `basis.ecutwfc` (own section) |
| Pseudopotentials | Nested under `[scf.pseudopotentials]` | Top-level `pseudopotentials:` section |

The TOML format also has ergonomic problems: TOML's `[[system.atoms]]` array-of-tables syntax is verbose and harder to read than YAML's list syntax for atom definitions.

Maintaining two parallel parsers doubles the surface area for bugs and means every new field must be added in two places. Only the weaker parser is actually used.

## Research: Full Inventory of Unexposed Settings

An exhaustive audit of the codebase found ~60 hardcoded numeric values controlling algorithm behavior. After filtering out pure physics constants (PZ LDA parameters, Coulomb constant, unit conversions) and values with performance implications that justify hardcoding (FFT-friendly prime factors `[2,3,5]`, GPU workgroup size), the following should be exposed.

### A. `ScfParams` fields not wired through `Settings`

The `to_scf_params()` conversion (`settings.rs:356-367`) silently drops these fields via `..Default::default()`:

| `ScfParams` field | Default | Belongs in | Notes |
|---|---|---|---|
| `mixing_mode` | `Plain` | `electrons` | `main.rs:124` overrides to `Kerker`; user can't choose |
| `energy_threshold` | `1e-5` | `scf` | Independent energy convergence criterion |
| `smearing_scheme` | `FermiDirac` | `electrons` | `SmearingType` enum exists in settings but never maps to `ScfParams` |
| `fft_grid` | `None` | `basis` | QE-match example needs `[20,20,20]`; no YAML path |
| `nspin` | `1` | `electrons` | 1 (unpolarized) or 2 (collinear spin-polarized) |
| `starting_magnetization` | `{}` | `electrons` | Per-element initial magnetization for spin-polarized |
| `tot_magnetization` | `None` | `electrons` | Fixed total magnetization constraint |

### B. Ewald summation parameters

All in `src/ewald.rs`:

| Value | Line | Default | What it controls |
|---|---|---|---|
| eta | 52 | `(N·π/Ω)^{1/3}` | Real/reciprocal space partition; auto is standard but suboptimal for some ionic systems |
| cutoff multiplier | 57, 89 | `10.0` | `g_max = 10η`, `r_max = 10/η`; QE uses 4–5x, 10x is safe but wasteful for large cells |
| self-interaction threshold | 106 | `1e-10` | Distance below which atom pair is treated as self-interaction |

### C. Fermi energy search parameters

All in `src/scf/smearing.rs`:

| Value | Line | Default | What it controls |
|---|---|---|---|
| bisection bounds factor | 66–67 | `10.0 * σ.max(0.1)` | Expansion of eigenvalue range for Fermi search; can fail for extreme smearing widths |
| bisection max iterations | 69 | `200` | Iteration cap for Fermi energy bisection |
| bisection convergence | 85 | `1e-14` | Tolerance for Fermi energy (eV) |
| zero-sigma threshold | 114+ | `1e-15` | Sigma below which smearing is treated as zero (step function); used by all 4 schemes |
| degenerate-state tolerance | 116+ | `1e-12` | Energy distance from E_F below which occupation is set to 0.5; used by all 4 schemes |
| entropy overflow cutoff | 215 | `30.0` | `|x| > 30` saturates FD entropy to avoid exp() overflow |
| entropy occupation floor | 219 | `1e-30` | Minimum occupation for `-f·ln(f)` entropy term |

### D. Density and occupation thresholds

| Value | File:Line | Default | What it controls |
|---|---|---|---|
| `RHO_FLOOR` | `consts.rs:19` | `1e-20` | Global density floor for XC evaluation (e/ų) |
| XC density floor | `potential/xc.rs:27,177,237,270` | `1e-30` | Density below which XC is skipped; **inconsistent with RHO_FLOOR** — should unify |
| `G2_ZERO_THRESHOLD` | `consts.rs:17` | `1e-12` | \|G\|² below which G-vector is treated as zero in Hartree/local potential |
| occupation skip threshold | `scf/density.rs:61` | `1e-15` | Occupation below which a band is excluded from density construction |

### E. Mixing internals

All in `src/scf/mixing.rs`:

| Value | Line | Default | What it controls |
|---|---|---|---|
| Kerker q_tf | 69–72 | auto from `(4/3 · π · ρ_avg)^{1/3}` | Thomas-Fermi screening vector; `MixingMode::Kerker { q_tf: Option }` exists but isn't in YAML |
| DIIS pivot tolerance | 231 | `1e-15` | Gauss elimination pivot in Anderson/Pulay solver |

### F. Initial density

In `src/scf/initial_density.rs`:

| Value | Line | Default | What it controls |
|---|---|---|---|
| Gaussian sigma | 28 | `1.0` Å | Width of Gaussian model charges for SAD initialization |
| normalization tolerance | 103 | `1e-15` | Threshold for density integral normalization check |
| small-gr Bessel threshold | 166–167 | `1e-10` | Below this `gr`, j₀ uses Taylor expansion instead of `sin(x)/x` |

### G. Non-local potential

In `src/potential/nonlocal.rs`:

| Value | Line | Default | What it controls |
|---|---|---|---|
| q-norm threshold | 175–178 | `1e-12` | \|q\| below which cos(θ) defaults to 1.0 for Legendre evaluation |
| projector sum cutoff | 187 | `1e-20` | Projector sum magnitude below which non-local term is skipped |
| Bessel small-x threshold | 250 | `1e-10` | Below this x, spherical Bessel j_l uses Taylor expansion |

### H. Symmetry

| Value | File:Line | Default | What it controls |
|---|---|---|---|
| search radius factor | `symmetry/detect.rs:100` | `1.5` | Multiplier for candidate lattice vector search radius |
| k-point rounding tolerance | `symmetry/kpoints.rs:100` | `1e-6` | Tolerance for snapping k-points to grid after symmetry rotation |

### I. GPU

In `src/gpu/mod.rs`:

| Value | Line | Default | What it controls |
|---|---|---|---|
| workgroup size | 65 | `256` | Compute shader thread count per workgroup; optimal is hardware-dependent (64–512) |
| complex buffer pool | 130 | `5` | Pre-allocated complex f32 buffers |
| real buffer pool | 142 | `3` | Pre-allocated real f32 buffers |

### Values NOT exposed (with justification)

| Value | Location | Why hardcoded |
|---|---|---|
| PZ LDA coefficients (γ, β₁, β₂, a, b, c, d) | `potential/xc.rs` | Published physical parameters — changing them produces a different functional, not a tuning knob |
| Smearing math (M-P coefficient 0.5, exchange 3/4) | `potential/xc.rs`, `smearing.rs` | Algorithm-defined constants from the original papers |
| Zeta clamp bounds (-1, 1) | `potential/xc.rs:275` | Mathematical constraint: spin polarization ∈ [-1, 1] |
| FFT-friendly factors [2, 3, 5] | `fft.rs:120-121` | **Performance:** ndrustfft is optimized for these; allowing arbitrary primes would degrade performance |
| Grid Nyquist formula `2*n_max + 1` | `fft.rs:105` | Physics requirement from sampling theorem, not tunable |

### Proposed full YAML schema

```yaml
system:
  lattice: [[ax,ay,az], [bx,by,bz], [cx,cy,cz]]
  atoms:
    - { symbol: Si, position: [0.0, 0.0, 0.0] }

basis:
  ecutwfc: 204.09               # eV, wavefunction cutoff
  ecutrho_ratio: 4              # charge density cutoff = ratio * ecutwfc
  fft_grid: null                # explicit [nx,ny,nz]; overrides ecutrho_ratio

kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]

scf:
  max_iter: 100
  conv_threshold: 1.0e-6        # density RMS (e/ų)
  energy_threshold: 1.0e-5      # energy change (eV)
  n_bands: null                  # null = auto

electrons:
  mixing_beta: 0.3
  mixing_ndim: 8
  mixing_mode: kerker            # plain | kerker
  kerker_q_tf: null              # null = auto from electron density
  diis_pivot_tolerance: 1.0e-15  # Gauss elimination pivot in Anderson/Pulay
  smearing: fermi_dirac          # fermi_dirac | gaussian | methfessel_paxton | cold | fixed
  smearing_width: 0.05           # eV
  occupations: smearing          # smearing | fixed
  nspin: 1                       # 1 = unpolarized, 2 = collinear
  starting_magnetization: {}     # e.g. { Fe: 0.5 }
  tot_magnetization: null        # null = self-consistent
  fermi_bisection_bounds: 10.0   # eigenvalue range expansion factor for Fermi search
  fermi_bisection_maxiter: 200
  fermi_bisection_tol: 1.0e-14   # eV
  occupation_threshold: 1.0e-15  # f < this → band excluded from density

xc:
  functional: pz                 # pz | pbe | pbe0 | hse06

symmetry:
  enabled: true
  time_reversal: true
  tolerance: 1.0e-5              # fractional coordinate tolerance
  search_radius_factor: 1.5      # candidate vector search multiplier
  kpoint_rounding_tol: 1.0e-6    # k-point grid snapping tolerance

pseudopotentials:
  Si: "Si.upf"

ewald:
  eta: null                      # null = auto (N·π/Ω)^{1/3}
  cutoff_multiplier: 10.0        # g_max = mult*eta, r_max = mult/eta
  self_interaction_threshold: 1.0e-10

numerics:
  rho_floor: 1.0e-20             # minimum density for XC (e/ų)
  g2_zero_threshold: 1.0e-12     # |G|² treated as zero in Coulomb sums
  zero_sigma_threshold: 1.0e-15  # σ below this → step-function occupations
  degenerate_state_tol: 1.0e-12  # |E - E_F| below this → occupation = 0.5
  entropy_overflow_cutoff: 30.0  # |x| > this saturates FD entropy
  entropy_occupation_floor: 1.0e-30
  projector_sum_cutoff: 1.0e-20  # non-local projector magnitude cutoff
  bessel_small_x_threshold: 1.0e-10  # spherical Bessel Taylor expansion cutoff
  q_norm_threshold: 1.0e-12      # |q| below this → cos(θ) = 1 for Legendre

initial_density:
  method: sad                    # sad | gaussian (future: random)
  gaussian_sigma: 1.0            # Å, width of model charges
  normalization_tol: 1.0e-15     # density integral check

gpu:
  workgroup_size: 256            # threads per compute workgroup (64–512)
  complex_buffers: 5             # pre-allocated complex buffer pool
  real_buffers: 3                # pre-allocated real buffer pool

output:
  verbosity: normal              # low | normal | high
  write_density: false
  write_bands: true
```

## Implementation

### Step 1: Extend `Settings` schema with all missing fields

Add every field from the research inventory to `settings.rs`. All new fields use `#[serde(default)]` with the same values currently hardcoded, so existing YAML files parse unchanged.

**Extend `ElectronSettings`:**

```rust
pub mixing_mode: MixingModeType,         // plain | kerker (default: kerker)
pub kerker_q_tf: Option<f64>,            // None = auto
pub diis_pivot_tolerance: f64,           // default: 1e-15
pub nspin: usize,                        // default: 1
pub starting_magnetization: HashMap<String, f64>,
pub tot_magnetization: Option<f64>,
pub fermi_bisection_bounds: f64,         // default: 10.0
pub fermi_bisection_maxiter: usize,      // default: 200
pub fermi_bisection_tol: f64,            // default: 1e-14
pub occupation_threshold: f64,           // default: 1e-15
```

**Extend `ScfSettings`:**

```rust
pub energy_threshold: f64,               // default: 1e-5
```

**Extend `BasisSettings`:**

```rust
pub fft_grid: Option<[usize; 3]>,
```

**Extend `SymmetrySettings`:**

```rust
pub search_radius_factor: f64,           // default: 1.5
pub kpoint_rounding_tol: f64,            // default: 1e-6
```

**New `EwaldSettings`:**

```rust
pub eta: Option<f64>,                    // None = auto
pub cutoff_multiplier: f64,              // default: 10.0
pub self_interaction_threshold: f64,     // default: 1e-10
```

**New `NumericsSettings`:**

```rust
pub rho_floor: f64,                      // default: 1e-20
pub g2_zero_threshold: f64,              // default: 1e-12
pub zero_sigma_threshold: f64,           // default: 1e-15
pub degenerate_state_tol: f64,           // default: 1e-12
pub entropy_overflow_cutoff: f64,        // default: 30.0
pub entropy_occupation_floor: f64,       // default: 1e-30
pub projector_sum_cutoff: f64,           // default: 1e-20
pub bessel_small_x_threshold: f64,       // default: 1e-10
pub q_norm_threshold: f64,               // default: 1e-12
```

**New `InitialDensitySettings`:**

```rust
pub method: InitialDensityMethod,        // sad | gaussian (default: sad)
pub gaussian_sigma: f64,                 // default: 1.0
pub normalization_tol: f64,              // default: 1e-15
```

**New `GpuSettings`:**

```rust
pub workgroup_size: u32,                 // default: 256
pub complex_buffers: usize,              // default: 5
pub real_buffers: usize,                 // default: 3
```

Also unify `potential/xc.rs`'s independent `1e-30` density floor with `NumericsSettings::rho_floor`.

Update `to_scf_params()` to map **every** field explicitly — remove `..Default::default()`.

**Files:** `src/settings.rs`

### Step 2: Wire `Settings` into `main.rs`

Replace `InputFile::from_file` with format detection based on file extension:

```rust
let config = match cli.input.extension().and_then(|e| e.to_str()) {
    Some("yaml" | "yml") => Settings::from_yaml_file(&cli.input)?,
    Some("toml") => Settings::from_toml_file(&cli.input)?,  // temporary bridge
    _ => return Err(PwdftError::Parse(
        "input file must have .yaml or .toml extension".into()
    )),
};
```

Refactor `main.rs` to use `Settings` as the single config type. The TOML branch would parse via `toml` then convert to `Settings` internally (temporary compatibility shim).

**Files:** `src/main.rs`, `src/settings.rs` (add `from_toml_file` bridge method)

### Step 3: Add TOML-to-Settings bridge

Add a `from_toml_str` / `from_toml_file` method on `Settings` that parses the old TOML format and maps it into the `Settings` struct. This provides backwards compatibility while we migrate example files:

```rust
impl Settings {
    pub fn from_toml_file(path: &Path) -> Result<Self> {
        let contents = std::fs::read_to_string(path)?;
        let old: InputFile = toml::from_str(&contents)
            .map_err(|e| PwdftError::Parse(e.to_string()))?;
        Ok(Self::from_legacy_input(old))
    }

    fn from_legacy_input(input: InputFile) -> Self {
        // Map InputFile fields → Settings sections
        // ...
    }
}
```

**Files:** `src/settings.rs`

### Step 4: Convert example files to YAML

Convert all 4 TOML examples to YAML equivalents:

| Old | New |
|---|---|
| `examples/si_scf.toml` | `examples/si_scf.yaml` |
| `examples/si_scf_converged.toml` | `examples/si_scf_converged.yaml` |
| `examples/si_scf_qe_match.toml` | `examples/si_scf_qe_match.yaml` |
| `examples/si_free_electron.toml` | `examples/si_free_electron.yaml` |

The existing `examples/si_scf_settings.yaml` serves as the template. Each new file should use the full sectioned layout with comments explaining each parameter.

**Files:** `examples/*.yaml` (4 new files), delete `examples/*.toml` (4 files)

### Step 5: Thread settings into call sites

This is the most invasive step. Every hardcoded value that moved into the YAML schema must be plumbed from `Settings` through to the code that uses it. The approach is to carry settings structs (or relevant fields) through the SCF context rather than reading global constants.

**Ewald:** Add `&EwaldSettings` parameter to `ewald_energy()`. Replace hardcoded `10.0` multiplier and auto-eta with values from settings. Pass `self_interaction_threshold` instead of `1e-10`.

**Numerics → potentials:** The SCF context (`src/scf/context.rs` or a new `NumericsContext`) carries `rho_floor`, `g2_zero_threshold`, etc. Pass to:

- `potential/xc.rs` — replace all 4 independent `1e-30` floors with `rho_floor`
- `scf/energy.rs` — replace `G2_ZERO_THRESHOLD` reads with context field
- `scf/mixing.rs` — pass `diis_pivot_tolerance` and `g2_zero_threshold`
- `potential/nonlocal.rs` — pass `projector_sum_cutoff`, `bessel_small_x_threshold`, `q_norm_threshold`

**Electrons → smearing:** Pass `fermi_bisection_bounds`, `fermi_bisection_maxiter`, `fermi_bisection_tol`, `zero_sigma_threshold`, `degenerate_state_tol`, `entropy_overflow_cutoff`, `entropy_occupation_floor` into `find_fermi_energy()` and the occupation/entropy functions.

**Electrons → density:** Pass `occupation_threshold` into `construct_density()`.

**Initial density:** Pass `InitialDensitySettings` into `initial_density()`.

**Symmetry:** Pass `search_radius_factor` into `SymmetryInfo::from_crystal()` and `kpoint_rounding_tol` into `reduce_kpoints()`.

**GPU:** Pass `GpuSettings` into `GpuAccelerator::new()` for workgroup size and buffer pool counts.

**`consts.rs` cleanup:** `G2_ZERO_THRESHOLD` and `RHO_FLOOR` can remain as **default values** (referenced by `NumericsSettings::default()`), but call sites should read from the settings struct rather than the global constant.

**Files:** `src/ewald.rs`, `src/potential/xc.rs`, `src/potential/nonlocal.rs`, `src/scf/mod.rs`, `src/scf/energy.rs`, `src/scf/mixing.rs`, `src/scf/smearing.rs`, `src/scf/density.rs`, `src/scf/initial_density.rs`, `src/symmetry/detect.rs`, `src/symmetry/kpoints.rs`, `src/gpu/mod.rs`, `src/consts.rs`

### Step 6: Remove `input.rs` and TOML dependency

Once all examples and docs reference YAML:

1. Delete `src/input.rs`
2. Remove `pub mod input;` from `src/lib.rs`
3. Remove the `toml` crate from `Cargo.toml`
4. Remove the TOML bridge code added in Step 2
5. Update `Cli` struct doc comment from "Path to TOML input file" to "Path to YAML input file"
6. Update CLAUDE.md references (run commands, conventions section)

**Files:** `src/input.rs` (delete), `src/lib.rs`, `Cargo.toml`, `src/main.rs`, `CLAUDE.md`

### Step 7: Update documentation and CLI help

- Update `CLAUDE.md`: change all TOML references to YAML, update example commands
- Update `Cli` struct `#[arg]` help text
- Verify no stale `.toml` references remain anywhere in `src/` or `tests/`

**Files:** `CLAUDE.md`, `src/main.rs`

## Verification

1. **Existing tests pass unchanged:** `cargo test` — the 18 existing `settings.rs` tests must still pass (minimal YAML with defaults)
2. **New section parsing tests:** Add tests for every new section (`ewald`, `numerics`, `initial_density`, `gpu`) and every new field on extended sections — both default-only and explicit-override variants
3. **Full roundtrip test:** Parse a YAML with every field explicitly set → `to_scf_params()` → assert **every** field matches (catches any `..Default::default()` silently dropping values)
4. **Numerics reach call sites:** Unit test that constructs an SCF with a non-default `rho_floor` (e.g. `1e-14`) and verifies XC evaluation uses it (not the old `1e-30` literal)
5. **Ewald override:** Test that explicit `ewald.eta` produces a different (correct) Ewald energy vs auto
6. **Example files parse:** `cargo run --release -- --input examples/si_scf.yaml` converges to the same energy as TOML version
7. **Band structure:** `cargo run --release -- --input examples/si_free_electron.yaml -o bands.tsv`
8. **No TOML remnants:** `grep -r '\.toml' src/ tests/` returns nothing (except Cargo.toml itself)
9. **Clippy clean:** `cargo clippy -q --all-targets`

## Estimated Effort

Three sessions.

- **Session 1:** Schema extension (step 1) + wire into main + TOML bridge + convert examples (steps 2–4). The `settings.rs` foundation is solid; this is additive.
- **Session 2:** Thread settings into call sites (step 5). This is the bulk of the work — ~13 files need their hardcoded values replaced with parameters. The changes are mechanical but touch deep call chains (SCF → smearing, SCF → potentials, SCF → density).
- **Session 3:** Remove TOML (step 6), update docs (step 7), write verification tests, clean up.
