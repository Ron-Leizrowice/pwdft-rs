---
id: CFGN
status: active
priority: low
complexity: medium
risk: low
depends_on: []
blocks: []
---

# CFGN: Expose Hardcoded Numerics as Configurable Settings

> **Re-scope 2026-04-19 (post-CFGN1).** CFGN1 (PR #114) landed knob #7
> (`initial_density.gaussian_sigma`) — the pattern is proven, the rest
> now breaks up into small independent follow-ups. **Read Section 11
> first** — it is the current scope. Sections 6–10 are still the fresh
> 2026-04-18 census but have been pared against what's landed and
> against what ECUT / ESPL (new proposals 2026-04-19) cover.
>
> **Re-scope 2026-04-18.** The original inventory below dates from before
> MODR (scf/mixing and symmetry/density folder splits, pseudopotential/upf
> folder split, driver.rs/driver_spin.rs/report.rs extraction), CAST
> (per-site cast invariants encoding physical bounds), DWGT (rustdoc
> `-D warnings` gate), FGRD (`MAX_FFT_DIM = 1024`), MXBA (`adaptive_beta`
> Settings landed), ITEV (iterative eigensolver with `DEFAULT_TOL` /
> `DEFAULT_MAX_RESTARTS` consts), and NCFX (the four `1e-30` XC density
> floors unified to `RHO_FLOOR`). Every file path and line number in the
> original inventory is stale or has been superseded. Read **Sections 6–10
> first** — they are the current census, phasing, and decisions. The
> original Sections 1–5 are kept verbatim for audit trail rather than
> struck so the deltas are visible.

## Section 1 — Problem (original; still accurate in spirit)

An audit of the codebase found ~60 hardcoded numeric values controlling algorithm behavior — convergence tolerances, cutoff radii, density floors, iteration limits, buffer sizes — scattered across 13 source files. These are all reasonable defaults, but an advanced researcher has no way to override them without editing source code. This blocks:

- **Sensitivity studies:** Can't sweep rho_floor or Ewald cutoff to assess numerical stability
- **Debugging:** Can't relax/tighten thresholds to isolate convergence issues
- **System-specific tuning:** Metals vs insulators, large vs small cells, extreme smearing regimes all benefit from different numerics
- **Hardware tuning:** GPU workgroup size is hardware-dependent (M2 vs A100)

Several values are also **inconsistent** — `potential/xc.rs` uses `1e-30` as a density floor in 4 independent locations while `consts.rs` defines `RHO_FLOOR = 1e-20`.

This proposal adds the settings structs and threads them through call sites. It depends on:
- **Proposal 32** (YAML migration, completed) having established `Settings` in `settings.rs` as the canonical config path
- **Proposal 33** (consolidate constants) to unify E2 and RHO_FLOOR into single canonical locations before we thread them — otherwise we'd be plumbing settings to duplicated definitions

The schema fields for the new sections (`ewald`, `numerics`, `initial_density`, `gpu`) need to be added to `Settings`, and the values plumbed through to call sites.

## Section 2 — Inventory (ORIGINAL; superseded by Section 6)

> **Stale-reference warnings inline.** File:line numbers in this section
> were valid at proposal time but most have drifted (MODR, NCFX, VNLM
> changed the contents). See Section 6 for the fresh census and Section 8
> for the full stale-reference list.

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

## Section 3 — Implementation (ORIGINAL; superseded by Section 10)

> **Superseded:** Steps 1–7 below reference file paths and APIs that have
> changed (no `src/potential/hartree.rs`; `src/scf/mod.rs` is now a thin
> dispatcher + `ScfParams`; mixing is a folder; densities are symmetrized
> in G-space). See Section 10 for the current phasing plan.

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

## Section 4 — Verification (ORIGINAL; still applicable)

1. **Existing tests pass unchanged:** All defaults match current hardcoded values, so behavior is identical
2. **Non-default override test:** Parse a YAML with `rho_floor: 1e-14`, run XC evaluation, verify it uses the override (not the old `1e-30` literal)
3. **Ewald override test:** Explicit `ewald.eta` produces a different Ewald energy vs auto; verify both are physically reasonable
4. **Smearing edge case test:** Set `fermi_bisection_bounds: 50.0` and `fermi_bisection_tol: 1e-16`, verify Fermi search still converges
5. **Consistency fix verified:** After unifying XC density floors, `grep -n '1e-30' src/potential/xc.rs` returns 0 matches
6. **Full SCF unchanged:** Si SCF energy with all-default settings matches the value before this change to within machine epsilon
7. **Clippy clean:** `cargo clippy -q --all-targets`

## Section 5 — Estimated Effort (ORIGINAL; superseded by Section 10)

Two sessions.

- **Session 1:** Steps 1–4 (Ewald, numerics/potentials, smearing/density, initial density). These are the highest-value changes and touch the deepest call chains. ~200 lines across 10 files.
- **Session 2:** Steps 5–7 (symmetry, GPU, wiring through `run_scf`). Lighter touch, plus verification tests. ~100 lines across 5 files.

---

# Re-scope — 2026-04-18

## Section 6 — Fresh census (supersedes Section 2)

A file-by-file walk of `src/` on `origin/main` @ `0be9290` (XCNI). Every
entry is cited by `file:line` against the current tree. Only entries
matching **all three** filters are listed: physically meaningful
(convergence/tolerance/cutoff — not a buffer size, unit conversion, or
published functional coefficient), currently hardcoded (literal in code,
or module-private `const`, not already a Settings field), and reasonably
user-facing (advanced researcher would override, not an internal
DIIS/LU pivot).

### 6.1 — Strong candidates (user-facing knobs, high value)

| # | file:line | current literal | purpose | proposed Settings field |
|---|-----------|-----------------|---------|-------------------------|
| 1 | `src/ewald.rs:65` | `10.0 * eta` (`g_max`) | Ewald reciprocal-space cutoff multiplier; QE's `alpha_mod.f90` uses 4–5× at convergence | `ewald.cutoff_multiplier: f64` (default `10.0`) |
| 2 | `src/ewald.rs:107` | `10.0 / eta` (`r_max`) | Ewald real-space cutoff multiplier; **must stay paired with #1** — the same eta-scaling produces balanced r/g truncation error | same knob as #1 |
| 3 | `src/ewald.rs:60` | `(N·π/Ω)^{1/3}` (auto `eta`) | Ewald partition parameter; auto-chosen to balance r/g work. Users rarely want to override but sensitivity studies do | `ewald.eta: Option<f64>` (default `None` = auto) |
| 4 | `src/scf/smearing.rs:73-74` | `10.0 * sigma.max(0.1)` | Fermi bisection bounds factor and sigma floor. For very cold systems (σ < 0.001 eV) the `.max(0.1)` is overly loose and can cause bisection to straddle unphysical bands | `electrons.fermi_search.bounds_factor: f64` (default `10.0`); `fermi_search.sigma_floor: f64` (default `0.1`) |
| 5 | `src/scf/smearing.rs:76` | `200` (max_iter) | Fermi bisection iteration cap | `electrons.fermi_search.max_iter: usize` (default `200`) |
| 6 | `src/scf/smearing.rs:92` | `1e-14` | Fermi bisection convergence tolerance (eV) | `electrons.fermi_search.tol: f64` (default `1e-14`) |
| 7 | `src/scf/initial_density.rs:28` | `DEFAULT_GAUSSIAN_SIGMA = 1.0` Å | SAD fallback Gaussian width when UPF has no PP_RHOATOM. `InitialDensityConfig.gaussian_sigma` is already a runtime `Option<f64>`, just unwired from Settings | `initial_density.gaussian_sigma: Option<f64>` (default `None` = 1.0 Å) |
| 8 | `src/eigensolver/iterative.rs:67` | `DEFAULT_TOL = f64::EPSILON * 128.0` (≈2.8e-14) | Iterative eigensolver residual-norm convergence threshold. Loosening to 1e-10 would trade SCF convergence rate for per-iter speed | `scf.iterative_eigensolver.tol: f64` (default `DEFAULT_TOL`) |
| 9 | `src/eigensolver/iterative.rs:74` | `DEFAULT_MAX_RESTARTS = 500` | Iterative eigensolver max Arnoldi restarts before falling back to Dense. Cheap failure case has this knob act as a "give up, use Dense" gate | `scf.iterative_eigensolver.max_restarts: usize` (default `500`) |
| 10 | `src/scf/smearing.rs:225` | `30.0` | Fermi-Dirac entropy reduced-variable cutoff (|(ε−E_F)/σ| > 30 → s=0). Safe at default σ but a user doing very-cold smearing (σ = 1e-4 eV) might notice | `electrons.fermi_search.entropy_cutoff: f64` (default `30.0`) — low priority |

### 6.2 — Medium candidates (borderline user-facing)

| # | file:line | current literal | purpose | proposed Settings field |
|---|-----------|-----------------|---------|-------------------------|
| 11 | `src/consts.rs:17` | `RHO_FLOOR = 1e-30` e/Å³ | XC density floor; NCFX unified the four `xc.rs` literals to this constant (post-CFGN-original) | `electrons.xc.rho_floor: f64` (default `RHO_FLOOR`) |
| 12 | `src/consts.rs:15` | `G2_ZERO_THRESHOLD = 1e-12` (Å⁻²) | `|G|²` treated as zero in Hartree and mixing. Consistently used via `crate::consts::G2_ZERO_THRESHOLD` in `scf/energy.rs:70,240`, `scf/mixing/{anderson,broyden}.rs`, `gpu/mod.rs:519` | `electrons.g2_zero_threshold: f64` (default `G2_ZERO_THRESHOLD`) |
| 13 | `src/scf/smearing.rs:{123,136,151,171}` | `1e-15` (zero-σ threshold) and `1e-12` (degenerate tolerance) — duplicated in all 4 occupation branches | σ below which → step function; (ε−E_F) below which → f=0.5 | bundle with fermi_search knobs (#4–6) |
| 14 | `src/scf/mixing/linalg.rs:39` | `1e-15` | Gauss elimination pivot safety in the mixer's DIIS linear solve. Borderline internal — a hand-tuned mixer debugger might want to see this, but normal users should never touch it | **Do not expose** |
| 15 | `src/scf/density.rs:61` | `1e-15` (occupation skip) | Bands with `f · w_k < 1e-15` skipped in density construction. Has an interaction with #4–6 | **Do not expose** initially — internal, paired with smearing |
| 16 | `src/scf/density.rs:92` | `1e-15` | Normalization integrand safety check | **Do not expose** |
| 17 | `src/pseudopotential/mod.rs:137` | `1e-12` | `|G|=0` check in `v_local_of_g`. Should use `G2_ZERO_THRESHOLD` const instead (or its sqrt) | **Code cleanup** — Code Reviewer; not a CFGN knob |

### 6.3 — Non-candidates (keep hardcoded, with justification)

| file:line | literal | why keep as-is |
|-----------|---------|----------------|
| `src/scf/grid.rs:26` | `MAX_FFT_DIM = 1024` | FGRD landed this as a `pub(crate) const` that gates CAST-documented integer-range assumptions on ~20 sites. Exposing via YAML would invalidate every `#[allow(...)] reason = "..."` that depends on it. **Not a CFGN candidate** — it's a compile-time safety invariant, not a tuning knob |
| `src/potential/xc.rs:33` | `XC_PARALLEL_THRESHOLD = 16_384` | Rayon dispatch-crossover point from calibrated M2 measurement. Performance tuning — Performance Engineer's domain |
| `src/gpu/mod.rs:72` | `WORKGROUP_SIZE: u32 = 256` | Baked into WGSL shader constants; changing it requires shader preprocessing. Not a Settings-accessible knob without a GPU-settings rewrite |
| `src/gpu/mod.rs:140` | `(0..5)` — complex buffer pool size | Implementation detail (number of concurrent dispatched pipelines, not a physical knob). **Original proposal's "real buffer pool = 3" no longer exists — that was removed when the pool was simplified** |
| `src/potential/nonlocal.rs:290` | `1e-20` (D_ij skip threshold) | Guard against pure-zero D_ij entries polluting GEMM output with NaN. Not user-facing |
| `src/potential/nonlocal.rs:377` | `1e-9` (`q-norm` threshold for Y_lm) | Branch between `Y_00=1/√(4π)` singular case and full Y_lm evaluation. Mathematical boundary, not tuning knob |
| `src/potential/nonlocal.rs:494`, `src/scf/potentials.rs:119`, `src/scf/initial_density.rs:171` | `1e-10` (Bessel/j₀ small-x thresholds) | Taylor expansion branch points for `j_l(x)/x` at x→0. Mathematical, not tunable |
| `src/ewald.rs:91`, `136` | `1e-12`, `1e-10` | Self-interaction and near-zero |G|² guards. Internal Ewald numerical safety |
| `src/symmetry/detect.rs:105` | `+ 1.5` | Half-unit search-radius padding in lattice-vector candidate generation. Pure geometry; not a physics knob |
| `src/symmetry/kpoints.rs:130` | `1e-6` | k-point grid snapping after symmetry rotation. Geometry-on-integer-grid threshold, closely tied to `SymmetrySettings::tolerance` but at a different scale |
| `src/consts.rs:11` | `E2_COULOMB = 14.399_645_351_950_548` | Published CODATA physical constant |

### 6.4 — Count

- **Strong candidates:** 10 (entries #1–10 above)
- **Medium candidates:** 2 that would actually be exposed (#11 RHO_FLOOR, #12 G2_ZERO_THRESHOLD); #13 bundles into #4–6
- **Non-candidates kept hardcoded:** 11 file-groups with explicit rationale

Total proposed Settings additions: **~12 fields**. The original proposal
listed ~25 distinct knobs (sections A–H); the fresh census finds only
10–12 that survive the "reasonably exposed" filter once we honor NCFX's
unification, MODR's refactors, and CAST's documented invariants.

## Section 7 — What CAST's audit implicitly documented

CAST landed ~51 `#[allow(clippy::cast_*, reason = "...")]` reasons
across 13 files. Each reason encodes an integer-range bound that the
cast relies on for safety. Reviewed against "could a CFGN user reasonably
want to override this?":

| CAST reason-fragment (file:line) | implied bound | CFGN candidate? |
|----------------------------------|---------------|-----------------|
| `src/scf/grid.rs:111,137` + 4 call sites in `src/symmetry/density/{mod,real_space,g_space}.rs` | `MAX_FFT_DIM ≤ 1024` | **No.** Exposing this would force `check-cast-safety` reviews on every linked site. A user with an 8k³ grid has bigger problems than Settings plumbing |
| `src/ewald.rs:68,73,78,110,115,120` | "n_i_max bounded by g_max = 10·eta" | **Yes.** But already captured by candidate #1 (`ewald.cutoff_multiplier`). If a user sets `cutoff_multiplier: 50.0`, the Ewald i32 casts still fit (g_max · |b_i| ≤ 50·η·|b_i| ≤ ~10⁵ for physical inputs) |
| `src/basis.rs:37,42,47` | "n_i_max bounded by ecut; exceeding i32::MAX would require ecut > 10¹⁸ eV" | **No.** `basis.ecutwfc` already in Settings; no further knob needed |
| `src/kpoints.rs:46,56,61,66`; `src/symmetry/kpoints.rs:107,135,142` | "MP mesh counts bounded by O(100)" | **No.** `kpoints.grid` already in Settings |
| `src/scf/grid.rs:59` | "ecutrho_ratio is a small input integer (typically 4)" | **No.** Already exposed |
| `src/potential/nonlocal.rs:142,231,249,295,359` | "each projector's `l` asserted non-negative" / "lmax is a non-negative angular-momentum bound" | **No.** `l_max` comes from the PP file — not a Settings concept. Max physical lmax is ~6 for any real pseudopotential; no user use-case for overriding |
| `src/symmetry/detect.rs:102` | "target_norm_sq from metric matrix" | **No.** Internal geometry |
| `src/gpu/mod.rs:198,290,354` | "n_grid ≤ 512³ fits in u32::MAX" | **No.** This is a free bound (32-bit grid index range), not a hardcoded cap |

**Net:** CAST's invariants reinforce the existing Settings shape; they
surface **zero new** CFGN candidates. The `MAX_FFT_DIM` question is
explicitly a "don't touch, it's load-bearing across 20 sites" answer.

## Section 8 — Stale references in Section 2 (original inventory)

Every file/symbol cited in the original Inventory (Section 2), verified
against `origin/main` @ `0be9290`:

| Original claim | Current reality | Note |
|----------------|-----------------|------|
| `src/ewald.rs:52` (eta) | Line **60** | MODR didn't touch; small line drift from doc expansion |
| `src/ewald.rs:57, 89` (cutoff multiplier) | Lines **65** and **107** | Same |
| `src/ewald.rs:106` (self-interaction threshold) | Line **136** | Same |
| `src/scf/smearing.rs:66-67,69,85,114+,116+,215,219` | Lines **73-74, 76, 92, 123+, 125+, 225, 229** | Largely unchanged, ~+10 line drift |
| `src/potential/xc.rs:27,177,237,270` (four `1e-30` density floors) | **GONE.** NCFX (2026-04-18) unified all four to `crate::consts::RHO_FLOOR` at lines 49, 212, 278, 288, 293, 321 | **Consistency bug fixed upstream of CFGN**. The "4 independent literals" claim is no longer true |
| `src/consts.rs:19` `RHO_FLOOR = 1e-20` | **`1e-30`** at line 17 | NCFX unified; original proposal's "inconsistent `1e-30` vs `1e-20`" complaint is resolved |
| `src/scf/density.rs:61` (occupation skip) | Line **61** | Unchanged |
| `src/scf/mixing.rs:69-72, 231` | **File does not exist.** MODR split into `src/scf/mixing/{mod,anderson,broyden,kerker,linalg}.rs`. Original "Kerker q_tf:69-72" ≈ `src/scf/mixing/kerker.rs:47-55`. "DIIS pivot tolerance:231" ≈ `src/scf/mixing/linalg.rs:39` | Stale path |
| `src/scf/initial_density.rs:28, 103, 166-167` | Lines **28 (DEFAULT_GAUSSIAN_SIGMA, unchanged), 107, 171** | Small drift |
| `src/potential/nonlocal.rs:175-178, 187, 250` | Lines **~377 (q-norm 1e-9, not 1e-12!), 290 (1e-20, unchanged), 494** | **Bug in original proposal:** claimed "q-norm threshold 1e-12" — actual value is **`1e-9`**, and the branch is on `q.norm() < eps` (real-space Y_lm), not `cos(θ)` as original claimed. VNLM rewrote this section to use Y_lm addition theorem |
| `src/symmetry/detect.rs:100` | Line **105**, factor `+ 1.5` unchanged |
| `src/symmetry/kpoints.rs:100` | Line **130** |
| `src/gpu/mod.rs:65` (workgroup size) | Line **72** (`WORKGROUP_SIZE = 256`) | Unchanged |
| `src/gpu/mod.rs:130` (complex buffer pool = 5) | Line **140** (`(0..5)`) | Still 5 |
| `src/gpu/mod.rs:142` (real buffer pool = 3) | **GONE.** There is no separate real buffer pool in the current `BufferPool`; the struct has `complex_bufs` + `complex_staging` + `g_squared_buf` only. Original proposal described a structure that was later simplified | **Stale; do not port** |
| Implementation Step 2 mentions `src/potential/hartree.rs` | **No such file.** Hartree is assembled inline in `scf::energy`, `scf::driver`, `scf::driver_spin` | Original proposal plumbing plan needs rewriting against the current scf/ layout |
| Implementation Step 7 mentions `run_scf()` in `src/scf/mod.rs` | `run_scf` now a thin dispatcher (lines ~222+); real work in `src/scf/driver.rs` and `src/scf/driver_spin.rs` | Threading plan updated in Section 10 |

**Summary of stale refs:** 100% of line numbers drifted (MODR + NCFX +
VNLM). Two **material errors**: (a) the "4 independent `1e-30` XC
literals" bug was fixed by NCFX and is not a CFGN item any more; (b) the
nonlocal q-norm threshold is `1e-9`, not `1e-12`, and the branch is on
the wrong quantity in the original proposal. One **phantom file**
(`potential/hartree.rs`) and one **phantom structure** (GPU "real buffer
pool = 3") cited in the original are fabrications against the current
codebase.

## Section 9 — MXBA tunables decision

MXBA's `AdaptiveBeta` carries four magic numbers seen at
`src/scf/mixing/mod.rs:120-123`:

```rust
growth_threshold: 1.2,
restore_threshold: 0.5,
damp_factor: 0.7,
restore_window: 3,
```

plus the implicit `beta_min = max(0.05·β_start, 0.01)` from the
constructor (line 119).

**Decision: Do NOT expose via CFGN.** Rationale (one sentence each):

- These are Eyert (1996, §3.3) paper-recommended defaults, not knobs a
  SCF user would tune from a YAML input — `tests/mxba_adaptive_beta_fe.rs`
  and the Fe-trajectory regression guard at
  `src/scf/mixing/mod.rs:437-526` encode an interaction with β_min and
  the restore-window that would break if a user slides one number.
- The **right way** to override these (if MXB2 or a follow-up mixer
  paper shows a better regime) is a new `mixing_mode` variant, not a
  user-facing parameter: the defaults embed an algorithmic contract.
- The feature is already opt-in (`adaptive_beta: false` by default), so
  users don't inherit the tuning unless they explicitly opt into the
  algorithm as documented.

MXB2's scope stays "defaults or a new mixer variant" — not "add four
Settings fields."

## Section 10 — Implementation phasing (supersedes Section 3)

With the census down from ~25 to ~12 knobs, a smaller two-PR split is
more tractable than the original 7-step plan. Each phase is an
independent PR that doesn't depend on the other.

### Phase 1 — Smearing / eigensolver / initial-density (highest value)

Knobs #4, #5, #6, #7, #8, #9 from Section 6.1 — these are the
user-facing ones a researcher doing a sensitivity study would actually
hit. Estimated scope: new `FermiSearchSettings` sub-struct under
`ElectronSettings`; a `scf.iterative_eigensolver` sub-struct with
`tol`/`max_restarts`; `initial_density.gaussian_sigma` plumbed through
`ScfParams` into `InitialDensityConfig.gaussian_sigma`. Threading into
the drivers:

- `src/scf/driver.rs` — call site for `find_fermi_energy` and
  `diagonalize_dispatch` (iterative path is keyed by `EigensolverKind`;
  add an `IterativeParams` struct carried on `ScfParams`).
- `src/scf/driver_spin.rs` — same call sites (spin driver).
- `src/scf/mod.rs::run_scf` — no change; the thin dispatcher already
  accepts `ScfParams`.

Changes to `ScfParams` + `Settings::to_scf_params`. No code physics
changes. Verification via:
- Existing `qe_validation.rs` Tier 1+2 systems must produce bit-identical
  results when YAML omits the new fields (defaults match hardcoded
  values).
- One new test per knob: override YAML, then check the override is
  observed (e.g., set `fermi_search.tol: 1e-8` → Fermi energy matches
  bisection with that tolerance, not 1e-14).

Estimated: ~150 lines across 6 files.

### Phase 2 — Ewald / RHO_FLOOR / G2_ZERO_THRESHOLD (nice-to-have)

Knobs #1, #2, #3, #11, #12 from Section 6. These are stress-test knobs
that advanced users wanting to debug a numerical residual would want,
but the defaults are already calibrated against QE and no user is
asking for them today.

- `ewald_energy` gains an `&EwaldSettings` argument.
- `Settings::rho_floor` and `Settings::g2_zero_threshold` become
  `NumericsSettings` fields threaded through `ScfParams`.
- XC call sites (`lda_xc`, `lda_xc_grid`, spin variants) gain a
  `rho_floor: f64` argument. The `crate::consts::RHO_FLOOR` default is
  kept as a fallback and referenced in `NumericsSettings::default()`.

Estimated: ~150 lines across 5 files.

### Phase 3 — Deferred (not in scope this release)

- Knob #10 (entropy cutoff). Low physical impact; ignore until a user
  reports a cold-smearing regression.
- Any MXBA `AdaptiveBeta` internals (Section 9: explicitly out of
  scope).
- GPU workgroup / buffer-pool sizes (Section 6.3: not a physics knob).

### Phase 1 → Phase 2 ordering

Phase 1 lands first because its knobs are the ones production users hit
in sensitivity studies. Phase 2 is a pure threading exercise over
already-cleaned constants and can slot in any time after Phase 1 merges.

---

## Change-log (re-scope session)

- **2026-04-19** — EM: re-scope pass 2 (post-CFGN1). See Section 11
  below for current state. Priority lowered from medium to low; the
  parent stays open as an umbrella, but each remaining knob should
  land as its own small proposal when a user asks for it. Frontmatter
  `depends_on` cleared (CNST/DDUP/SIMP all landed earlier this week);
  complexity downgraded from large to medium (10 knobs left, none
  large individually).
- **2026-04-18** — Researcher: full re-scope against `origin/main`
  @ `0be9290`. Original inventory (Section 2) kept verbatim. Added
  Sections 6–10 as the authoritative census, phasing, and decisions.
  Census shrank from ~25 to ~12 exposed knobs after honoring NCFX
  (consolidated XC floors), MODR (file-path renames), CAST (invariant
  documentation), and MXBA (four-constant don't-expose call).

---

## Section 11 — Re-scope 2026-04-19 (supersedes Sections 6–10)

### What landed

- **CFGN1 — `initial_density.gaussian_sigma`** (PR #114, 2026-04-19).
  Knob #7 from Section 6.1. Proved the pattern: sub-struct under
  `Settings`, `Default` points at the canonical constant,
  `Settings::to_scf_params` threads it, validator rejects
  non-positive/non-finite, bit-identical when omitted from YAML.
  Four unit tests in `src/settings.rs`.

### What adjacent proposals cover (NOT CFGN's scope)

- **ECUT** — per-PP recommended `ecutwfc` from a PseudoDojo `.standard`
  table. ECUT is about **replacing** the hardcoded 204.09 eV default
  with a PP-aware default, not exposing a hardcoded knob. `ecutwfc`
  itself is already a Settings field (never hardcoded); ECUT only
  changes the default-value policy. See `proposals/ECUT-pp-recommended-ecut.md`.
- **ESPL** — split `ElectronSettings` + drop default `scf.max_iter`
  from 100 to 50. `max_iter` is already a Settings field (never
  hardcoded); ESPL changes its default. See
  `proposals/ESPL-electrons-settings-split.md`.

These two adjacent proposals re-tune defaults for *already-exposed*
settings; CFGN's remaining work is genuinely about exposing new knobs.

### What's actually left in CFGN's scope

Ten knobs, unchanged from Section 6 except for #7 (landed). Originally
split into two phases; post-CFGN1 they are better landed as individual
proposals when a user asks for one, because they're independent:

**Phase-1 territory (Fermi + iterative eigensolver; 5 knobs):**
- #4 `electrons.fermi_search.bounds_factor` (default 10.0)
- #5 `electrons.fermi_search.max_iter` (default 200)
- #6 `electrons.fermi_search.tol` (default 1e-14)
- #8 `scf.iterative_eigensolver.tol` (default ~2.8e-14)
- #9 `scf.iterative_eigensolver.max_restarts` (default 500)

**Phase-2 territory (Ewald + numerics floors; 5 knobs):**
- #1/#2 `ewald.cutoff_multiplier` (default 10.0)
- #3 `ewald.eta: Option<f64>` (default None = auto)
- #11 `electrons.xc.rho_floor` (default RHO_FLOOR = 1e-30)
- #12 `electrons.g2_zero_threshold` (default G2_ZERO_THRESHOLD = 1e-12)

**Out of scope** (moved from old Phase 3 into a firm "no"):
- MXBA `AdaptiveBeta` internals — Section 9 decision stands; algorithmic
  contract, not a user knob.
- GPU workgroup / buffer-pool sizes — Section 6.3; WGSL-constant
  coupling makes them not-trivially-exposable.
- Knob #10 (`entropy_cutoff`) — low physical impact; defer until a user
  reports a cold-smearing regression.

### Recommendation

Keep CFGN open as an umbrella, priority **low**. Each remaining knob is
small enough (≤ ~40 lines) to land as a standalone proposal when a real
user workflow demands it. Do not pre-emptively plumb all 10 — the YAML
surface gets bloated with defaults nobody overrides. The CFGN1 pattern
is the template: one knob per PR, validator-enforced, bit-identical
when omitted.

### Status of the original Phase plan (Section 10)

Sections 6–10 remain accurate as a census. The phasing in Section 10 is
now a **menu**, not a sequence: `Phase 1 → Phase 2` ordering no longer
applies because the knobs are individually tiny. Pull from either list
on demand.
