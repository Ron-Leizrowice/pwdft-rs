---
id: HD5I
status: archived
priority: low
complexity: large
risk: medium
depends_on: []
blocks: []
archived_on: 2026-04-19
archived_reason: Deferred since initial backlog (2026-04-16) with no user demand. Checkpoint/restart and structured output are real eventual needs, but not what research workflows are blocked on today. Archived to keep the active backlog focused; re-open when a concrete user workflow requires it (e.g., long-running geometry relaxation that must survive a crash).
---

# HD5I: HDF5 Restart and Structured Output — ARCHIVED

> **ARCHIVED 2026-04-19.** Deferred for 3 days with no user demand and no blocking dependency. Archived to keep the active backlog focused. Re-open as a fresh proposal (HDF5-rooted naming) when a real workflow needs restart or structured binary output — likely when MD or geometry-optimization lands.

> **Note:** Line numbers reference the pre-ScfContext codebase (src/scf/mod.rs was ~1127 lines, now ~709). Verify locations before implementing.

## Motivation

As calculations grow larger (more atoms, higher cutoffs, more k-points), two needs emerge:

1. **Restart capability:** If an SCF calculation is interrupted or needs to be continued with tighter convergence, the density `rho_g`, wavefunctions, and Fermi energy must be saved to disk. Currently there is no checkpoint mechanism.

2. **Structured output:** Band structure data is written as TSV (via `src/bandstructure.rs`). Eigenvalues, densities, and potentials for post-processing (band plots, DOS, charge density visualization) are better served by a self-describing format.

HDF5 is the standard binary format in computational physics. It provides:
- Efficient storage of large multi-dimensional arrays (densities on 64^3+ grids)
- Self-describing metadata (units, grid dimensions, crystal structure)
- Portable across platforms and languages (Python h5py, Julia HDF5.jl, C/Fortran)
- Compression for large datasets

## Dependencies

Add:
```toml
hdf5 = { version = ">=0.9", optional = true }
```

Behind a feature flag since HDF5 requires the system library:
```toml
[features]
hdf5 = ["dep:hdf5"]
```

System requirement: `brew install hdf5` (macOS) or `apt install libhdf5-dev` (Linux).

## Scope of Changes

### New file: `src/checkpoint.rs`

Handles saving and loading SCF state:

```rust
use hdf5::File;

pub struct Checkpoint {
    pub rho_g: Vec<Complex64>,
    pub eigenvalues: Vec<Vec<f64>>,
    pub fermi_energy: f64,
    pub iteration: usize,
}

impl Checkpoint {
    pub fn save(&self, path: &str, crystal: &Crystal, params: &ScfParams) -> Result<()> {
        let file = File::create(path)?;

        // Crystal structure metadata
        let crystal_group = file.create_group("crystal")?;
        crystal_group.new_dataset::<f64>()
            .shape((3, 3))
            .create("lattice")?
            .write(&lattice_array)?;

        // SCF state
        let scf_group = file.create_group("scf")?;
        // Store complex as interleaved [re, im, re, im, ...]
        scf_group.new_dataset::<f64>()
            .shape((self.rho_g.len(), 2))
            .create("rho_g")?
            .write(&rho_g_interleaved)?;

        scf_group.new_attr::<f64>().create("fermi_energy")?
            .write_scalar(&self.fermi_energy)?;
        scf_group.new_attr::<u64>().create("iteration")?
            .write_scalar(&(self.iteration as u64))?;

        Ok(())
    }

    pub fn load(path: &str) -> Result<Self> {
        let file = File::open(path)?;
        let scf = file.group("scf")?;
        // ... read datasets back ...
    }
}
```

### File: `src/scf/mod.rs`

Add checkpoint save inside the SCF loop (e.g., every 10 iterations or on convergence):

```rust
// After line 287 (convergence check):
#[cfg(feature = "hdf5")]
if iter % 10 == 9 || delta < params.conv_threshold {
    let ckpt = checkpoint::Checkpoint {
        rho_g: rho_g.clone(),
        eigenvalues: eigenvalues_all.clone(),
        fermi_energy,
        iteration: iter + 1,
    };
    ckpt.save("pwdft_checkpoint.h5", crystal, params)?;
}
```

Add checkpoint load at the start of `run_scf`:

```rust
// After line 196 (initial density):
#[cfg(feature = "hdf5")]
if let Some(ref ckpt_path) = params.checkpoint_path {
    if let Ok(ckpt) = checkpoint::Checkpoint::load(ckpt_path) {
        info!("Resuming from checkpoint at iteration {}", ckpt.iteration);
        rho_g = ckpt.rho_g;
        // reconstruct rho_r from rho_g via inverse FFT
    }
}
```

### File: `src/bandstructure.rs`

Add HDF5 output alongside TSV:

```rust
#[cfg(feature = "hdf5")]
pub fn write_hdf5(path: &str, kpoints: &[KPoint], eigenvalues: &[Vec<f64>]) -> Result<()> {
    let file = hdf5::File::create(path)?;
    // k-point coordinates: (n_kpts, 3)
    // eigenvalues: (n_kpts, n_bands)
    // labels: string dataset
}
```

### File: `src/input.rs`

Add optional `checkpoint` field to input TOML:

```toml
[scf]
checkpoint = "previous_run.h5"  # resume from this file
checkpoint_interval = 10        # save every N iterations
```

## HDF5 File Structure

```
pwdft_output.h5
├── crystal/
│   ├── lattice          (3, 3) float64 — lattice vectors in Angstrom
│   ├── positions        (n_atoms, 3) float64 — fractional coordinates
│   ├── species          (n_atoms,) string — element symbols
│   └── attrs: volume, n_atoms
├── scf/
│   ├── rho_g            (n_pw, 2) float64 — density Fourier coefficients [re, im]
│   ├── eigenvalues      (n_kpts, n_bands) float64
│   ├── occupations      (n_kpts, n_bands) float64
│   └── attrs: fermi_energy, total_energy, n_iterations, converged
├── kpoints/
│   ├── coordinates      (n_kpts, 3) float64
│   ├── weights          (n_kpts,) float64
│   └── labels           (n_kpts,) string (optional)
└── bands/ (if band structure)
    ├── k_distances      (n_kpts,) float64
    └── eigenvalues      (n_kpts, n_bands) float64
```

## Risks

- System dependency (libhdf5). Mitigated by making it optional behind a feature flag.
- HDF5 Rust bindings (`hdf5-rs`) can be tricky to compile. The crate auto-detects the system library via `pkg-config`.
- Checkpoint consistency: if the crystal structure or basis set changes between runs, a checkpoint from a previous run is invalid. The checkpoint file should store enough metadata to detect this.

## Expected Impact

- **Capability:** Restart interrupted calculations. Post-process results in Python without re-running.
- **Interoperability:** Standard HDF5 files readable by any scientific computing tool.
- **Development:** Makes it easy to save intermediate states for debugging convergence issues.
