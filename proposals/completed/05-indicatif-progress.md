# Proposal 05: Add indicatif for SCF progress display

> **Note:** Line numbers reference the pre-ScfContext codebase (src/scf/mod.rs was ~1127 lines, now ~709). Verify locations before implementing.

## Motivation

The SCF loop (`src/scf/mod.rs`, line 287) currently logs convergence via `log::info!`:

```rust
info!("SCF iter {}: E_fermi = {:.6} eV, delta_rho = {:.2e}", iter + 1, fermi_energy, delta);
```

This produces a wall of text for long calculations (50+ iterations). Users must scroll through log output to find the convergence trend, and there's no wall-time information.

`indicatif` provides progress bars and spinners that would give at-a-glance status:

```
SCF [████████░░░░░░░░] 12/50  E_f=-5.284 eV  Δρ=2.3e-04  [18s elapsed, ~25s remaining]
```

## Dependencies

Add:
```toml
indicatif = ">=0.17"
```

Lightweight, no transitive dependencies beyond `console` and `unicode-width`.

## Scope of Changes

### File: `src/scf/mod.rs` — `run_scf` (lines 144-329)

Add a progress bar before the SCF loop:

```rust
use indicatif::{ProgressBar, ProgressStyle};

let pb = ProgressBar::new(params.max_iter as u64);
pb.set_style(
    ProgressStyle::with_template(
        "SCF [{bar:30}] {pos}/{len}  E_f={msg}  [{elapsed} elapsed, {eta} remaining]"
    ).unwrap()
    .progress_chars("█░")
);
```

Inside the loop (after line 287), replace or supplement the `info!` call:

```rust
pb.set_position((iter + 1) as u64);
pb.set_message(format!("{fermi_energy:.4} eV  Δρ={delta:.1e}"));
```

On convergence (line 289):
```rust
pb.finish_with_message(format!("converged at {:.6} eV (Δρ={delta:.1e})", fermi_energy));
```

On failure (line 325):
```rust
pb.abandon_with_message("did not converge");
```

### Preserving log output

Keep the existing `info!` calls for when users pipe output to a file (`RUST_LOG=info cargo run ... > log.txt`). The progress bar is a visual supplement, not a replacement. `indicatif` detects non-terminal output and automatically falls back to plain text.

### Optional: Local potential precomputation progress

`compute_v_local_on_fft_grid` (lines 332-366) can take several seconds for large grids. A spinner during this phase would be useful:

```rust
let sp = ProgressBar::new_spinner();
sp.set_message("Computing V_local on FFT grid...");
// ... computation ...
sp.finish_with_message("V_local computed");
```

### Optional: Band structure progress

If band structure calculations iterate over many k-points, a progress bar over k-points would also be natural.

## Risks

- Minimal. `indicatif` is widely used and well-maintained.
- Progress bars can interfere with `env_logger` output. Fix by using `indicatif`'s `ProgressBar::println()` for log messages that should appear above the bar, or by using `indicatif::MultiProgress` if needed.

## Expected Impact

- **UX:** Users can see at a glance whether the calculation is converging, how many iterations remain, and wall-clock timing.
- **Performance:** Zero. Progress bar updates are ~microseconds.
- **Code:** ~15 lines of new code.
