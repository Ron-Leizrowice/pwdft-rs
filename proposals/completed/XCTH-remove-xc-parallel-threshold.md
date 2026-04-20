---
id: XCTH
status: active
priority: low
complexity: trivial
risk: low
depends_on: []
blocks: []
---

# XCTH: Remove `XC_PARALLEL_THRESHOLD`; always use rayon in XC grid kernels

## Problem

`src/potential/xc.rs:63` defines `const XC_PARALLEL_THRESHOLD: usize = 16_384;`,
the only size-gated rayon dispatch in the entire `src/` tree. Both
`lda_xc_grid` (line 105) and `lda_xc_spin_grid` (line 263) branch on grid
size: below 16 384 points they use a hand-rolled sequential `for` loop with
manual `Vec::with_capacity`, above the threshold they fall through to a
`par_iter().unzip()` path. Each function therefore carries two parallel-but-
independent code paths, two docstrings explaining the threshold, and a
hard-coded constant calibrated to one specific machine (Apple M3 Max,
performance-engineer logbook entry 2026-04-17).

Two problems with the status quo:

1. **Local inconsistency.** Every other `par_iter` call site in the engine
   parallelizes unconditionally (verified by grep — `scf/density.rs:49`,
   `scf/energy.rs:329/382/427/518`, `scf/potentials.rs:41`,
   `scf/context.rs:122/136`, `scf/driver.rs:381`,
   `scf/driver_spin.rs:270/284`). Only XC has a size gate. New
   contributors reading these files have to reason about why XC is
   special; the answer is "it isn't, we just happened to micro-benchmark
   this one."

2. **The constant is hardware-specific.** `XC_PARALLEL_THRESHOLD = 16_384`
   was picked on Apple M3 Max where rayon's fork/join costs ~70 µs per
   region. A different machine (Linux x86 with a different `RAYON_NUM_THREADS`,
   a future M-series chip, or any CI runner) has a different crossover.
   The constant masquerades as a portable tuning parameter; it isn't one.

The original justification is in the docstring (`src/potential/xc.rs:48-62`)
and is honest about the tradeoff:

| n       | sequential | parallel | net           |
|---------|------------|----------|---------------|
| 4 096   | 43 µs      | 113 µs   | +161% (worse) |
| 32 768  | 407 µs     | 195 µs   | −52%          |
| 262 144 | 3 035 µs   | 459 µs   | −85%          |

The crossover is real — at n=4 096 the unconditional rayon path is ~70 µs
slower than sequential. The question is whether 70 µs per call matters at
SCF scale. It doesn't.

## Research

### How often is the small-grid path actually hit in production?

`lda_xc_grid` is called once per SCF iteration from `src/scf/driver.rs`
(non-spin) and `lda_xc_spin_grid` is called once per iteration from
`src/scf/driver_spin.rs`. A typical SCF runs 15–30 iterations.

Production FFT grids are 24³ = 13 824 to 48³ = 110 592 points (CLAUDE.md
§ Architecture quotes "typical production FFT grids are 24³–48³"). The
24³ case sits *just below* the 16 384 threshold and currently takes the
sequential path; 32³ and above take the parallel path.

If we always parallelize, the 24³ case pays the ~70 µs fork/join overhead
on every XC call. Worst case impact on a 30-iteration spin SCF:

```text
2 calls/iter × 30 iter × ~70 µs/call ≈ 4.2 ms per SCF
```

At 24³, a single SCF run is on the order of seconds (eigensolve dominates
at >100× the XC cost — see PERF-2026-04-18 headline numbers: V_NL apply +
eigensolve = 75 ms per k-point per iter at n_pw=725). 4.2 ms over 30
iterations is below measurement noise.

### What about benchmarks that explicitly target small grids?

`benches/scf_benchmarks.rs:298-322` exercises `lda_xc_grid` and
`lda_xc_spin_grid` at n ∈ {256, 512, 4 096, 16 384, 32 768, 262 144}. The
n=256, n=512, and n=4 096 cases will regress visibly:

| Bench (current → after XCTH)            | Expected change          |
|------------------------------------------|--------------------------|
| `lda_xc_grid_n256`                       | 2.7 µs → ~70 µs          |
| `lda_xc_grid_n512`                       | 5.4 µs → ~70 µs          |
| `lda_xc_grid_n4096`                      | 42.7 µs → ~110 µs        |
| `lda_xc_grid_n16384` and above           | unchanged (already parallel) |

This is a real regression *in the bench*, but it does not correspond to a
real regression in any user-facing SCF — no production input file produces
a 256-point or 512-point FFT grid (those would be ~2 Å on a side, smaller
than a single atom). The n=4 096 case corresponds to ~16³, which would
require ecut < 30 eV on a small unit cell — well below the SCF defaults.

Two ways to handle the bench regression:

a. **Drop the n=256, n=512, n=4 096 bench cases.** They exist to track
   sequential XC performance at sizes that don't appear in production.
   Keep n=16 384 and above.
b. **Keep them as-is and accept the regression.** They would still be
   useful as a "rayon overhead floor" reference but would no longer be
   labeled as production-relevant.

I recommend (a) — they serve no validation purpose and would otherwise
trigger a benchmark-pass investigation every quarter.

### Alternative considered: keep the threshold, drop the calibration

We could keep a threshold but set it conservatively (e.g., 1 024 — below
any production grid) just to shield pathological micro-cases. Rejected: it
still leaves two code paths, two docstrings, and a magic number to
maintain. The whole point of the proposal is to delete code, not to
re-tune a constant.

### Alternative considered: rayon's own adaptive scheduling

rayon has `with_min_len(N)` on parallel iterators that sets a per-chunk
minimum. We could call `par_iter().with_min_len(8192)` and let rayon decide
whether to spawn workers. Rejected: it doesn't actually skip the par_iter
machinery for tiny inputs (the per-region overhead is the bulk of the
fork/join cost, not the per-chunk dispatch), and it introduces a
different magic number with worse semantics. Cleaner to just commit to
one parallel code path.

## Implementation

### Step 1 — `src/potential/xc.rs`

Delete the threshold constant and the sequential branch in both functions.
Trim the docstrings.

```rust
// DELETE lines 41-63 (the threshold doc-block + const).

// REPLACE lda_xc_grid (current lines 95-124):
/// Compute ε_xc(r) and V_xc(r) on a real-space grid.
///
/// `rho_r`: electron density on real-space grid (e/ų).
///
/// Returns (exc_r, vxc_r): energy density and potential on the grid (eV).
///
/// Each point is an independent evaluation, so this is parallelized
/// unconditionally via rayon — consistent with every other grid kernel
/// in the engine.
pub fn lda_xc_grid(rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    rho_r
        .par_iter()
        .map(|&rho| {
            let xc = lda_xc(rho);
            (xc.exc, xc.vxc)
        })
        .unzip()
}

// REPLACE lda_xc_spin_grid (current lines 256-293):
/// Spin-polarized XC on a real-space grid.
///
/// Returns (exc_r, vxc_up_r, vxc_down_r) in eV.
///
/// rayon's `unzip` handles only 2-tuples, so we unzip into an
/// ((exc, vxc_up), vxc_down) shape and then flatten.
pub fn lda_xc_spin_grid(
    rho_up_r: &[f64],
    rho_down_r: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    debug_assert_eq!(
        rho_up_r.len(),
        rho_down_r.len(),
        "spin channels must share grid size",
    );

    let ((exc, vxc_up), vxc_down): ((Vec<f64>, Vec<f64>), Vec<f64>) = rho_up_r
        .par_iter()
        .zip(rho_down_r.par_iter())
        .map(|(&ru, &rd)| {
            let xc = lda_xc_spin(ru, rd);
            ((xc.exc, xc.vxc_up), xc.vxc_down)
        })
        .unzip();

    (exc, vxc_up, vxc_down)
}
```

### Step 2 — `benches/scf_benchmarks.rs`

Trim the small-grid bench cases that no longer exercise distinct code paths:

```rust
// scf_benchmarks.rs:301
- for &n in &[256_usize, 512, 4_096, 16_384, 32_768, 262_144] {
+ for &n in &[16_384_usize, 32_768, 262_144] {
```

### Step 3 — `proposals/CFGN-configurable-numerics.md`

Remove the `XC_PARALLEL_THRESHOLD` row from the constants table at line 314
(the constant no longer exists, so it can't be configurable).

### Step 4 — Logbook note

Performance Engineer logbook (`.claude/logbooks/performance-engineer.md`)
gets a one-line entry recording the XCPR-era calibration's removal so the
2026-04-17 measurement is not re-discovered as a "regression" later.

## Verification

```bash
.claude/bin/machine-lock acquire "Core Engineer" "XCTH validation"
cargo test                                                # all 265 pass
cargo test --features gpu                                 # GPU path unchanged (lda_xc_grid is CPU-side)
cargo clippy -q --all-targets                             # no new warnings
cargo clippy -q --all-targets --features gpu              # ditto
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps            # docstring trim is clean
cargo bench --bench scf_benchmarks -- xc_grid             # confirm n>=16k unchanged
.claude/bin/machine-lock release
```

Acceptance:

- All tests pass with no physics drift (the kernel is bit-identical;
  `par_iter().unzip()` preserves order in rayon).
- Existing `xc_grid/lda_xc_grid_n16384` and larger benches stay within
  noise of their PERF-2026-04-18 numbers (162 µs / 193 µs / 433 µs).
- A representative SCF (e.g., `examples/si_scf.yaml`) shows no measurable
  wall-time change — predicted impact is sub-millisecond on top of a
  multi-second SCF.

## Out of scope

- No changes to GPU paths. `gpu/shaders/xc.wgsl` already runs every grid
  point in parallel; no threshold there to remove.
- No changes to non-XC parallel call sites. They already parallelize
  unconditionally — this proposal brings XC in line with them, not the
  other way around.
