---
id: LOGH
status: active
priority: low
complexity: trivial
risk: low
depends_on: []
blocks: []
---

# LOGH: Logging hygiene — eprintln cleanup (MELG + TXEP + PCEP)

Bundle of three trivial cleanups around `eprintln!` use in production and
test code. None of these add or remove behaviour; they just convert
unstructured stderr writes to either `log` macros (production) or
deletion (test diagnostics that the surrounding assertions already pin).
Bundled because each is one-touch and the failure mode is the same:
unstructured stderr that shouldn't be there.

## Problem

### MELG — `main.rs` SCF result reporting bypasses `log`

`src/main.rs:133-143` writes the SCF summary via five raw `eprintln!`
calls:

```rust
eprintln!("SCF converged in {} iterations", result.n_iterations);
eprintln!("Total energy: {:.6} eV", result.total_energy);
eprintln!("Fermi energy: {:.6} eV", result.fermi_energy);
for (ik, evs) in result.eigenvalues.iter().enumerate() {
    if ik < 3 || ik == result.eigenvalues.len() - 1 {
        eprintln!(
            "  k-point {ik}: bands = {:?}",
            evs.iter().map(|e| format!("{e:.4}")).collect::<Vec<_>>()
        );
    }
}
```

The same binary already calls `env_logger::init()` (line 26) and uses
`log::info!` for every other status message (lines 32-122). The SCF
summary is the most important post-run output but it is the *only*
status emitted unfiltered. Two consequences:

1. Users running `RUST_LOG=warn` (or filtering noise) still get the SCF
   summary spam. There's no way to silence it.
2. Users running `RUST_LOG=info,pwdft_rs::scf=debug` to capture detailed
   SCF traces lose the summary's stream consistency — `info!` lines go
   to stderr through the env_logger formatter, the SCF summary skips
   the formatter.

### TXEP — Debug `eprintln!` in `kpoints.rs` test sticks around forever

`src/symmetry/kpoints.rs:177,182-186,188` has three `eprintln!` calls
inside `test_si_4x4x4_reduces_to_8`:

```rust
eprintln!("n_ops: {}", symmetry.n_ops);
// ... loop ...
eprintln!(
    "IBZ k{i}: w={:.6} k=({:.4},{:.4},{:.4})",
    kp.weight, kp.k.x, kp.k.y, kp.k.z
);
// ...
eprintln!("total weight: {total_w}");
```

These are diagnostic prints from the SYKP audit (closed). The test now
asserts `ibz.len() <= 10 && ibz.len() >= 8` (line 204), which already
captures the conclusion the prints were guiding the developer toward.
Every `cargo test` run dumps ~12 lines of unfiltered stderr that nobody
reads.

### PCEP — Debug `eprintln!` residue across `pseudopotential/` tests

13 `eprintln!` macro calls in pseudopotential test code:

| File                                          | Lines (one entry = one `eprintln!`)              |
|-----------------------------------------------|--------------------------------------------------|
| `src/pseudopotential/mod.rs`                  | 247, 256, 272, 281                               |
| `src/pseudopotential/upf/convert.rs`          | 239, 313, 341, 377, 404, 450, 476, 514, 540      |

Two patterns:

a. **"Skipping" diagnostics** (mod.rs lines 247, 272): printed when a PP
   has no `rho_atom` and the test early-returns. The early `return;`
   already documents the no-op path; the print is dead documentation.
b. **"Computed value" diagnostics** (the rest): print a quantity that
   the next 3-5 lines of the test assert against a reference value with
   `assert_relative_eq!` or `assert!((x - ref).abs() < tol)`. The
   assertion already pins the number; the print just spams stderr on
   green runs and adds no information on red runs (the assertion's
   panic message already includes the value).

Examples:

```rust
// upf/convert.rs:239 — print value that line 242-245 asserts against
eprintln!("Si partial core charge = {q_core:.6} e");
assert!(
    (q_core - 0.74).abs() < 0.05,
    "Si partial core charge = {q_core:.6} e (expected ≈ 0.74 e)"
);

// upf/convert.rs:313 — print value, assert is in the next 5 lines
eprintln!("Si ρ_core(G=0) = {rho_g0:.6e} e/Å³  (ref 1.8476e-2)");
```

These came in during the NCFX/NLCC validation push when the developer
was iterating on numerical tolerances and wanted to *see* the numbers
move. Once the tolerances were fixed and the asserts pinned, the
prints became noise.

## Implementation

### Step 1 — MELG: replace `eprintln!` with `log::info!` in `src/main.rs`

```rust
// src/main.rs:133-143 — replace eprintln! with info!
info!("SCF converged in {} iterations", result.n_iterations);
info!("Total energy: {:.6} eV", result.total_energy);
info!("Fermi energy: {:.6} eV", result.fermi_energy);
for (ik, evs) in result.eigenvalues.iter().enumerate() {
    if ik < 3 || ik == result.eigenvalues.len() - 1 {
        info!(
            "  k-point {ik}: bands = {:?}",
            evs.iter().map(|e| format!("{e:.4}")).collect::<Vec<_>>(),
        );
    }
}
```

The `use log::info;` import on line 4 already covers this — no Cargo or
import changes.

Decision point: should the SCF summary be `info` or stay highly visible
(e.g., `warn`)? Recommend `info` for consistency with the surrounding
status output. Users running with the default env_logger config
(`RUST_LOG` unset) currently see only `error` — but the existing
crystal/basis/symmetry status lines on lines 32-103 are already `info`,
so the user experience pre/post change is consistent: either the user
asks for `info` (and gets everything) or they don't (and get nothing
besides errors). The pre-LOGH inconsistency is the bug.

### Step 2 — TXEP: delete debug prints in `kpoints.rs::test_si_4x4x4_reduces_to_8`

```rust
// src/symmetry/kpoints.rs — delete these three blocks:
// line 177:  eprintln!("n_ops: {}", symmetry.n_ops);
// lines 181-186: the for-loop with eprintln!("IBZ k{i}: ...")
// line 188: eprintln!("total weight: {total_w}");
```

The assertion on line 204-207 still pins the IBZ count to [8, 10]. If a
future change breaks the count, the assertion's panic message already
contains `ibz.len()` — no diagnostic loss.

### Step 3 — PCEP: delete or fold debug prints in pseudopotential tests

For each of the 13 sites:

- **`src/pseudopotential/mod.rs:247, 272`** ("skipping" prints):
  delete. The early `return;` is self-documenting.
- **`src/pseudopotential/mod.rs:256, 281`** ("integral / z_valence",
  multi-line `eprintln!` blocks): delete. The next assert (lines
  260-264, 286-290) panics with the same data on failure.
- **`src/pseudopotential/upf/convert.rs:239, 313, 341, 377, 404, 450, 476, 514, 540`**:
  delete. Each is followed within 5 lines by an `assert!` whose panic
  message includes the same values (most use the `eprintln!` format
  string verbatim, e.g. line 244 already echoes `q_core` in its panic).

If any of the assert messages is *missing* the printed value (rather
than just printing a less-detailed version), fold the value into the
assert message before deleting the print — don't lose information just
to delete a line. Spot-check on a per-site basis during the PR.

## Verification

```bash
.claude/bin/machine-lock acquire "Core Engineer" "LOGH validation"
cargo test                                              # all 265 pass; no behaviour change
cargo test -- --nocapture 2>&1 | grep -c '^IBZ k\|^Si \|^Fe \|^Cu \|^Mn \|^total weight'
                                                         # → 0 (was ~30 lines pre-LOGH)
cargo run --release --quiet -- --input examples/si_scf.yaml 2>/dev/null | wc -l
                                                         # SCF output gone from stderr-default
RUST_LOG=info cargo run --release --quiet -- --input examples/si_scf.yaml 2>&1 | grep -c 'SCF converged\|Total energy\|Fermi energy'
                                                         # → 3 (visible at info level)
cargo clippy -q --all-targets                           # no new warnings
.claude/bin/machine-lock release
```

Acceptance:

- `grep -rn 'eprintln!' src/` returns no matches in `pseudopotential/`,
  `symmetry/kpoints.rs`, or `main.rs` (other matches in `src/` may
  remain — those are out of scope; this proposal addresses 21
  enumerated sites: MELG=5, TXEP=3, PCEP=13).
- A `cargo run` with `RUST_LOG` unset produces no `SCF converged` /
  `Total energy` lines on stderr.
- The same run with `RUST_LOG=info` produces exactly the same SCF
  summary content as pre-LOGH.

## Out of scope

- Other `eprintln!` sites in `src/` (if any survive — not all of `src/`
  was audited for this proposal). A general "no `eprintln!` in `src/`"
  lint is a separate item.
- The SCF summary's content or format. Only the *transport* changes.

## Note on stack direction

LOGH is a terminal `eprintln!` → `log` conversion. Per PROF,
`log` + `env_logger` is the long-term observability layer for the
project.
