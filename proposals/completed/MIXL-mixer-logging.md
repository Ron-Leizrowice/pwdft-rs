---
id: MIXL
status: active
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# MIXL: Mixer initialization & event logging

## Problem

The SCF density mixers (`Anderson`, `Broyden`, `PeriodicPulay`, plus the
`Kerker` preconditioner) make several runtime decisions that today leave
no trace in the SCF log. When an SCF run misbehaves (slow convergence,
limit cycle, surprise restart) the only observable about the mixer is
its iteration-by-iteration β line in `scf::report::log_iteration`. That
β line is the *output* of the mixer's decisions; it does not record the
*inputs* (mode, initial parameters) or *events* (DIIS history truncation,
adaptive-β triggers, Kerker auto-q_TF estimation).

Concrete gaps:

1. **No mixer-init log.** When SCF starts, nothing logs which mixer is
   active, its β, history depth, or whether Kerker is on. Compare to
   the existing `info!` lines in `src/main.rs:32-103` for crystal,
   basis, symmetry, and pseudopotentials — every other major SCF
   ingredient self-announces. The mixer is silent.

2. **Auto-q_TF chosen but not logged.** When the user writes
   `mixing_mode: kerker` in YAML without specifying `q_tf`, the code at
   `src/scf/mixing/anderson.rs:72` and `src/scf/mixing/broyden.rs:67`
   calls `auto_q_tf_squared(n_electrons, omega)` and uses the result
   silently. There is no way to see what q_TF the auto-estimate picked
   without running a debugger or adding a print and recompiling.

3. **DIIS history truncation is silent.** Anderson and Broyden bound
   history at `max_history`. When the buffer overflows, the oldest
   entry is dropped — invisibly. A pathological case where the mixer
   restarts its DIIS history every iteration looks indistinguishable
   from one that never restarts.

4. **Adaptive-β trigger reasons are silent.** `AdaptiveBeta::update`
   (`src/scf/mixing/mod.rs:134-168`) damps β when the residual grew
   (`ratio > growth_threshold`) and restores when it dropped over a
   `restore_window`. Only the *resulting* β is visible via the per-iter
   log. The MXBA Fe failure mode (β floors at β_min, restore never
   fires) is exactly the kind of thing a `debug!` line at each trigger
   would have made obvious — and is the kind of thing the next mixer
   tuning effort will need.

The single existing log statement in any mixer file is the
singular-DIIS-overlap fallback at `src/scf/mixing/linalg.rs:40`
(`log::warn!`). Keep it; the four gaps above are everything else.

The driver already has the plumbing for adaptive β
(`IterationReport.beta`, `SpinIterationFields.mag_beta`); the mixers
just need to surface the remaining state. Per PROF, the project stays
on `log` + `env_logger` for observability — this proposal adds content
into that layer.

## Implementation

Single PR, four small additions. None change behaviour — only log output.

### Step 1 — `Mixer::log_init` summary called once at SCF start

Add a method on `Mixer` (in `src/scf/mixing/mod.rs`) that emits a single
`info!` line describing the active configuration. Called from
`scf::driver::run_scf_unpolarized` and `scf::driver_spin::run_scf_spin`
after the mixer is constructed. Target output (one line):

```
Mixer: Anderson + Kerker  β=0.300  history=8  adaptive_β=off  kerker=on (q_TF=auto, est. 1.42 Å⁻¹)
```

Implementation sketch (illustrative — exact return-type plumbing left
to the implementer):

```rust
// in src/scf/mixing/mod.rs
impl Mixer {
    /// Emit a one-line `info!` describing mixer mode, β, history, and Kerker.
    pub(crate) fn log_init(&self, mode: &MixingMode, max_history: usize, adaptive: bool) {
        log::info!(
            "Mixer: {}  β={:.3}  history={}  adaptive_β={}  kerker={}",
            mode_label(mode),
            self.current_beta(),
            max_history,
            if adaptive { "on" } else { "off" },
            self.kerker_summary(),
        );
    }
}
```

`mode_label` is a free helper returning `&'static str` for the four
non-PeriodicPulay variants and an owned `String` (or `Cow<str>`) for
PeriodicPulay (which carries `period`). `kerker_summary` returns "off"
when no Kerker, "on (q_TF=user, X.YZ Å⁻¹)" when user-specified, and
"on (q_TF=auto, est. X.YZ Å⁻¹)" when auto-estimated.

`kerker_summary` needs access to the chosen q_TF. Today neither
`AndersonMixer` (`src/scf/mixing/anderson.rs:34-48`) nor `BroydenMixer`
(`src/scf/mixing/broyden.rs:30-46`) retains q_TF — they store the
precomputed `kerker_weights` Vec but discard the scalar. Add a
`kerker_q_tf: Option<f64>` field to both, populated alongside
`kerker_weights` in each constructor (~3 lines per mixer).

### Step 2 — `debug!` on DIIS history truncation

In `src/scf/mixing/anderson.rs::AndersonMixer::mix` and
`src/scf/mixing/broyden.rs::BroydenMixer::mix`, when `history_*.len() >
max_history` and the oldest entry is removed, emit one debug line:

```rust
log::debug!(
    "DIIS history at capacity ({max_history}); dropping oldest residual"
);
```

Since this is `debug!`, it's silent at default log level — but available
when `RUST_LOG=pwdft_rs::scf::mixing=debug`. Adds ~2 lines to each
mixer.

### Step 3 — `debug!` on adaptive-β triggers

In `src/scf/mixing/mod.rs::AdaptiveBeta::update`, when β changes, log
the trigger and direction:

```rust
// inside the `ratio > growth_threshold` arm:
log::debug!(
    "AdaptiveBeta: damp β {current_beta:.3} → {new_beta:.3} (residual ratio {ratio:.2} > {growth:.2})"
);
// inside the `restore_streak >= restore_window` arm:
log::debug!(
    "AdaptiveBeta: restore β {current_beta:.3} → {new_beta:.3} (after {window} consecutive ratios < {restore:.2})"
);
```

Adds ~6 lines. `debug!` not `info!` because in adverse conditions these
fire every iteration — info-level would flood the log.

### Step 4 — Wire `log_init` into both drivers

In `src/scf/driver.rs::run_scf_unpolarized` and
`src/scf/driver_spin.rs::run_scf_spin`, after `Mixer::new(...)`, call
`mixer.log_init(...)`. Spin driver calls it twice (once per channel) or
once with a "(spin)" annotation — either is fine; pick whichever yields
the cleaner line in the output.

## Verification

```bash
.claude/bin/machine-lock acquire "Core Engineer" "MIXL validation"
cargo test                                              # 265 still pass; no behaviour change
cargo clippy -q --all-targets                           # no new warnings
RUST_LOG=info cargo run --release --quiet -- --input examples/si_scf.yaml 2>&1 | grep '^.*Mixer:'
                                                         # → exactly 1 line
RUST_LOG=pwdft_rs::scf::mixing=debug cargo run --release --quiet -- --input examples/si_scf.yaml 2>&1 \
    | grep -c 'DIIS history\|AdaptiveBeta:'
                                                         # → 0 on a healthy SCF (no truncation, no β trigger);
                                                         # > 0 on a long or oscillating SCF
.claude/bin/machine-lock release
```

Acceptance:

- A single `Mixer:` line appears at SCF start showing mode, β, history,
  adaptive flag, Kerker status (and q_TF if Kerker is on).
- Auto-q_TF: feed `mixing_mode: kerker` with no `q_tf` and confirm the
  log line shows the estimated value.
- Adaptive-β triggers: run `tests/mxba_adaptive_beta_fe.rs` (the
  `#[ignore]`d Fe trajectory test) with `RUST_LOG=...=debug`; expect to
  see damp triggers fire repeatedly and no restore triggers — the
  failure-mode signal that drove MXBA's `default = false`.
- Behaviour unchanged: `cargo test` still passes bit-identical SCF
  outputs (every change is a log emission, no state mutation).

## Out of scope

- Per-iteration β lines (`scf::report::log_iteration` already does this
  when adaptive is on; MIXL does not duplicate or alter that path).
- Mixer-internal numerical diagnostics like DIIS condition numbers or
  Broyden's Jacobian rank. Useful but a different audit (`MAUD`-style).
- Surfacing mixer state in the final SCF summary
  (`scf::report::log_convergence_summary`). That summary is about
  energies, not mixers. A "mixer epilog" line could be added later if
  diagnostically useful.
