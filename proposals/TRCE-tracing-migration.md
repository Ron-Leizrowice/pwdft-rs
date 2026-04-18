---
id: TRCE
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# TRCE: Migrate logging from `log` + `env_logger` to `tracing`

## Problem

The project uses `log` + `env_logger` for emit/consume and `indicatif` for SCF progress bars. That stack works but has three limitations that the roadmap will bump into:

1. **No span-based timing.** To profile per-iteration cost (eigensolve vs FFT vs mixing vs XC) we either instrument by hand (touching hot loops) or run `cargo flamegraph` externally. `tracing` offers `#[tracing::instrument]` + `tracing-timing` / `tracing-flame`, which produce span timings from annotations without surgery at every call site. This matters directly for ongoing perf work (`ALOC F-7/F-12`, `GOPT`, `WFRX` Technique 2).
2. **No structured fields.** Current iteration log lines in `src/scf/report.rs` are hand-formatted strings. Anything downstream that wants to parse them (a validation dashboard, a regression-threshold check) re-parses text. `tracing` events carry typed fields (`delta = 1.2e-6`, `iter = 12`), so consumers read them structurally.
3. **Progress bar + log coexistence is ad-hoc.** `src/scf/driver.rs` and `src/scf/report.rs` mix `indicatif::ProgressBar` and `log::info!` without coordination — log lines can disrupt the bar. `tracing-indicatif` integrates spans with the progress bar so per-iter timing and the bar share one consistent output layer.

Current usage is modest — this is the right time to migrate, before any of the larger perf/validation efforts scale the log footprint:

| File | `log::info!`/`warn!`/`error!` sites | Direct `log::X` refs | `indicatif` refs |
|---|---|---|---|
| `src/main.rs` | 8 | 1 | 0 |
| `src/scf/driver.rs` | 8 | 5 | 3 |
| `src/scf/driver_spin.rs` | 6 | 2 | 0 |
| `src/scf/report.rs` | 21 | 1 | 2 |
| `src/scf/context.rs` | 6 | 1 | 0 |
| `src/scf/mixing/linalg.rs` | 1 | 1 | 0 |
| `src/gpu/mod.rs` | 3 | 1 | 0 |
| **Total** | **53** | **12** | **5** |

7 files, 53 emit sites, 5 indicatif sites. Not a rewrite of the codebase — a targeted dep swap plus annotations.

## Research

**Alternatives considered:**

| Approach | Verdict |
|---|---|
| `tracing` + `tracing-subscriber` (proposed) | Pure-Rust, de-facto standard for Rust observability. Keeps `log` compat layer so transitive deps (`faer`, `wgpu`, `ndrustfft` use `log`) continue to emit. Adds `#[instrument]` + span timing. |
| Stay on `log` + `env_logger`, add manual `Instant::now()` timers | Low cost, but every new hot path needs hand-rolled timing + a place to store it. Does not scale with perf work. |
| `slog` | Predates `tracing`; smaller ecosystem; no longer the idiomatic choice for new Rust code. |
| `defmt` | For embedded; not applicable. |

**Dependency footprint of the migration:**

Add:
- `tracing` (~50 kloc, de-facto standard, tokio-rs ownership)
- `tracing-subscriber` (with `env-filter`, `fmt` features — `RUST_LOG`-compatible filtering for drop-in replacement of `env_logger`)
- `tracing-indicatif` (~500 LoC, maintained; wraps spans with `indicatif::ProgressBar` per-span)

Remove: `env_logger`. Keep `log` in transitive dep graph via `tracing-log` feature of `tracing-subscriber` (zero-cost bridge).

Build-time cost: `tracing` + subscriber adds ~2–3 s to cold clean builds. Acceptable.

**`RUST_LOG`-compat keeps the user interface identical:**

`tracing-subscriber`'s `EnvFilter` reads `RUST_LOG=info,pwdft_rs::scf=debug` exactly like `env_logger`. No user-facing change at the CLI.

## Implementation

Four phases, each independently shippable.

### Phase A — dep swap (mechanical, 1 hr)

1. `Cargo.toml`: replace `log = ">=0.4"` and `env_logger = ">=0.11"` with:
   ```toml
   tracing = ">=0.1"
   tracing-subscriber = { version = ">=0.3", features = ["env-filter", "fmt", "tracing-log"] }
   tracing-indicatif = ">=0.3"
   ```
2. `src/main.rs`: replace `env_logger::init()` with a `tracing_subscriber::fmt()` builder that uses `EnvFilter::from_default_env()` and a `tracing_indicatif::IndicatifLayer`.
3. Touch every `use log::{info, warn, error, debug, trace};` and replace with `use tracing::{...}`. Macro call sites (`info!`, `warn!`, etc.) are already source-compatible — `tracing`'s macros are drop-in replacements for `log`'s at the same level.
4. Delete `use log;` and `log::X` qualified refs (12 sites — mechanical).

### Phase B — span instrumentation (medium, 2–3 hrs)

Add `#[tracing::instrument(skip_all, fields(...))]` on the hot-path entry points. This is where the profiling win lives.

Priority targets (in order of cost in a typical Si 4×4×4 SCF):

| Function | File | Fields to emit |
|---|---|---|
| `run_scf_unpolarized` | `src/scf/driver.rs` | `n_atoms`, `n_pw`, `n_kpts` |
| `run_scf_spin` | `src/scf/driver_spin.rs` | same + `nspin=2` |
| SCF iteration body (inline `scf_iter` span) | `src/scf/driver.rs`, `driver_spin.rs` | `iter`, `delta`, `de` |
| `diagonalize_lowest` | `src/eigensolver/dense.rs` | `n_pw`, `n_bands` |
| `build_hamiltonian_with_v_eff` | `src/scf/potentials.rs` | `n_pw` |
| `NonlocalPotential::add_to_hamiltonian` | `src/potential/nonlocal.rs` | `n_proj`, `n_pw` |
| `symmetrize_density_g` | `src/symmetry/density/g_space.rs` | `n_ops`, `n_grid` |
| `lda_xc` grid loop entry | `src/potential/xc.rs` | `n_grid`, `nspin` |
| `Mixer::mix` | `src/scf/mixing/mod.rs` | `history_size`, `mode` |

Rule of thumb: any function whose wall-time is > 1% of an SCF iteration gets an `#[instrument]`. Anything below, skip.

### Phase C — `IterationReport` → structured events (small, 1 hr)

`src/scf/report.rs:46` `log_iteration` currently builds a formatted string. Emit a `tracing::info!` event with structured fields (`iter`, `e_total`, `delta`, `hf_diff`, etc.) and move the human-readable formatting into a custom `fmt::Layer` for the default subscriber. Keep the `tracing-indicatif`-driven progress bar for CLI users; the event stream becomes machine-readable for any future validation-dashboard consumer.

### Phase D — update docs + env (small, 15 min)

1. CLAUDE.md: replace the implicit `env_logger` references with `tracing-subscriber`; note `RUST_LOG=...` still works.
2. `.claude/agents/*.md`: nothing to change (agent protocol doesn't reference log crate).
3. README / example: `RUST_LOG=info cargo run ...` continues to work — document nothing new.

## Verification

1. `cargo test` — all 265+ tests pass. Tests don't assert on log output, so the dep swap is invisible to them.
2. `cargo clippy -q --all-targets` and `cargo clippy -q --all-targets --features gpu` — clean.
3. `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` — clean.
4. **Runtime equivalence at the CLI:**
   ```bash
   RUST_LOG=info cargo run --release -- --input examples/si_scf.yaml 2>pre.log    # before rebase
   # apply TRCE, then
   RUST_LOG=info cargo run --release -- --input examples/si_scf.yaml 2>post.log
   diff pre.log post.log  # should be empty modulo line-number/timestamp formatting
   ```
   Acceptance criterion: same set of log lines emitted at `info` level (ordering and exact text may differ by timestamp format — check the content semantically, not byte-for-byte).
5. **New span-timing capability demonstrated:**
   ```bash
   RUST_LOG=info,pwdft_rs::scf=debug cargo run --release -- --input examples/si_scf.yaml
   # span_close events now show per-function wall-time
   ```
   Before TRCE this requires `cargo flamegraph` + sudo; after TRCE it requires an env var.
6. No performance regression: re-run `cargo bench --bench scf_benchmarks` before/after; the `#[instrument]` macro expansion costs ~30 ns per call at the default subscriber level, so a function called < 10k times per SCF adds < 0.3 ms total — well under benchmark noise.

## Open questions

- **Scope of Phase B:** the proposal lists 9 `#[instrument]` targets, but the threshold ("> 1% of an SCF iteration") means the final list should come from a profiler run, not from guesswork. Before Phase B lands, run `samply` or `cargo flamegraph` once and annotate only the functions that actually show up hot. Keeps the annotation footprint minimal.
- **Should Phase C split to a separate proposal?** Phase A+B are pure mechanical + observability. Phase C changes the iteration log format, which could affect any off-line tool that greps the log. No such tool exists today, but if one appears during review, split C out.

## Non-goals

- Not introducing OpenTelemetry / remote tracing export. Local-only spans.
- Not rewriting the progress bar. `indicatif` stays; `tracing-indicatif` just wraps it for span integration.
- Not changing the semantic level of existing log lines (info stays info, warn stays warn).
