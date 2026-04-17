---
id: QLN2
status: active
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# QLN2: Quality Lints — GPU + Benches Follow-up

## Origin

QLNT (PR #25, 2026-04-17) Code Reviewer flagged 4 lint hits outside its scope:

1. **`criterion::black_box` deprecation** — 19 sites across `benches/scf_benchmarks.rs` and `benches/gpu_benchmarks.rs`. Migrate to `std::hint::black_box`.
2. **`src/gpu/mod.rs:79`** — manually reimplements `div_ceil`. Candidate for `clippy::manual_div_ceil` enable + replace.
3. **`src/gpu/mod.rs:234`** — collapsible `if` (clippy `collapsible_if`).
4. **`tests/gpu_consistency.rs:312`** — `eprintln!("Gamma eigenvalues: {:?}", evs);` violates the `uninlined_format_args` lint enabled in CLIP. Slipped in because the GPU test binary wasn't part of CLIP's validation. Either fix here OR enable `--features gpu` in the clippy CI gate.

## Implementation

1. Replace all `criterion::black_box(x)` with `std::hint::black_box(x)` in `benches/`.
2. Enable `clippy::manual_div_ceil` in `Cargo.toml`. Replace `(a + b - 1) / b` pattern at `src/gpu/mod.rs:79` with `a.div_ceil(b)`.
3. Collapse the nested `if` at `src/gpu/mod.rs:234`.
4. Inline the format args at `tests/gpu_consistency.rs:312`.
5. Validate with both `cargo clippy -q --all-targets` and `cargo clippy -q --all-targets --features gpu`.

## Verification

- `cargo clippy -q --all-targets --features gpu` — clean
- `cargo bench --no-run` — benches still compile
- `cargo test --features gpu` — 231 pass, 9 ignored (matches QLNT baseline)

## Notes

- Trivially mergeable; <30 minutes Code Reviewer work.
- The CLIP-CI-gap is a process issue worth flagging separately to the EM. For now, fixing the offender locally is enough.
