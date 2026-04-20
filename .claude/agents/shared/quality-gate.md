# Quality gate

Every PR must pass `/quality-gate` before the EM will merge. That skill runs, in order, under the machine lock:

1. `cargo clippy -q --fix --allow-dirty --allow-staged --all-targets`
2. `cargo clippy -q --all-targets`
3. `cargo clippy -q --all-targets --features gpu`
4. `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps`
5. `cargo test` (Tier-1)

Both clippy invocations are required. Without `--features gpu`, the `pwdft/pwdft-core/src/gpu/` source tree and GPU-only test binaries are not linted and warnings accumulate silently.

## Tier-2 tests

If the PR touches `pwdft/pwdft-core/src/{scf,potential,symmetry,pseudopotential,eigensolver,gpu}/`, `basis.rs`, `fft.rs`, `ewald.rs`, `crystal.rs`, `kpoints.rs`, or bumps a numerics dep (faer / ndrustfft / nalgebra / ndarray) in `Cargo.toml`, also run `/test --tier2` and include the outcome in the PR body's Test Plan.

The Tier-2 run hits the pre-TSPL physics-blocker ignores (VGCH heavy-atom cells, Al ecut, C mixer, MXBA) that fail by design. The authoritative skip list lives in the `#[ignore]` reason strings in `pwdft/pwdft-core/tests/qe_validation.rs` and `pwdft/pwdft-core/tests/mxba_adaptive_beta_fe.rs` — report outcomes honestly, do not hide expected failures.

Doc-only, proposal-only, and lint-only PRs skip Tier 2.

## Don't

- Suppress warnings instead of fixing them. Refactor when `too_many_arguments` fires.
- `#[allow]` rustdoc warnings. Fix the prose.
- Skip a step because "it's not related." Every merge runs the full gate.
- Introduce new `unwrap()` or `panic!()` in production code paths.
