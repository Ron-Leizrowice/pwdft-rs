# Quality gate

Every PR must pass these checks before the EM will merge. Run them locally under the machine lock:

```bash
.claude/bin/machine-lock run "<role>" "quality gate" -- bash -c '
  cargo clippy -q --fix --allow-dirty --allow-staged --all-targets &&
  cargo clippy -q --all-targets &&
  cargo clippy -q --all-targets --features gpu &&
  RUSTDOCFLAGS="-D warnings" cargo doc --no-deps &&
  cargo test'
```

## Why both clippy invocations

Without `--features gpu`, the `pwdft/pwdft-core/src/gpu/` tree and the GPU-only test binaries are not linted and warnings accumulate silently. Both runs are required.

## Rustdoc

`RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` is part of the gate. Fix broken intra-doc links, unescaped brackets, and links at private items by editing the prose — do not `#[allow]` rustdoc warnings. The `-D warnings` flag must travel through `RUSTDOCFLAGS`; current cargo rejects `cargo doc -- -D warnings`.

## Tier-2 tests

Any PR that touches `pwdft/pwdft-core/src/scf/`, `potential/`, `symmetry/`, `pseudopotential/`, `eigensolver/`, `basis.rs`, `fft.rs`, `ewald.rs`, `gpu/`, `crystal.rs`, `kpoints.rs`, or bumps a numerics dep (faer / ndrustfft / nalgebra / ndarray) in `Cargo.toml` must also run:

```bash
.claude/bin/machine-lock run "<role>" "tier-2" -- cargo test -- --ignored
```

Report the outcome in the PR body. Doc-only, proposal-only, and lint-only PRs skip Tier 2. Known pre-TSPL physics-blocker ignores (VGCH heavy-atom cells, Al ecut, C mixer, MXBA) fail by design; the authoritative skip list lives in the `#[ignore]` reason strings in `pwdft/pwdft-core/tests/qe_validation.rs` and `pwdft/pwdft-core/tests/mxba_adaptive_beta_fe.rs`.

## Don't

- Suppress warnings instead of fixing them. Refactor when `too_many_arguments` fires.
- Skip a quality-gate step because "it's not related." Every merge runs the full gate.
- Introduce new `unwrap()` or `panic!()` in production code paths.
