---
id: TSPL
status: active
priority: high
complexity: low
risk: low
depends_on: [TPRF]
blocks: []
---

# TSPL: Bifurcate the test suite — fast-default tier vs. heavy-opt-in tier

## Problem

The current `cargo test` runs every test, including heavy integration SCF calculations that dominate wall time. Observed on the M3 Max (2026-04-19, default features, `[profile.test] opt-level=0` baseline before TPRF):

| Test binary | Wall time |
|---|---|
| `tests/vgc5_per_component_si.rs` | ~6 min 38 s |
| `tests/qe_validation.rs` | ~2 min 33 s |
| `tests/spin_polarization.rs` | ~1 min 10 s |
| `tests/nonlocal_sym.rs` | ~39 s |
| all other binaries combined | < 30 s |

Total: **~11 minutes** for the full gate, of which **~10 minutes** is four heavy SCF suites. TPRF (`[profile.test] opt-level=3`) is expected to shrink each of those by 10–50× — but even post-TPRF, the heaviest SCF suites will dominate, and every new SCF regression test will be paid for on every PR gate.

User directive (2026-04-19): "we should definitely bifurcate the test suite between quick unit tests for routine code-work, and heavy full SCF calculations that only run when genuinely relevant to the changes made."

## Proposal

Split the test suite into two tiers, enforced by the `#[ignore]` attribute with a clear reason string:

**Tier 1 — routine (`cargo test`):**

- All in-`src/` unit tests.
- Lightweight integration tests: free-electron bands, KB projector, parallel consistency, ScfParams validation, symmetry detection.
- Any SCF-containing test that runs under ~5 s wall post-TPRF.
- Target budget: **under 60 s wall for the full tier**.

**Tier 2 — heavy (`cargo test -- --ignored`):**

- Full SCF validation against QE (`tests/qe_validation.rs`).
- `vgc5` per-component + MADOC band-sum identity on Si and Fe (SCF-heavy, conv≤1e-8).
- Fe spin-polarization integration tests with tight convergence.
- Anything that runs > 5 s wall post-TPRF.
- Gate-label: `#[ignore = "TSPL tier-2: heavy SCF, run with --ignored when touching SCF/density/mixing/XC/NLCC/symmetry/eigensolver/GPU code"]`

A Tier-2 test's ignore reason must name the code paths whose changes warrant running it. Example: an ecut-convergence sweep would name `basis.rs`, `fft.rs`, `scf/grid.rs`.

## Running each tier

- `cargo test` — Tier 1 only (default-fast).
- `cargo test -- --ignored` — Tier 2 only.
- `cargo test -- --include-ignored` — both tiers.

## When Tier 2 is "genuinely relevant"

| Touched path | Tier-2 run? |
|---|---|
| `src/scf/**`, `src/potential/**`, `src/symmetry/**`, `src/pseudopotential/**`, `src/eigensolver/**`, `src/basis.rs`, `src/fft.rs`, `src/ewald.rs`, `src/gpu/**` | **Yes, mandatory** |
| `src/crystal.rs`, `src/kpoints.rs` (geometry/lattice changes) | Yes, mandatory |
| `Cargo.toml` dep bumps on faer/ndrustfft/nalgebra/ndarray | Yes |
| Dev-only files: `.claude/`, `proposals/`, `README.md`, bench code, QE reference data (`qe_validation/*.in|*.out|*.toml`) | No |
| Docs-only (docstring edits, `CLAUDE.md`) | No |
| Cargo.toml lint-only changes | No |

Plan: embed this table verbatim in `CLAUDE.md § Testing`. Reviewer is responsible for confirming Tier-2 ran when the PR's diff touches a qualifying path.

## Future: CI-side enforcement

Out of scope for this proposal, but the natural next step is a GitHub Actions nightly that runs `cargo test -- --include-ignored`, plus a per-PR workflow that runs `--ignored` conditionally on path globs. Tracked as a follow-up.

## Implementation plan

1. Land TPRF first (so the Tier-1 budget of "under 60 s" is realistic).
2. Measure each integration-test binary's post-TPRF wall time on the M3 Max.
3. Tag Tier-2 tests with `#[ignore = "TSPL tier-2: ..."]` + the triggering paths list.
4. Update `CLAUDE.md § Testing` with the two-tier policy + the triggering-paths table.
5. Update the PR template (if one exists) to include a "Tier-2 run?" checkbox.
6. Verify default `cargo test` wall time is under the 60 s budget on the M3 Max.

## Non-goals

- Not proposing a custom test runner, test-harness rewrite, or CI-side nightly today.
- Not proposing per-test timing assertions — runtime regressions are noisy.
- Not moving ignored tests to a separate binary crate (`tests/heavy/*.rs`) — that restricts `--include-ignored` ergonomics for no clear gain.

## Open questions

- **Where does `tests/qe_validation.rs` sit?** It already has per-test `#[ignore]`s for the unblocked-but-failing subset (Si E_F, C, Al, Fe). The whole binary should move to Tier 2 — we want QE match on every SCF-path change, not every commit.
- **Does `tests/gpu_consistency.rs` belong in Tier 2?** It's only ever compiled with `--features gpu`, and today's default gate skips it. For consistency, yes — `#[ignore]` inside the gpu feature gate.
- **Does `tests/wfrx_subspace_consistency.rs` belong in Tier 2?** WFRX Technique 1 runs two SCFs back-to-back (with and without subspace projection). Almost certainly Tier 2.

## Risk

- **Tier 2 silently breaks between runs.** Mitigation: the nightly CI workflow (tracked as follow-up). Absent that, EM discipline is the backstop.
- **Developer ergonomic regression:** contributors forget to run Tier 2 when they should. Mitigation: `CLAUDE.md` documents the table; PR template prompts; EM catches on review.
