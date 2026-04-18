---
id: TRV2
status: proposed
priority: medium
complexity: medium
risk: low
depends_on: []
blocks: []
---

# TRV2: Fresh Test-Suite Review — Post-PCFX/CCMX/NCFX/GGAP Coverage Pass

## 1. Scope

### What TACC and TAUD already closed

- **TAUD (2026-04-17, completed)** — 5 silent-pass patterns, 7
  tolerance tightenings, 1 stale `#[ignore]`, convergence-variant
  specificity. All 5 fix-PRs landed.
- **TACC (2026-04-18, completed)** — 2 sibling silent-pass arms
  (`anderson.rs:548`, `broyden.rs:348`), 6 stale `#[ignore]` reasons in
  `qe_validation.rs` retargeted to SYKP/MPSH/VGCMP, 1 zip-code in
  `lapack_smoke`, deletion of `tests/fe_debug.rs`, ITEV `#[ignore]`
  resolution, `final_delta` comment hygiene.

### What TRV2 adds (post-TACC landings in scope)

GGAP Phase A (XcEvaluator dispatcher), MADOC Phase A (SCF energy
docstrings), TYPE-A (i32→i8/i16 narrowings), ERR2 P0 (lints + test
banners), XCNI (safety trap), plus CCMX/PCFX/NCFX areas TACC only
spot-checked.

Five directions not covered by TAUD/TACC:

1. Physics-path coverage for recent landings (NLCC beyond Si/Fe, GGAP
   driver integration, CCMX basis-change identity).
2. The MADOC Phase A band-sum identity
   `E_band = e_kin + e_loc + e_nl + 2·e_H + e_vxc` — documented but
   unpinned.
3. Silent-pass regressions since TACC-I (spot check — none found).
4. Brittle patterns TACC didn't audit (hardcoded `fft_grid` across 6
   call sites, watchdog iteration-count assertions).
5. Bench gaps: spin mixer, heavy-atom V_NL, eigensolver at ecut=600.

Methodology: surface-level sweep of `tests/` + in-src `#[cfg(test)] mod
tests` + `benches/`. All file:line references verified against
`origin/main` at commit `1f568be`.

## 2. Findings by category

| Category | Count | Most severe |
|----------|-------|-------------|
| 1 — True gap (physics not tested) | **5** | §1.1 MADOC band-sum identity |
| 2 — Silent pass (can't fail on claimed bug) | **1** | §2.1 VNLM GEMM path |
| 3 — Brittle (breaks on benign change) | **3** | §3.2 iteration-count watchdog |
| 4 — Organizational (duplication) | **2** | §4.1 12× `si_crystal()` sites |
| 5 — Bench gap (production-scale / new-physics) | **4** | §5.1 eigensolver caps at ecut=400 |
| **Total** | **15** | |

---

### Category 1 — True gaps

#### 1.1 MADOC band-sum identity undocumented by any assertion [highest severity]

MADOC Phase A added a 10-line derivation in
`src/scf/energy.rs:596-608` documenting:

```text
E_band = e_kinetic + e_local + e_nonlocal + 2·e_hartree + e_vxc
```

This identity is the invariant that breaks if any energy component gets
sign-flipped or scaled wrong. The existing self-check in
`tests/vgc5_per_component_si.rs:205-212` only pins
`Σ_components = total_energy` — trivially true because `total_energy`
*is* that sum. It does **not** pin the band-sum route.

A half-×`e_H` factor bug or sign-flipped `e_vxc` at the Hamiltonian
assembly site would change both `e_band` and `total_energy`
consistently under the current sum-identity check, evading VGC5. Only a
test that independently computes
`e_band_expected = e_kin + e_loc + e_nl + 2·e_H + e_vxc` and compares to
the reported `e_band` can catch this class.

**Fix:** new `tests/madoc_band_sum_identity.rs` with one test per
spin case (nspin=1 Si Γ-only, nspin=2 Fe BCC 4×4×4). Tolerance
~`1e-6·|E_band|` at `conv_threshold = 1e-8`. ~80 LOC.

#### 1.2 NLCC tested only on Si and Fe; 7 other LDA PPs with NLCC untouched

NCFX closed the NLCC unit/radial-weight bug universally in
`src/pseudopotential/upf/convert.rs:182-260`. NLCC Part A/B/C pinned
`ρ_core(G)` against Python/SciPy for Si and Fe only (4 test sites):

| File:line | Element |
|-----------|---------|
| `convert.rs:305` | Si G=0 |
| `convert.rs:332` | Si shell 1 |
| `convert.rs:369` | Fe G=0 |
| `convert.rs:395` | Fe shell 1 |

Every LDA PP with `core_correction="T"` is now blind-tested:

```
pseudopotentials/nc/lda/{Al, Co, Cr, Cu, Fe, Mn, Ni, Ti}.upf  all have NLCC
```

Cu (3s/3p/3d semicore edge case) and Mn (magnetic reference) would be
the highest-value additions. NCFX's bug class is unit/weight — it
could regress on one PP layout without touching Si/Fe.

**Fix:** extend `scripts/validate/rho_core_g_reference.py` to Cu + Mn
rows; add 4 test pins mirroring Si shell-1 pattern. ~25 LOC Rust +
~15 Python.

#### 1.3 CCMX basis-change identity (forward/inverse round-trip) is not unit-tested

The CCMX inline hot-loop at `src/scf/driver_spin.rs:537-577`:

```text
ρ_total = ρ↑ + ρ↓      m = ρ↑ − ρ↓
ρ↑ = (ρ_total + m)/2   ρ↓ = (ρ_total − m)/2
```

The round-trip invariant (`recombine(split(u, d)) == (u, d)` to machine
precision) is untested. `driver_spin.rs` has **zero** `#[test]`
(verified). The integration test
`tests/spin_polarization.rs:233::test_ccmx_fe_free_magnetization_converges`
catches only errors large enough to break convergence. An off-by-one
sign swap on the ρ↓ branch that converges to a wrong alternative state
would evade every current assertion.

**Fix:** extract `fn split_rho_m(rho_up, rho_down) -> (rho_total, m)`
and `fn recombine(rho_total, m) -> (rho_up, rho_down)` helpers into
`src/scf/mixing/mod.rs`; add `#[test]` with randomised 64-point input at
`relative_eq!(epsilon=1e-15)`. ~40 LOC.

#### 1.4 GGAP Phase A `XcEvaluator::Pz` driver-integration invariant is structural-only

GGAP's 6 unit tests (`src/potential/xc.rs:727-866`) pin the dispatcher
in isolation. `tests/spin_polarization.rs:446-504` pins driver-entry
dispatch (`run_scf` rejects Pbe/Pbe0/Hse06). But no test proves the
constructed `XcEvaluator::Pz` is **actually used** by the driver — a
refactor that silently reverts to a direct `lda_xc_grid` call would
produce identical LDA results today and only break when Phase B lands.

**Severity: low.** This becomes testable naturally in Phase B (any
bypass produces wrong PBE energies). Phase B tests cover it for free;
no Phase A action needed. Flagged here so the Phase B reviewer
remembers to verify the call path.

#### 1.5 GPU nspin=2 path has zero test coverage

`tests/gpu_consistency.rs` — 6 tests, zero `nspin: 2`. The GPU fast-path
in `scf::driver::eval_xc_with_gpu` is only wired for nspin=1; the spin
driver `scf::driver_spin::run_scf_spin` does not take a `&GpuAccelerator`
parameter at all (verified). Current state is either "GPU silently does
nothing for nspin=2" or "GPU path falls through to CPU" — the tests
cannot distinguish.

**Fix:** `tests/gpu_consistency.rs::test_gpu_vs_cpu_scf_nspin2` mirroring
`test_gpu_vs_cpu_scf_direct_comparison` (line 406) with `nspin = 2` +
`starting_magnetization = 0.1` on Si. Asserts totals match within 0.1 eV
and magnetizations within 1e-3 μB. ~70 LOC. Same graceful
CPU-only-host fallback as every other GPU test.

---

### Category 2 — Silent-pass patterns TACC missed

#### 2.1 VNLM `add_to_hamiltonian` GEMM path has no direct unit test

`src/potential/nonlocal.rs:317-340` assembles `H += B·D·B^H` via a
single `faer::linalg::matmul::matmul`. Correctness is pinned via
integration tests (`nonlocal_symmetry.rs`, `kb_projector_validation.rs`)
and indirectly via VGC5. A `matmul` argument-ordering swap (e.g. `B^H B`
instead of `B B^H`) would show up in the integration tests as "SCF
result weird" without pointing at the VNLM routine.

**Severity: medium.** Not a pure silent-pass — integration tests do
catch the bug, just slowly. A direct unit test would collapse the
root-cause chain.

**Fix:** `#[test] fn gemm_vs_nested_loop_reference()` in
`src/potential/nonlocal.rs::tests` with n_pw=12 toy basis, hand-built B
and D, bit-exact match. ~40 LOC.

---

### Category 3 — Brittle patterns

#### 3.1 Hardcoded `fft_grid: Some([16, 16, 16])` in 6 test sites

```
tests/parallel_consistency.rs:140 (test_scf_serial_vs_parallel)
tests/parallel_consistency.rs:288 (test_scf_kerker_serial_vs_parallel)
tests/spin_polarization.rs:70      (test_si_nspin2_matches_nspin1)
tests/spin_polarization.rs:187     (test_fe_spin_xc_consistency_regression)
tests/gpu_consistency.rs:58        (test_gpu_hartree_on_realistic_density)
tests/gpu_consistency.rs:563       (test_gpu_scf_kerker_converges)
```

**Severity: low.** All pair with `BasisSet::new(..., 100.0)` on Si
where automatic grid selection would also pick 16³. Defensive
documentation, not drift risk. Fold into Category 4.1 dedup.

#### 3.2 CCMX iteration-count watchdog (`n_iterations < max_iter`) is loose

`tests/spin_polarization.rs:329-333` asserts
`result.n_iterations < params.max_iter` where `max_iter = 80`.
Post-CCMX empirical is 14 iters. A future mixer change that slows
Fe CCMX convergence to 60-70 iters (still below 80) would pass
silently, masking a meaningful regression.

**Severity: low-medium.** The other three assertions in the same test
(|HF-KS| < 1e-3 eV, magnetization < 0.05 μB, `final_delta < 1e-2`)
catch the physics regressions. The iter count is efficiency only.

**Fix:** tighten to `r.n_iterations < 30` (≈2× the empirical 14). ~1
line. Alternative: document as watchdog with comment.

#### 3.3 GPU energy zip-codes — already in TAUD backlog

`tests/gpu_consistency.rs:262` (`fermi_energy > -5.0 && < 10.0`) and
`:302` (`5σ Fermi tail`) are TAUD findings 2.6 and 2.7. Verify those
PR-Cs landed; no new TRV2 work needed.

---

### Category 4 — Organizational

#### 4.1 `si_crystal()` in 12 sites; `fe_bcc()` in 3

```
src/potential/local.rs:87              src/symmetry/kpoints.rs:158
src/symmetry/detect.rs:224             src/symmetry/density/real_space.rs:130
src/symmetry/density/g_space.rs:291    src/scf/initial_density.rs:195
benches/scf_benchmarks.rs:35
tests/parallel_consistency.rs:20       tests/spin_polarization.rs:21
tests/nonlocal_symmetry.rs:21          tests/kb_projector_validation.rs:51
tests/gpu_consistency.rs:34
```

12 copies of the same `a = 5.431`, Si@(0,0,0), Si@(¼,¼,¼) fixture.
TACC finding 6 already flagged; post-TACC PRs (`mxba_adaptive_beta_fe.rs`,
`vgc5_per_component_si.rs`) inherited the pattern. Count ticks up.

**Fix:** shared `tests/common/fixtures.rs` + `src/testing.rs` behind
`#[cfg(any(test, feature = "testing"))]`. ~700 LOC net deletion.

#### 4.2 CCMX logic inline in driver; no `mixing/` test module covers it

Resolved as a byproduct of 1.3 — extracting `split_rho_m` / `recombine`
helpers into `src/scf/mixing/mod.rs` places CCMX logic where the rest
of the mixer API + tests live.

---

### Category 5 — Bench coverage

#### 5.1 `eigensolver` benches cap at ecut=400; production is ecut≈600 (n≈1363)

`benches/scf_benchmarks.rs:70`:
```rust
for &ecut in &[100.0, 200.0, 400.0] {
```

Line 69 comment says `ecut=600 crashes Accelerate due to libc++ TMO
bug`. **Stale** — we moved to faer, no Accelerate in stack. Faer's
`SelfAdjointEigen` handles 1363×1363 cleanly on M2.

Without the 1363 data point, any eigensolver regression visible only at
production scale (cache-blocking, bandwidth saturation) goes unnoticed,
blocking accurate ITEV-vs-Dense characterization when ITEV unblocks.

**Fix:** delete stale comment, add `600.0` to the ecut iteration. ~5 LOC.

#### 5.2 `hamiltonian` benches use only Si; heavy atoms (Fe, Cu) unbenched

VNLM GEMM-size scaling is 8–16× larger for Fe (s,p,d) or Cu
(s,p,d + semicore). VNLM's "single GEMM wins at heavy-atom scale"
claim has no measured numbers in the bench suite.

**Fix:** add `fe_pp()` / `fe_bcc()` fixtures and one
`vnl_apply_fe_n{n}` bench at ecut=400. ~25 LOC. Reuses fixture helper
from §4.1.

#### 5.3 CCMX coupled-channel mixer unbenched

Per iter at 72³ (`n_grid = 373_248`) the coupled-channel path allocates
4× `Vec<f64>` of size `n_grid` + 2 independent mixer invocations, ~12
MB of transient allocation. Allocator-pressure regressions invisible.

**Fix:** extract
`mix_coupled_channel(mixer_total, mixer_mag, up_in, down_in, up_out, down_out)`
helper (same helper §1.3 recommends for unit-testability); bench at
dims {18³, 36³, 72³}. ~80 LOC including helper + test + bench.

#### 5.4 NLCC `compute_core_density` unbenched

One-time setup per `ScfContext::new`, not per-iter. Matters only for
post-SCF band-structure scans with many k-points. Low priority.

**Fix:** defer to post-MODR Perf Eng profiling pass.

---

## 3. Recommended fix-PR sequence

Sequenced "biggest regression-exposure win first, cheapest second".
Each PR is independently landable.

### PR 1 — MADOC band-sum identity pin (§1.1)

New `tests/madoc_band_sum_identity.rs` (~80 LOC). Catches factor-2
Hartree bugs, `e_vxc` sign flips, and any compensating regression VGC5
misses. **Highest severity; no dependencies.**

### PR 2 — CCMX basis-change round-trip test (§1.3 + §4.2)

Extract `split_rho_m` / `recombine` helpers into
`src/scf/mixing/mod.rs`; add round-trip `#[test]`. ~40 LOC. Closes the
organizational gap (§4.2) as a byproduct and provides the helper PR 6
needs for the mixer bench.

### PR 3 — NLCC coverage extension (§1.2)

Extend `scripts/validate/rho_core_g_reference.py` to Cu + Mn; add 4
test pins. ~40 LOC total. One-shot, closes element-coverage concern
for VGCH/MPSH downstream.

### PR 4 — GPU nspin=2 consistency test (§1.5)

`test_gpu_vs_cpu_scf_nspin2` in `tests/gpu_consistency.rs`. ~70 LOC.
Graceful CPU-only-host fallback. Lower priority (GPU + nspin=2 is a
small slice).

### PR 5 — Fixture deduplication (§3.1 + §4.1)

Shared `tests/common/fixtures.rs` + `src/testing.rs`. Migrates all 12
`si_crystal()` sites + centralizes `fft_grid = [16, 16, 16]` rationale.
~700 LOC net deletion. Hygiene, no physics regression guarded — lands
after 1-4.

### PR 6 — Bench expansion (§5.1 + §5.2 + §5.3)

Eigensolver at ecut=600, heavy-atom Fe in hamiltonian bench, CCMX
mixer bench. ~100 LOC added. Depends on PR 2 helper + PR 5 fixtures
for ergonomics; can land standalone with inline fixtures if rushed.

---

## 4. Acceptance criteria for "test suite is complete post-audit"

1. Every post-2026-04-18 landed physics proposal has at least one
   regression guard:
   - MADOC identity: PR 1 ✓
   - CCMX basis change: PR 2 ✓
   - NLCC ≥ 3 elements: PR 3 ✓
   - GGAP dispatcher (6 unit + 1 driver-entry): already landed ✓
   - PCFX (5 unit tests incl. 18³ non-symmorphic): already landed ✓
   - TYPE-A narrowings (90°/3-fold/compose): already landed ✓
   - GPU nspin=2: PR 4 ✓
2. No hardcoded `[16, 16, 16]` fft_grid literal outside
   `src/testing.rs` / `tests/common/fixtures.rs`. (PR 5.)
3. Every bench group has ≥1 data point at n_pw ≥ 1000 or grid ≥ 64³.
   (PR 6 closes eigensolver; symmetry already covers 72³; xc_grid
   already covers 262k.)
4. Zero silent-pass `(Err, Err) => {}` or `if let Ok(..)` patterns in
   `tests/` or `src/`. (TAUD + TACC closed; TRV2 reconfirms no new
   ones landed post-TACC.)
5. Every `#[ignore]` reason cites a specific proposal ID or external
   blocker. (TACC refreshed 2026-04-18; verified.)
6. Both the `Σ_components = total_energy` identity (VGC5) and the
   band-sum double-counting identity (PR 1) hold to machine precision.
7. Fixture duplication ≤ 3 sites per crystal post-PR-5.

No acceptance criterion requires chasing VGCMP-blocked `#[ignore]`d QE
validation tests — tracked by VGCH/MPSH, not by TRV2.

---

## What this is NOT

- Not re-auditing TAUD/TACC findings.
- Not proposing a new test framework (proptest, snapshot, custom
  macros). Uses existing `#[test]` + `approx::relative_eq!` only.
- Not touching `qe_validation.rs` `#[ignore]` reasons (TACC refreshed).
- Not a re-review of XCNI (coverage verified adequate).
- Not proposing tests for unlanded work (HYBR, Phase-1 GPU spin).

## Open questions

- §1.4 GGAP driver integration — defer to Phase B, or add a
  `debug_assert!` in the driver now? Recommend defer.
- §5.1 eigensolver ecut=600 — verify on M2 that faer handles n=1363
  cleanly. If it crashes, update comment rather than remove.
- §3.2 watchdog — tighten to `< 30` or keep with a comment?
  Recommend tighten.

## Flagged for follow-up

- `benches/scf_benchmarks.rs:69` — stale Accelerate TMO comment
  (library no longer used). **Technical Writer** or folded into PR 6.
- `src/scf/driver_spin.rs:537-577` — basis change inlined; should be a
  testable helper. **Core Engineer** (PR 2).
- `src/potential/nonlocal.rs:317-340` — GEMM has no direct unit-test.
  **Core Engineer** (§2.1, low priority).
- `scripts/validate/rho_core_g_reference.py` — extend to Cu/Mn.
  **Researcher** (PR 3 Python half).

## References

- `proposals/completed/TACC-test-accuracy-audit.md` (closed 2026-04-18).
- `proposals/completed/TAUD-test-quality-audit.md` (closed 2026-04-17).
- `proposals/completed/MADOC-mathematical-documentation.md` Phase A —
  documented identity for §1.1.
- `proposals/completed/CCMX-coupled-channel-mixer.md` — basis change
  for §1.3.
- `proposals/completed/NCFX-nlcc-unit-and-radial-weight-fix.md` —
  universal fix tested only Si/Fe for §1.2.
- `proposals/completed/PCFX-symmetrize-rho-g-space.md` — well-tested
  (including 18³ non-symmorphic case); no gap.
- `proposals/completed/GGAP-gga-pbe-functional.md` Phase A — 6 unit
  tests adequate; §1.4 deferred.
- `proposals/completed/TYPE-A-integer-narrowings.md` — test-covered by
  existing compose/inverse tests.
- `.claude/logbooks/code-reviewer.md` — 2026-04-17/18 context.
