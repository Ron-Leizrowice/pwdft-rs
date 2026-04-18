---
id: TYPE
status: completed
priority: low
complexity: small-medium
risk: low
depends_on: []
blocks: []
---

## Status (2026-04-18)

**Phase A landed** as PR `TYPE-A: i32→i8 SpaceGroupOp rotations + i32→i16 BasisSet Miller + dead index_map`. All three items in the ranked-implementable table were implemented as three commits.

**Measured wins (Apple M3 Max, machine lock held):**

| Bench | Before | After | Δ |
|-------|--------|-------|---|
| `symmetrize_density_g_n18_ops48` | 629 µs | 624 µs | −0.8% |
| `symmetrize_density_g_n36_ops48` | 2.580 ms | 2.590 ms | +0.4% (noise) |
| `symmetrize_density_g_n72_ops48` | 18.58 ms | **17.70 ms** | **−4.7%** |
| `symmetrize_density_g_n72_ops8` | 11.24 ms | 10.63 ms | −5.4% |

Consistent with the proposal's predicted 5–15% headline; smaller-n cases are dominated by FFT wall time, not the symmetrizer inner loop.

**Memory reclaimed:** ~11.6 kB per `BasisSet` (30 kB from the deleted `index_map` + 4.4 kB from `miller` narrowing at n_pw = 725) and 27 B per `SpaceGroupOp` (~1.3 kB per `SymmetryInfo` at 48 ops).

**`index_map` decision:** option (a) — deleted entirely, replaced with an O(n_pw) linear-scan `index_of`. The 12 test-only call sites already used `Option<usize>` so no migration pattern was needed; linear-scan cost is invisible next to SCF wall time.

**KB form-factor f32 storage escalation (§2.3) still pending EM decision post-ITEV.** Not escalated as part of Phase A.

---

# TYPE: Numeric-Type Efficiency Audit (integer narrowings + targeted f32 survey)

## Why This Is A Proposal, Not An Edit

This is a **read-only audit** of every numeric-type choice in `src/` that could plausibly be narrowed (integer width reduction) or precision-reduced (f32). The DFT codebase is held to QE-equivalent precision, and QE uses `REAL(DP) = f64` everywhere in the CPU path. **Any f32 recommendation here is justified site-by-site against that baseline.** The one pattern we have already committed to is the GPU `f32` / CPU `f64` boundary for grid-point kernels in `src/gpu/*` — that is the template, not a precedent to generalize.

Three type-efficiency axes, ordered by expected ROI:
1. **Integer-width reductions** — free wins, measurable storage/cache savings.
2. **f32 at transient storage sites** — case-by-case, most rejected.
3. **Complex<f32> for reduced-to-real temporaries** — all rejected below.

No implementation here. Per-change implementation PRs follow approval.

## Baseline (from existing logbook measurements, n_pw = 725, ecut = 200 Ry)

| Quantity | Type | Elements | Bytes | Hot-path? |
|---|---|---|---|---|
| `BasisSet::miller` | `Vec<[i32; 3]>` | 725 | 8 700 | cold (built once) |
| `BasisSet::pw` | `Vec<Vector3<f64>>` | 725 | 17 400 | warm (kinetic, V_NL q-vec) |
| `BasisSet::index_map` | `HashMap<(i32,i32,i32), usize>` | 725 | ~30 000 | very cold (construction only) |
| `SpaceGroupOp::rotation` | `[[i32; 3]; 3]` | per op (≤ 48) | 36 | warm in `symmetrize_density_g` |
| Wavefunction coeffs `Mat<Complex64>` | `f64+f64i` | 725 × n_bands | 11 600/band | hot — eigensolve sink |
| `fft::FFT3D` buf_a/buf_b | `Array3<Complex64>` | dnx·dny·dnz | ~100k @ 24³ | hot |

Ranking below is by `(memory or cache win) / (correctness risk)`.

---

## 1. Integer narrowings — RECOMMEND (low risk, measurable)

### 1.1 `SpaceGroupOp::rotation: [[i32; 3]; 3]` → `[[i8; 3]; 3]`

**Site:** `src/symmetry/operations.rs:18` (`SymmOp`) and `:27` (`SpaceGroupOp`).

**Current cost:** 36 bytes per op, 1 cache line can hold ~1.7 ops. For Fd-3m Si with 48 ops that is 1 728 B (27 cache lines).

**Narrowed cost:** 9 B per op → under 2 cache lines for all 48 ops.

**Justification:** Rotation entries in the fractional basis are strictly in `{-1, 0, 1}` for cubic groups and in `{-2, -1, 0, 1, 2}` for rare hexagonal settings (those arise from c/a ≠ √3 axis conventions; empirically pwdft-rs' `detect::find_symmetry_operations` never produces |R_ij| > 3). `i8` range `[-128, 127]` is > 40× the empirical bound.

**Where the savings show up:** `symmetry::density::g_space::symmetrize_density_g` reads `op.rotation` in the inner loop (SYMP already parallelized this path). For Si 48³ with 48 ops, the op table fits comfortably in L1 with i8; at i32 it still fits in L1 but burns 4× more line traffic than needed. **Expected win at 48³·48ops: 5–15% in `symmetrize_density_g` inner loop.**

**Gotchas:**
- `SymmOp::det`, `SymmOp::compose` do 3×3 integer arithmetic; with i8, `r[1][1] * r[2][2]` can overflow (max 2·2 = 4, fine) but the triple product `r[0][0] * (r[1][1]*r[2][2] − ...)` can reach 27 for a worst-case compose. Use `i16` for temporaries or leave these arithmetic paths promoting to `i32`. Trivial — just `let a = i32::from(r[i][j])` at method entry.
- `spglib`-style external code would expect `i32`; we don't have that dependency. `detect.rs` produces these rotations internally, so the type is fully owned.
- `apply()` already casts to `f64` at the call site — no change there.

**What would break:** Nothing at the test level. The only test touching `rotation` directly is `test_90_degree_rotation_order_4` and friends, all of which use `from_flat(&[i32; 9])` — that constructor becomes `from_flat(&[i8; 9])` or stays `i32` and narrows internally.

**Validation:** run `cargo test --release` (symmetry + PCFX regression tests pin G-space invariants); Si SCF energy must match bit-for-bit because the computation path doesn't change numerically — same integer values, narrower storage.

**ROI:** high. Single struct field change, ~6 method-internal widening casts, and a search-and-replace in tests.

### 1.2 Miller indices `BasisSet::miller: Vec<[i32; 3]>` → `Vec<[i16; 3]>`

**Site:** `src/basis.rs:11` and the `HashMap<(i32,i32,i32), usize>` on line 13.

**Current cost at n_pw = 725:** 725 × 12 B = 8 700 B.

**Narrowed cost:** 725 × 6 B = 4 350 B.

**Justification:** `n_i_max ≤ ceil(g_max / |b_i|)`. For Si ecut = 100 Ry, a = 5.43 Å, |b_i| ≈ 2π/2.72 ≈ 2.31 Å⁻¹; g_max = sqrt(ecut/HBAR2_OVER_2M) ≈ sqrt(100 / 3.81) ≈ 5.12 Å⁻¹; n_max ≈ 3. At ecut = 400 Ry on a 20 Å supercell, n_max could reach ~40. **i16 range [-32768, 32767] covers every physically reasonable ecut·Ω combination by 3 orders of magnitude.** Reaching i16 overflow would need ecut > 10⁷ Ry, physically meaningless.

**Where the savings show up:**
- `symmetry::density::g_space::flat_to_miller` and `miller_to_flat` are called once per destination G per op inside `symmetrize_density_g`. Narrower Miller triples → shorter rotate_miller input, possibly shorter intermediate (if the compiler keeps it narrow).
- `scf::grid::miller_to_idx` is called for every (atom, G-vector) in V_local/V_NL setup.
- `basis.miller_indices()` is iterated in `FftGrid::basis_to_fft` — one pass per k-point setup. Cold path, but the return type propagates everywhere downstream.

**Honest caveat:** most of the math in `miller_to_idx` promotes to `i32` anyway (the `(ni % d) + d` pattern does `i32` arithmetic explicitly). Storage narrows, compute doesn't. **Expected win: 5–10% in the setup-phase cache traffic, negligible in steady-state SCF.**

**Gotchas:**
- `HashMap<(i32,i32,i32), usize>` on line 13 — if we narrow to i16, the key becomes `(i16,i16,i16)`. The map is only used in `index_of` (a lookup during debugging/tests); its internal hashing cost is unchanged. Narrowing here saves ~8 B/entry × 725 = ~5.8 kB.
- `index_of(n1: i32, n2: i32, n3: i32)` API takes i32; it'd need to either cast internally (`if let Some(&i) = map.get(&(n1.try_into().ok()?, ...))`) or change signature to i16. Keep the `i32` signature; internally `try_into()` with `ok()?` returning `None` for out-of-range queries preserves the existing semantics (out-of-range → not in basis).
- **Miller storage outside `BasisSet`:** `symmetry::density::g_space::rotate_miller` takes `[i32; 3]` and returns `[i32; 3]`. If we narrow the BasisSet field but widen back for symmetry rotation (R·n can produce values up to 3·n_max at cubic R, still well within i16 but we'd want `let n = [i16::from(m[0]), ...]`). Cleanest: keep the rotation math in i32 and widen on read.

**What would break:** Any consumer that indexes `basis.miller_indices()[ig][0]` with `i32` inference needs explicit widening. The compiler will flag these — not a silent bug.

**ROI:** medium. More surface area than #1.1, but the n_pw × 3 array appears in basis construction, FFT grid mapping, and symmetry inner loops. Savings are storage-dominated (~4.4 kB at n_pw = 725, scales linearly with n_pw).

### 1.3 `BasisSet::index_map: HashMap<(i32,i32,i32), usize>` — convert to lookup table or drop

**Site:** `src/basis.rs:13`, used only by `BasisSet::index_of`.

**Honest surprise:** The `index_map` is allocated on every `BasisSet::new` call but only used by `index_of`, which `grep` shows is called from **tests only** (none of the production code paths call it — `FftGrid::basis_to_fft` works the other direction, Miller→FFT via `miller_to_idx`). **This is dead memory on the production path.**

**Proposed change:** Either
- (a) lazily construct `index_map` on first `index_of` call (`OnceLock<HashMap<...>>`), or
- (b) replace with a sorted-array binary search over `miller_indices()` + `partial_cmp` on the triple.

Savings at n_pw = 725: ~30 kB per `BasisSet`. One `BasisSet` per run, so this is small absolutely — but it's also free: nobody needs it.

**Caveat:** flag for **Code Reviewer** — this is dead-code territory, not perf. I recommend lazy construction since `index_of` is genuinely public API. Still worth calling out in this audit because it surfaced during the type sweep.

### 1.4 `kpoints::monkhorst_pack` signature `(n1: u32, n2: u32, n3: u32)` is already tight

Not an action item — just noting that `u32` for mesh sizes is already appropriately narrow. `Vec<KPoint>` itself contains a `Vector3<f64>` and `f64 weight` which are correct (k-points enter phase factors and need full precision).

---

## 2. f32 sites — MOSTLY REJECT (correctness risk)

### 2.1 [REJECTED] Wavefunction coefficients `Mat<Complex64>` → `Mat<Complex32>`

**Why rejected:** ψ_{n,k}(G) enters ρ(r) = Σ |ψ|² which flows back into V_eff → eigensolve. At ecut = 200 Ry and conv_threshold = 1e-6, we pin 10⁻¹⁰ fractional changes in total energy (see QE cross-validation suite). f32's 7 decimal digits gives ~10⁻⁷ relative error per coefficient; squared and accumulated across n_bands · n_k · n_pw ≈ 20 · 10 · 725 = 145 000 terms, the density error floor rises to ~10⁻⁶ — right at the convergence threshold. Tests that pin |E_HF - E_KS| < 10⁻⁴ eV (`test_fe_bcc_fm_vs_qe` et al.) would likely fail.

**Storage savings that were on offer:** 11.6 kB/band → 5.8 kB/band at n_pw = 725. For n_bands = 8 and 10 k-points this is ~464 kB peak. Eigensolver dense matrix at n_pw² complex = 4.2 MB — wavefunction storage isn't the bottleneck.

**Verdict:** do not touch. The memory win is dwarfed by the O(n_pw²) Hamiltonian.

### 2.2 [REJECTED] V_local(r), V_xc(r), V_H(r) on the FFT grid → `f32`

**Why rejected:** These ARE transient-after-construction on the CPU path, and the GPU path already does this (`gpu/shaders/*.wgsl`). But on the CPU path, the final consumer is `build_hamiltonian` via `v_eff_fft_to_g`, which needs `Complex64` because H is `Mat<Complex64>` for the eigensolver. Any f32 intermediate requires f32→f64 conversion at the handoff — not free. **On Apple M3 Max the XC path is bound by `cbrt`/`ln` scalar latency (see logbook 2026-04-17 XCPR entry), not SIMD width; M1/M2 NEON is 128-bit so 2 × f64 or 4 × f32 — no vector-width win from f32 on transcendentals.** The cbrt/ln latency is identical f32 vs f64 on Apple Silicon.

**Savings that were on offer:** halves one or two Vec<f64> of size n_grid. At 32³ grid = 32 768 elements that's 262 kB → 131 kB, one-time peak.

**Verdict:** do not touch on CPU. The existing GPU f32 path already captures this win where it matters (memory bandwidth to VRAM).

### 2.3 [REJECTED, escalate to EM] KB projector form factors `form_factors: Vec<Vec<f64>>`

**Site:** `src/potential/nonlocal.rs:86`.

**Analysis:** F_i(q) is a radial Bessel transform — smooth function of |k+G|. At ecut = 200 Ry, F_i values span ~5 orders of magnitude (projector values are sharply peaked at small q). f32 would give ~10⁻⁷ relative precision, which multiplied into D_ij · F^H leg of V_NL could perturb eigenvalues by ~10⁻⁵ eV.

**Savings:** at n_pw = 725 with ~4 projectors per Si atom and 2 atoms = 8 projectors × 725 × 8 B = 46 kB. Narrowing halves this.

**Why this is an "EM decision":** F_i is neither a pure accumulator nor a pure storage value — it's a Bessel integral result that then feeds a BLAS-3 multiplication. The logbook flags `vnl_new = 47 ms` (n_pw = 725); much of that is Bessel transform cost and cache traffic. **Potentially** f32 storage with f64 accumulation in the GEMM could win 10–20% on `vnl_new`, but the risk of a few-meV total-energy drift is real and would need QE cross-validation at 1e-8 convergence thresholds to catch.

**Recommendation:** escalate as a separate proposal contingent on ITEV landing (which removes V_NL build from the dominant cost), NOT as part of this audit's implementable set.

**Validation protocol if someone does this:** full test suite + tighten QE cross-check thresholds (`test_si_scf_vs_qe`) from 1e-4 eV to 1e-6 eV, run 10 SCF cases, check residuals.

### 2.4 [REJECTED] FFT scratch buffers `ndrustfft<f64>` → `ndrustfft<f32>`

**Why rejected:** Per-iteration SCF does ~6 FFTs on the full density grid (forward, V_H, inverse, inverse-normalized). At 32³ = 32 768 complex × 16 B = 524 kB, the buffer is too big for L2 on Apple M3 Max (4 MB shared). Narrowing to f32 halves to 262 kB — still L2-resident. Bandwidth-bound transforms might see 15–25% gain, but ρ(r) needs to return to f64 for XC, so every FFT incurs an f64→f32→f32→f64 conversion chain. **Existing FFTB (completed) already captured the buffer-reuse win; further narrowing is a functional change requiring a second ndrustfft handler path.**

**Verdict:** defer indefinitely. Revisit only if `criterion --bench fft` shows FFT dominant post-ITEV (currently 2.1% of user CPU per the logbook).

### 2.5 [REJECT — with note] `basis.kinetic_energy(): Vec<f64>`

Already `f64`. No action. Noted because it could naively look like a "transient" target, but it's the diagonal of H — must stay f64.

---

## 3. Complex<f32> for reduced-to-real temporaries — REJECT

The candidate pattern is `Σ_G ρ(G) · exp(iG·R)` (the structure factor paths). These are computed with `Complex64::cis` and summed into ≤ n_atoms accumulators. Savings are negligible (n_pw × 16 B = 11.6 kB at n_pw = 725, one-shot per call). The accumulator must be f64 regardless, so narrowing the per-G multiply doesn't actually save memory (they're loop-local values in SIMD registers) and introduces a conversion dance. Do not touch.

---

## Adjacent findings (mention, don't prescribe)

### A. `PseudopotentialData::r_grid, rab, v_local, rho_atom, core_charge` are `Vec<f64>`

Radial grids are typically 500–2 000 points. `Vec<f64>` is the correct choice (radial integrals need full precision); `Box<[f64]>` would save `Vec`'s 8-byte capacity field but not change compile semantics. **Not worth touching.**

### B. `SpaceGroupOp::translation: [f64; 3]` is already right

Fractional translations need full precision for symmetry tolerance checks (`tolerance = 1e-5` at `symmetry::detect`). Do not narrow.

### C. `SymmetryInfo::{has_inversion, has_time_reversal}: bool`

Already one byte each — no `bitflags` packing motivated at n_groups = 1.

### D. SoA vs AoS for G-vector data

`BasisSet` keeps `pw: Vec<Vector3<f64>>` (AoS) and `miller: Vec<[i32; 3]>` (AoS). The two arrays are iterated together in `NonlocalPotential::new` (line 120 creates `q_vecs` by summing k + g). Converting to SoA (three `Vec<f64>` for Gx, Gy, Gz) would enable auto-vectorized q-norm computation — but we don't have that hot loop today; Bessel transforms dominate. **Not worth proposing** under the "no baseline, no opt" rule.

### E. Surprise: `BasisSet::index_map` is only used by tests

See §1.3. Flag for **Code Reviewer** as dead-memory dead-code.

### F. `miller_to_idx` `as i32` / `as usize` dance

`src/scf/grid.rs:108` does `(n1 % dims[0] as i32) + dims[0] as i32) as usize % dims[0]`. This is correct and clippy-allowed. **If Miller narrows to i16 per §1.2**, that cast becomes `n1 as i32` which the compiler will handle for free. No action independently.

---

## Ranked implementable items

| Rank | Change | Site | Savings | Risk | Verification |
|---|---|---|---|---|---|
| 1 | `SpaceGroupOp::rotation` i32 → i8 | `symmetry/operations.rs` | ~30 B/op × ≤48 ops + cache-line reduction in `symmetrize_density_g` inner loop | very low | `cargo test` (symmetry suite + PCFX regression); Si SCF bit-identical total energy |
| 2 | `BasisSet::miller` `[i32;3]` → `[i16;3]` + map key `(i32,i32,i32)` → `(i16,i16,i16)` | `basis.rs`, downstream `miller_to_idx` call sites | ~4.4 kB at n_pw = 725 + cache-line reduction in G-loop setup | low | `cargo test`; no numeric impact — integer values unchanged |
| 3 | `BasisSet::index_map` → lazy / sorted-array (**flag to Code Reviewer**) | `basis.rs` | ~30 kB per `BasisSet`, semi-dead on production | low | `cargo test` |

## Validation protocol (for any change above)

1. `cargo test --release` — full suite including QE cross-validation.
2. Run `examples/si_scf.yaml` and `examples/si_scf_converged.yaml`, compare `total_energy`, `fermi_energy`, first 8 eigenvalues at Γ. Bit-identical expected for §1.1 and §1.2 because no arithmetic changes occur — just narrower storage for values that already fit.
3. `cargo clippy -q --all-targets` AND `cargo clippy -q --all-targets --features gpu`.
4. For §1.1 specifically: verify `SymmOp::{compose, inverse, det}` correctness on the hexagonal P6/mmm op set (|R_ij| = 2 entries) if we keep any `i8` arithmetic internal. A targeted unit test with R = 2·identity (synthetic) asserts no overflow.

## Out of scope for this audit

- GPU f32 path itself: already correctly scoped; leave `gpu/*` alone.
- `Complex<f64>` in wavefunctions: see §2.1.
- Anything touching accumulators, `eigenvalues: Vec<f64>`, `total_energy: f64`, or `rho_g: Vec<Complex64>` in `ScfResult`. These are the published quantities — immutable by precision policy.

## Future-bench items (record, do NOT implement in this proposal)

- After §1.1+§1.2: `criterion --bench scf_benchmarks -- symm` for `symmetrize_density_g` at 48³·48ops; expect 5–15% speedup.
- After §1.3 if Code Reviewer adopts: measure `BasisSet::new` wall time at n_pw = 725 — expect ~10% reduction (HashMap construction is O(n_pw)).
- Post-ITEV: re-profile to see if `vnl_new` (currently ~47 ms at n_pw = 725) still warrants §2.3 reconsideration.

## References

- Logbook: `.claude/logbooks/performance-engineer.md` — 2026-04-16 baseline, 2026-04-17 XCPR/FMAD/ITEV entries for bottleneck ranking.
- `src/gpu/mod.rs:10-20` — established GPU f32 / CPU f64 precision boundary and documented error budget (reference for "where f32 is justified").
- QE 7.5: `PW/src/symme.f90` uses `INTEGER` (4-byte) for `s(:,:,ns)` rotation tables — we can do better with i8 because we don't interop.
- Proposal `CFGN` — orthogonal; this audit does not overlap (CFGN exposes hardcoded *values*; this exposes inefficient *types*).
- Proposal `CLSS` — similarly orthogonal; CLSS is about `cast_lossless` lint hygiene, not about changing storage types.
