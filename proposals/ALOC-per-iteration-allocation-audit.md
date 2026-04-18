---
id: ALOC
status: active
priority: medium
complexity: medium
risk: low
depends_on: []
blocks: []
---

# ALOC: Per-iteration allocation audit for the SCF hot loop

## §1 Scope

**In scope.** This is a static read-only audit of heap allocations that
occur **inside the SCF iteration body** — the body of the
`for iter in 0..ctx.params.max_iter` loops in
[`src/scf/driver.rs:161-393`](../src/scf/driver.rs) (non-spin) and
[`src/scf/driver_spin.rs:135-507`](../src/scf/driver_spin.rs) (spin).
Callees of that body (`lda_xc_grid`, `assemble_v_eff`, `hartree_on_fft_grid`,
`compute_density`, `symmetrize_density_g`, `hartree_energy`, `xc_energy_*`,
`add_core_density`, `density_r_to_g`, `real_to_g_space`, `build_hamiltonian_with_v_eff`,
`NonlocalPotential::add_to_hamiltonian`, `mixer.mix`) are audited as well
whenever they allocate on behalf of the loop.

**Out of scope.**

- **One-time-setup allocations** in `ScfContext::new`
  ([`src/scf/context.rs`](../src/scf/context.rs)), `NonlocalPotential::new`
  ([`src/potential/nonlocal.rs:106`](../src/potential/nonlocal.rs)),
  initial density (SAD), `compute_v_local`, `compute_core_density`. These
  run once before the SCF loop, so even a large `Vec` there is irrelevant
  to per-iteration overhead.
- **Mixer history buffers** (`history_in`, `history_res` in Anderson /
  Broyden / PRPL). They grow monotonically to `max_history` (≤ 8 by
  default), then recycle. Not a per-iteration allocation; audited below
  only where mixers allocate *in addition to* history.
- **The FFT scratch path.** Already addressed by FFTB (completed).
  `FFT3D::forward` / `FFT3D::inverse` reuse `buf_a`/`buf_b` and do not
  allocate per call.
- **GPU transfer buffers** (behind `feature = "gpu"`). `GpuAccelerator`
  already uses a `BufferPool` of pre-allocated f32 staging buffers; the
  f32→f64 conversions on the CPU side happen inside the GPU helper
  methods, not in the driver body.
- **The convergence branch** (the `if rho_converged && energy_converged`
  block at [driver.rs:314-388](../src/scf/driver.rs)). Runs once at
  convergence, not per iteration.
- **Eigensolve internals** (faer dense Hermitian). Proposals WFRX + ITEV
  address eigensolver cost directly. The driver-side handling of
  `kpoint_results` is audited.
- **The VGC5 per-component decomposition** (only computed at convergence).

**Convention.** Throughout the tables below:

| Symbol     | Meaning                                                                 | Typical range                          |
|------------|-------------------------------------------------------------------------|----------------------------------------|
| `n_pw`     | Number of plane waves at a k-point                                      | 89 (ecut=100) – 725 (ecut=400)         |
| `n_bands`  | Bands computed per k-point                                              | 4–16                                   |
| `n_grid`   | FFT grid points `nx·ny·nz`                                              | 5 832 (18³) – 373 248 (72³)            |
| `n_k`      | k-points in the IBZ                                                     | 1 (Γ) – ~64 (4×4×4 with symmetry)      |
| `n_atoms`  | Atoms in the unit cell                                                  | 2 (Si) – 8+ (Fe, supercells)           |
| `n_sym`    | Space-group operations                                                  | 1 (identity) – 48 (Fd-3m)              |

**Hardware reference.** All savings estimates target Apple M2 (8-core, 8
GB unified LPDDR5), release profile, where `malloc`/`free` for a
medium-sized `Vec` typically costs 300 ns–1 µs and a 64³ `Vec<f64>`
allocation+touch is ~30 µs. Grid sizes that clear the 2 MB L2
per-core boundary (64³ `Vec<Complex64>` = 4 MB) force a fresh OS page
walk on first write, roughly doubling the cost.

## §2 Findings

The hot loop body is analyzed in execution order. For each finding:

- `file:line` — primary source location.
- `size` — per-allocation byte count as a function of the size symbols above.
- `per iter` — number of allocations of that kind per SCF iteration.
- `fix` — suggested remediation.
- `µs/iter (M2 estimate)` — rough static estimate: product of
  (alloc count × touch/memset cost at that size). "touch cost" uses ~1 ns/byte
  at ≤ L2, ~2 ns/byte above. A ×1 multiplier is assumed unless stated.

### 2.1 Non-spin driver (`src/scf/driver.rs`)

#### F-1. Hartree potential — `hartree_on_fft_grid` allocates full Complex64 grid

- **Site:** [`energy.rs:232-247`](../src/scf/energy.rs) (called from `driver.rs:175, 199`).
- **What:** `par_iter().zip(...).collect()` → `Vec<Complex64>` of length `n_grid`.
- **Size:** `16 · n_grid` bytes (16 bytes per Complex64). 64³ = 4 MB.
- **Per iter:** 1 (non-spin) / 1 (spin — shared Hartree).
- **Fix:** Pass in a preallocated `&mut [Complex64]` of length `n_grid`
  from a new `ScfWorkspace` struct stored alongside `rho_r`/`rho_g`.
  Convert from returning `Vec` to filling a buffer: pattern mirrors
  `density_r_to_g(&mut fft, &rho_r, &mut rho_g)`.
- **Est. saving:** **~8–40 µs/iter** at 24³–48³ grids.

#### F-2. XC potential — two `Vec<f64>` from `lda_xc_grid` + G-space FFT of V_xc

- **Sites:**
  - [`xc.rs:75-94`](../src/potential/xc.rs): `lda_xc_grid` returns
    `(Vec<f64>, Vec<f64>)` of length `n_grid` each.
  - [`energy.rs:186-190`](../src/scf/energy.rs): `real_to_g_space` wraps
    `density_r_to_g` and allocates a fresh `Vec<Complex64>` of length
    `n_grid`.
  - [`driver.rs:187,189`](../src/scf/driver.rs): both are consumed
    immediately.
  - [`driver.rs:266`](../src/scf/driver.rs): a **second** pair of
    `lda_xc_grid` calls on the OUTPUT density for the Harris-Foulkes
    / E_KS split.
- **Size:** `8 · n_grid` bytes each for exc / vxc; `16 · n_grid` for vxc_g.
- **Per iter:** 4 `Vec<f64>` (2 calls × 2 vecs) + 1 `Vec<Complex64>` (vxc_g).
- **Fix:** Extend `lda_xc_grid` to an in-place variant
  `lda_xc_grid_into(rho_r, exc_r_out, vxc_r_out)`. Store `exc_r`,
  `vxc_r`, `vxc_g` in the workspace. Reuse across the two calls (input
  density and output density) since the first pair is no longer live
  when the second fires.
- **Est. saving:** **~50–200 µs/iter** at 24³–48³. At 64³ non-NLCC
  materials this is one of the top three savings.

#### F-3. `add_core_density` — `Vec<f64>` even when there is no NLCC

- **Site:** [`energy.rs:219-229`](../src/scf/energy.rs), called at
  `driver.rs:179, 265` and twice per iter in the spin driver.
- **What:** Branch A (NLCC off) does `rho_val.to_vec()` → a fresh
  `Vec<f64>` of length `n_grid`, wasted — the caller could use
  `&rho_val` directly. Branch B (NLCC on) allocates via `.collect()`.
- **Size:** `8 · n_grid` bytes.
- **Per iter:** 2 non-spin / 4 spin.
- **Fix:** Two remedies, landing in order:
  1. Change the signature to return `Cow<'_, [f64]>` so the non-NLCC
     path borrows.
  2. When NLCC *is* active, thread a workspace buffer through.
- **Est. saving:** **~5–20 µs/iter** (NLCC off, Si LDA) /
  **~10–40 µs/iter** (NLCC on, Fe). The NLCC-off path is a pure win
  because the allocation is currently completely wasted.

#### F-4. `assemble_v_eff` — per-iteration `Vec<Complex64>` of length `n_grid`

- **Site:** [`energy.rs:193-205`](../src/scf/energy.rs), called at
  `driver.rs:199, 196`.
- **What:** Parallel `map + collect` over three zipped `[Complex64]`
  inputs → new `Vec<Complex64>`.
- **Size:** `16 · n_grid` bytes. 64³ = 4 MB.
- **Per iter:** 1 non-spin / 2 spin (V_eff for up and down channels).
- **Fix:** `assemble_v_eff_into(v_local, v_h, vxc_g, out: &mut [Complex64])`
  stored in workspace. Safe to do with rayon
  (`par_iter_mut().zip(par_iter).zip(par_iter).zip(par_iter)`).
- **Est. saving:** **~10–45 µs/iter** non-spin; ~2× for spin.

#### F-5. Per-k-point Hamiltonian matrix — `faer::Mat::<Complex64>::zeros(n_pw, n_pw)`

- **Site:** [`scf/potentials.rs:147`](../src/scf/potentials.rs).
- **What:** Fresh zeroed n_pw × n_pw complex matrix every iteration, per k-point.
- **Size:** `16 · n_pw²` bytes. **725² × 16 = 8.4 MB per k-point.**
- **Per iter:** `n_k` non-spin / `2 · n_k` spin. For Si 4×4×4 with 10 IBZ
  k-points and n_pw=725: **84 MB/iter** of short-lived heap. This is the
  **single largest per-iteration allocation in the loop**.
- **Fix:** Two tiers:
  1. Cache one `faer::Mat<Complex64>` per k-point in `ScfContext`
     (total `n_k · 16 · n_pw²` bytes of resident state — same as already
     done for `vnl_cache` entries, which each hold a `B` of size
     `n_pw × n_channels`). Mutate in place each iteration via
     `h.fill(Complex64::zero())` followed by the existing kinetic-fill
     loop, then add V_eff, then VNLM GEMM.
  2. Alternatively, a `ThreadLocal<faer::Mat<Complex64>>` pool keyed
     by (n_pw, n_pw), allocated lazily, reused across iterations.
     Avoids sizing the pool to `n_k`; good for systems where
     `n_k > n_threads`.
- **Est. saving:** **~200–2 000 µs/iter** for Si-class
  (n_pw=89 → ~0.13 MB, cheap; n_pw=725 → 8.4 MB, above L2 → page-walk
  on first touch). Above ~n_pw=500 this becomes the headline allocation
  saving in the driver body.
- **Correctness watch:** Cached matrices must be zeroed before reuse.
  `build_hamiltonian_with_v_eff` currently relies on `Mat::zeros` semantics,
  then writes every diagonal and every off-diagonal. A one-liner
  `h.fill(Complex64::zero())` before the kinetic loop preserves that
  invariant.

#### F-6. Eigenvalue extraction — `Vec<Vec<f64>>` clone per k-point

- **Site:** [`driver.rs:214`](../src/scf/driver.rs) and
  [`driver_spin.rs:184-185`](../src/scf/driver_spin.rs).
- **What:** `kpoint_results.iter().map(|r| r.eigenvalues.clone()).collect()` —
  clones the eigenvalue vector out of each `EigenResult` so the
  subsequent `all_kpoint_wavefns` move doesn't borrow-check-fail.
- **Size:** `8 · n_bands · n_k` bytes total (tiny, typically < 1 KB). Not
  a size concern; the concern is the double `Vec<Vec<f64>>` (outer +
  inner per-k allocations).
- **Per iter:** 1 outer + `n_k` inner non-spin; 2× spin.
- **Fix:** Refactor `EigenResult` to `(eigenvalues: Vec<f64>, eigenvectors: Mat)`
  and use `kpoint_results.into_iter().map(|r| (r.eigenvalues,
  r.eigenvectors)).unzip()`. This consumes `kpoint_results` once and
  extracts both fields without cloning. Micro-fix; listed for
  completeness.
- **Est. saving:** **< 5 µs/iter.** Mostly a code-clarity win, not a
  performance win.

#### F-7. `compute_density` — per-worker per-band `psi_g` allocation

- **Site:** [`scf/density.rs:65`](../src/scf/density.rs) inside the rayon
  `fold` closure:
  `let mut psi_g = vec![Complex64::new(0.0, 0.0); n_grid];`
- **What:** Fresh `Vec<Complex64>` of length `n_grid` **for every band of
  every k-point**. Inner loop `for ib in 0..n_bands`.
- **Size:** `16 · n_grid` bytes. Over one iter: `n_k · n_bands` allocs,
  each 16·n_grid bytes.
- **Per iter:** `n_k · n_bands`. Si 4×4×4 / 4 bands / 32³ grid = 40 allocs
  × 524 KB = **21 MB/iter of Complex64 grid allocations inside this
  function alone**.
- **Fix:** Move `psi_g` outside the `for ib in 0..n_bands` loop —
  allocate once per rayon worker (inside the `fold`'s init closure),
  zero it before each band, and reuse. `psi_g` is written into
  completely during the `for ig in 0..n_pw` sparse-scatter + the
  inverse FFT, so a `.iter_mut().for_each(|v| *v = Complex64::zero())`
  memset before the scatter is adequate. The per-worker allocation
  pattern is already in use for `rho_acc` via `fold`'s init closure —
  extend it to `psi_g` alongside `fft_local`.
- **Est. saving:** **~200–1 000 µs/iter** at 32³–48³. Likely the single
  largest fix outside F-5.

  (There is a secondary related issue in the same function:
  `let mut fft_local = FFT3D::new(dnx, dny, dnz);` is created per
  *k-point*. Each `FFT3D::new` itself allocates two `Array3<Complex64>`
  scratch buffers and builds per-axis FFT handlers — it's a
  non-trivial allocation. The rayon `fold` init closure handles this:
  it's built once per worker, reused per k-point. But ndrustfft's
  `FftHandler::new` also allocates twiddle factors per call; see
  **F-12** below.)

#### F-8. `density_r_to_g` caller-side staging — `rho_g`, `rho_g_new`, `rho_total_new_g`

- **Sites:**
  - [`driver.rs:142`](../src/scf/driver.rs): `rho_g` outside the loop
    (one-time, fine).
  - [`driver.rs:262`](../src/scf/driver.rs): `rho_g_new` allocated inside
    the loop.
  - [`driver_spin.rs:139, 294`](../src/scf/driver_spin.rs): two more
    per-iter allocations.
- **Size:** `16 · n_grid` bytes per Vec.
- **Per iter:** 1 non-spin / 2 spin. Plus the per-iter
  `real_to_g_space(vxc_r, ...)` at `driver.rs:189` = +1 allocation
  (already covered under F-2 as the `vxc_g` buffer).
- **Fix:** Hoist `rho_g_new` into `ScfWorkspace` alongside `rho_g`. The
  function already takes `&mut [Complex64]`; the caller just needs to
  stop allocating the buffer fresh each iter.
- **Est. saving:** **~8–40 µs/iter** at 24³–48³.

#### F-9. `density_diff` — no allocation, fine

- **Site:** [`energy.rs:163-171`](../src/scf/energy.rs).
- **What:** Pure `iter().zip().map().sum()`. Zero allocations.
- **Note:** Called 1× non-spin / 2× spin. No action.

#### F-10. `symmetrize_density_g` — several per-call allocations

- **Site:** [`symmetry/density/g_space.rs:156-271`](../src/symmetry/density/g_space.rs).
- **What:**
  - Line 173: `let mut rho_g: Vec<Complex64> = rho_r.iter().map(...).collect();`
    — `16 · n_grid` bytes.
  - Line 211: `let mut rho_g_sym: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n_grid];`
    — another `16 · n_grid` bytes.
  - Line 213-217: `let r_ts: Vec<[[i32; 3]; 3]>` — `36 · n_sym` bytes
    (tiny; ≤ 1.7 KB at n_sym=48).
- **Size:** 2 × 16 · n_grid bytes (~8 MB per call at 64³).
- **Per iter:** 1 non-spin / 2 spin.
- **Fix:** Accept two workspace buffers `(rho_g_tmp, rho_g_sym_tmp)`
  from the caller (workspace-owned). `r_ts` can be precomputed once at
  SCF entry and stored in `ScfContext` (rotations don't change during
  the SCF loop — they're fixed by crystal symmetry). That cache avoids
  the per-iter `.iter().map(transpose_rotation).collect()`.
- **Est. saving:** **~20–100 µs/iter** at 24³–64³. Possibly larger
  because both allocations are written to via long outer loops (page-walk
  cost is paid in full).

#### F-11. Mixer — per-iter allocations in `Mixer::mix`

- **Sites (Anderson, [anderson.rs:117-121](../src/scf/mixing/anderson.rs)):**
  The `mix` path calls `push_history` → `diis_step`.
  - `push_history` (line 130): builds `raw_residual` via
    `&rho_out_arr - &rho_in_arr` → `Array1::zeros`-equivalent. Allocates.
  - `push_history` (line 135): `precondition_residual(&raw_residual.to_vec(), ...)`
    — calls `raw_residual.to_vec()`, a fresh `Vec<f64>` clone of the
    residual Array1 (length n_grid).
  - `precondition_residual` ([kerker.rs:24-42](../src/scf/mixing/kerker.rs)):
    allocates `res_g: Vec<Complex64>` and the final `Vec<f64>` output.
    Two further `n_grid` allocations.
  - `push_history` (line 153): `rho_in_arr.to_owned()` — new `Array1<f64>`
    of length n_grid (kept in history → amortized, not a *per-iter*
    cost past the first few iters; fine).
  - `diis_step` (line 194-196): builds `dr: Vec<Array1<f64>>` of length
    `mm = history_len - 1`. Per-entry allocation: `Array1::from(&a - &b)`.
  - `diis_step` (line 198): `a_mat`, `b_vec` — `mm × mm` and `mm` float
    vectors. Tiny (mm ≤ 8).
  - `diis_step` (line 212-219): the mixed density is built via
    `alpha_last * (&hist + &(β · R))` — ndarray expression templates
    that materialize a fresh `Array1<f64>` per term in the DIIS
    combination (≈ mm intermediate arrays of length n_grid).
- **Sites (Broyden, [broyden.rs:108-213](../src/scf/mixing/broyden.rs)):**
  - Line 112: `raw_residual` `Vec<f64>` of length n_grid — allocated.
  - Line 115: `residual` = either preconditioned (3 allocations inside
    `precondition_residual`, same as Anderson) or `raw_residual` moved.
  - Line 139-140: `dv` and `df` differences — two `Vec<f64>` of length
    n_grid. Appended to `history_dv/df`; amortized after
    `max_history` warmup.
  - Line 153-154: `prev_rho_in = Some(rho_in.to_vec())` +
    `prev_residual = Some(residual.clone())` — **two full `Vec<f64>`
    clones of length n_grid every iter**. Freed at the start of the
    *next* iter (when the `Option` is overwritten). Net: 2 × n_grid
    transient allocations *held* across iterations but `mem::swap`-able.
  - Line 164, 198-199: `beta_mat`, `work`, `corr_dv`, `corr_df` — two
    `m × m` vectors (tiny) + two full-grid `Vec<f64>` (2 × 8 · n_grid).
    The latter two are per-iter n_grid allocations (not amortized).
- **Per iter (Anderson + Kerker, n_grid = 32³ ≈ 33k):** at least 5
  full-grid allocations per `mix` (raw_residual Array1 + to_vec +
  residual_g Complex64 + precond out Vec + dr-chain in diis_step).
- **Per iter (Broyden + Kerker):** 6+ full-grid allocations.
- **Fix:** The mixer is a good candidate for a cross-cutting
  `MixerWorkspace`: preallocate `residual_r`, `residual_g`
  (for Kerker), `dv_scratch`, `df_scratch`, plus a small pool for the
  DIIS `a_mat` / `b_vec`. For Broyden specifically, the
  `prev_rho_in.clone()` and `prev_residual.clone()` at line 153-154 can
  be replaced with `std::mem::swap`s into pre-owned buffers (the logic
  is "store current, reuse next iter" — perfect for a 2-slot ring).
- **Est. saving:** **~30–100 µs/iter** per mixer (grid-dependent).
  Larger for spin (mixer called twice per iter: ρ_total and m).

#### F-12. FFT creation inside `compute_density`'s per-k-point closure

- **Site:** [`scf/density.rs:57`](../src/scf/density.rs).
- **What:** `let mut fft_local = FFT3D::new(dnx, dny, dnz);` inside the
  rayon `fold` init closure. The closure's init runs **once per worker
  thread** (rayon convention), not per k-point, so in steady state
  this is amortized across iterations… **if** the rayon pool reuses
  workers. If the pool respins (e.g. shared global pool under pressure
  from other rayon work in XCPR), each respawn rebuilds the
  FftHandlers (allocates twiddle tables per axis).
- **Size:** `FFT3D::new` allocates 2 × `Array3<Complex64>` of size
  `n_grid` (2 × 16 · n_grid bytes), plus small twiddle tables per axis
  (O(N log N) complex entries — hundreds of KB at 64³).
- **Per iter:** 0 in steady state (amortized via `fold` init). But
  debug: `scf::density::compute_density` is called **every SCF
  iteration**, and each call re-runs `fold(|| { FFT3D::new(...) }, …)`.
  Rayon's init closure fires once per worker *per `fold` call*, not
  once per entire program. So this is in fact **~n_threads FFT3D
  allocations per SCF iter**, each building twiddle tables fresh.
- **Fix:** Thread a `&[ThreadLocal<FFT3D>]` (or a `&[FFT3D]` sized to
  `rayon::current_num_threads()`) through `DensityGrid` from
  `ScfContext`, preallocated once. Reset with a cheap `.iter_mut()`
  loop if FFT3D state ever mutates (our implementation is stateless
  between calls — buffers are overwritten at each call start).
- **Est. saving:** **~100–500 µs/iter** at 32³–48³ (twiddle-table
  reallocation dominates at larger grids). This could be the dark-horse
  fix — the cost scales with `n_threads · n_grid · log n_grid`.

### 2.2 Spin driver (`src/scf/driver_spin.rs`) — additional findings

Most of the non-spin findings apply at 2× rate (two Hamiltonians, two
V_xc channels, two `compute_density` calls, two mixer calls). New
spin-specific findings:

#### F-13. CCMX basis change — six `Vec<f64>` allocations per iter

- **Site:** [`driver_spin.rs:473-506`](../src/scf/driver_spin.rs).
- **What:** At the end of each iter, the CCMX mixer requires
  `(ρ↑, ρ↓) → (ρ_total, m) → mix → (ρ_total_new, m_new) → (ρ↑_new, ρ↓_new)`.
  This is expressed as six `.iter().zip().map().collect()` calls:
  `rho_total_in`, `m_in`, `rho_total_out`, `m_out`,
  `rho_up_r` (new), `rho_down_r` (new). Each `Vec<f64>` is `8 · n_grid`
  bytes.
- **Per iter:** 6 `Vec<f64>` of length n_grid = 48 · n_grid bytes. At
  32³: ~1.6 MB. At 48³: ~5.3 MB.
- **Fix:** Add four workspace buffers (`basis_tmp_a`, `basis_tmp_b`,
  `rho_up_new`, `rho_down_new`) to the spin workspace. The mixer's
  return type remains `Vec<f64>` but can become `&mut [f64]`-in-place
  at the same time as F-11.
- **Est. saving:** **~15–50 µs/iter** at 32³–48³.

#### F-14. Spin driver: `rho_total_r`, `occ_all`, `rho_xc_total`, `rho_xc_total_in`

- **Sites:** [`driver_spin.rs:137-138, 280-281, 297-298, 336`](../src/scf/driver_spin.rs).
- **What:**
  - `rho_total_r` (line 137): `rho_up_r + rho_down_r` — `Vec<f64>` of
    length n_grid.
  - `rho_total_new` (line 280): same shape, built from the symmetrized
    up+down.
  - `occ_all` (line 298): `Vec<Vec<f64>>` concat of `occ_up` and
    `occ_down`. Small (`n_k · n_bands` f64), but allocates
    `2 · n_k` inner Vecs.
  - `rho_xc_total` (line 297): output of `add_core_density` — covered
    by **F-3**.
  - `rho_xc_total_in` (line 336): second `add_core_density`, INPUT
    density — covered by **F-3**.
  - `rho_up_xc`, `rho_down_xc`, `rho_up_xc_out`, `rho_down_xc_out` (lines
    146-147, 307-308): four more `add_core_density` calls per iter.
- **Per iter:** 6 full-grid `Vec<f64>` in the spin driver, in addition
  to the 4 under F-3. Total `add_core_density`-style allocations in
  spin = 4 per iter; other full-grid sums = 2 per iter.
- **Fix:** Fold into the same workspace buffers as F-3 (one shared
  `Cow<'_, [f64]>` or write-through buffer) and F-13 (the total/m
  basis buffers already cover `rho_total_r` and `rho_total_new`).
- **Est. saving:** **~30–80 µs/iter** spin.

#### F-15. Spin driver: per-spin V_eff_up, V_eff_down

- **Site:** [`driver_spin.rs:154-155`](../src/scf/driver_spin.rs). Two
  `assemble_v_eff` calls → **2×** the non-spin F-4.
- **Fix:** Two workspace buffers.
- **Est. saving:** already counted in F-4; spin factor 2×.

### 2.3 Nonlocal (`src/potential/nonlocal.rs`)

#### F-16. `NonlocalPotential::add_to_hamiltonian` is allocation-free

- Already a single `matmul` into the caller's `h` matrix. Good.
- **No action.**

The allocating code (`NonlocalPotential::new`) runs once per k-point at
SCF entry (stored in `vnl_cache`) — out of scope.

### 2.4 Energy-convergence dual evaluation

#### F-17. Double Hartree+XC evaluation for Harris-Foulkes

- **Site:** [`driver.rs:266-285`](../src/scf/driver.rs).
- **What:** Every iter recomputes XC on the OUTPUT density
  (`xc::lda_xc_grid(&rho_new_for_xc)` → 2 fresh Vec<f64>) **and** calls
  `hartree_energy` a second time on `rho_g` (no allocations, just a
  scalar reduction). The spin driver has a similar structure but also
  re-computes `lda_xc_spin_grid` a second time
  ([driver_spin.rs:309-310](../src/scf/driver_spin.rs)).
- **Size:** 2 × `8 · n_grid` bytes (non-spin) / 3 × `8 · n_grid` bytes
  (spin) per iter in addition to the first XC call.
- **Fix:** Already partially covered by F-2. The second XC call would
  use the same `exc_r_out`/`vxc_r_out` workspace slots once F-2 lands.
- **Est. saving:** bundled into F-2.

### 2.5 Summary table

| # | Site (file:line) | Size/alloc | Allocs/iter | Fix sketch | Est. µs/iter (M2) |
|---|---|---|---|---|---|
| F-1  | energy.rs:232 (Hartree) | 16·n_grid | 1 | workspace `v_h` buffer | 8–40 |
| F-2  | xc.rs:75 + driver.rs:189, 266 (XC + vxc_g) | 8·n_grid × 2, 16·n_grid | 4 vec<f64> + 1 complex | `lda_xc_grid_into` + workspace | 50–200 |
| F-3  | energy.rs:219 (add_core_density) | 8·n_grid | 2 non-spin / 4 spin | `Cow<'_, [f64]>` + workspace | 5–20 / 10–40 |
| F-4  | energy.rs:193 (assemble_v_eff) | 16·n_grid | 1 / 2 spin | `_into` variant | 10–45 |
| F-5  | scf/potentials.rs:147 (build_hamiltonian Mat) | 16·n_pw² | n_k / 2·n_k | cache per-k Mat or thread-local pool | **200–2 000** |
| F-6  | driver.rs:214 (eigenvalue clone) | 8·n_bands·n_k | 1 outer + n_k inner | `.into_iter().unzip()` | <5 |
| F-7  | density.rs:65 (psi_g per band) | 16·n_grid | n_k·n_bands | move outside band loop, reuse per worker | **200–1 000** |
| F-8  | driver.rs:262 (rho_g_new) | 16·n_grid | 1 / 2 spin | workspace | 8–40 |
| F-10 | g_space.rs:173, 211 (symmetrize_density_g) | 16·n_grid × 2 | 2 per call × 1/2 | workspace + precomputed R^T | 20–100 |
| F-11 | anderson/broyden mixer | 8/16·n_grid × 5-6 | 5–6 | MixerWorkspace | 30–100 |
| F-12 | density.rs:57 (FFT3D::new in fold init) | ~32·n_grid + twiddles | n_threads | thread-local FFT3D pool | **100–500** |
| F-13 | driver_spin.rs:473 (CCMX basis change) | 8·n_grid × 6 | 6 | 4 workspace buffers | 15–50 |
| F-14 | driver_spin.rs:137,280 (rho_total_r, _new) | 8·n_grid × 2 | 2 | workspace | 10–30 |
| F-17 | driver.rs:266 (second XC call) | absorbed into F-2 | — | — | — |

**Aggregate estimate (non-spin, Si 4×4×4 @ n_pw=725, 32³ grid, n_k=10):**

| Bucket                   | Low (µs/iter) | High (µs/iter) |
|--------------------------|---------------|-----------------|
| Grid-level Vecs (F-1/2/3/4/8/10) | 100 | 450 |
| Per-k Mat allocation (F-5)       | 200 | 2 000 |
| psi_g in density (F-7)           | 200 | 1 000 |
| Mixer (F-11)                     | 30  | 100 |
| FFT3D pool (F-12)                | 100 | 500 |
| Micro (F-6, F-9)                 | <10 | <10 |
| **Total non-spin**               | **~640 µs** | **~4 060 µs** |

**Spin 2×** on the grid-level and eigenpath findings; add F-13/F-14 →
**~1.2–5.5 ms/iter saved** at production scale.

**Context.** Post-FFTB eigensolve on Si n_pw=725 is ~72 ms/iter (see
performance-engineer logbook, 2026-04-18). Allocation savings of
1–4 ms/iter are **1.5–5% of current iter wall time**, enough to show
up as ~10–30% of the non-eigen portion. Post-WFRX (when eigensolve
drops), the relative share of allocations grows proportionally —
this is the motivation: ALOC is strategic cleanup *ahead of* WFRX, not
a replacement for it.

## §3 Fix sequence (highest savings per complexity first)

Ordered by **estimated savings / implementation complexity**, with the
constraint that correctness-adjacent changes land behind a workspace
abstraction introduced in step 1.

### Step 1. Introduce `ScfWorkspace` skeleton (prerequisite, no net savings)

**Scope.** A new module `src/scf/workspace.rs` with a `ScfWorkspace`
struct (and a sibling `ScfWorkspaceSpin` for nspin=2). Fields are the
buffers needed by the driver body, all preallocated in
`ScfWorkspace::new(n_grid, n_k, n_pw, n_bands)`:

```rust
pub(crate) struct ScfWorkspace {
    // Grid-level Complex64 buffers (n_grid each):
    pub rho_g: Vec<Complex64>,
    pub rho_g_new: Vec<Complex64>,
    pub v_h: Vec<Complex64>,
    pub vxc_g: Vec<Complex64>,
    pub v_eff: Vec<Complex64>,
    // Grid-level f64 buffers (n_grid each):
    pub exc_r: Vec<f64>,
    pub vxc_r: Vec<f64>,
    pub exc_r_out: Vec<f64>,
    pub vxc_r_out: Vec<f64>,
    pub rho_for_xc: Vec<f64>,      // used only when NLCC active
    // Symmetry workspace (2 × n_grid Complex64):
    pub sym_rho_g: Vec<Complex64>,
    pub sym_rho_g_out: Vec<Complex64>,
    // Per-k Hamiltonian matrices (n_k × n_pw²):
    pub h_cache: Vec<faer::Mat<Complex64>>,
    // Per-worker psi_g (n_grid) and FFT3D — stored via thread_local
    // or a rayon thread-id indexed Vec:
    pub psi_workers: Vec<WorkerLocal>,
}
```

Pass `&mut ScfWorkspace` into the driver as a second argument alongside
`&mut ScfContext`. Initial PR just threads the struct through; no
functions change behavior. CI should be bit-identical.

**Why first.** Every subsequent step plugs into this struct. Landing it
standalone keeps review burden low.

### Step 2. F-5 + F-7 (the big two)

**F-5 (Hamiltonian matrix cache)** and **F-7 (`psi_g` outside band
loop)** together deliver the bulk of the projected savings (400–3 000
µs/iter). Both are structurally local changes:

- F-5: extend `build_hamiltonian_with_v_eff` to a `fill` variant taking
  `&mut faer::Mat<Complex64>`. Zero the matrix, then run the existing
  kinetic and off-diagonal loops. Add `h_cache` to `ScfWorkspace`.
- F-7: hoist `psi_g` to the `fold` init closure in `compute_density`;
  reset with `iter_mut().for_each(|v| *v = Complex64::zero())` each
  band.

**Verification:** `tests/kb_projector_consistency.rs`,
`tests/qe_validation.rs`, and `tests/free_electron_bands.rs` all exercise
the SCF loop end-to-end. A workspace that reuses a Mat must produce
bit-identical eigenvalues and densities. Run the full
`cargo test --features gpu` before/after.

### Step 3. F-2 (XC grid in-place) + F-1 + F-4 + F-8

These four are structurally similar: change a function from
"allocate and return Vec" to "fill caller-provided slice". The
patterns:

- `lda_xc_grid_into(rho_r: &[f64], exc_out: &mut [f64], vxc_out: &mut [f64])`
- `hartree_on_fft_grid_into(rho_g: &[Complex64], g_squared: &[f64], out: &mut [Complex64])`
- `assemble_v_eff_into(v_local, v_h, v_xc, out: &mut [Complex64])`
- Make `density_r_to_g` the canonical XXX-r-to-g entry point (already
  `&mut` in signature); remove `real_to_g_space` wrapper that allocates
  fresh.

Land as a single PR (they touch the same workspace fields and share a
test: density+energy equivalence on Si Γ-only).

### Step 4. F-10 (symmetrize_density_g workspace + rot-transpose cache)

Move the two per-call allocations in `symmetrize_density_g` to the
workspace. Pre-compute the transposed rotation list in `ScfContext`
(it's a property of `symmetry`, not of iteration state).

**Complexity:** slightly higher because the function's public API
currently takes only `&mut [f64]` + a symmetry ref. Adding `workspace:
&mut SymmetryWorkspace` (or inlining the scratch into the `FFT3D`
struct — they're same-sized) requires touching the spin driver's two
call sites as well.

### Step 5. F-11 (Mixer workspace)

Introduce `MixerWorkspace { residual_r, residual_g, dv_scratch, df_scratch }`.
Rework Anderson/Broyden/PRPL mix paths to use slice-based intermediates.
Rewrite `precondition_residual` as a slice-in/slice-out transform.

**Complexity:** Anderson uses `ndarray::Array1` for its DIIS linear
algebra (dot products etc.), which allocates on subtraction. Either
rewrite those expressions to scalar accumulations or keep `Array1` for
the history slots but use raw slices for the transient residual. Prefer
the latter — less churn in the DIIS math, which is well-tested.

### Step 6. F-12 (FFT3D thread-local pool)

Introduce `FftPool` sized to `rayon::current_num_threads()` in
`ScfContext`, indexed in the `compute_density` closure.

**Complexity:** the indexing scheme needs to match rayon's thread
assignment. Either use `rayon::current_thread_index()` or
`thread_local!` macros. The simplest path: use a `Mutex<Vec<FFT3D>>`
pop/push pool (borrow, use, return). The locking overhead is cheaper
than allocating a `FftHandler` per call.

### Step 7. F-13 + F-14 (spin basis-change allocations)

Small PR; lands after F-11 because it reuses the same workspace
abstraction.

### Step 8. F-3 (Cow for add_core_density), F-6 (EigenResult unzip)

Cleanup; can land any time.

## §4 Benchmarking plan

### 4.1 Extend `benches/scf_benchmarks.rs`

Add a new bench group `scf_iter_end_to_end` that measures one full SCF
iteration body (not just the eigensolve) at production scale:

- `si_gamma_ecut200_grid24` — n_pw ≈ 283, n_grid = 13 824 (micro)
- `si_gamma_ecut400_grid32` — n_pw ≈ 893, n_grid = 32 768 (typical)
- `si_2x2x2_ecut400_grid32` — n_k = 4 (multi-k)
- `si_4x4x4_ecut400_grid32` — n_k = 10 (production)
- `fe_2x2x2_spin_ecut400_grid48` — spin driver, larger grid (stress)

For each: set up `ScfContext` and initial density outside `Criterion::iter`,
then time **one** iteration body. The setup cost (vnl_cache, v_local,
Ewald) is excluded.

**Acceptance metric:** µs/iter saved at each scale. Headline number is
the n_pw=725, n_grid=32 768 case.

### 4.2 Allocation-count verification (post-fix)

Once WFRX releases the machine lock, run under Instruments.app's
Allocations instrument for 5-iteration SCF trajectories on the
`si_4x4x4` config. Expected pre-fix count: ≥ 30 distinct allocation
sites per iter (mostly `Vec<_>::allocate`, `faer::Mat::alloc`,
`Array1::zeros`). Expected post-Step 2: all `faer::Mat::alloc` calls
inside the iteration body disappear; `Vec<_>::allocate` calls tied to
grid-level buffers drop by ~2/3.

The specific Instruments query: filter `Count` for allocations tagged
with `pwdft_rs::scf::driver` in the callstack, normalized by iteration
count.

### 4.3 Numerical equivalence gate (each fix)

Every PR in the ALOC sequence must pass:

- `cargo test` — all 265 tests green (tolerance unchanged).
- Si Γ-only SCF (`examples/si_scf.yaml`): final `total_energy` and
  `e_band` must match pre-PR to **machine epsilon** (< 1e-10 eV
  absolute). Workspace reuse is a refactor, not a physics change —
  any bit-drift is a bug.
- For F-7 (`psi_g` reuse): verify the density is bit-identical by
  comparing `rho_r` on the first, fifth, and converged iterations vs.
  a reference run.
- For F-11 (mixer workspace): verify the mixer's `final_delta`
  trajectory is identical iteration-by-iteration.

### 4.4 Watch items (not headline savings, but monitor)

- **Hot grid size that's larger than L2.** Above ~64³, allocations
  become *more* expensive (OS page walk) but the relative win of reuse
  is also higher. The bench should include at least one config with
  n_grid ≥ 48³ so the allocator regression mode is exercised.
- **`Mat` zeroing cost for F-5.** A 725² zeroed matrix is 8.4 MB of
  writes. If that's equivalent to the allocation cost, the fix is
  net-neutral (just moves bytes from the allocator into the clear).
  Prefer `fill` after explicit reset-by-touch; verify by adding a
  skip-the-zero microbench (should show ~15–30% speedup over
  allocate-and-fill).

## Appendix: non-findings explicitly ruled in-scope but clean

- **`basis.g_vectors()` / `basis.miller_indices()`** — both return
  `&[...]`, no per-call allocation. Fine.
- **`hartree_energy`, `lda_xc_energy`, `density_diff`** — pure
  reductions. Zero allocations. Fine.
- **`add_to_hamiltonian`** — post-VNLM, a single `matmul` into the
  caller's matrix. Allocation-free. Fine.
- **`NonlocalPotential::new`** — heavy allocator, but runs once at SCF
  entry (stored in `ScfContext.vnl_cache`). Out of scope.
- **`dense::diagonalize_lowest`** — faer internals are opaque to us;
  WFRX/ITEV address them directly.
- **String allocations in `info!` / `log::warn!` / `format!`.** Grep
  shows the loop body has 2 formatting sites: the per-iter progress log
  (`log_iteration`, emitted every iter — inherently allocating) and the
  rare post-convergence HF warning. The progress log uses `format!`
  once per iter; its cost is <5 µs and cannot easily be avoided without
  a log-line-skip mode. **No action; not worth the API churn.**
