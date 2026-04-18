---
id: VNLM
title: V_NL Hamiltonian assembly via single GEMM
priority: medium
complexity: small-medium
risk: low
depends_on: []
blocks: []
status: active
author: Performance Engineer
date: 2026-04-18
---

# VNLM — V_NL Hamiltonian assembly via single GEMM

## Problem

`NonlocalPotential::add_to_hamiltonian` (`src/potential/nonlocal.rs:118-208`) is
the #2 SCF per-iteration cost after the dense eigensolve (post-FFTB/FMAD
profiling, Apple M2, machine-locked):

| n_pw | vnl_apply (current) |
|------|---------------------|
| 89   | 398 µs              |
| 259  | 3.18 ms             |
| 725  | 27.4 ms             |

At n_pw = 725 this is ~29 ms per k-point per SCF iteration, i.e. ~3 % of the
eigensolve itself (836 ms) but the #2 line item. On Si/ecut = 200 with
n_kpt = 10 and ~15 SCF iterations, this is ~4 s of wall time — worth claiming.

## Root cause

The current implementation expands V_NL using the compact KB identity
`Σ_m Y_lm(q̂) Y*_lm(q̂') = (2l+1)/(4π) P_l(cos θ)`, which folds the angular
dependence into a Legendre polynomial of `cos θ_{G,G'}`. This makes the inner
summand *non-separable* in `(G, G')`, forcing a doubly-nested O(n_pw²) loop
with an O(n_proj²) kernel body:

```
for ig in 0..n_pw {
    for jg in 0..n_pw {
        // O(n_proj²) projector sum with P_l(q̂_i · q̂_j)
        for i,j in projectors { vnl += F_i(|q_i|) · D_ij · F_j(|q_j|) · angular(ig,jg) }
        h[(ig,jg)] += sf_sum · (vnl / Ω)
    }
}
```

## Fix

Restore the factorized form that QE uses
(`qe-7.5/upflib/ylmr2_gpu.f90` and `init_us_2_acc.f90`). Define the KB
projector matrix indexed by an expanded channel `(atom, radial_projector, m)`:

```
B[G, α] = (1/√Ω) · exp(−iG·τ_α) · F_{i_α}(|k+G|) · Y_{l_α m_α}(q̂_{k+G})
```

Then

```
H_NL = B · D · B^H
```

where `D` is block-diagonal: nonzero entry `D[(a,i,m), (a',j,m')]` requires
`a == a'`, `l_i == l_j`, `m == m'`. Two GEMMs lift the cost from
`O(n_pw² · n_proj²)` to `O(n_pw · n_channels + n_pw² · n_channels)` (still
n_pw² because the *output* matrix has n_pw² entries — but the inner kernel is
a BLAS-3 GEMM, not nested scalar loops).

Implementation sketch:

```rust
// One-time per call: build B (n_pw × n_channels).
let mut b = Mat::<Complex64>::zeros(n_pw, n_channels);
// fill with exp(-iG·τ) · F_α(|k+G|) · Y_lm(q̂) / √Ω

// Build per-atom/ℓ/m block of D_α · B^H (n_channels × n_pw).
let mut db_h = Mat::<Complex64>::zeros(n_channels, n_pw);
for each (atom, l, m) block { db_h_block = D_block · b_block.adjoint() }

// Single GEMM into H (n_pw × n_pw).
faer::linalg::matmul::matmul(h.as_mut(), Accum::Add, b.as_ref(), db_h.as_ref(),
                              Complex64::ONE, par);
```

The identity `(2l+1)/(4π) P_l(cos θ) = Σ_m Y_lm Y*_lm` is exact in exact
arithmetic; in f64 it matches to ~1e-14, far below every existing tolerance
(V_NL hermiticity 1e-10, diagonal shell degeneracy 1e-8, V_NL cross-check
against Python reference 1e-6).

## Baseline (Apple M2, machine-locked, commit origin/main@9a5e9e9)

criterion `hamiltonian/vnl_apply_*` on Si FCC:

| n_pw | current |
|------|---------|
| 89   | 398 µs  |
| 259  | 3.18 ms |
| 725  | 27.4 ms |

## Projected speedup

GEMM-vs-scalar-loop speedups on Apple M2 (NEON via faer/gemm):
- ~5-10× at n_pw = 89 (GEMM setup overhead limits the low end);
- ~10-20× at n_pw = 725 (loop cost fully dominates).

Target: n_pw = 725: 27 ms → ~3 ms.

## Measured result (Apple M2, machine-locked, branch `VNLM/blocked-matmul`)

| n_pw | vnl_apply before | vnl_apply after | speedup | vnl_new before | vnl_new after |
|------|------------------|-----------------|---------|----------------|---------------|
| 89   | 398 µs           | 80 µs           | **5.0×** | 5.41 ms        | 7.08 ms       |
| 259  | 3.18 ms          | 776 µs          | **4.1×** | 15.5 ms        | 20.95 ms      |
| 725  | 27.4 ms          | 5.24 ms         | **5.2×** | 43.2 ms        | 78.5 ms       |

`vnl_new` regressed by ~1.4-1.8× because we now build the full KB projector
matrix `B` and pre-apply `D` at construction. This is a one-time cost per
k-point (cached in `ScfContext.vnl_cache`); `vnl_apply` runs once per SCF
iteration. Break-even is at ~2 SCF iterations — any real SCF run is a net win.

Net V_NL cost per k-point over 15 SCF iterations at n_pw = 725:

- Before: `43 + 15 × 27.4 = 454 ms`
- After : `78 + 15 × 5.24 = 157 ms`
- **Overall: 2.9× faster**

## Validation plan

1. All existing tests pass with original `relative_eq!` tolerances. Critical
   ones: `tests/nonlocal_symmetry.rs` (hermiticity ≤ 1e-10, shell degeneracy
   ≤ 1e-8), `tests/kb_projector_validation.rs`, `tests/qe_validation.rs` (Si
   total energy pinned).
2. Bit-level cross-check test: run old vs new on a fixed Si H at Γ and an
   off-Γ k-point, assert per-element max |ΔH| < 1e-10. Added in the PR.
3. `cargo test` on both feature flags; both `cargo clippy` invocations.
4. `cargo bench --bench scf_benchmarks -- hamiltonian/vnl_` before + after.

## Scope

- `src/potential/nonlocal.rs` only. Adds a small `ylmr.rs` helper submodule
  or inlined real-spherical-harmonic routine (QE-style recurrence).
- Does NOT touch GPU V_NL (if any future GPU V_NL lands, it follows the same
  math; this proposal is CPU-only).
- Does NOT touch `scf/mod.rs`, mixing modules, or symmetry (concurrent MODR
  / MXBA work).

## Risks & mitigations

- **Numerical drift in pinned tests.** Mitigation: the addition theorem
  is exact; f64 errors from `sin/cos/sqrt` stay at ULP. If any pinned
  number shifts meaningfully (> 1e-8 relative), that's a bug — stop and
  investigate.
- **Memory.** B is n_pw × n_channels complex (~8 n_pw × n_channels × 16 B).
  For Si (n_channels = 2 atoms × 4 projectors = 8) at n_pw = 725, that is
  ~90 KB — negligible. DB^H is the same size.
- **Complex Y_lm via recurrence** (per-G point, l=0..lmax): O((lmax+1)²) work,
  cheap. For Si PP `lmax = 1` so it's 4 values per G; for Fe GBRV it's
  `lmax = 2` (9 values per G). Well within bandwidth.
