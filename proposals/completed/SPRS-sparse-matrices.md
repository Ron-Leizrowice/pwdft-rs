---
id: SPRS
status: archived
priority: low
complexity: large
risk: medium
depends_on: []
blocks: []
archived_on: 2026-04-19
archived_reason: Sparsity assumption is wrong for this codebase. The plane-wave Hamiltonian has O(n_pw²) dense non-zeros — kinetic is diagonal but V_eff is a convolution in G-space, which fills the matrix. Sparse storage only pays off in real-space-grid or atomic-orbital bases, neither of which this project uses or plans to use. Archived without implementation.
---

# SPRS: Sparse Matrix Support — ARCHIVED

> **ARCHIVED 2026-04-19.** Sparsity assumption is wrong for plane-wave DFT: the KS Hamiltonian in a plane-wave basis is structurally dense because V_eff enters as a G-space convolution. Sparse storage would pay off only in real-space-grid or atomic-orbital bases; this project uses neither. Keep this file for historical context only.
>
> **Note:** Line numbers reference the pre-ScfContext codebase. Proposal 01 (faer) is now completed — faer's own `SparseColMat` may be preferable to adding sprs.

## Motivation

The Hamiltonian is currently a dense `DMatrix<Complex64>` (or `faer::Mat<c64>` if Proposal 01 is adopted). For the current test case (Si with ecut=100, ~50-200 basis functions), dense storage is optimal. But as the code targets larger systems:

| System | Atoms | ecut (eV) | n_pw | Dense H (MB) | Eigensolve time |
|--------|-------|-----------|------|---------------|-----------------|
| Si (2 atoms) | 2 | 100 | ~60 | 0.06 | <1ms |
| Si (2 atoms) | 2 | 400 | ~300 | 1.4 | ~10ms |
| Si (8 atoms) | 8 | 400 | ~1200 | 22 | ~1s |
| MgO slab | 32 | 500 | ~5000 | 380 | ~30s |

At n_pw > 1000, the dense eigensolve dominates and memory becomes significant. However, the Hamiltonian is often sparse because:

- **Kinetic energy** is purely diagonal (always sparse)
- **V_eff(G-G')** decays for large |G-G'| (potential is smooth in real space)
- **Non-local potential** has rank at most n_atoms * n_projectors (typically 4-20), so V_NL is a low-rank update

This structure is exploitable by iterative eigensolvers, which need only the matrix-vector product `H|v>` rather than the full matrix.

## Dependencies

Add:

```toml
sprs = { version = ">=0.11", optional = true }
```

Behind a feature flag:

```toml
[features]
sparse = ["dep:sprs"]
```

## Scope of Changes

### Phase 1: Sparse Hamiltonian storage (optional path)

**New file: `src/eigensolver/sparse.rs`**

Define a sparse Hamiltonian representation:

```rust
use sprs::CsMat;
use num_complex::Complex64;

/// Sparse Hamiltonian for iterative eigensolvers.
///
/// H = T (diagonal) + V_eff (sparse) + V_NL (low-rank)
pub struct SparseHamiltonian {
    /// Kinetic energy: diagonal, stored as Vec<f64>
    kinetic: Vec<f64>,
    /// V_eff: sparse matrix (only entries with |V_eff(G-G')| > threshold)
    v_eff: CsMat<Complex64>,
    /// Non-local: stored as projector vectors for matrix-free application
    /// V_NL |v> = Σ_i |β_i> D_ij <β_j|v>
    projectors: Vec<Vec<Complex64>>,  // n_proj vectors of length n_pw
    dij: Vec<Vec<f64>>,               // D_ij coupling matrices
}

impl SparseHamiltonian {
    /// Matrix-vector product H|v> without forming the full matrix.
    pub fn apply(&self, v: &[Complex64], result: &mut [Complex64]) {
        let n = v.len();

        // Kinetic (diagonal)
        for i in 0..n {
            result[i] = Complex64::new(self.kinetic[i], 0.0) * v[i];
        }

        // V_eff (sparse matrix-vector multiply)
        // sprs provides efficient SpMV
        let v_eff_v = &self.v_eff * v;
        for i in 0..n {
            result[i] += v_eff_v[i];
        }

        // Non-local (low-rank: project, multiply D, unproject)
        // Cost: O(n_pw * n_proj) instead of O(n_pw^2)
        for (projectors_for_type, dij) in ... {
            let projections: Vec<Complex64> = projectors_for_type.iter()
                .map(|beta| beta.iter().zip(v).map(|(b, v)| b.conj() * v).sum())
                .collect();
            // Apply D_ij and accumulate
            for (i, beta_i) in projectors_for_type.iter().enumerate() {
                let coeff: Complex64 = dij.iter()...;
                for (k, &b) in beta_i.iter().enumerate() {
                    result[k] += b * coeff;
                }
            }
        }
    }
}
```

### Phase 2: V_eff sparsification

**File: `src/scf/mod.rs` — `build_hamiltonian_with_v_eff` (lines 411-440)**

Currently iterates over all n^2 pairs (lines 429-436). Add a threshold to skip small elements:

```rust
let v = v_eff_fft[fft_idx];
if v.norm() > 1e-12 {
    // insert into sparse matrix builder
    triplets.push((i, j, v));
}
```

For smooth potentials, this typically retains only 10-30% of entries.

### Phase 3: Integration with iterative eigensolver

This proposal naturally pairs with the iterative eigensolver path in Proposal 01 (faer Phase 4). The sparse Hamiltonian's `apply` method provides the matrix-vector product needed by LOBPCG or Davidson solvers. The key advantage: cost per iteration drops from O(n_pw^2) (dense matrix-vector) to O(nnz + n_pw * n_proj) where nnz is the number of nonzero V_eff entries.

## When This Matters

This proposal is **premature for the current codebase** (n_pw < 300). It becomes valuable when:

- n_pw > 1000 (systems with > 8 atoms at moderate cutoff)
- An iterative eigensolver is implemented (otherwise sparse storage doesn't help — dense zheev needs the full matrix regardless)

The natural sequencing is: **Proposal 01 (faer) -> iterative eigensolver -> this proposal**.

## Risks

- Sparse matrix overhead (index storage, cache-unfriendly access) can make small-system performance worse. Only use for n_pw above a crossover threshold (~500-800).
- V_eff sparsification with a threshold introduces a controllable approximation. Must validate that converged energies are unchanged at the chosen threshold.
- `sprs` is less actively maintained than `faer`. If faer (Proposal 01) is adopted, faer's own `SparseColMat` type may be preferable to adding a separate sparse matrix crate.

## Expected Impact

- **Scaling:** Enables calculations on systems with 1000+ plane waves where dense storage is prohibitive.
- **Performance:** Combined with an iterative eigensolver, reduces eigensolve scaling from O(n^3) to O(n * n_bands * n_iter_krylov).
- **Memory:** Reduces Hamiltonian storage from O(n^2) to O(nnz).
- **Priority:** Low for current development. File for when larger systems are targeted.
