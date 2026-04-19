# Vendored faer — `iterate_lanczos` MAX_REORTH fix

This repo vendors `faer` v0.24.0 from upstream
[`codeberg.org/sarah-quinones/faer`](https://codeberg.org/sarah-quinones/faer.git)
with **one** local edit: a bounded-retry cap on the full-reorthogonalization
loop inside `partial_self_adjoint_eigen_imp`. Everything else matches
upstream tag `v0.24.0`.

## Why we vendor

`faer::matrix_free::eigen::partial_self_adjoint_eigen` (the
implicitly-restarted Arnoldi partial eigensolver that `ITEV` calls) hangs on
near-null Krylov vectors. The upstream reorthogonalization loop has a
relative-tolerance convergence test that can never fire when `Vnext`
collapses into the span of the existing basis:

```text
converged[i] = r.abs() < f * Vnext.norm_l2();
```

If `Vnext.norm_l2()` → 0 (which happens when successive SCF iterations
produce near-identical Krylov vectors — e.g. late in a converging run on an
ill-conditioned Hamiltonian), both sides shrink together and the test never
trips. The loop spins forever.

The fix caps the retry count at `MAX_REORTH = 3` (Parlett & Kahan "Twice Is
Enough", 1966; Stewart *Matrix Algorithms Vol II* §5.2; Golub & Van Loan
§9.2.4). On breakdown, the existing post-loop norm check at the surrounding
scope breaks the outer `for j in ..` loop with a shorter Krylov basis, and
the caller's restart logic takes over.

## The diff

Single-file, 18-line addition inside
`faer/faer/src/operator/self_adjoint_eigen/mod.rs`:

```rust
let ref f = from_f64::<T::Real>(Ord::max(j, 8) as f64) * eps::<T::Real>();
// Parlett & Kahan, "Twice Is Enough" (1966); Stewart, "Matrix
// Algorithms Vol II" §5.2; Golub & Van Loan §9.2.4. Full
// reorthogonalization needs a retry bound — without it, a
// Krylov vector that collapses into the existing basis
// (norm_l2 → 0) leaves the relative-tolerance convergence test
// stuck because |r| and |Vnext| shrink together. Two passes
// are sufficient for non-degenerate inputs; three gives margin.
const MAX_REORTH: usize = 3;
let mut reorth = 0;
loop {
    reorth += 1;
    let mut all_true = true;
    for i in 0..j {
        // ... existing body unchanged ...
    }
    if all_true {
        break;
    }
    if reorth >= MAX_REORTH {
        // Breakdown: Vnext is (numerically) in the span of the
        // existing basis. The norm check immediately below will
        // see Vnext.norm_l2() near zero and break the outer
        // for-j loop; the caller restarts with the partial
        // subspace.
        break;
    }
}
```

The full patch is preserved at the end of this file for reference.

## Cargo wiring

`Cargo.toml` at repo root uses `[patch.crates-io]` so the vendored copy
takes over wherever the dep graph would otherwise pull `faer 0.24.0` /
`faer-traits 0.24.0` from crates.io:

```toml
[patch.crates-io]
faer = { path = "./faer/faer" }
faer-traits = { path = "./faer/faer-traits" }
```

Both crates need patching — our dep tree pulls `faer-traits` in
transitively, and cargo refuses to have two versions of a shared trait
crate (the vendored copy and the crates.io 0.24.0 one).

`Cargo.lock` records the patched versions. Running `cargo build` or
`cargo check` in any worktree just works: the path is relative to the
`Cargo.toml` file, and every worktree carries a checked-out copy of
`faer/` via normal git tracking.

## Editing the vendored faer

Treat `faer/` as first-class source. A Proposal that needs to change the
vendored solver (e.g. fix defect 1 of `ITEV`'s adaptive `n_request`
padding) can edit files under `faer/faer/src/...` directly inside its
worktree, run the pwdft-rs quality gate (`cargo test`, both clippy
invocations, `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps`), and commit
the faer-side + pwdft-rs-side changes together on its feature branch.

The vendored faer carries a separate `Cargo.lock` at `faer/Cargo.lock`.
Running `cargo test` inside `faer/` (to exercise faer's own test suite)
writes build artifacts to `faer/target/`, which is gitignored.

## Upstream-submit procedure

The fix applied above is not in upstream faer. When the user is ready to
submit it:

1. Clone upstream fresh:
   ```bash
   git clone https://codeberg.org/sarah-quinones/faer.git /tmp/faer-upstream
   cd /tmp/faer-upstream
   git checkout v0.24.0    # or `main` if more recent upstream work landed
   ```
2. Apply the diff at the end of this file. It applies cleanly on top of
   v0.24.0; for newer tags, inspect the context lines around
   `fn iterate_lanczos` first.
3. Commit with the original message (preserved below, authored by
   `RonLeizrowice-Pelanor <ron@pelanor.io>` on 2026-04-19).
4. Push to a fork, open a PR to `sarah-quinones/faer`.
5. When upstream merges + cuts a new release, bump `faer = ">=0.xx"` in
   this repo's `Cargo.toml`, delete the `[patch.crates-io]` block, and
   delete `faer/` from the tree. The vendored copy becomes redundant.

## Original commit metadata

```
commit 6a5edcd5744a5291f4e65d19385f3a4551f9e694
Author: RonLeizrowice-Pelanor <ron@pelanor.io>
Date:   Sun Apr 19 10:55:03 2026 +0800

    iterate_lanczos: bound reorthogonalization retries to avoid infinite loop

    The full-reorthogonalization loop in iterate_lanczos uses a relative
    convergence test `|r| < eps*sqrt(j) * |Vnext|`. When Vnext collapses
    into the span of the existing basis (norm_l2 → 0, which happens for
    near-degenerate Krylov subspaces), |r| and |Vnext| shrink together
    and the test never fires — the loop spins forever.

    Add a MAX_REORTH = 3 cap (Parlett & Kahan "Twice Is Enough"; Stewart
    Vol II §5.2; Golub & Van Loan §9.2.4). On breakdown the post-loop
    `if norm > zero() { ... } else { break; }` check at lines 60-67 already
    handles exit correctly — breaking the outer for-j loop with a shorter
    Krylov basis — so the caller's restart logic takes over.

    Reported by a downstream user (pwdft-rs) where SCF iterations on
    ill-conditioned Hamiltonians hung because successive SCF steps built
    near-identical Krylov vectors. Repro: call
    partial_self_adjoint_eigen inside a loop that re-uses the same
    operator (simulating converged SCF); the third or fourth call hangs.
```

## Full patch (as applied)

```diff
diff --git a/faer/src/operator/self_adjoint_eigen/mod.rs b/faer/src/operator/self_adjoint_eigen/mod.rs
index 48f4acf..c3bad88 100644
--- a/faer/src/operator/self_adjoint_eigen/mod.rs
+++ b/faer/src/operator/self_adjoint_eigen/mod.rs
@@ -39,7 +39,17 @@ fn iterate_lanczos<T: ComplexField>(
 		}
 		let ref f =
 			from_f64::<T::Real>(Ord::max(j, 8) as f64) * eps::<T::Real>();
+		// Parlett & Kahan, "Twice Is Enough" (1966); Stewart, "Matrix
+		// Algorithms Vol II" §5.2; Golub & Van Loan §9.2.4. Full
+		// reorthogonalization needs a retry bound — without it, a
+		// Krylov vector that collapses into the existing basis
+		// (norm_l2 → 0) leaves the relative-tolerance convergence test
+		// stuck because |r| and |Vnext| shrink together. Two passes
+		// are sufficient for non-degenerate inputs; three gives margin.
+		const MAX_REORTH: usize = 3;
+		let mut reorth = 0;
 		loop {
+			reorth += 1;
 			let mut all_true = true;
 			for i in 0..j {
 				if !converged[i] {
@@ -56,6 +66,14 @@ fn iterate_lanczos<T: ComplexField>(
 			if all_true {
 				break;
 			}
+			if reorth >= MAX_REORTH {
+				// Breakdown: Vnext is (numerically) in the span of the
+				// existing basis. The norm check immediately below will
+				// see Vnext.norm_l2() near zero and break the outer
+				// for-j loop; the caller restarts with the partial
+				// subspace.
+				break;
+			}
 		}
 		let norm = Vnext.norm_l2();
 		if norm > zero() {
```
