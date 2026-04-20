---
id: DLTB
status: active
priority: high
complexity: trivial
risk: low
depends_on: []
blocks: []
---

# DLTB: Fix Convergence Failure Delta Bug

## Problem

When the SCF loop fails to converge, the reported `delta` in the error is always wrong:

**Unpolarized** (`src/scf/mod.rs:340-343`):

```rust
Err(PwdftError::ConvergenceFailure {
    iterations: ctx.params.max_iter,
    delta: density_diff(&rho_r, &rho_r, ctx.omega, ctx.n_grid), // rho_r vs rho_r = 0!
})
```

This compares `rho_r` with itself, which is always 0.

**Spin-polarized** (`src/scf/mod.rs:597-600`):

```rust
Err(PwdftError::ConvergenceFailure {
    iterations: ctx.params.max_iter,
    delta: 0.0, // hardcoded!
})
```

In both cases, users cannot tell how far from convergence the calculation was.

## Implementation

### Step 1: Track last delta in `run_scf`

Add `let mut last_delta = f64::INFINITY;` before the loop (near line 174). Inside the loop, after computing `delta` (line 262), add `last_delta = delta;`. Replace line 342:

```rust
// Before:
delta: density_diff(&rho_r, &rho_r, ctx.omega, ctx.n_grid),

// After:
delta: last_delta,
```

### Step 2: Track last delta in `run_scf_spin`

Add `let mut last_delta = f64::INFINITY;` before the loop (near line 403). Inside the loop, after computing `delta` (line 526), add `last_delta = delta;`. Replace line 599:

```rust
// Before:
delta: 0.0,

// After:
delta: last_delta,
```

## Verification

```bash
cargo test
```

Set `max_iter = 1` in a test to force convergence failure and verify the delta is non-zero and matches the last iteration's printed value.

## Estimated Effort

Under 15 minutes. Four lines changed.
