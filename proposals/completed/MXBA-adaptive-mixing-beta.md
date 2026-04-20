---
id: MXBA
status: completed
priority: medium
complexity: medium
risk: medium
depends_on: []
blocks: []
outcome: landed-default-off
---

# MXBA: Adaptive Mixing Beta

## Completion note (2026-04-18)

Landed as opt-in (default `adaptive_beta = false`). Implementation uses the
Eyert 1996 §3.3 residual-norm monitor in `src/scf/mixing/mod.rs`
(`AdaptiveBeta` helper) plugged into both `AndersonMixer::push_history`
and `BroydenMixer::mix`. `PeriodicPulayMixer` inherits via its inner
Anderson. `Mixer::current_beta()` exposes the effective β for logging;
the SCF iteration line in both `run_scf` and `run_scf_spin` now prints
β (and `(tot, mag)` for nspin=2).

**Defaults chosen:**

- `growth_threshold = 1.2`, `damp_factor = 0.7`
- `restore_threshold = 0.5`, `restore_window = 3`
- `β_min = max(0.05·β_start, 0.01)`, clamp ceiling = β_start

**Default `adaptive_beta = false`.** Empirical finding: on Fe BCC CCMX
(the existing `test_ccmx_fe_free_magnetization_converges` regression)
adaptive β *hurts*. The residual plateaus at Δρ ≈ 0.34 for the first
~5 iters while Anderson accumulates DIIS history; the monitor reads
this as "slow convergence" and damps β all the way to β_min ≈ 0.017,
at which point the DIIS extrapolation is starved and the SCF cannot
escape. Fixed-β 0.3 + Anderson DIIS converges the same system in 14
iters. Documented failure: `tests/mxba_adaptive_beta_fe.rs`
(`#[ignore]`d, asserts the failure so a change to the monitor that
fixes this case will fire loudly). The user opt-in path leaves the
existing integration test suite byte-identical.

**Followup ideas (not in this PR):**

- Tune growth_threshold and restore_threshold on a metallic test case
  where β adaptation actually helps (e.g. C diamond 30 Ry plain mixing,
  per the proposal's original motivation — blocked on that system
  converging reproducibly).
- Consider skipping adaptive updates during a "DIIS warm-up" window
  (first `max_history` iterations) to avoid damping before the
  Anderson extrapolator has enough history to do useful work.
- Per-channel adaptive tuning for nspin=2 (total vs magnetization).

## Problem

The mixing parameter `mixing_beta` (`src/scf/mixing.rs:61`, `src/scf/mixing.rs:209`) is
a single user-specified scalar that is applied unchanged to every SCF iteration
of every mixer (Anderson, Kerker-Anderson, Broyden). A fixed β is never
simultaneously correct for the whole run:

- **Early iterations** (large residual, far from the fixed point): β should be
  small to prevent the residual norm from growing. A large β can cause charge
  sloshing — especially for metals and for transition-metal oxides where the
  Jacobian eigenvalues of the SCF map span many orders of magnitude.
- **Late iterations** (small residual, near fixed point): β should be close to 1
  because the approximate inverse Jacobian (Anderson coefficients, Broyden rank-m
  update) is already doing most of the work. A small β throws away that work and
  slows final convergence.

Two concrete symptoms in this codebase:

1. **C diamond at ecut=30 Ry does not converge in 80 iterations with
   `MixingMode::Plain`** (researcher logbook, 2026-04-17 — QEVL notes). A smaller
   β in the first 5 iterations is one standard remedy.
2. **Fe BCC (nspin=2) at ecut=15 Ry collapses to NM under PseudoDojo** (same
   logbook). Metals with exchange splitting benefit from per-channel β scaling.

Adaptive mixing β is a cheap ($<10$ flop/iter) addition that handles both.

## Background and references

No published adaptive-β scheme is universally standard — this is a place where
different codes diverge.

- **VASP** exposes `AMIX`, `BMIX`, `AMIX_MAG`, `BMIX_MAG` and the Kerker model
  but holds them fixed during the run; adaptivity is manual. Reference:
  [VASP wiki / IMIX](https://www.vasp.at/wiki/index.php/IMIX).
- **ABINIT** `iscf=7` (Anderson) and `iscf=17` (RMM-DIIS) can use a variable
  `diemix` via `diemixmag`, again manually tuned.
- **QE** `mix_rho.f90:15` receives `alphamix` as `INTENT(IN)` and does **not**
  adapt it across iterations. Adaptivity in QE is done *outside* `mix_rho`
  in user-scripted restart loops.

The implementable schemes are from the numerical-analysis literature:

1. **Residual-monitor rule** (Eyert, J. Comp. Phys. 124, 271 (1996), §5):
   multiply β by $\gamma_+ > 1$ when $\|R_n\| < c_+ \|R_{n-1}\|$ with
   $c_+ \approx 0.8$, and by $\gamma_- < 1$ when $\|R_n\| > c_- \|R_{n-1}\|$
   with $c_- \approx 1.0$. Clamp β to $[\beta_{\min}, \beta_{\max}]$. Simple,
   robust, 4 flops/iter.
2. **Optimal-damping for Broyden** (Marks & Luke, PRB 78, 075114 (2008), Eq. 21):
   minimize $\|R_{n+1}\|^2$ along the search direction by one line-search step.
   Requires one extra Hamiltonian build; more expensive.
3. **Periodic reset with β_0** (Banerjee et al., JCTC 12, 3053 (2016)):
   every $k$-th iteration reset β to a small "safe" value, otherwise use
   β_large. Related to the PRPL (periodic Pulay) proposal — see below.

This proposal targets **scheme 1** because it is O(1) per iteration, has no
physics assumptions, and plugs cleanly into the existing `Mixer` enum. Scheme 2
is deferred as a follow-up (adds a Hamiltonian build per iteration; needs
benchmarking). Scheme 3 is proposal PRPL.

## Proposed Rust API

### Option A — wrapper struct (recommended)

Add an `AdaptiveBeta` wrapper that holds the previous residual norm and mutates
the inner mixer's β each call:

```rust
// src/scf/mixing.rs, new type

/// Residual-monitor adaptive β (Eyert, J. Comp. Phys. 124, 271 (1996), §5).
///
/// Adjusts the mixing parameter between iterations based on residual norm:
/// - If ‖R_n‖ ≤ c_down · ‖R_{n-1}‖: multiply β by γ_up (converging well)
/// - If ‖R_n‖ ≥ c_up   · ‖R_{n-1}‖: multiply β by γ_down (diverging or flat)
/// - Otherwise: keep β unchanged
///
/// Clamps β to [β_min, β_max]. Typical values:
///   γ_up=1.2, γ_down=0.5, c_down=0.8, c_up=1.0, β_min=0.02, β_max=0.8.
pub struct AdaptiveBeta {
    beta_min: f64,
    beta_max: f64,
    gamma_up: f64,
    gamma_down: f64,
    c_down: f64,
    c_up: f64,
    prev_residual_norm: Option<f64>,
}

impl AdaptiveBeta {
    /// Construct with literature defaults from Eyert 1996.
    #[must_use]
    pub fn eyert_defaults() -> Self { /* ... */ }

    /// Update and return the new β based on the current residual norm.
    pub fn step(&mut self, residual: &[f64], current_beta: f64) -> f64 { /* ... */ }
}
```

Inject into both `AndersonMixer::mix` and `BroydenMixer::mix` at the
point where the raw residual is already computed — `src/scf/mixing.rs:129`
(Anderson) and `src/scf/mixing.rs:278` (Broyden). Replace `self.beta` with
`self.adaptive_beta.step(&raw_residual, self.beta)` for both call sites.

### Option B — dispatch via MixingMode

Reject. Would duplicate `MixingMode::Plain { adaptive: bool }` and
`MixingMode::Broyden { kerker: bool, adaptive: bool }`, inflating the cartesian
product. Option A keeps the mixer algorithm and the β schedule orthogonal.

### Settings

`src/settings.rs` `ScfSettings`:

```yaml
scf:
  mixing_mode: broyden   # unchanged
  mixing_beta: 0.3       # initial β
  adaptive_beta:         # optional table; if absent, β stays fixed
    enabled: true
    beta_min: 0.02
    beta_max: 0.8
    gamma_up: 1.2
    gamma_down: 0.5
    c_down: 0.8
    c_up: 1.0
```

When `adaptive_beta.enabled = false` (or the table is absent), the current
fixed-β behavior is preserved byte-for-byte. This is the backward-compatibility
invariant.

## Risk assessment

- **Behavioral change when enabled.** Running existing test cases with
  `adaptive_beta: true` will produce different iteration counts and
  intermediate densities. Final converged energy should be unchanged.
  Mitigated by: default OFF; regression tests compare against fixed-β golden
  numbers only.
- **Stability edge cases.** Oscillating residual norms (e.g., when Pulay
  coefficients go through a near-singular region) can cause β to bounce
  between $\beta_{\min}$ and $\beta_{\max}$. Mitigation: hard clamp + a
  hysteresis window ($c_{\text{down}} = 0.8$, $c_{\text{up}} = 1.0$: strictly
  non-overlapping so no oscillation near the identity case).
- **Harris-Foulkes sanity check.** With adaptive β the SCF still converges to
  the same fixed point, so $|E_{HF} - E_{KS}| \to 0$ at the same rate in the
  convergent tail. Near-divergence runs where current code fails outright will
  now converge with $\Delta$β-dependent $|E_{HF} - E_{KS}|$ trajectories — fine,
  because they converge.
- **Interaction with Kerker / Broyden.** Both preconditioners operate on
  residuals; adaptive β changes the *step size* along the preconditioned
  direction, not the direction itself. No mathematical coupling. Ship the
  adaptive β as a pure wrapper; the existing tests for Kerker and Broyden are
  still valid.

## Verification plan

1. **Unit test — monotone convergent residuals increase β.**
   Feed a sequence $\|R_n\| = 0.5^n$ to `AdaptiveBeta::step` and verify β
   grows monotonically until clamped at `beta_max`.
2. **Unit test — diverging residuals reduce β.**
   Feed $\|R_n\| = 1.2^n$ and verify β monotonically drops to `beta_min`.
3. **Unit test — hysteresis window.**
   Feed $\|R_n\| = \|R_{n-1}\| \cdot 0.9$ (in between `c_down` and `c_up`)
   and verify β is unchanged.
4. **Regression — backward compatibility.**
   `test_kerker_vs_plain_scf_convergence` and
   `test_broyden_vs_plain_scf_convergence` (`src/scf/mixing.rs:715`,
   `src/scf/mixing.rs:904`) must produce the same converged energies when
   `adaptive_beta` is disabled (default). Add a companion test that turns
   adaptive β on and asserts the same converged energy within 1e-6 eV.
5. **Convergence test — C diamond.**
   Add `tests/adaptive_beta_c_diamond.rs` reproducing the researcher logbook
   case (ecut=30 Ry, `MixingMode::Plain`, `mixing_beta: 0.3`, 80 iters) and
   check that (a) plain mixing fails to converge, (b) plain + adaptive β
   converges in ≤ 80 iterations. This makes the concrete "adaptive β helps"
   claim falsifiable and rerunnable.
6. **Convergence test — Si insulator.**
   Run the existing `tests/qe_validation.rs` Si case with adaptive β ON
   and verify total energy matches QE within the existing 0.05 eV tolerance
   (once VERF blockers lift — currently `#[ignore]`'d).
7. **No QE comparison needed** for the schedule itself — it is a black-box
   numerical trick, not a physics change. QE does not adapt β so there is no
   reference value to match.

## Implementation sketch

- `src/scf/mixing.rs`:
  - New `pub struct AdaptiveBeta` with `Default` for Eyert defaults.
  - New field `adaptive_beta: Option<AdaptiveBeta>` on both `AndersonMixer`
    and `BroydenMixer`. Wire `Mixer::new` at `src/scf/mixing.rs:377` to
    accept an `Option<AdaptiveBeta>` and propagate.
  - Two-line change at `src/scf/mixing.rs:129` (Anderson) and
    `src/scf/mixing.rs:278` (Broyden): call `.step(&raw_residual, self.beta)`
    and overwrite `self.beta` before the rest of the routine uses it.
- `src/settings.rs`: new optional `adaptive_beta: Option<AdaptiveBetaSettings>`
  in `ScfSettings`.
- `src/scf/context.rs:~140`: forward the setting to `Mixer::new`.
- Three unit tests in `src/scf/mixing.rs` `#[cfg(test)] mod tests`.
- One integration test: `tests/adaptive_beta_c_diamond.rs`.

## Estimated effort

4–6 hours:

- 1 h — `AdaptiveBeta` struct + unit tests.
- 1 h — wire into `Mixer::new`, `AndersonMixer::mix`, `BroydenMixer::mix`.
- 1 h — settings plumbing (`settings.rs`, `ScfParams`, `ScfContext`).
- 1 h — integration test (C diamond).
- 1 h — documentation pass (mixing.rs module docstring, CLAUDE.md mention).
- 1 h — clippy / test-suite stabilization.

## Success criteria

1. `adaptive_beta: false` (default) reproduces existing energies bit-for-bit.
2. `adaptive_beta: true` converges C diamond at ecut=30 Ry inside 80 iters.
3. All existing mixing tests pass unchanged when adaptive β is off.
4. Documentation: mixing.rs module docstring explains the Eyert rule and cites
   J. Comp. Phys. 124, 271 (1996).

## Note (2026-04-18, MXB1 verification)

The coded recipe in `src/scf/mixing/mod.rs` (hysteresis band `[restore_threshold,
growth_threshold] = [0.5, 1.2]`, multiplicative `damp_factor = 0.7` applied on
`ratio > 1.2`, symmetric restore `β ← β/damp_factor` only after a
`restore_window = 3` streak of `ratio < 0.5`) does **not** match the §5
formula quoted above (γ_up=1.2, γ_down=0.5, c_down=0.8, c_up=1.0, damp on
`ratio > 1.0`). It matches the §3.3 "residual-norm monitor" sketch that the
implementation's own docstrings (`src/scf/mixing/mod.rs:31-42, 80-87`;
`anderson.rs:45`; `broyden.rs:42, 121`) consistently cite. The §5 block in
this archived proposal is a draft residue — the constants were superseded
during implementation and the citation was not updated.

**Paper full text was not accessible** (ScienceDirect preview only; no
open preprint; secondary sources confirm the paper exists but don't
reprint §3.3/§5). So the section label "§3.3" is carried over from the
code's own comments, and the *exact* equation number remains unverified.
The constants themselves appear empirically tuned (no sweep script
currently recorded) — recommend MXB2 author add
`scripts/validate/mxba_eyert_sweep.py` when re-tuning. Escalate to EM
if a verbatim paper quote is required.
