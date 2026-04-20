---
id: PRPL
status: completed
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# PRPL: Periodic Pulay Mixing

## Problem

The BROY proposal (`proposals/completed/BROY-broyden-mixing.md`, Step 4)
sketched a "periodic Pulay" mixer — plain linear mixing on every iteration
*except* every $k$-th, where a Pulay/Anderson extrapolation is applied — and
deferred it without implementation. The current code (`src/scf/mixing.rs`) has
three mixers (plain, Anderson with optional Kerker, modified Broyden) but no
periodic variant.

Motivation (what periodic Pulay fixes that the current mixers don't):

1. **Anderson/Broyden can diverge on early iterations** when the density
   residual is dominated by long-wavelength charge sloshing and the history
   window is still small (< 3 entries). Plain linear mixing with a conservative
   β is more robust at the start; switching to Pulay only after the residual
   has stabilized gives the best of both worlds.
2. **DIIS stagnation.** When the overlap matrix $\langle R_i | R_j \rangle$
   (`src/scf/mixing.rs:167-175`) becomes near-singular, the current Anderson
   falls back to uniform coefficients (`src/scf/mixing.rs:478-482`), which is
   equivalent to throwing away all Pulay information. Periodic Pulay resets
   cleanly instead of silently degrading to plain mixing with a bad β.
3. **Banerjee et al. (JCTC 12, 3053, 2016)** demonstrated 30–50% iteration-count
   reduction on transition-metal-oxide test cases vs. continuous Anderson. The
   paper is the canonical reference.

The hook for this in the current codebase is tiny — the mixer already owns
its history; we just gate when to use it.

## Background and references

- **Banerjee et al., JCTC 12, 3053 (2016), §2.2, Eq. 5.** Periodic Pulay
  with period $k$ and relaxation parameter β. The paper reports $k=3$–$5$
  works best for LDA/GGA on insulators and semiconductors, and $k=5$–$8$ for
  metals.
- **Kresse & Furthmüller, Comp. Mater. Sci. 6, 15 (1996), §4.2.** The
  original VASP mixer uses a "delay + DIIS" hybrid equivalent to periodic
  Pulay with $k=1$ after $n_0$ preparatory iterations. PRPL can be seen as
  a generalization where the delay repeats.
- **Fang & Saad, Numer. Linear Algebra Appl. 16, 197 (2009), §4.** Anderson
  acceleration is equivalent to multi-secant Broyden; periodic application
  preserves that equivalence on the Pulay steps while giving plain mixing
  the bookkeeping-free steps in between.

## Proposed Rust API

Add a new variant to `MixingMode` (`src/scf/mixing.rs:26-43`):

```rust
pub enum MixingMode {
    Plain,
    Kerker { q_tf: Option<f64> },
    Broyden { kerker: bool },
    /// Periodic Pulay (Banerjee et al., JCTC 12, 3053 (2016)).
    ///
    /// Plain linear mixing with β on every iteration except every `period`-th,
    /// where Anderson/Pulay extrapolation is performed using the accumulated
    /// history. Optional Kerker preconditioning applied to the linear step.
    PeriodicPulay {
        period: usize,   // k in the paper; default 3
        kerker: bool,
    },
}
```

And a new mixer struct next to `AndersonMixer` and `BroydenMixer`:

```rust
pub struct PeriodicPulayMixer {
    anderson: AndersonMixer,  // holds history + does the Pulay step
    period: usize,             // k
    iteration: usize,          // 1-based counter
    beta_linear: f64,          // β for the linear-mixing fallback
    kerker_weights: Option<Vec<f64>>,  // shared with the Anderson inner
}

impl PeriodicPulayMixer {
    pub fn mix(&mut self, rho_in: &[f64], rho_out: &[f64], fft: &mut FFT3D)
        -> Vec<f64>
    {
        self.iteration += 1;
        let raw_residual: Vec<f64> = rho_out.iter().zip(rho_in)
            .map(|(&o, &i)| o - i).collect();

        // Always accumulate history — cheap, and lets the Pulay step use it
        self.anderson.push_history(rho_in, &raw_residual);

        if self.iteration.is_multiple_of(self.period) && self.anderson.history_len() >= 2 {
            // Pulay step — solve DIIS and do Anderson mix
            self.anderson.diis_step(rho_in, rho_out, fft)
        } else {
            // Linear step — β·R, optionally Kerker-preconditioned
            let residual = if let Some(w) = &self.kerker_weights {
                precondition_residual(&raw_residual, w, fft)
            } else {
                raw_residual
            };
            rho_in.iter().zip(&residual)
                .map(|(&r, &rr)| self.beta_linear.mul_add(rr, r))
                .collect()
        }
    }
}
```

The only non-trivial extraction is splitting `AndersonMixer::mix`
(`src/scf/mixing.rs:126-190`) into two helpers:

- `push_history(rho_in, residual)` — appends and trims, no DIIS solve.
- `diis_step(rho_in, rho_out, fft)` — runs the existing DIIS linear-system
  solve using the already-accumulated history.

The current `AndersonMixer::mix` can then be rewritten as
`{ push_history(); diis_step(); }` — one-for-one equivalent.

Dispatch in `Mixer::new` (`src/scf/mixing.rs:377-397`):

```rust
MixingMode::PeriodicPulay { period, kerker } => Mixer::PeriodicPulay(
    PeriodicPulayMixer::new(beta, max_history, *period, *kerker, g_squared, n_electrons, omega),
),
```

### Settings

YAML (`src/settings.rs`):

```yaml
scf:
  mixing_mode: periodic_pulay
  mixing_beta: 0.3
  mixing_ndim: 8        # history depth for the Pulay step
  pulay_period: 3       # k; only read when mixing_mode == periodic_pulay
```

## Risk assessment

- **Low risk.** The algorithm is a strict superset of plain + Anderson. With
  `period = 1`, it reduces to plain Anderson; with `period = ∞`, it is plain
  linear mixing. Both limits are already in the test suite.
- **Settings shape.** Adding `pulay_period: usize` to `ScfParams` is a small
  additive change. If `mixing_mode != PeriodicPulay`, the field is ignored;
  the validator should not reject its presence.
- **History accumulation cost.** Accumulating history every iteration (not
  only on Pulay steps) uses $\sim 2 N_{\text{grid}}$ extra memory per slot;
  the default `max_history = 4` makes this negligible ($\sim 8 N_{\text{grid}}$
  doubles ≈ 2 MB at $128^3$).
- **No physics risk.** Periodic Pulay changes only the SCF convergence path,
  not the converged density. All energy-correctness tests pass unchanged at
  convergence.

## Verification plan

1. **Unit test — reduces to plain for `period = ∞`.**
   Construct with `period = usize::MAX` and verify results match
   `AndersonMixer::new(… Plain …)` with zero history (i.e., pure linear
   mixing) for 10 iterations.
2. **Unit test — reduces to Anderson for `period = 1`.**
   Same inputs on both a `PeriodicPulayMixer` with `period=1` and a plain
   Anderson mixer; assert bit-wise equality of outputs over 10 iterations.
3. **Unit test — history accumulates between Pulay steps.**
   Run 6 iterations with `period = 3`. After iteration 5, assert
   `anderson.history_len() == 5` (accumulated, not reset).
4. **Unit test — `period = 3` triggers Pulay exactly on iter 3, 6, 9.**
   Instrument a counter; assert the DIIS solve is called on those iterations
   and nowhere else.
5. **Convergence test — Si insulator.**
   Add to `tests/` or extend `src/scf/mixing.rs` test module: run Si ecut=100
   Γ-only for 30 iterations with `MixingMode::PeriodicPulay { period: 3, kerker: false }`
   and `MixingMode::Plain`, assert (a) converged energy matches within 1e-6 eV,
   (b) iter count for Pulay ≤ iter count for plain.
6. **Convergence test — Fe BCC (metallic).**
   Same pattern against Broyden: should converge in similar or fewer
   iterations (the paper's strongest claim is for metals / TMOs).
7. **No QE reference.** Periodic Pulay is a convergence-path optimization,
   not a physics change. Once converged, QE and pwdft-rs total energies
   remain the same apples-to-apples at the existing 0.05 eV tolerance.

## Implementation sketch

- `src/scf/mixing.rs`:
  - Refactor `AndersonMixer::mix` into `push_history` + `diis_step` (no
    behavior change; covered by existing tests).
  - Add `PeriodicPulayMixer` struct + `impl` block.
  - Extend `MixingMode` with `PeriodicPulay { period, kerker }`.
  - Extend `Mixer` enum + `Mixer::new` + `Mixer::mix` dispatch.
  - Four new unit tests (periods 1, 3, ∞; history accumulation).
- `src/settings.rs`:
  - Extend `MixingMode`-deserializer to accept `periodic_pulay`.
  - Add `pulay_period: Option<usize>` with default `Some(3)` when
    `mixing_mode == periodic_pulay`.
- `src/scf/context.rs:~140`: forward to `Mixer::new`.
- Optionally update `examples/*.yaml` with a commented-out
  `# pulay_period: 3` line so users can discover it.

## Estimated effort

3–4 hours:

- 0.5 h — refactor `AndersonMixer::mix` into two methods.
- 1 h — `PeriodicPulayMixer` struct + dispatch.
- 0.5 h — unit tests.
- 0.5 h — settings plumbing.
- 0.5 h — integration test.
- 0.5 h — docstrings + mixing.rs module doc update + clippy pass.

## Success criteria

1. `period = 1` matches Anderson bit-for-bit on 10 iterations (fixed seed).
2. `period = usize::MAX` matches plain linear mixing bit-for-bit.
3. On Si insulator Γ-only test, periodic-Pulay iter count ≤ plain Anderson
   iter count for `conv_threshold = 1e-6`.
4. Module docstring cites Banerjee et al. JCTC 12, 3053 (2016).
5. Default behavior unchanged — `mixing_mode: plain` users see no diff.
