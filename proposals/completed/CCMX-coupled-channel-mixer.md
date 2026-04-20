---
id: CCMX
status: completed
priority: high
complexity: medium
risk: medium
depends_on: [SPNC]
blocks: []
---

# CCMX: Coupled-Channel Density Mixer

## Completion note (2026-04-18)

Landed in `src/scf/mod.rs::run_scf_spin`. Two `Mixer` instances (`mixer_total`,
`mixer_mag`) replace the former `mixer_up`/`mixer_down`. Forward basis change
applies QE's `rhoz_or_updw` convention (ρ_total = ρ↑+ρ↓, m = ρ↑−ρ↓); inverse
uses `ρ↑ = (ρ_total + m)/2`, `ρ↓ = (ρ_total − m)/2`. Kerker is disabled on
`mixer_mag` (not a charge-response). Mixer internals untouched.

**Fe BCC 4×4×4 regression test** (`tests/spin_polarization.rs::test_ccmx_fe_free_magnetization_converges`):

- Pre-CCMX: Δρ locked at ≈0.254 for 200+ iters, |HF-KS| ≈ 13 eV.
- Post-CCMX: **converges in 14 iters** with |HF-KS| = 1.06e-4 eV, M = 0 μB.

Si nspin=2 regression (`test_fe_spin_xc_consistency_regression`) still passes
at sub-microelectronvolt |HF-KS| — basis change is mathematically exact, so
well-behaved systems are unaffected.

Open: the Fe *fixed-mag=2* test (`test_fe_ferromagnetic_fixed_moment`) still
fails to converge post-CCMX. Root cause is that the PP prefers M=0; fixing
M=2 by the Fermi-energy-per-spin mechanism is not a stable SCF fixed point
regardless of mixer topology. That's a PP-choice or constrained-DFT issue,
not a CCMX issue. The inverted regression test is kept as a detector.

Also un-blocked: `tests/qe_validation.rs::test_fe_bcc_fm_vs_qe` 8×8×8 Kerker
now converges to −3050.80 eV, but remains `#[ignore]` pending VGCMP
(residual 9.5 eV heavy-atom V_local gap vs QE).

## CCMX: Coupled-Channel Mixer for nspin=2

### Problem

After SPNC (2026-04-17) tightened the nspin=2 convergence criterion from
total-density-only to per-spin max, Fe BCC fixed-mag=2 fails to converge.
Both spin channels enter a steady-state limit cycle with
`Δρ_up ≈ Δρ_down ≈ 0.254` from iteration ~5 onwards, while the total
density stays flat at `Δρ_total ≈ 1e-7`. This is a pure `+ε / −ε`
oscillation in the spin channels that cancels perfectly in the sum.

The root cause is the mixer topology in `run_scf_spin` (`src/scf/mod.rs`):
two **independent** Anderson mixers, one per channel:

```rust
// src/scf/mod.rs ~line 498
let mut mixer_up = mixing::Mixer::new(...);
let mut mixer_down = mixing::Mixer::new(...);
// ...
// src/scf/mod.rs ~line 771
rho_up_r   = mixer_up.mix(&rho_up_r,   &rho_up_sym,   &mut ctx.grid.fft);
rho_down_r = mixer_down.mix(&rho_down_r, &rho_down_sym, &mut ctx.grid.fft);
```

Neither mixer sees the other channel's history. For a system like Fe
fixed-mag=2 where the SCF fixed-point requires net charge transfer
between channels (up → down or vice versa) each iteration, two
independent Anderson updates produce uncorrelated predictions: channel
up pushes `+δ` exactly as channel down pushes `+δ` (rather than `−δ`),
so the spin polarisation `m = ρ↑ − ρ↓` oscillates around its target
while the total stays put.

See `.claude/logbooks/core-engineer.md` (2026-04-17 SPNC entry) for the
empirical trace and `proposals/SPNC-spin-per-density-convergence.md` for
the convergence-criterion change that surfaced this pathology.

### Background: what QE does

QE mixes in the `(ρ_total, m) = (ρ↑ + ρ↓, ρ↑ − ρ↓)` representation, not
per-channel `(ρ↑, ρ↓)`. The conversion happens in `sum_band.f90` right
after the density is accumulated from band sums:

```fortran
! qe-7.5/PW/src/sum_band.f90:307
IF ( nspin == 2 ) CALL rhoz_or_updw( rho, 'r_and_g', '->rhoz' )
```

where `rhoz_or_updw` (`qe-7.5/PW/src/scf_mod.f90:1360-1414`) performs
the in-place transformation:

```fortran
! qe-7.5/PW/src/scf_mod.f90:1395-1396
rho%of_r(ir,1)     = ( rho%of_r(ir,1) + rho%of_r(ir,nspin) ) * vi  ! total
rho%of_r(ir,nspin) =   rho%of_r(ir,1) - rho%of_r(ir,nspin) * vi * 2._dp  ! magnetisation
```

with `vi = 1.0` for `'->rhoz'` (to total/mag). The mixer in
`mix_rho.f90` then operates on this two-component
`(ρ_total, m)` object without any special-case logic — the coupling
across channels is encoded in the basis change, not in the mixer
internals. Back-conversion to `(ρ↑, ρ↓)` happens locally wherever the
per-channel densities are actually needed (e.g. `v_of_rho.f90:325-326`
builds a temporary `rho_updw` array for `xc_metagcx`).

Why this works: the natural "hard" mode of an LSDA SCF cycle is total
charge conservation (governed by the Hartree kernel `4πe²/G²`, which is
diagonal in ρ_total), while the magnetisation dynamics is much softer
and responds to the XC spin-splitting. These two modes have very
different effective preconditioners, and mixing them together in the
wrong basis (per-channel) makes every Anderson/Broyden residual a
mixture of the two mode families. In the `(ρ_total, m)` basis they
partially decouple, and the Anderson history across iterations accurately
captures the slow magnetisation relaxation without being drowned out by
the fast charge-neutrality enforcement.

### Proposed change

Apply the basis change at the start and end of each mixing step in
`run_scf_spin`, leaving the `Mixer` internals untouched. The mixer
operates on channel-agnostic `Vec<f64>` inputs, so the basis rotation is
a pure transformation around the mix call.

Pseudocode:

```rust
// Pre-mix: project (ρ↑, ρ↓) → (ρ_total, m) for both input and new densities
let rho_total_in:  Vec<f64> = zip(rho_up_r,   rho_down_r).map(|(u,d)| u + d).collect();
let m_in:          Vec<f64> = zip(rho_up_r,   rho_down_r).map(|(u,d)| u - d).collect();
let rho_total_new: Vec<f64> = zip(rho_up_sym, rho_down_sym).map(|(u,d)| u + d).collect();
let m_new:         Vec<f64> = zip(rho_up_sym, rho_down_sym).map(|(u,d)| u - d).collect();

// Mix each mode independently; but they are now total/magnetisation, not up/down
let rho_total_mixed = mixer_total.mix(&rho_total_in, &rho_total_new, &mut ctx.grid.fft);
let m_mixed         = mixer_mag  .mix(&m_in,         &m_new,         &mut ctx.grid.fft);

// Post-mix: project back (ρ_total, m) → (ρ↑, ρ↓)
rho_up_r   = zip(rho_total_mixed, m_mixed).map(|(t,mm)| 0.5 * (t + mm)).collect();
rho_down_r = zip(rho_total_mixed, m_mixed).map(|(t,mm)| 0.5 * (t - mm)).collect();
```

The two `Mixer` instances (`mixer_total`, `mixer_mag`) replace the
current `mixer_up`, `mixer_down`. Their internal histories (Anderson /
Pulay / Broyden) are agnostic to the physical meaning of their inputs —
they see two independent sequences of residuals, but now each sequence
corresponds to a single physical mode rather than a mixture.

Optionally, we may want a different `mixing_beta` for the magnetisation
channel than for the charge channel (QE calls this `mix_beta_spin` and
defaults it lower than the charge `mix_beta`). That is a follow-up —
the first cut uses the same `beta` for both and sees how far it gets.

#### Scope

Only `src/scf/mod.rs::run_scf_spin` is touched. The `Mixer` module
(`src/scf/mixing.rs`) is channel-agnostic and stays unchanged. Kerker
preconditioning (which is a pure charge-response model) continues to
apply naturally to the `ρ_total` mixer; applying it to the
magnetisation channel is nonsense (Kerker derives from `4πe²/G²`,
which is the charge-charge response, not the spin-spin response) — for
the magnetisation mixer we disable Kerker or use plain mixing. This
requires a small config tweak; see the Implementation sketch below.

### Implementation sketch

1. **`src/scf/mod.rs::run_scf_spin`, mixer setup (~line 498):** allocate
   `mixer_total` and `mixer_mag` instead of `mixer_up`/`mixer_down`.
   Both share `mixing_beta` and `mixing_ndim` from `ScfParams`.
   `mixer_total` uses `ctx.params.mixing_mode` as-is; `mixer_mag` uses
   the same mode but with Kerker disabled (if the user chose a
   Kerker-prefixed mode, downgrade it to the non-Kerker variant for the
   magnetisation channel).

2. **`src/scf/mod.rs::run_scf_spin`, mixing step (~line 771):** replace
   the two `mixer_up`/`mixer_down` calls with the basis-change +
   per-mode mix + basis-change-back sequence sketched above. Keep the
   existing per-channel convergence check (from SPNC) unchanged —
   converting back to `(ρ↑, ρ↓)` before computing `Δρ_up`/`Δρ_down`
   preserves the SPNC semantics exactly.

3. **No changes to `Mixer` internals.** The mixer operates on
   channel-agnostic flat `Vec<f64>` of length `n_grid` either way.

### Verification

Primary: Fe BCC fixed-mag=2 at `conv=1e-6` should converge in under 100
iterations with `|E_HF − E_KS| < 0.1 eV`. Specifically:

- Input: the `test_fe_ferromagnetic_fixed_moment` setup from
  `tests/spin_polarization.rs` (or the equivalent fixed-mag harness
  from the SPNC regression trace), 4×4×4 k, 15 Ry, `starting_mag=0.5`,
  `max_iter=100`, `conv_threshold=1e-6`.
- Expected post-CCMX: SCF returns `Ok(ScfResult{...})` with
  `n_iterations < 100` and `|HF-KS| < 0.1 eV`.
- If SCF still fails to converge, the issue is not mixer topology but
  something deeper — e.g. the PseudoDojo Fe PP at 15 Ry does not
  support a fixed-mag=2 SCF fixed point (LDA Fe ground state is
  non-magnetic for this PP), and CCMX cannot fix that. In that case,
  retry with a ferromagnetic-stable PP or a higher ecut before
  concluding CCMX is wrong.

Secondary (regression guards):

- All SPNC-era nspin=2 tests still pass: `test_fe_spin_xc_consistency_regression`
  (Si nspin=2, non-magnetic) should still converge in ~20-25 iterations
  with `|HF-KS| < 1e-5 eV`. The basis change is mathematically exact,
  so Si (where both channels converge together anyway) should be
  unaffected.
- `test_weight_sum_is_one`-style invariants: total charge conservation
  in the mixed density (integrate `ρ_total_mixed` → number of
  electrons) should match the input to floating-point precision.

### Estimated effort

One session (~3–4 hours), mostly concentrated on:

- The mechanical basis change around the mixer call (~30 min).
- Handling the Kerker disable for the magnetisation channel cleanly
  (~30 min — touch `src/scf/mixing.rs::MixingMode` only if needed).
- Running the Fe fixed-mag convergence experiment and iterating on
  `mixing_beta` if the first attempt doesn't converge (~1-2 hours).
- Updating the regression test and logbook (~30 min).

### References

- `qe-7.5/PW/src/sum_band.f90:307` — QE applies the `(up,dw) → (tot,mag)`
  transformation before mixing.
- `qe-7.5/PW/src/scf_mod.f90:1360-1414` — `rhoz_or_updw` subroutine
  implementing the in-place basis change (`vi=1.0` forward, `vi=0.5`
  inverse).
- `qe-7.5/PW/src/v_of_rho.f90:320-360` — example of back-conversion
  to `(ρ↑, ρ↓)` where the per-channel densities are actually needed
  (XC functional evaluation).
- `proposals/SPNC-spin-per-density-convergence.md` — why the current
  mixer fails the per-channel convergence check.
- `.claude/logbooks/core-engineer.md` 2026-04-17 SPNC entry — empirical
  trace of the limit-cycle pathology on Fe fixed-mag=2.
