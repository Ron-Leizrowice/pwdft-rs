# Proposal 09: Kerker Preconditioning

**Status: COMPLETED**

## Result

Implemented Kerker preconditioning as a G-space filter on the density residual before Anderson mixing. P(G) = |G|²/(|G|² + q_TF²) suppresses the dangerous long-wavelength charge sloshing that plagues metals and large cells.

Key design decision: kept the mixer in real space (simpler, compatible with existing code) and applied the preconditioner as an FFT→filter→IFFT step on the residual. This avoids rewriting the entire Anderson algorithm in complex arithmetic.

## Changes

- `src/scf/mixing.rs`: `MixingMode` enum (Plain/Kerker), `AndersonMixer` takes g_squared and computes Kerker weights at construction, `precondition_residual()` applies filter via FFT
- `src/scf/mod.rs`: `ScfParams.mixing_mode` field, mixer construction passes g_squared
- `src/main.rs`: defaults to Kerker with auto q_TF
- Auto q_TF from Thomas-Fermi formula: q_TF² = 4(3π²ρ)^{1/3}/π

## Tests

- `test_kerker_suppresses_g0`: uniform residual (pure G=0) heavily suppressed by Kerker
- `test_auto_q_tf_reasonable`: auto q_TF ∈ [0.5, 5.0] Å⁻¹ for Si density
- Plain mixing backward compatible (all existing tests pass unchanged)
