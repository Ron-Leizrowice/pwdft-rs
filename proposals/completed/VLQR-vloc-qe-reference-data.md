---
id: VLQR
status: completed
priority: medium
complexity: small
risk: low
depends_on: []
blocks: []
---

# VLQR: V_local QE Reference Data Re-extraction

## Origin

TAUD PR D (PR #26, 2026-04-17) un-ignored `tests/kb_projector_validation.rs::test_vloc_comparison_with_qe` and observed:

- V_local(G=0): ours **+1.343 eV** vs QE Cube reference **−1.003 eV** (sign flip; |Δ| = 2.35 eV)
- |V_local(G=(1,0,0))|: ours 5.468 vs QE 6.968 (|Δ| = 1.50 eV)
- |V_local(G=(1,1,1))|: ours 5.468 vs QE 6.968 (|Δ| = 1.50 eV)

This contradicts **VGCMP Phase 1**, which independently proved our `V_local(G)` matches a Python Simpson-rule reference computed directly from the UPF file to **machine precision** (max |Δ| = 2.78e-9 Ry across Si's first 20 G-shells).

**Conclusion:** the QE reference data inside this test is wrong, not our Rust code. The test currently loads the QE side via a Cube file extracted with `pp.x plot_num=2`, then FFTs it to G-space and reads off shells (`tests/kb_projector_validation.rs:976-978`). The Cube → FFT → G-shell pipeline likely introduces a convention mismatch (sign, normalization, or grid alignment).

## Proposed work

Pick one of:

1. **Re-extract reference via QE's `vlocal_mod` dump.** Modify the qe-runner skill (or write a one-off Fortran patch to QE) to dump `vlocal_of_G` from `qe-7.5/PW/src/init_run.f90` directly. This is the same array our code targets, in QE's own convention, no Cube round-trip.

2. **Delete the test** since VGCMP Phase 1 already covers `V_local(G)` validation against a known-good Python reference. Removes a maintenance burden.

3. **Diagnose the Cube → G-shell mismatch.** Cheaper if the Cube format's sign convention or normalization is the only issue (then patch the test, not QE). But it's a one-off — VGCMP Phase 1 is the durable validation.

**Recommendation:** option 2 (delete) plus a brief comment in `kb_projector_validation.rs` pointing readers at `tests/vgcmp_vloc_cross_check.rs` and `proposals/completed/VGCMP-vloc-g-cross-check.md`. The Cube round-trip is intrinsically lossy (Cube stores real-space samples, not G-space), so a tight cross-check via Cube is fragile by design.

## Verification

After option 2:

- `cargo test test_vloc_comparison_with_qe` — gone.
- VGCMP cross-check tests still pass (already do).

## References

- `proposals/completed/VGCMP-vloc-g-cross-check.md` — Phase 1 result.
- `tests/vgcmp_vloc_cross_check.rs` — the durable replacement.
- `tests/kb_projector_validation.rs:967-1037` — the un-ignored failing test plus its diagnostic comment block.
