# Researcher Logbook

Entries: date, what was validated, discrepancies found (with numbers), references used. Physics findings only — not code quality or docs.

## 2026-04-16 — Handoff and orientation

### QE discrepancy status

| System | QE (eV) | Ours (eV) | ΔE (eV) |
|--------|---------|-----------|---------|
| Si (2 atoms, 15 Ry, 4x4x4) | -231.61 | -218.28 | 13.3 |
| Fe BCC (1 atom, 16 Ry, 4x4x4) | -3059.46 | -3104.83 | 45.4 |
| C diamond | — | — | Does not converge |

Eigenvalue degeneracy breaking at Gamma confirms the bug is in V_local/V_NL form factors, not energy accounting. Error scales with Z (1.66 eV/el for Si, 2.84 eV/el for Fe).

**Root cause:** Radial quadrature — O(h²) trapezoidal vs QE's O(h⁴) Simpson. Compounded by V_local 1/r singularity (QE uses bounded erf/r subtraction). Proposals SIMP + VERF.

### Formula audit: all 23 items verified correct

No formula-level bugs found. Full comparison against QE 7.5 source. The only issue is numerical (quadrature quality).

### Open physics questions (low priority)

1. IBZ reduction gives 10 k-points for Si 4×4×4 vs QE's 8 — BZ boundary tolerance issue
2. Spin exchange formula at `xc.rs:249` uses equivalent but non-standard weighted-average form — needs documenting
3. `total_energy()` docstring omits V_local(G=0)·N_el correction — formula is correct, doc is incomplete

### Key references

PZ: PRB 23, 5048 (1981). KB: PRL 48, 1425 (1982). NLCC: PRB 26, 1738 (1982). QE source: `vloc_mod.f90`, `simpsn.f90`, `setlocal.f90`.
