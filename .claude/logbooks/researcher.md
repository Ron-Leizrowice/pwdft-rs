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

## 2026-04-16 — KBTF: KB projector test failure investigation

Investigated 3 failing tests in `tests/kb_projector_validation.rs`. PR #1 on branch `KBTF/kb-test-failures`.

| Test | Classification | Root cause |
|------|---------------|------------|
| `test_09` (D_ij vs HGH h^l_ij) | **Test bug** | Test assumed D_ij = raw HGH h^l_ij (3×3). Si.upf has 6 projectors (2 per l=0,1,2); QE diagonalizes h^l and absorbs eigenvector rotation into projectors. UPF D_ij is diagonal 6×6, not raw h^l_ij. |
| `test_07` (form factor decay) | **Known limitation** | l=1 projector |F(24.5)|/|F_max| = 0.106 > 0.1 threshold. Trapezoidal quadrature artifact at high q. SIMP would fix. |
| `test_vloc` (V_local vs QE) | **Known limitation** | `v_local_of_g` uses bare Coulomb subtraction (V+Z/r); QE uses erf(r)/r subtraction (numerically superior, avoids cancellation at large r). VERF would fix. |

**Key finding on HGH→UPF mapping:** QE's UPF conversion diagonalizes each l-block of h^l_ij, stores eigenvalues as diagonal D_ij, and rotates projectors accordingly. This means UPF D_ij values (e.g., D[0,0]=11.13 Ry) cannot be compared against published HGH h^l_ij (e.g., h^0_11=2.95 Ry). Tests 05/06/08/10 passing confirms D_ij + projectors are self-consistent.

**Confirms root cause from orientation:** Both test_07 and test_vloc failures trace to trapezoidal quadrature + bare Coulomb subtraction — same root cause as the 13.3 eV Si discrepancy. SIMP + VERF remain the correct fix path.

### Key references

PZ: PRB 23, 5048 (1981). KB: PRL 48, 1425 (1982). NLCC: PRB 26, 1738 (1982). QE source: `vloc_mod.f90`, `simpsn.f90`, `setlocal.f90`.
