#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy>=1.26",
#     "scipy>=1.11",
# ]
# ///
"""
VGCMP Phase 4 — Independent Python reference for assembled Hamiltonian
diagonal H[G,G] at k=Γ for Si FCC (2 atoms).

Purpose
-------
Phases 1-3 closed the per-term form-factor checks (V_local(G), β_l(q), D_ij
all match QE to machine precision). This Phase 4 script assembles the
*diagonal* Kohn-Sham Hamiltonian matrix element at k=Γ for a handful of
G-vectors and provides a reference value that the Rust test
(`tests/vgcmp_assembled_h_cross_check.rs`) must reproduce using
`build_hamiltonian_with_v_eff` + `NonlocalPotential::add_to_hamiltonian`.

We decompose the diagonal into two terms that have independent assembly
logic and cross-check each separately:

1. **Kinetic** (trivial):
        T(G) = (ħ²/2m) · |G|²     (eV·Å² times Å⁻² = eV)
   In QE-native Hartree units,
        T(G) = 0.5 · |G|²_Bohr⁻² Ha = |G|²_Bohr⁻² Ry
   The conversion from Rust's HBAR2_OVER_2M (SI-derived ~3.80998 eV·Å²)
   to HA_TO_EV × 0.5 × BOHR_TO_ANG² (CODATA-consistent) is the first
   check: they should agree to ~1e-7 relative.

2. **Non-local (KB separable form)**, diagonal contribution:
        V_NL(G, G) = (1/Ω) · Σ_atoms S_atom(G-G'=0) · [sum over projectors]
                   = (N_atoms_of_type / Ω) · Σ_{i,j: l_i=l_j} F_i(|G|) · D_ij · F_j(|G|)
                                              × (2l+1)/(4π) · P_l(cosθ=1)

   where cosθ = (G·G')/(|G||G'|) = 1 on the diagonal. Units: Ry.

We deliberately **exclude V_local and V_eff(G=0)** from the Phase 4
reference for two reasons:

(a) In pwdft-rs (`src/scf/context.rs:93-94`), the code sets
    `v_local_fft[0] = 0` — the V_local(G=0) term is absorbed into the
    Ewald/background. This matches QE's convention; see
    `qe-7.5/PW/src/setlocal.f90`.
(b) V_H(G=0) = 0 by the Coulomb-G²-divergence fix, so the *only*
    non-zero V_eff(G=0) contribution is V_xc(G=0) — a constant scalar
    that shifts the entire diagonal uniformly and cannot contribute to
    degeneracy-breaking or state-by-state errors.

Phase 4 therefore isolates the G-dependent part of the diagonal,
which is where any assembly-stage bug (structure factor, (2l+1)/(4π)
angular factor, 1/Ω prefactor, D_ij summation pattern) would show up.
To exercise the V_eff pathway without recomputing the SCF pipeline in
Python, the Rust test loads V_eff_fft = 0 — strictly isolating the
kinetic + V_NL diagonal assembly.

Output CSV columns
------------------
    shell_index, miller_n1, miller_n2, miller_n3,
    g_bohr_inv, g2_int_tpba2, q_bohr_inv,
    kinetic_ry, v_nl_diag_ry, h_diag_ry

(Several G-vectors may belong to the same shell; we pick one
representative per shell with Miller indices (n1, n2, n3) given by
the lexicographic first vector of each shell.)

Si FCC geometry (matches Phase 1)
---------------------------------
- lattice constant a = 5.431 Å = 10.2638... Bohr
- primitive volume Ω = a³/4 Bohr³
- reciprocal lattice vectors (BCC):
    b1 = (2π/a)·(-1, 1, 1), b2 = (2π/a)·( 1,-1, 1), b3 = (2π/a)·( 1, 1,-1)
- k = Γ = (0, 0, 0)
- 2 atoms/cell at fractional (0, 0, 0) and (1/4, 1/4, 1/4)
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.special import spherical_jn

# ---------------------------------------------------------------------------
# Constants (QE native units: Ry, Bohr)
# ---------------------------------------------------------------------------

BOHR_TO_ANG = 0.529_177_210_903
RY_TO_EV = 13.605_693_122_994

# ---------------------------------------------------------------------------
# Manual UPF v2 parsing (lifted from Phase 1/2 scripts, no external deps)
# ---------------------------------------------------------------------------


def _extract_attr(text: str, name: str) -> str:
    m = re.search(rf'{re.escape(name)}\s*=\s*"([^"]*)"', text)
    if m is None:
        raise ValueError(f"attribute {name!r} not found")
    return m.group(1).strip()


def _extract_block(text: str, tag: str) -> np.ndarray:
    m = re.search(rf"<{re.escape(tag)}\b[^>]*>(.*?)</{re.escape(tag)}>", text, re.DOTALL)
    if m is None:
        raise ValueError(f"tag <{tag}> not found")
    body = m.group(1)
    values = [float(tok) for tok in body.split() if tok.strip()]
    return np.asarray(values, dtype=np.float64)


def _extract_beta_projector(text: str, idx: int, mesh: int) -> tuple[int, np.ndarray]:
    tag = f"PP_BETA.{idx}"
    m = re.search(
        rf"<{re.escape(tag)}\b([^>]*)>(.*?)</{re.escape(tag)}>",
        text,
        re.DOTALL,
    )
    if m is None:
        raise ValueError(f"tag <{tag}> not found")
    attrs = m.group(1)
    body = m.group(2)
    l = int(_extract_attr(attrs, "angular_momentum"))
    values = np.asarray(
        [float(tok) for tok in body.split() if tok.strip()],
        dtype=np.float64,
    )
    if values.size != mesh:
        raise ValueError(f"{tag}: expected {mesh} values, got {values.size}")
    return l, values


def parse_upf_si(path: Path) -> dict:
    text = path.read_text()
    m_hdr = re.search(r"<PP_HEADER\b([^>]*)/?>", text)
    if m_hdr is None:
        raise ValueError("PP_HEADER not found")
    hdr = m_hdr.group(1)
    z_val = float(_extract_attr(hdr, "z_valence"))
    mesh = int(_extract_attr(hdr, "mesh_size"))
    n_proj = int(_extract_attr(hdr, "number_of_proj"))

    r = _extract_block(text, "PP_R")
    rab = _extract_block(text, "PP_RAB")

    for name, arr in (("PP_R", r), ("PP_RAB", rab)):
        if arr.size != mesh:
            raise ValueError(f"{name}: expected {mesh} values, got {arr.size}")

    projectors: list[tuple[int, np.ndarray]] = []
    for i in range(1, n_proj + 1):
        l, chi = _extract_beta_projector(text, i, mesh)
        projectors.append((l, chi))

    # D_ij: PP_DIJ is a flat row-major matrix of size n_proj^2, values in Ry.
    dij_flat = _extract_block(text, "PP_DIJ")
    if dij_flat.size != n_proj * n_proj:
        raise ValueError(
            f"PP_DIJ: expected {n_proj * n_proj} values, got {dij_flat.size}"
        )
    dij = dij_flat.reshape(n_proj, n_proj)

    return {
        "z_valence": z_val,
        "mesh_size": mesh,
        "n_proj": n_proj,
        "r_bohr": r,
        "rab_bohr": rab,
        "projectors": projectors,
        "dij_ry": dij,
    }


# ---------------------------------------------------------------------------
# F_l(q) in QE native units (Bohr), same as Phase 2 reference.
# ---------------------------------------------------------------------------


def f_l_of_q(
    r_bohr: np.ndarray,
    chi: np.ndarray,
    l: int,
    q_bohr_inv: float,
) -> float:
    """F_l(q) = 4π · ∫₀^∞ χ(r) · j_l(q·r) · r dr   [Bohr^(3/2)]."""
    qr = q_bohr_inv * r_bohr
    jl = spherical_jn(l, qr)
    integrand = chi * jl * r_bohr
    integral = simpson(integrand, x=r_bohr)
    return 4.0 * math.pi * integral


# ---------------------------------------------------------------------------
# Non-local diagonal contribution V_NL(G, G) per atom type.
#
# For a single type with N_atoms atoms, block-diagonal D_ij in l:
#
#   V_NL(G, G) = (N_atoms / Ω) · Σ_{i,j: l_i=l_j}
#                    F_i(|G|) · D_{ij} · F_j(|G|) · (2 l_i + 1) / (4π)
#
# P_l(1) = 1 so the angular factor reduces to (2l+1)/(4π). Units: Ry
# (because D_ij is in Ry; F in Bohr^(3/2); Ω in Bohr³; F·F·/Ω is Bohr^(3/2)
# × Bohr^(3/2) / Bohr³ = 1 (dimensionless) so V_NL = Ry).
# ---------------------------------------------------------------------------


def v_nl_diag_ry(
    pp: dict,
    q_bohr_inv: float,
    n_atoms: int,
    omega_bohr3: float,
) -> float:
    """Diagonal non-local matrix element V_NL(G, G) in Ry at |k+G| = q.

    Since we evaluate at the diagonal (G' = G), cosθ = 1 exactly and the
    angular factor collapses to (2l+1)/(4π).
    """
    projectors = pp["projectors"]
    dij = pp["dij_ry"]
    four_pi = 4.0 * math.pi

    # Precompute F_i(q) for all projectors.
    f_vals = np.array(
        [f_l_of_q(pp["r_bohr"], chi, l, q_bohr_inv) for (l, chi) in projectors]
    )
    l_arr = np.array([l for (l, _) in projectors], dtype=np.int64)
    n_proj = len(projectors)

    total = 0.0
    for i in range(n_proj):
        for j in range(n_proj):
            if l_arr[i] != l_arr[j]:
                continue
            li = l_arr[i]
            angular = (2 * li + 1) / four_pi  # P_l(1) = 1
            total += f_vals[i] * dij[i, j] * f_vals[j] * angular

    return (n_atoms / omega_bohr3) * total


# ---------------------------------------------------------------------------
# Si FCC G-vector enumeration at k = Γ.
#
# We pick one representative G per shell, choosing the Miller triple that
# appears first in the enumeration. The reciprocal lattice vectors (in
# Bohr⁻¹) are b_i = (2π/a_Bohr) · B_i where B_i are integer vectors:
#     B_1 = (-1, 1, 1), B_2 = (1, -1, 1), B_3 = (1, 1, -1)
# ---------------------------------------------------------------------------


def si_fcc_shells(n_shells: int, a_bohr: float) -> list[tuple[int, tuple[int, int, int], float]]:
    """Return [(g2_int, (n1, n2, n3), |G| in Bohr⁻¹)] for the first n_shells shells.

    g2_int is |G|² in units of (2π/a)².
    The Miller triple (n1, n2, n3) is one representative G with integer
    coefficients on the BCC reciprocal basis (b_1, b_2, b_3). The
    Cartesian G = n1·b1 + n2·b2 + n3·b3.
    """
    tpba = 2.0 * math.pi / a_bohr
    # Enumerate (n1, n2, n3) with small |n_i| and accumulate one rep per shell.
    reps: dict[int, tuple[int, int, int]] = {}
    n_max = 4
    for n1 in range(-n_max, n_max + 1):
        for n2 in range(-n_max, n_max + 1):
            for n3 in range(-n_max, n_max + 1):
                if (n1, n2, n3) == (0, 0, 0):
                    continue
                x = -n1 + n2 + n3
                y = n1 - n2 + n3
                z = n1 + n2 - n3
                g2_int = x * x + y * y + z * z
                if g2_int not in reps:
                    reps[g2_int] = (n1, n2, n3)

    # Sort by g2_int ascending, include |G|=0 at the head with Miller (0,0,0).
    sorted_g2 = sorted(reps.keys())[: n_shells - 1]
    out: list[tuple[int, tuple[int, int, int], float]] = [(0, (0, 0, 0), 0.0)]
    for g2 in sorted_g2:
        out.append((g2, reps[g2], tpba * math.sqrt(g2)))
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    upf_path = repo_root / "pseudopotentials" / "nc" / "lda" / "Si.upf"
    out_csv = repo_root / "scripts" / "validate" / "vgcmp_phase4_reference.csv"

    if not upf_path.exists():
        print(f"ERROR: Si UPF not found at {upf_path}", file=sys.stderr)
        return 1

    pp = parse_upf_si(upf_path)
    print(
        f"Loaded Si UPF: Z_val={pp['z_valence']}, mesh={pp['mesh_size']}, "
        f"n_proj={pp['n_proj']}"
    )
    l_list = [l for (l, _) in pp["projectors"]]
    print(f"projector l's: {l_list}")
    print(f"D_ij (Ry) diagonal: {[float(pp['dij_ry'][i, i]) for i in range(pp['n_proj'])]}")

    # Si FCC geometry (matches Phase 1 choices; a = 5.431 Å).
    a_ang = 5.431
    a_bohr = a_ang / BOHR_TO_ANG
    omega_bohr3 = a_bohr**3 / 4.0
    tpba = 2.0 * math.pi / a_bohr
    n_atoms = 2
    print(
        f"Si FCC: a = {a_ang} Å = {a_bohr:.6f} Bohr;  Ω = {omega_bohr3:.4f} Bohr³;"
        f"  2π/a = {tpba:.6f} Bohr⁻¹"
    )

    # First 5 shells: (|G|=0, shell 1 = (2π/a)·√3, shell 2 = ·√8, etc.)
    shells = si_fcc_shells(n_shells=5, a_bohr=a_bohr)

    rows: list[tuple[int, int, int, int, float, int, float, float, float, float]] = []
    print()
    header = (
        f"{'shell':>5}  {'n1,n2,n3':>10}  {'|G|² (int)':>10}  "
        f"{'|G| (Bohr⁻¹)':>14}  {'T (Ry)':>14}  {'V_NL (Ry)':>14}  {'H_diag (Ry)':>14}"
    )
    print(header)
    print("-" * len(header))
    for idx, (g2_int, miller, g_bohr) in enumerate(shells):
        # Kinetic in QE native units: T = 0.5 * |G|² in Ha = |G|²_Bohr⁻² Ry
        # because 1 Ha = 2 Ry, so 0.5 * |G|² Ha = |G|² Ry in Bohr⁻¹ input.
        kinetic_ry = g_bohr * g_bohr  # Ry = Ha × 2 × 0.5 = Ha. In Bohr units, T = |G|² Ry.

        # Non-local diagonal.
        v_nl_ry = v_nl_diag_ry(pp, g_bohr, n_atoms, omega_bohr3)

        # Phase 4 "H_diag" definition: kinetic + V_NL only.
        # V_eff(G=0) is handled separately (see module docstring).
        h_diag_ry = kinetic_ry + v_nl_ry

        rows.append(
            (
                idx,
                miller[0], miller[1], miller[2],
                g_bohr,
                g2_int,
                g_bohr,  # q = |k+G| at k=Γ equals |G|
                kinetic_ry,
                v_nl_ry,
                h_diag_ry,
            )
        )
        print(
            f"{idx:>5d}  ({miller[0]:>2d},{miller[1]:>2d},{miller[2]:>2d})  {g2_int:>10d}  "
            f"{g_bohr:>14.6f}  {kinetic_ry:>14.6e}  {v_nl_ry:>14.6e}  {h_diag_ry:>14.6e}"
        )

    # Write CSV with per-term and summed reference values.
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "shell_index", "miller_n1", "miller_n2", "miller_n3",
                "g_bohr_inv", "g2_int_tpba2", "q_bohr_inv",
                "kinetic_ry", "v_nl_diag_ry", "h_diag_ry",
            ]
        )
        for row in rows:
            w.writerow(
                [
                    row[0], row[1], row[2], row[3],
                    f"{row[4]:.12e}", row[5], f"{row[6]:.12e}",
                    f"{row[7]:.12e}", f"{row[8]:.12e}", f"{row[9]:.12e}",
                ]
            )

    # Pretty-printed H_diag in eV for log.
    print("\nH_diag (kinetic + V_NL) in eV:")
    for row in rows:
        print(f"  shell={row[0]}  H = {row[9] * RY_TO_EV:+.6f} eV")

    print(f"\nWrote {len(rows)} shells to {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
