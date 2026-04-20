"""Assembled Hamiltonian diagonal reference generator (VGCMP Phase 4).

Computes H[G,G] at k=Γ for Si FCC (kinetic + V_NL; V_eff excluded by design).
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

from pwdft_validation.reference.beta import f_l_of_q_bohr
from pwdft_validation.units import BOHR_TO_ANG, RY_TO_EV
from pwdft_validation.upf import UpfData, parse_upf


def _v_nl_diag_ry(pp: UpfData, q_bohr_inv: float, n_atoms: int, omega_bohr3: float) -> float:
    """V_NL(G,G) [Ry] at |k+G|=q; cosθ=1 so angular factor = (2l+1)/(4π)."""
    assert pp.dij_ry is not None
    four_pi = 4.0 * math.pi
    f_vals = np.array([f_l_of_q_bohr(pp, l, chi, q_bohr_inv) for (l, chi) in pp.projectors])
    l_arr = np.array([l for (l, _) in pp.projectors], dtype=np.int64)
    total = 0.0
    for i in range(pp.n_proj):
        for j in range(pp.n_proj):
            if l_arr[i] != l_arr[j]:
                continue
            angular = (2 * int(l_arr[i]) + 1) / four_pi
            total += f_vals[i] * pp.dij_ry[i, j] * f_vals[j] * angular
    return (n_atoms / omega_bohr3) * total


def _si_fcc_shells(n_shells: int, a_bohr: float) -> list[tuple[int, tuple[int, int, int], float]]:
    """Shell list for Si FCC including G=0: [(g2_int, (n1,n2,n3), |G| Bohr⁻¹)]."""
    tpba = 2.0 * math.pi / a_bohr
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
                g2 = x * x + y * y + z * z
                if g2 not in reps:
                    reps[g2] = (n1, n2, n3)
    sorted_g2 = sorted(reps.keys())[: n_shells - 1]
    tail: list[tuple[int, tuple[int, int, int], float]] = [(g2, reps[g2], tpba * math.sqrt(g2)) for g2 in sorted_g2]
    return [(0, (0, 0, 0), 0.0), *tail]


def generate(pseudo_dir: Path, out_csv: Path) -> int:
    """Generate ``vgcmp_phase4_reference.csv`` (VGCMP Phase 4)."""
    upf_path = pseudo_dir / "Si.upf"
    if not upf_path.exists():
        print(f"ERROR: Si UPF not found at {upf_path}", file=sys.stderr)
        return 1

    pp = parse_upf(upf_path)
    if pp.dij_ry is None:
        print("ERROR: Si UPF has no PP_DIJ block", file=sys.stderr)
        return 1

    a_ang = 5.431
    a_bohr = a_ang / BOHR_TO_ANG
    omega_bohr3 = a_bohr**3 / 4.0
    n_atoms = 2
    print(f"Si FCC: a={a_ang} Å = {a_bohr:.6f} Bohr;  Ω={omega_bohr3:.4f} Bohr³")

    shells = _si_fcc_shells(5, a_bohr)
    rows: list[tuple] = []
    header = f"{'shell':>5}  {'n1,n2,n3':>10}  {'|G|² (int)':>10}  {'|G| (Bohr⁻¹)':>14}  {'T (Ry)':>14}  {'V_NL (Ry)':>14}  {'H_diag (Ry)':>14}"
    print(header)
    print("-" * len(header))
    for idx, (g2_int, miller, g_bohr) in enumerate(shells):
        kinetic_ry = g_bohr * g_bohr
        v_nl_ry = _v_nl_diag_ry(pp, g_bohr, n_atoms, omega_bohr3)
        h_diag_ry = kinetic_ry + v_nl_ry
        rows.append((idx, miller[0], miller[1], miller[2], g_bohr, g2_int, g_bohr, kinetic_ry, v_nl_ry, h_diag_ry))
        print(
            f"{idx:>5d}  ({miller[0]:>2d},{miller[1]:>2d},{miller[2]:>2d})  {g2_int:>10d}  "
            f"{g_bohr:>14.6f}  {kinetic_ry:>14.6e}  {v_nl_ry:>14.6e}  {h_diag_ry:>14.6e}"
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "shell_index",
                "miller_n1",
                "miller_n2",
                "miller_n3",
                "g_bohr_inv",
                "g2_int_tpba2",
                "q_bohr_inv",
                "kinetic_ry",
                "v_nl_diag_ry",
                "h_diag_ry",
            ]
        )
        for row in rows:
            w.writerow(
                [
                    row[0],
                    row[1],
                    row[2],
                    row[3],
                    f"{row[4]:.12e}",
                    row[5],
                    f"{row[6]:.12e}",
                    f"{row[7]:.12e}",
                    f"{row[8]:.12e}",
                    f"{row[9]:.12e}",
                ]
            )

    print("\nH_diag (kinetic + V_NL) in eV:")
    for row in rows:
        print(f"  shell={row[0]}  H = {row[9] * RY_TO_EV:+.6f} eV")
    print(f"\nWrote {len(rows)} shells to {out_csv}")
    return 0
