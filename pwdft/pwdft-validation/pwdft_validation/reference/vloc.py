"""V_local(G) reference generators (VGCMP Phase 1 and VGCH Phase 1c).

Covers:
- ``generate_si``:     V_local(G) for the first 20 |G| shells of Si FCC.
- ``generate_heavy``:  V_local(G=0) for 11 VGCH-scope elements.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.special import erf

from pwdft_validation.units import BOHR_TO_ANG, E2_RY_BOHR, RY_TO_EV
from pwdft_validation.upf import UpfData, parse_upf

# ---------------------------------------------------------------------------
# Core math
# ---------------------------------------------------------------------------


def v_local_of_g_qe_units(
    pp: UpfData,
    g_bohr_inv: float,
    omega_bohr3: float,
) -> float:
    """V_local(G) [Ry] using the QE erf-subtracted convention.

    G=0 uses the bare-Coulomb form; G≠0 uses erf subtraction matching
    ``qe-7.5/upflib/vloc_mod.f90`` lines 136-148.
    """
    assert pp.v_local_ry is not None, "UPF has no PP_LOCAL block"
    r = pp.r_bohr
    v = pp.v_local_ry
    z = pp.z_valence
    four_pi = 4.0 * math.pi

    if g_bohr_inv < 1e-12:
        with np.errstate(divide="ignore", invalid="ignore"):
            coulomb = np.where(r > 0.0, z * E2_RY_BOHR / r, 0.0)
        integral = simpson(r**2 * (v + coulomb), x=r)
        return four_pi / omega_bohr3 * integral

    g = g_bohr_inv
    gr = g * r
    with np.errstate(divide="ignore", invalid="ignore"):
        sin_over_g = np.where(gr > 1e-10, np.sin(gr) / g, r * (1.0 - gr * gr / 6.0))
    short = r * v + z * E2_RY_BOHR * erf(r)
    integral = simpson(short * sin_over_g, x=r)
    tail = four_pi * z * E2_RY_BOHR * math.exp(-g * g / 4.0) / (omega_bohr3 * g * g)
    return four_pi / omega_bohr3 * integral - tail


def v_local_g0_qe_units(pp: UpfData, omega_bohr3: float) -> float:
    """V_local(G=0) [Ry] via (4π/Ω) ∫ r²[V+Ze²/r] dr."""
    return v_local_of_g_qe_units(pp, 0.0, omega_bohr3)


def fcc_g_shells(n_shells: int, a_bohr: float) -> list[tuple[int, float, float]]:
    """First *n_shells* non-zero |G| shells for Si FCC (BCC reciprocal lattice).

    Returns list of (|G|² integer in (2π/a)² units, |G| in Bohr⁻¹, |G|² in Bohr⁻²).
    """
    tpba = 2.0 * math.pi / a_bohr
    seen: dict[int, None] = {}
    n_max = 8
    for n1 in range(-n_max, n_max + 1):
        for n2 in range(-n_max, n_max + 1):
            for n3 in range(-n_max, n_max + 1):
                x = -n1 + n2 + n3
                y = n1 - n2 + n3
                z = n1 + n2 - n3
                g2 = x * x + y * y + z * z
                if g2 > 0:
                    seen.setdefault(g2, None)
    shells_int = sorted(seen.keys())[:n_shells]
    return [(g2, tpba * math.sqrt(g2), (tpba**2) * g2) for g2 in shells_int]


def primitive_volume_bohr3(bravais: str, a_bohr: float) -> float:
    if bravais == "fcc":
        return a_bohr**3 / 4.0
    if bravais == "bcc":
        return a_bohr**3 / 2.0
    raise ValueError(f"unknown bravais {bravais!r}")


# ---------------------------------------------------------------------------
# Generator functions
# ---------------------------------------------------------------------------

_HEAVY_CELLS = [
    # (elem, bravais, a_ang, n_atoms_this_species)
    ("Si", "fcc", 5.431, 2),
    ("C", "fcc", 3.567, 2),
    ("Al", "fcc", 4.05, 1),
    ("Fe", "bcc", 2.87, 1),
    ("Cu", "fcc", 3.61, 1),
    ("Ga", "fcc", 5.653, 1),
    ("As", "fcc", 5.653, 1),
    ("Na", "fcc", 5.614, 1),
    ("Cl", "fcc", 5.614, 1),
    ("Mg", "fcc", 4.212, 1),
    ("O", "fcc", 4.212, 1),
]


def generate_si(pseudo_dir: Path, out_csv: Path) -> int:
    """Generate ``vloc_g_si_reference.csv`` (VGCMP Phase 1)."""
    upf_path = pseudo_dir / "Si.upf"
    if not upf_path.exists():
        print(f"ERROR: Si UPF not found at {upf_path}", file=sys.stderr)
        return 1

    pp = parse_upf(upf_path)
    print(f"Loaded Si UPF: Z_val={pp.z_valence}, mesh={pp.mesh_size}")

    a_ang = 5.431
    a_bohr = a_ang / BOHR_TO_ANG
    omega_bohr3 = a_bohr**3 / 4.0
    tpba = 2.0 * math.pi / a_bohr
    print(f"Si FCC: a={a_ang} Å = {a_bohr:.6f} Bohr;  Ω={omega_bohr3:.4f} Bohr³;  2π/a={tpba:.6f} Bohr⁻¹")

    shells = fcc_g_shells(20, a_bohr)
    rows: list[tuple] = []
    print(f"\n{'shell':>5}  {'|G|² (int)':>10}  {'|G| (Bohr⁻¹)':>14}  {'V_loc(G) (Ry)':>16}")
    print("-" * 56)
    for idx, (g2_int, g_bohr_inv, _) in enumerate(shells):
        v_g = v_local_of_g_qe_units(pp, g_bohr_inv, omega_bohr3)
        rows.append((idx, g_bohr_inv, g2_int, v_g))
        print(f"{idx:>5d}  {g2_int:>10d}  {g_bohr_inv:>14.6f}  {v_g:>16.8e}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["shell_index", "g_bohr_inv", "g2_units_of_tpba2", "v_local_g_ry"])
        for row in rows:
            w.writerow([row[0], f"{row[1]:.12e}", row[2], f"{row[3]:.12e}"])
    print(f"\nWrote {len(rows)} shells to {out_csv}")
    return 0


def generate_heavy(pseudo_dir: Path, out_csv: Path) -> int:
    """Generate ``vgch_vloc_heavy.csv`` (VGCH Phase 1c — 11 heavy-atom PPs)."""
    rows = []
    print(f"{'elem':>4}  {'Z_val':>5}  {'Ω (Bohr³)':>12}  {'V(G=0) Ry':>14}  {'V(G=0) eV':>14}  {'N_atoms':>7}")
    print("-" * 70)
    for elem, bravais, a_ang, n_atoms in _HEAVY_CELLS:
        upf_path = pseudo_dir / f"{elem}.upf"
        if not upf_path.exists():
            print(f"  SKIP {elem}: {upf_path} missing", file=sys.stderr)
            continue
        pp = parse_upf(upf_path)
        a_bohr = a_ang / BOHR_TO_ANG
        omega_bohr3 = primitive_volume_bohr3(bravais, a_bohr)
        v_g0_ry = v_local_g0_qe_units(pp, omega_bohr3)
        v_g0_ev = v_g0_ry * RY_TO_EV
        rows.append((elem, pp.z_valence, omega_bohr3, v_g0_ry, v_g0_ev, n_atoms))
        print(
            f"{elem:>4}  {pp.z_valence:>5.1f}  {omega_bohr3:>12.4f}  {v_g0_ry:>14.8f}  {v_g0_ev:>14.8f}  {n_atoms:>7}"
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["element", "z_valence", "omega_bohr3", "v_local_g0_ry", "v_local_g0_ev", "n_atoms_this_species"])
        for row in rows:
            w.writerow([row[0], row[1], f"{row[2]:.6f}", f"{row[3]:.10e}", f"{row[4]:.10e}", row[5]])
    print(f"\nWrote {out_csv}")
    return 0
