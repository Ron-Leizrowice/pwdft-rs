#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy>=1.26",
#     "scipy>=1.11",
# ]
# ///
"""
VGCH diagnostic — Reference V_local(G=0) for every heavy-atom PP,
independently in Python, using the QE bare-Coulomb form

    V(0) = (4π/Ω) ∫ r² [V_loc(r) + Z·e²/r] dr  (Ry units)

For each element (Si/Fe/Cu/Ga/As/Na/Cl/Mg/O/Al/C), report V_local(G=0)
in Ry, eV, and — for the 1-atom BCC/FCC cells — the expected
E_local_g0_shift = N_el · V_local_cell(G=0).

This is fed into the pwdft-rs side's diagnostic: any per-species
V_local(G=0) mismatch vs this reference would indicate a bug in
`v_local_of_g` at G=0 (the bare-Coulomb integrand with 1/r) that
scales with Z.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.special import erf  # noqa: F401 (kept for parity w/ Si script)

BOHR_TO_ANG = 0.529_177_210_903
RY_TO_EV = 13.605_693_122_994
E2_RY_BOHR = 2.0


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
    return np.asarray([float(t) for t in body.split() if t.strip()], dtype=np.float64)


def parse_upf(path: Path) -> dict:
    text = path.read_text()
    m_hdr = re.search(r"<PP_HEADER\b([^>]*)/?>", text)
    if m_hdr is None:
        raise ValueError("PP_HEADER not found")
    hdr = m_hdr.group(1)
    z_val = float(_extract_attr(hdr, "z_valence"))
    mesh = int(_extract_attr(hdr, "mesh_size"))
    r = _extract_block(text, "PP_R")
    rab = _extract_block(text, "PP_RAB")
    v_loc = _extract_block(text, "PP_LOCAL")
    for name, arr in (("PP_R", r), ("PP_RAB", rab), ("PP_LOCAL", v_loc)):
        if arr.size != mesh:
            raise ValueError(f"{name}: expected {mesh} values, got {arr.size}")
    return {
        "z_valence": z_val,
        "mesh_size": mesh,
        "r_bohr": r,
        "rab_bohr": rab,
        "v_local_ry": v_loc,
    }


def v_local_g0_qe_units(r_bohr: np.ndarray, v_ry: np.ndarray, z_val: float, omega_bohr3: float) -> float:
    """V_local(G=0) in Ry via (4π/Ω) ∫ r²[V+Ze²/r] dr."""
    four_pi = 4.0 * math.pi
    with np.errstate(divide="ignore", invalid="ignore"):
        coulomb = np.where(r_bohr > 0.0, z_val * E2_RY_BOHR / r_bohr, 0.0)
    integrand = r_bohr**2 * (v_ry + coulomb)
    integral = simpson(integrand, x=r_bohr)
    return four_pi / omega_bohr3 * integral


# ---------------------------------------------------------------------------
# Cells (match pwdft/pwdft-core/tests/qe_validation.rs).
# ---------------------------------------------------------------------------
CELLS = [
    # (elem, bravais, a_ang, n_atoms_of_THIS_species, n_atoms_primitive_cell)
    ("Si", "fcc", 5.431, 2, 2),
    ("C", "fcc", 3.567, 2, 2),
    ("Al", "fcc", 4.05, 1, 1),
    ("Fe", "bcc", 2.87, 1, 1),
    ("Cu", "fcc", 3.61, 1, 1),
    ("Ga", "fcc", 5.653, 1, 2),
    ("As", "fcc", 5.653, 1, 2),
    ("Na", "fcc", 5.614, 1, 2),
    ("Cl", "fcc", 5.614, 1, 2),
    ("Mg", "fcc", 4.212, 1, 2),
    ("O", "fcc", 4.212, 1, 2),
]


def primitive_volume_bohr3(bravais: str, a_bohr: float) -> float:
    if bravais == "fcc":
        return a_bohr**3 / 4.0
    if bravais == "bcc":
        return a_bohr**3 / 2.0
    raise ValueError(bravais)


def main() -> int:
    from pwdft_validation import CSV_REF_DIR, PSEUDO_DIR

    pp_dir = PSEUDO_DIR / "nc" / "lda"
    out_csv = CSV_REF_DIR / "vgch_vloc_heavy.csv"

    rows = []
    print(f"{'elem':>4}  {'Z_val':>5}  {'Ω (Bohr³)':>12}  {'V(G=0) Ry':>14}  {'V(G=0) eV':>14}  {'N_atoms':>7}")
    print("-" * 70)
    for elem, bravais, a_ang, n_atoms_species, _n_atoms_prim in CELLS:
        upf_path = pp_dir / f"{elem}.upf"
        if not upf_path.exists():
            print(f"  SKIP {elem}: {upf_path} missing", file=sys.stderr)
            continue
        pp = parse_upf(upf_path)
        a_bohr = a_ang / BOHR_TO_ANG
        omega_bohr3 = primitive_volume_bohr3(bravais, a_bohr)
        v_g0_ry = v_local_g0_qe_units(pp["r_bohr"], pp["v_local_ry"], pp["z_valence"], omega_bohr3)
        v_g0_ev = v_g0_ry * RY_TO_EV
        rows.append((elem, pp["z_valence"], omega_bohr3, v_g0_ry, v_g0_ev, n_atoms_species))
        print(
            f"{elem:>4}  {pp['z_valence']:>5.1f}  {omega_bohr3:>12.4f}  "
            f"{v_g0_ry:>14.8f}  {v_g0_ev:>14.8f}  {n_atoms_species:>7}"
        )

    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["element", "z_valence", "omega_bohr3", "v_local_g0_ry", "v_local_g0_ev", "n_atoms_this_species"])
        for row in rows:
            w.writerow([row[0], row[1], f"{row[2]:.6f}", f"{row[3]:.10e}", f"{row[4]:.10e}", row[5]])
    print(f"\nWrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
