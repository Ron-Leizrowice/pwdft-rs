"""ρ_core(G) NLCC core-density Fourier coefficients (VGCH NLCC audit)."""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

from pwdft_validation.units import BOHR3_TO_ANG3, BOHR_TO_ANG
from pwdft_validation.upf import parse_upf


@dataclass(frozen=True, slots=True)
class _NlccSystem:
    name: str
    element: str
    a_ang: float
    lattice: str  # "fcc" or "bcc"


def rho_core_of_g_qe_units(r_bohr: np.ndarray, rho_core: np.ndarray, g_bohr_inv: float, omega_bohr3: float) -> float:
    """ρ_core(G) [e/Bohr³] via (4π/Ω) ∫ ρ_core(r) · j₀(|G|r) · r² dr."""
    gr = g_bohr_inv * r_bohr
    if g_bohr_inv < 1e-12:
        j0 = np.ones_like(r_bohr)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            j0 = np.where(gr > 1e-10, np.sin(gr) / gr, 1.0 - gr * gr / 6.0)
    integral = simpson(rho_core * r_bohr**2 * j0, x=r_bohr)
    return 4.0 * math.pi / omega_bohr3 * integral


def _fcc_shells(n: int, a_bohr: float) -> list[tuple[int, float]]:
    tpba = 2.0 * math.pi / a_bohr
    seen: dict[int, None] = {}
    for n1 in range(-6, 7):
        for n2 in range(-6, 7):
            for n3 in range(-6, 7):
                x, y, z = -n1 + n2 + n3, n1 - n2 + n3, n1 + n2 - n3
                g2 = x * x + y * y + z * z
                if g2 > 0:
                    seen.setdefault(g2, None)
    return [(g2, tpba * math.sqrt(g2)) for g2 in sorted(seen)[:n]]


def _bcc_shells(n: int, a_bohr: float) -> list[tuple[int, float]]:
    tpba = 2.0 * math.pi / a_bohr
    seen: dict[int, None] = {}
    for n1 in range(-6, 7):
        for n2 in range(-6, 7):
            for n3 in range(-6, 7):
                x, y, z = n2 + n3, n1 + n3, n1 + n2
                g2 = x * x + y * y + z * z
                if g2 > 0:
                    seen.setdefault(g2, None)
    return [(g2, tpba * math.sqrt(g2)) for g2 in sorted(seen)[:n]]


_SYSTEMS: list[_NlccSystem] = [
    _NlccSystem("si", "Si", 5.431, "fcc"),
    _NlccSystem("fe", "Fe", 2.87, "bcc"),
    # TRV2 Finding #3 — Cu FCC covers the 3s/3p/3d semicore edge case (Z_val=19).
    # Lattice constant from data/qe/cu_fcc_scf.in: celldm(1) = 6.8219 Bohr.
    _NlccSystem("cu", "Cu", 6.8219 * BOHR_TO_ANG, "fcc"),
    # TRV2 Finding #3 — Mn (Z_val=15, magnetic reference). α-Mn has a complex
    # 58-atom cubic ground state; for NLCC regression only the cell volume
    # matters, so we use a simple BCC container with a = 2.89 Å (close to Fe's
    # a = 2.87 Å — puts Mn's ρ_core(G) in the same |G|-shell range as Fe's).
    _NlccSystem("mn", "Mn", 2.89, "bcc"),
    # VGCH-2F Part C session-2 — H-C4: extend NLCC ρ_core(G) regression
    # coverage to every NLCC-active PP in a Class A heavy-atom QE cell. Cu
    # was pinned by TRV2; Ga/As close GaAs zinc-blende, O closes MgO, Cl
    # closes NaCl. Mg and Na PPs have core_correction=F — no pin needed.
    _NlccSystem("ga", "Ga", 10.6829 * BOHR_TO_ANG, "fcc"),  # GaAs: 5.6530 Å
    _NlccSystem("as", "As", 10.6829 * BOHR_TO_ANG, "fcc"),  # GaAs: 5.6530 Å
    _NlccSystem("o", "O", 7.9586 * BOHR_TO_ANG, "fcc"),  # MgO: 4.2115 Å
    _NlccSystem("cl", "Cl", 10.6078 * BOHR_TO_ANG, "fcc"),  # NaCl: 5.6133 Å
]


def generate(pseudo_dir: Path, out_csv: Path) -> int:
    """Generate ``rho_core_g_reference.csv`` (NLCC audit — Si, Fe, Cu, Mn, Ga, As, O, Cl)."""
    rows: list[tuple] = []

    for cfg in _SYSTEMS:
        upf_path = pseudo_dir / f"{cfg.element}.upf"
        if not upf_path.exists():
            print(f"SKIP {cfg.name}: {upf_path} missing", file=sys.stderr)
            continue
        pp = parse_upf(upf_path)
        if not pp.has_nlcc or pp.rho_core_ebohr3 is None:
            print(f"SKIP {cfg.name}: no NLCC block", file=sys.stderr)
            continue

        a_bohr = cfg.a_ang / BOHR_TO_ANG
        if cfg.lattice == "fcc":
            omega_bohr3 = a_bohr**3 / 4.0
            shells = [(0, 0.0), *_fcc_shells(5, a_bohr)]
        else:
            omega_bohr3 = a_bohr**3 / 2.0
            shells = [(0, 0.0), *_bcc_shells(5, a_bohr)]

        q_core = simpson(4.0 * math.pi * pp.r_bohr**2 * pp.rho_core_ebohr3, x=pp.r_bohr)
        print(f"=== {cfg.name.upper()} (a={cfg.a_ang:.4f} Å, Ω={omega_bohr3:.4f} Bohr³) ===")
        print(f"  Q_core = {q_core:.6f} e")
        for idx, (g2_int, g_bohr_inv) in enumerate(shells):
            rho_g_bohr = rho_core_of_g_qe_units(pp.r_bohr, pp.rho_core_ebohr3, g_bohr_inv, omega_bohr3)
            rho_g_ang = rho_g_bohr / BOHR3_TO_ANG3
            rows.append((cfg.name, idx, g_bohr_inv, g2_int, rho_g_bohr, rho_g_ang))
            print(f"  {idx:>5d}  {g2_int:>10d}  {g_bohr_inv:>14.6f}  {rho_g_bohr:>22.10e}  {rho_g_ang:>20.10e}")
        print()

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "system",
                "shell_index",
                "g_bohr_inv",
                "g2_units_of_tpba2",
                "rho_core_g_e_per_bohr3",
                "rho_core_g_e_per_ang3",
            ]
        )
        for row in rows:
            w.writerow([row[0], row[1], f"{row[2]:.12e}", row[3], f"{row[4]:.12e}", f"{row[5]:.12e}"])
    print(f"Wrote {len(rows)} rows to {out_csv}")
    return 0
