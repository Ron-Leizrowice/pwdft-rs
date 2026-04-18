#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy>=1.26",
#     "scipy>=1.11",
# ]
# ///
"""
NLCC audit — Independent Python reference for ρ_core(G).

Purpose
-------
Compute the NLCC core-density Fourier coefficient
    ρ_core(G) = (4π/Ω) ∫₀^∞ ρ_core(r) · j₀(|G|r) · r² dr
for the first few |G| shells of Si FCC, Fe BCC, Cu FCC, and Mn BCC,
using the bare ρ_core(r) stored in `PP_NLCC` of the UPF file.  This
matches QE's `upflib/rhoc_mod.f90:107-115` (`init_tab_rhc`):

    aux(ir)     = upf%rho_atc(ir) * rgrid%r2(ir) * sin(qr)/(qr)
    tab_rhc(iq) = fpi * simpson(aux, rab) / omega

The CSV output is consumed by the Rust unit tests in
`src/pseudopotential/upf/convert.rs` that pin ρ_core(G=0) and
ρ_core(|G|>0) for Si, Fe, Cu, and Mn as regression guards against a
future NLCC regression.  Si and Fe cover the original NCFX scope
(PR #40); Cu and Mn extend coverage per TRV2 Finding #3: Cu exercises
the 3s/3p/3d semicore edge case (Z_val=19), Mn the magnetic reference
(Z_val=15).

Convention (QE native units)
----------------------------
- PP_NLCC stores bare ρ_core(r) in e/Bohr³ (NOT 4πr²·ρ — that's PP_RHOATOM).
- `r`, `rab` in Bohr; ρ in e/Bohr³; G in Bohr⁻¹; Ω in Bohr³.
- Result ρ_core(G) has units of e/Bohr³ (volumetric density in G-space).
- Multiply by 1/BOHR_TO_ANG³ ≈ 6.748 to convert to e/Å³ (pwdft-rs internal).

Si FCC geometry
---------------
- lattice constant a = 5.431 Å = 10.2638... Bohr
- primitive volume Ω = a³/4 Bohr³
- reciprocal lattice: BCC with |G|² in (2π/a)² = 3, 4, 8, 11, 12, 16, ...

Fe BCC geometry
---------------
- lattice constant a = 2.87 Å = 5.4237... Bohr
- primitive volume Ω = a³/2 Bohr³
- reciprocal lattice: FCC with |G|² in (2π/a)² = 2, 4, 6, 8, ...

Output CSV columns
------------------
    system, shell_index, g_bohr_inv, g2_units_of_tpba2,
    rho_core_g_e_per_bohr3, rho_core_g_e_per_ang3
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BOHR_TO_ANG = 0.529_177_210_903
BOHR3_TO_ANG3 = BOHR_TO_ANG**3

# ---------------------------------------------------------------------------
# Manual UPF v2 parsing (subset sufficient for PP_NLCC)
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


def parse_upf_nlcc(path: Path) -> dict:
    """Parse the minimum subset of UPF v2 we need for ρ_core(G)."""
    text = path.read_text()
    m_hdr = re.search(r"<PP_HEADER\b([^>]*)/?>", text)
    if m_hdr is None:
        raise ValueError("PP_HEADER not found")
    hdr = m_hdr.group(1)
    core_flag = _extract_attr(hdr, "core_correction").lower()
    if core_flag != "t":
        raise ValueError(f"PP has no NLCC (core_correction={core_flag!r})")
    mesh = int(_extract_attr(hdr, "mesh_size"))

    r = _extract_block(text, "PP_R")
    rab = _extract_block(text, "PP_RAB")
    rho_core = _extract_block(text, "PP_NLCC")

    for name, arr in (("PP_R", r), ("PP_RAB", rab), ("PP_NLCC", rho_core)):
        if arr.size != mesh:
            raise ValueError(f"{name}: expected {mesh} values, got {arr.size}")

    return {
        "mesh_size": mesh,
        # Keep native UPF units: r in Bohr, ρ in e/Bohr³.
        "r_bohr": r,
        "rab_bohr": rab,
        "rho_core_e_per_bohr3": rho_core,
    }


# ---------------------------------------------------------------------------
# ρ_core(G) in QE native units (e/Bohr³)
# ---------------------------------------------------------------------------


def rho_core_of_g_qe_units(
    r_bohr: np.ndarray,
    rho_core: np.ndarray,
    g_bohr_inv: float,
    omega_bohr3: float,
) -> float:
    """
    Compute ρ_core(G) [e/Bohr³] via the spherical Bessel transform
        ρ_core(G) = (4π/Ω) ∫ ρ_core(r) · j₀(|G|r) · r² dr.

    j₀(x) = sin(x)/x with the x→0 limit j₀(0)=1.
    """
    four_pi = 4.0 * math.pi
    gr = g_bohr_inv * r_bohr
    if g_bohr_inv < 1e-12:
        j0 = np.ones_like(r_bohr)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            j0 = np.where(gr > 1e-10, np.sin(gr) / gr, 1.0 - (gr * gr) / 6.0)
    integrand = rho_core * r_bohr**2 * j0
    integral = simpson(integrand, x=r_bohr)
    return four_pi / omega_bohr3 * integral


# ---------------------------------------------------------------------------
# G shells
# ---------------------------------------------------------------------------


def fcc_g_shells(n_shells: int, a_bohr: float) -> list[tuple[int, float]]:
    """First n_shells non-zero |G| shells for FCC real-space lattice (BCC rec)."""
    tpba = 2.0 * math.pi / a_bohr
    seen: dict[int, None] = {}
    n_max = 6
    for n1 in range(-n_max, n_max + 1):
        for n2 in range(-n_max, n_max + 1):
            for n3 in range(-n_max, n_max + 1):
                x = -n1 + n2 + n3
                y = n1 - n2 + n3
                z = n1 + n2 - n3
                g2_int = x * x + y * y + z * z
                if g2_int == 0:
                    continue
                seen.setdefault(g2_int, None)
    shells_int = sorted(seen.keys())[:n_shells]
    return [(g2_int, tpba * math.sqrt(g2_int)) for g2_int in shells_int]


def bcc_g_shells(n_shells: int, a_bohr: float) -> list[tuple[int, float]]:
    """First n_shells non-zero |G| shells for BCC real-space lattice (FCC rec)."""
    tpba = 2.0 * math.pi / a_bohr
    seen: dict[int, None] = {}
    n_max = 6
    for n1 in range(-n_max, n_max + 1):
        for n2 in range(-n_max, n_max + 1):
            for n3 in range(-n_max, n_max + 1):
                # Reciprocal of BCC is FCC with b_i = (2π/a) * (0,1,1),(1,0,1),(1,1,0)
                x = n2 + n3
                y = n1 + n3
                z = n1 + n2
                g2_int = x * x + y * y + z * z
                if g2_int == 0:
                    continue
                seen.setdefault(g2_int, None)
    shells_int = sorted(seen.keys())[:n_shells]
    return [(g2_int, tpba * math.sqrt(g2_int)) for g2_int in shells_int]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    out_csv = repo_root / "scripts" / "validate" / "rho_core_g_reference.csv"

    systems = [
        {
            "name": "si",
            "upf": repo_root / "pseudopotentials" / "nc" / "lda" / "Si.upf",
            "a_ang": 5.431,
            "lattice": "fcc",
        },
        {
            "name": "fe",
            "upf": repo_root / "pseudopotentials" / "nc" / "lda" / "Fe.upf",
            "a_ang": 2.87,
            "lattice": "bcc",
        },
        # TRV2 Finding #3 — Cu FCC covers the 3s/3p/3d semicore edge case
        # (Z_val=19).  Lattice constant matches `qe_validation/cu_fcc_scf.in`:
        # celldm(1) = 6.8219 Bohr = 3.6100 Å.
        {
            "name": "cu",
            "upf": repo_root / "pseudopotentials" / "nc" / "lda" / "Cu.upf",
            "a_ang": 6.8219 * BOHR_TO_ANG,  # 3.610017 Å
            "lattice": "fcc",
        },
        # TRV2 Finding #3 — Mn (Z_val=15, magnetic reference).  α-Mn has a
        # complex 58-atom cubic ground state; for NLCC regression only the
        # cell volume matters, so we use a simple BCC container with a =
        # 2.89 Å (close to Fe's a = 2.87 Å — puts Mn's ρ_core(G) in the
        # same |G|-shell range as Fe's for comparable sensitivity).
        {
            "name": "mn",
            "upf": repo_root / "pseudopotentials" / "nc" / "lda" / "Mn.upf",
            "a_ang": 2.89,
            "lattice": "bcc",
        },
    ]

    rows: list[tuple] = []

    for sys_cfg in systems:
        name = sys_cfg["name"]
        pp = parse_upf_nlcc(sys_cfg["upf"])
        a_bohr = sys_cfg["a_ang"] / BOHR_TO_ANG
        if sys_cfg["lattice"] == "fcc":
            omega_bohr3 = a_bohr**3 / 4.0
            shells = [(0, 0.0)] + fcc_g_shells(5, a_bohr)
        elif sys_cfg["lattice"] == "bcc":
            omega_bohr3 = a_bohr**3 / 2.0
            shells = [(0, 0.0)] + bcc_g_shells(5, a_bohr)
        else:
            raise ValueError(f"unknown lattice {sys_cfg['lattice']!r}")

        q_core = simpson(
            4.0 * math.pi * pp["r_bohr"] ** 2 * pp["rho_core_e_per_bohr3"],
            x=pp["r_bohr"],
        )
        print(
            f"=== {name.upper()} (a = {sys_cfg['a_ang']} Å = {a_bohr:.4f} Bohr, "
            f"Ω = {omega_bohr3:.4f} Bohr³) ==="
        )
        print(f"  Q_core = ∫ 4πr²ρ_core(r) dr = {q_core:.6f} e")
        print(f"  {'shell':>5}  {'|G|² (int)':>10}  {'|G| (Bohr⁻¹)':>14}  "
              f"{'ρ_core(G) (e/Bohr³)':>22}  {'ρ_core(G) (e/Å³)':>20}")
        print("-" * 84)

        for idx, (g2_int, g_bohr_inv) in enumerate(shells):
            rho_g_bohr = rho_core_of_g_qe_units(
                pp["r_bohr"], pp["rho_core_e_per_bohr3"], g_bohr_inv, omega_bohr3
            )
            rho_g_ang = rho_g_bohr / BOHR3_TO_ANG3
            rows.append(
                (name, idx, g_bohr_inv, g2_int, rho_g_bohr, rho_g_ang)
            )
            print(
                f"  {idx:>5d}  {g2_int:>10d}  {g_bohr_inv:>14.6f}  "
                f"{rho_g_bohr:>22.10e}  {rho_g_ang:>20.10e}"
            )
        print()

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "system",
            "shell_index",
            "g_bohr_inv",
            "g2_units_of_tpba2",
            "rho_core_g_e_per_bohr3",
            "rho_core_g_e_per_ang3",
        ])
        for row in rows:
            w.writerow([
                row[0], row[1], f"{row[2]:.12e}", row[3],
                f"{row[4]:.12e}", f"{row[5]:.12e}",
            ])
    print(f"Wrote {len(rows)} rows to {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
