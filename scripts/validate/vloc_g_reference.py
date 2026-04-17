#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy>=1.26",
#     "scipy>=1.11",
# ]
# ///
"""
VGCMP Phase 1 — Independent Python reference for V_local(G).

Purpose
-------
Compute V_local(G) for the first 20 distinct |G| shells of Si FCC
(a = 5.431 Å) directly from `pseudopotentials/nc/lda/Si.upf`, using the
QE erf-subtracted form (QE 7.5 `upflib/vloc_mod.f90`, lines 136-148).

This is an *independent* second witness of the integral — we parse the
UPF by hand (no external UPF library) and use `scipy.integrate.simpson`
for the radial integral. The output CSV is consumed by
`tests/vgcmp_vloc_cross_check.rs` which asserts shell-by-shell agreement
with pwdft-rs's `PseudopotentialData::v_local_of_g`.

Convention (QE native units)
----------------------------
- lengths in Bohr, energies in Ry
- e² = 2 (Gaussian Rydberg units: e² Ry·Bohr = 2)
- for G = 0:
      V(0) = (4π/Ω) ∫₀^∞ r² · [V(r) + Z·e²/r] dr        (Ry·Bohr³)
    with the bracketed term short-ranged because V(r) → −Z·e²/r as r→∞.
- for G ≠ 0 (erf subtraction, matches QE):
      V(G) = (4π/Ω) ∫₀^∞ [r·V(r) + Z·e²·erf(r)] · sin(Gr)/G dr
             − 4π·Z·e² · exp(−G²/4) / (Ω·G²)                     (Ry)
    r, G in Bohr / Bohr⁻¹; the erf argument is r in Bohr (treated as
    dimensionless, matching QE's convention — the Gaussian width is
    1 Bohr).

Si FCC geometry
---------------
- lattice constant a = 5.431 Å = 10.2638... Bohr
- primitive volume Ω = a³/4 Bohr³
- reciprocal lattice: BCC with |G|² in units of (2π/a)² = 3, 4, 8, 11,
  12, 16, 19, 20, 24, 27, 32, 35, 36, 40, 43, 44, 48, 51, 52, 56 for
  the first 20 non-zero shells.

Output CSV columns
------------------
    shell_index, |G| (Bohr⁻¹), |G|² ((2π/a)²), V_local_G (Ry)
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.special import erf

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BOHR_TO_ANG = 0.529_177_210_903
# e² in Rydberg atomic units: e²/r gives potential in Ry when r in Bohr.
#     1 Ry = e²/(2·a₀); so e² = 2 Ry·Bohr in these units.
E2_RY_BOHR = 2.0

# ---------------------------------------------------------------------------
# Manual UPF v2 parsing (no external dependencies)
# ---------------------------------------------------------------------------


def _extract_attr(text: str, name: str) -> str:
    """Extract XML attribute `name="value"` from *text*."""
    m = re.search(rf'{re.escape(name)}\s*=\s*"([^"]*)"', text)
    if m is None:
        raise ValueError(f"attribute {name!r} not found")
    return m.group(1).strip()


def _extract_block(text: str, tag: str) -> np.ndarray:
    """Extract numeric data between `<TAG ...>` and `</TAG>`."""
    m = re.search(rf"<{re.escape(tag)}\b[^>]*>(.*?)</{re.escape(tag)}>", text, re.DOTALL)
    if m is None:
        raise ValueError(f"tag <{tag}> not found")
    body = m.group(1)
    # Split on whitespace, parse floats; ignore empty fields.
    values = [float(tok) for tok in body.split() if tok.strip()]
    return np.asarray(values, dtype=np.float64)


def parse_upf(path: Path) -> dict:
    """Parse the minimum subset of UPF v2 we need for V_local(G)."""
    text = path.read_text()
    # Attributes live inside <PP_HEADER ...>. Limit attr search to that tag
    # to avoid catching unrelated occurrences in PP_INFO (e.g. "z=14.00").
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
        # Keep native UPF units: r in Bohr, v in Ry.
        "r_bohr": r,
        "rab_bohr": rab,
        "v_local_ry": v_loc,
    }


# ---------------------------------------------------------------------------
# V_local(G) in QE native units (Ry, Bohr)
# ---------------------------------------------------------------------------


def v_local_of_g_qe_units(
    r_bohr: np.ndarray,
    v_ry: np.ndarray,
    z_val: float,
    g_bohr_inv: float,
    omega_bohr3: float,
) -> float:
    """
    Compute V_local(G) [Ry] using the QE erf-subtracted convention.

    At G = 0 the erf form degenerates; we use the bare-Coulomb form which is
    algebraically identical in the limit and numerically clean for the
    ONCVPSP linear mesh (r₀ = 0 handled by setting 1/r contribution to 0 at
    the origin).
    """
    four_pi = 4.0 * math.pi

    if g_bohr_inv < 1e-12:
        # G = 0: (4π/Ω) ∫ r² [V(r) + Z·e²/r] dr
        with np.errstate(divide="ignore", invalid="ignore"):
            coulomb = np.where(r_bohr > 0.0, z_val * E2_RY_BOHR / r_bohr, 0.0)
        integrand = r_bohr**2 * (v_ry + coulomb)
        # Use the UPF-provided r grid directly; simpson's adaptive handling
        # takes care of the non-uniform case (this mesh is uniform dr=0.01).
        integral = simpson(integrand, x=r_bohr)
        return four_pi / omega_bohr3 * integral

    # G ≠ 0: erf-subtracted short-range part plus analytic Coulomb tail.
    g = g_bohr_inv
    gr = g * r_bohr
    with np.errstate(divide="ignore", invalid="ignore"):
        sin_gr_over_g = np.where(gr > 1e-10, np.sin(gr) / g, r_bohr * (1.0 - (gr * gr) / 6.0))
    # Short-range integrand: [r·V(r) + Z·e²·erf(r)] · sin(Gr)/G
    short = r_bohr * v_ry + z_val * E2_RY_BOHR * erf(r_bohr)
    integrand = short * sin_gr_over_g
    integral = simpson(integrand, x=r_bohr)

    # Analytic Coulomb tail contribution (Fourier transform of Z·e²·erf(r)/r):
    #   FT[Z·e²·erf(r)/r] = 4π·Z·e²·exp(−G²/4)/G²
    tail = four_pi * z_val * E2_RY_BOHR * math.exp(-g * g / 4.0) / (omega_bohr3 * g * g)

    return four_pi / omega_bohr3 * integral - tail


# ---------------------------------------------------------------------------
# Si FCC |G| shells (in Bohr⁻¹)
# ---------------------------------------------------------------------------


def si_fcc_g_shells(n_shells: int, a_bohr: float) -> list[tuple[int, float, float]]:
    """
    Return list of (|G|² in (2π/a)² units, |G| in Bohr⁻¹, |G|² in Bohr⁻²) for
    the first *n_shells* distinct non-zero |G| shells of FCC with lattice
    constant *a_bohr* (Bohr). We include |G|=0 as shell index 0 at the caller.

    The reciprocal lattice of FCC is BCC with primitive vectors
        b1 = (2π/a)·(-1, 1, 1)
        b2 = (2π/a)·( 1,-1, 1)
        b3 = (2π/a)·( 1, 1,-1)
    """
    tpba = 2.0 * math.pi / a_bohr
    seen: dict[int, None] = {}
    # Enumerate G = n1 b1 + n2 b2 + n3 b3 with |n_i| ≤ n_max.
    # The integer |G|² is |n1(-1,1,1) + n2(1,-1,1) + n3(1,1,-1)|² in (2π/a)².
    n_max = 8
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
    return [
        (g2_int, tpba * math.sqrt(g2_int), (tpba**2) * g2_int)
        for g2_int in shells_int
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    upf_path = repo_root / "pseudopotentials" / "nc" / "lda" / "Si.upf"
    out_csv = repo_root / "scripts" / "validate" / "vloc_g_si_reference.csv"

    if not upf_path.exists():
        print(f"ERROR: Si UPF not found at {upf_path}", file=sys.stderr)
        return 1

    pp = parse_upf(upf_path)
    print(f"Loaded Si UPF: Z_val={pp['z_valence']}, mesh={pp['mesh_size']}")
    print(
        f"  r_bohr: [{pp['r_bohr'][0]:.4e}, {pp['r_bohr'][-1]:.4f}]"
        f"  rab[0]={pp['rab_bohr'][0]:.4e}"
    )
    print(
        f"  v_local_ry: [{pp['v_local_ry'][0]:.4e}, {pp['v_local_ry'][-1]:.4e}]"
    )

    # Si FCC geometry
    a_ang = 5.431
    a_bohr = a_ang / BOHR_TO_ANG
    omega_bohr3 = a_bohr**3 / 4.0
    tpba = 2.0 * math.pi / a_bohr
    print(f"Si FCC: a = {a_ang} Å = {a_bohr:.6f} Bohr;  Ω = {omega_bohr3:.4f} Bohr³")
    print(f"  2π/a = {tpba:.6f} Bohr⁻¹")

    shells = si_fcc_g_shells(20, a_bohr)

    # Compute V_local(G) for each shell.
    rows: list[tuple[int, float, int, float]] = []
    print()
    print(
        f"{'shell':>5}  {'|G|² (int)':>10}  {'|G| (Bohr⁻¹)':>14}  {'V_loc(G) (Ry)':>16}"
    )
    print("-" * 56)
    for idx, (g2_int, g_bohr_inv, _g2_bohr2) in enumerate(shells):
        v_g = v_local_of_g_qe_units(
            pp["r_bohr"], pp["v_local_ry"], pp["z_valence"], g_bohr_inv, omega_bohr3
        )
        rows.append((idx, g_bohr_inv, g2_int, v_g))
        print(f"{idx:>5d}  {g2_int:>10d}  {g_bohr_inv:>14.6f}  {v_g:>16.8e}")

    # Write CSV
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["shell_index", "g_bohr_inv", "g2_units_of_tpba2", "v_local_g_ry"])
        for row in rows:
            w.writerow([row[0], f"{row[1]:.12e}", row[2], f"{row[3]:.12e}"])

    print(f"\nWrote {len(rows)} shells to {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
