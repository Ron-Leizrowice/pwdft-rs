"""KB projector form-factor β_l(q) generators (VGCMP Phase 2 and VGCH Phase 1b).

Covers:
- ``generate_si``:     β_l(q) for Si using scipy Simpson on a q-grid.
- ``generate_heavy``:  β_l(q) for 11 heavy-atom PPs using QE-style quadrature.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import spherical_jn

from pwdft_validation.integrate import simpson_qe
from pwdft_validation.units import BOHR_TO_ANG
from pwdft_validation.upf import UpfData, parse_upf

# ---------------------------------------------------------------------------
# Core math
# ---------------------------------------------------------------------------

_Q_GRID_SI = np.linspace(0.1, 7.0, 20)

_Q_GRID_HEAVY = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

_HEAVY_ELEMENTS = ["Si", "C", "Al", "Fe", "Cu", "Ga", "As", "Na", "Cl", "Mg", "O"]


def f_l_of_q_bohr(pp: UpfData, l: int, chi: np.ndarray, q_bohr_inv: float) -> float:
    """F_l(q) = 4π · Simpson[χ · j_l(qr) · r ; rab]   [Bohr^(3/2)].

    Uses QE-style quadrature (``simpson_qe``) with the log-mesh Jacobian rab.
    ``chi`` = r·β(r) in Bohr^(-1/2) as stored in PP_BETA.
    """
    qr = q_bohr_inv * pp.r_bohr
    jl = spherical_jn(l, qr)
    integrand = chi * jl * pp.r_bohr
    return 4.0 * math.pi * simpson_qe(integrand, pp.rab_bohr)


def f_l_of_q_ang(pp: UpfData, l: int, chi: np.ndarray, q_per_ang: float) -> float:
    """F_l(q) in Å^(3/2) — pwdft-rs internal convention.

    Converts r/rab/chi to Å before integrating so the result matches
    ``src/potential/nonlocal.rs::bessel_transform_projector``.
    The dimensionless j_l argument q·r is invariant under unit change.
    """
    r_ang = pp.r_bohr * BOHR_TO_ANG
    rab_ang = pp.rab_bohr * BOHR_TO_ANG
    chi_ang = chi / math.sqrt(BOHR_TO_ANG)
    qr = q_per_ang * r_ang
    jl = spherical_jn(l, qr)
    integrand = chi_ang * jl * r_ang
    return 4.0 * math.pi * simpson_qe(integrand, rab_ang)


# ---------------------------------------------------------------------------
# Generators
# ---------------------------------------------------------------------------


def generate_si(pseudo_dir: Path, out_csv: Path) -> int:
    """Generate ``beta_q_si_reference.csv`` (VGCMP Phase 2)."""
    upf_path = pseudo_dir / "Si.upf"
    if not upf_path.exists():
        print(f"ERROR: Si UPF not found at {upf_path}", file=sys.stderr)
        return 1

    pp = parse_upf(upf_path)
    print(f"Loaded Si UPF: Z_val={pp.z_valence}, mesh={pp.mesh_size}, n_proj={pp.n_proj}")
    angular_momenta = [l for (l, _) in pp.projectors]
    print(f"projectors: l = {angular_momenta}")
    print(f"q-grid: {len(_Q_GRID_SI)} values in [{_Q_GRID_SI[0]}, {_Q_GRID_SI[-1]}] Bohr⁻¹")

    rows: list[tuple] = []
    print(f"\n{'proj':>4}  {'l':>2}  {'q (Bohr⁻¹)':>11}  {'F_l(q) (Bohr^3/2)':>22}")
    print("-" * 50)
    for pi, (l, chi) in enumerate(pp.projectors):
        for q in _Q_GRID_SI:
            f_val = f_l_of_q_bohr(pp, l, chi, float(q))
            rows.append((pi, l, float(q), f_val))
            print(f"{pi:>4d}  {l:>2d}  {q:>11.6f}  {f_val:>22.12e}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["projector_index", "l", "q_bohr_inv", "F_l_q_bohr_3halves"])
        for pi, l, q, f_val in rows:
            w.writerow([pi, l, f"{q:.12e}", f"{f_val:.12e}"])
    print(f"\nWrote {len(rows)} rows to {out_csv}")
    return 0


def generate_heavy(pseudo_dir: Path, out_csv: Path) -> int:
    """Generate ``vgch_beta_l_heavy.csv`` (VGCH Phase 1b — 11 heavy-atom PPs)."""
    rows = []
    hdr = f"{'elem':>4} {'i':>3} {'l':>3} {'q (Bohr⁻¹)':>11} {'F_Bohr (Bohr^3/2)':>22} {'F_Å (Å^3/2)':>22}"
    print(hdr)
    print("-" * len(hdr))

    for elem in _HEAVY_ELEMENTS:
        upf_path = pseudo_dir / f"{elem}.upf"
        if not upf_path.exists():
            print(f"# SKIP {elem}: {upf_path} missing", file=sys.stderr)
            continue
        pp = parse_upf(upf_path)
        print(f"# {elem}: Z_val={pp.z_valence:.2f}, n_proj={pp.n_proj}, ls={[l for (l, _) in pp.projectors]}")
        for i, (l, chi) in enumerate(pp.projectors):
            for q_bohr in _Q_GRID_HEAVY:
                f_bohr = f_l_of_q_bohr(pp, l, chi, q_bohr)
                q_ang = q_bohr / BOHR_TO_ANG
                f_ang = f_l_of_q_ang(pp, l, chi, q_ang)
                rows.append((elem, i, l, q_bohr, q_ang, f_bohr, f_ang))
                print(f"{elem:>4} {i:>3} {l:>3} {q_bohr:>11.6f} {f_bohr:>22.12e} {f_ang:>22.12e}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["element", "proj_index", "l", "q_bohr_inv", "q_per_ang", "F_bohr_3halves", "F_ang_3halves"])
        for row in rows:
            w.writerow([row[0], row[1], row[2], f"{row[3]:.12e}", f"{row[4]:.12e}", f"{row[5]:.12e}", f"{row[6]:.12e}"])
    print(f"\nWrote {len(rows)} rows to {out_csv}")
    return 0
