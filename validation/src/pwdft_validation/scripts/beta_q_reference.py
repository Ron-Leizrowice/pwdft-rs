#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy>=1.26",
#     "scipy>=1.11",
# ]
# ///
"""
VGCMP Phase 2 — Independent Python reference for KB projector form factors β_l(q).

Purpose
-------
For each β projector in `pseudopotentials/nc/lda/Si.upf`, compute the
Bessel-transformed form factor

    F_l(q) = 4π ∫₀^∞ χ(r) · j_l(q·r) · r dr

where χ(r) is the UPF-stored radial function (QE convention: `upf%beta(r,nb)`
**is already** r·β(r), i.e. the stored array has the factor of r baked in).
This matches QE 7.5 `upflib/beta_mod.f90:111-116`:

    call sph_bes (kkbeta, r, qi, l, besr)
    aux(ir) = upf(nt)%beta(ir, nb) * besr(ir) * r(ir)
    call simpson (kkbeta, aux, rab, vqint)
    tab_beta(iq, nb, nt) = vqint * (4π / √Ω)

We strip the `(4π / √Ω)` normalization; it enters through the Hamiltonian
assembly, not the form factor itself. Concretely this script computes the
bare `4π · Simpson[χ · j_l · r]` integral — *exactly* what pwdft-rs's
`bessel_transform_projector` (`src/potential/nonlocal.rs:223`) computes.

Output CSV columns
------------------
    projector_index, l, q_bohr_inv, F_l_q_bohr_3halves

Units throughout the CSV: QE native (Bohr, Ry).
- q in Bohr⁻¹
- F_l(q) in Bohr^(3/2)
  (since χ is in Bohr^(-1/2), r in Bohr: F = 4π · [Bohr^(-1/2)] · 1 · [Bohr] · [Bohr] = Bohr^(3/2))

Run
---
    uv run validation/src/pwdft_validation/scripts/beta_q_reference.py
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
# Constants
# ---------------------------------------------------------------------------

BOHR_TO_ANG = 0.529_177_210_903

# q-grid: 20 equally-spaced values in [0.1, 7.0] Bohr⁻¹.
# q_max ≈ 7 Bohr⁻¹ corresponds to ecut ≈ 25 Ry = q_max²/2, covering the
# standard Si plane-wave envelope. q_min = 0.1 avoids the q=0 edge where
# F_l(q=0) is identically 0 for all l > 0 and the comparison degenerates.
Q_MIN_BOHR = 0.1
Q_MAX_BOHR = 7.0
N_Q = 20

# ---------------------------------------------------------------------------
# Manual UPF v2 parsing (no external dependencies)
# ---------------------------------------------------------------------------


def _extract_attr(text: str, name: str) -> str:
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
    values = [float(tok) for tok in body.split() if tok.strip()]
    return np.asarray(values, dtype=np.float64)


def _extract_beta_projector(text: str, idx: int, mesh: int) -> tuple[int, np.ndarray]:
    """Extract (l, χ(r) values) for PP_BETA.idx.

    The PP_BETA tag header contains `angular_momentum="L"`; the body contains
    `mesh` floats (the stored array is r·β(r), a.k.a. χ, in Bohr^(-1/2)).
    """
    tag = f"PP_BETA.{idx}"
    # Extract the tag with its attributes.
    m = re.search(
        rf"<{re.escape(tag)}\b([^>]*)>(.*?)</{re.escape(tag)}>",
        text,
        re.DOTALL,
    )
    if m is None:
        raise ValueError(f"tag <{tag}> not found")
    attrs = m.group(1)
    body = m.group(2)
    l_str = _extract_attr(attrs, "angular_momentum")
    l = int(l_str)
    values = np.asarray(
        [float(tok) for tok in body.split() if tok.strip()],
        dtype=np.float64,
    )
    if values.size != mesh:
        raise ValueError(f"{tag}: expected {mesh} values, got {values.size}")
    return l, values


def parse_upf(path: Path) -> dict:
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

    return {
        "z_valence": z_val,
        "mesh_size": mesh,
        "n_proj": n_proj,
        "r_bohr": r,
        "rab_bohr": rab,
        "projectors": projectors,
    }


# ---------------------------------------------------------------------------
# F_l(q) in QE native units (Bohr)
# ---------------------------------------------------------------------------


def f_l_of_q(
    r_bohr: np.ndarray,
    chi: np.ndarray,
    l: int,
    q_bohr_inv: float,
) -> float:
    """
    F_l(q) = 4π ∫₀^∞ χ(r) · j_l(q·r) · r dr

    `chi` is the UPF-stored quantity (= r·β(r) in Bohr^(-1/2)). `r` in Bohr.
    `q` in Bohr⁻¹. Result in Bohr^(3/2).

    Uses `scipy.integrate.simpson(y, x=r)` on the UPF mesh (non-uniform allowed).
    Uses `scipy.special.spherical_jn(l, x)` which implements j_l via the standard
    recurrences and is numerically stable for l ≤ ~10 and moderate x.
    """
    qr = q_bohr_inv * r_bohr
    jl = spherical_jn(l, qr)
    integrand = chi * jl * r_bohr
    integral = simpson(integrand, x=r_bohr)
    return 4.0 * math.pi * integral


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    repo_root = Path(__file__).resolve().parents[4]
    upf_path = repo_root / "pseudopotentials" / "nc" / "lda" / "Si.upf"
    out_csv = repo_root / "validation" / "reference" / "csv" / "beta_q_si_reference.csv"

    if not upf_path.exists():
        print(f"ERROR: Si UPF not found at {upf_path}", file=sys.stderr)
        return 1

    pp = parse_upf(upf_path)
    print(f"Loaded Si UPF: Z_val={pp['z_valence']}, mesh={pp['mesh_size']}, n_proj={pp['n_proj']}")
    print(f"  r_bohr: [{pp['r_bohr'][0]:.4e}, {pp['r_bohr'][-1]:.4f}]  rab[0]={pp['rab_bohr'][0]:.4e}")

    # Build q-grid (inclusive bounds).
    q_grid = np.linspace(Q_MIN_BOHR, Q_MAX_BOHR, N_Q)
    print(f"q-grid: {N_Q} values in [{Q_MIN_BOHR}, {Q_MAX_BOHR}] Bohr⁻¹")

    # Print a header table to stdout; project angular momenta.
    angular_momenta = [l for (l, _) in pp["projectors"]]
    print(f"projectors: l = {angular_momenta}")

    # Compute F_l(q) for each projector and q.
    rows: list[tuple[int, int, float, float]] = []
    print()
    print(f"{'proj':>4}  {'l':>2}  {'q (Bohr⁻¹)':>11}  {'F_l(q) (Bohr^3/2)':>22}")
    print("-" * 50)
    for pi, (l, chi) in enumerate(pp["projectors"]):
        for q in q_grid:
            f_val = f_l_of_q(pp["r_bohr"], chi, l, float(q))
            rows.append((pi, l, float(q), f_val))
            print(f"{pi:>4d}  {l:>2d}  {q:>11.6f}  {f_val:>22.12e}")

    # Write CSV
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["projector_index", "l", "q_bohr_inv", "F_l_q_bohr_3halves"])
        for pi, l, q, f_val in rows:
            w.writerow([pi, l, f"{q:.12e}", f"{f_val:.12e}"])

    print(f"\nWrote {len(rows)} rows to {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
