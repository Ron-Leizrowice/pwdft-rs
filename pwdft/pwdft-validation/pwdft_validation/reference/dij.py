"""D_ij KB coupling matrix reference generator (VGCMP Phase 3)."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

from pwdft_validation.upf import parse_upf


def _print_dij_grid(l_per_proj: list[int], dij_ry: np.ndarray) -> None:
    n = dij_ry.shape[0]
    header = "         " + "".join(f"  j={j} (l={l_per_proj[j]})   " for j in range(n))
    print(header)
    print("        " + "-" * (15 * n))
    for i in range(n):
        row_label = f"i={i} (l={l_per_proj[i]})"
        cells = "".join(f"{dij_ry[i, j]:+14.6e} " for j in range(n))
        print(f"{row_label:>8}  {cells}")


def generate(pseudo_dir: Path, out_csv: Path) -> int:
    """Generate ``dij_si_reference.csv`` (VGCMP Phase 3)."""
    upf_path = pseudo_dir / "Si.upf"
    if not upf_path.exists():
        print(f"ERROR: Si UPF not found at {upf_path}", file=sys.stderr)
        return 1

    pp = parse_upf(upf_path)
    if pp.dij_ry is None:
        print("ERROR: Si UPF has no PP_DIJ block", file=sys.stderr)
        return 1

    l_per_proj = [l for (l, _) in pp.projectors]
    print(f"Loaded Si UPF: n_proj={pp.n_proj}, l_per_proj={l_per_proj}")
    print("D_ij matrix (Ry), row-major from UPF v2 <PP_DIJ> block:")
    print()
    _print_dij_grid(l_per_proj, pp.dij_ry)
    print()

    mismatches = 0
    max_off = 0.0
    for i in range(pp.n_proj):
        for j in range(pp.n_proj):
            if l_per_proj[i] != l_per_proj[j] and abs(pp.dij_ry[i, j]) > 1e-20:
                mismatches += 1
                max_off = max(max_off, abs(pp.dij_ry[i, j]))
    if mismatches == 0:
        print("Block-diagonal in l: OK (all off-block entries are exactly zero).")
    else:
        print(f"WARN: {mismatches} off-block entries have |D| > 1e-20; max = {max_off:.3e} Ry")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["i", "j", "dij_ry"])
        for i in range(pp.n_proj):
            for j in range(pp.n_proj):
                w.writerow([i, j, f"{pp.dij_ry[i, j]:.16e}"])
    print(f"\nWrote {pp.n_proj * pp.n_proj} rows to {out_csv}")
    return 0
