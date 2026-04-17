#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy>=1.26",
# ]
# ///
"""
VGCMP Phase 3 — Independent Python reference for the D_ij KB coupling matrix.

Purpose
-------
Parse `pseudopotentials/nc/lda/Si.upf` by hand (no external UPF library) and
extract the `<PP_DIJ>` block. UPF v2 stores `dij` in Ry·e (i.e. the KB
coupling coefficient has energy units of Ry; the "·e" is formal since
projectors are dimensionless in Ry atomic units). The matrix is stored
row-major with `n_proj × n_proj` entries.

This is the third of three "form factor" cross-checks (Phase 1: V_local(G);
Phase 2: β_l(q); Phase 3: D_ij). Si's ONCVPSP LDA pseudopotential has
6 projectors (l = 0, 0, 1, 1, 2, 2), so D_ij is 6×6. QE diagonalizes each
l-block, so off-diagonal entries within a block may be non-zero while
different-l blocks are zero.

Output
------
- Prints the full 6×6 matrix in Ry to stdout in a readable grid.
- Writes `scripts/validate/dij_si_reference.csv` with columns:
      i, j, dij_ry
  in row-major order (i is the outer loop).

Pair with `tests/vgcmp_dij_cross_check.rs` which loads the same UPF via
`pseudopotential::load` and asserts element-wise agreement to 1e-12 eV
absolute (pure float round-off modulo the Ry→eV conversion).
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Manual UPF v2 parsing (no external dependencies beyond numpy)
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
    values = [float(tok) for tok in body.split() if tok.strip()]
    return np.asarray(values, dtype=np.float64)


def parse_upf_dij(path: Path) -> tuple[int, list[int], np.ndarray]:
    """
    Parse the minimum subset of UPF v2 we need for D_ij.

    Returns
    -------
    n_proj : int
        Number of KB projectors (from `<PP_HEADER number_of_proj="…">`).
    l_per_proj : list[int]
        Angular momentum of each projector (from `<PP_BETA.i angular_momentum=…>`).
    dij_ry : np.ndarray, shape (n_proj, n_proj)
        D_ij coupling matrix in Ry, reshaped row-major from the flat UPF block.
    """
    text = path.read_text()

    # Attributes live inside <PP_HEADER ...>. Limit search to that tag to
    # avoid catching unrelated occurrences in PP_INFO.
    m_hdr = re.search(r"<PP_HEADER\b([^>]*)/?>", text)
    if m_hdr is None:
        raise ValueError("PP_HEADER not found")
    hdr = m_hdr.group(1)
    n_proj = int(_extract_attr(hdr, "number_of_proj"))

    # Per-projector angular momentum from <PP_BETA.i ...>.
    l_per_proj: list[int] = []
    for i in range(1, n_proj + 1):
        m = re.search(rf"<PP_BETA\.{i}\b([^>]*)>", text)
        if m is None:
            raise ValueError(f"PP_BETA.{i} not found")
        l_per_proj.append(int(_extract_attr(m.group(1), "angular_momentum")))

    # D_ij is a flat list of n_proj² Ry values, row-major.
    dij_flat = _extract_block(text, "PP_DIJ")
    if dij_flat.size != n_proj * n_proj:
        raise ValueError(
            f"PP_DIJ: expected {n_proj * n_proj} values, got {dij_flat.size}"
        )

    dij = dij_flat.reshape((n_proj, n_proj))
    return n_proj, l_per_proj, dij


# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------


def print_dij_grid(l_per_proj: list[int], dij_ry: np.ndarray) -> None:
    """Print the D_ij matrix in a readable grid."""
    n = dij_ry.shape[0]
    # Column header (projector index + l).
    header = "         " + "".join(f"  j={j} (l={l_per_proj[j]})   " for j in range(n))
    print(header)
    print("        " + "-" * (15 * n))
    for i in range(n):
        row_label = f"i={i} (l={l_per_proj[i]})"
        cells = "".join(f"{dij_ry[i, j]:+14.6e} " for j in range(n))
        print(f"{row_label:>8}  {cells}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    upf_path = repo_root / "pseudopotentials" / "nc" / "lda" / "Si.upf"
    out_csv = repo_root / "scripts" / "validate" / "dij_si_reference.csv"

    if not upf_path.exists():
        print(f"ERROR: Si UPF not found at {upf_path}", file=sys.stderr)
        return 1

    n_proj, l_per_proj, dij_ry = parse_upf_dij(upf_path)
    print(f"Loaded Si UPF: n_proj={n_proj}, l_per_proj={l_per_proj}")
    print(f"D_ij matrix (Ry), row-major from UPF v2 <PP_DIJ> block:")
    print()
    print_dij_grid(l_per_proj, dij_ry)
    print()

    # Block-diagonal sanity check: zeros between different-l blocks.
    mismatches = 0
    max_off_block = 0.0
    for i in range(n_proj):
        for j in range(n_proj):
            if l_per_proj[i] != l_per_proj[j] and abs(dij_ry[i, j]) > 1e-20:
                mismatches += 1
                max_off_block = max(max_off_block, abs(dij_ry[i, j]))
    if mismatches == 0:
        print("Block-diagonal in l: OK (all off-block entries are exactly zero).")
    else:
        print(
            f"WARN: {mismatches} off-block entries have |D| > 1e-20; "
            f"max |D_off_block| = {max_off_block:.3e} Ry"
        )

    # Write CSV in row-major order.
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["i", "j", "dij_ry"])
        for i in range(n_proj):
            for j in range(n_proj):
                w.writerow([i, j, f"{dij_ry[i, j]:.16e}"])

    print(f"\nWrote {n_proj * n_proj} rows to {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
