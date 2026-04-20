"""Filesystem layout for the pwdft-rs workspace.

One source of truth for every path used by a validation script. Scripts import
these constants instead of re-deriving paths with ``Path(__file__).parents[…]``
chains, which break whenever the tree is reorganized.
"""

from __future__ import annotations

from pathlib import Path

from pyprojroot import find_root, has_file

PROJECT_ROOT = find_root(has_file("uv.lock"))

DATA_DIR: Path = PROJECT_ROOT / "validation"
QE_REF_DIR: Path = DATA_DIR / "reference" / "qe"
CSV_REF_DIR: Path = DATA_DIR / "reference" / "csv"

PSEUDO_DIR: Path = PROJECT_ROOT / "pseudopotentials"
INPUTS_DIR: Path = PROJECT_ROOT / "inputs"
