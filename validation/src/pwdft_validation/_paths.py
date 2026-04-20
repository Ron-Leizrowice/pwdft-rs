"""Filesystem layout for the pwdft-rs workspace.

One source of truth for every path used by a validation script. Scripts import
these constants instead of re-deriving paths with ``Path(__file__).parents[…]``
chains, which break whenever the tree is reorganized.
"""

from __future__ import annotations

from pathlib import Path

# This file lives at ``<repo>/validation/src/pwdft_validation/_paths.py``, so
# ``parents[3]`` is the workspace root.
REPO_ROOT: Path = Path(__file__).resolve().parents[3]

VALIDATION_DIR: Path = REPO_ROOT / "validation"
QE_REF_DIR: Path = VALIDATION_DIR / "reference" / "qe"
CSV_REF_DIR: Path = VALIDATION_DIR / "reference" / "csv"

PSEUDO_DIR: Path = REPO_ROOT / "pseudopotentials"
INPUTS_DIR: Path = REPO_ROOT / "inputs"
