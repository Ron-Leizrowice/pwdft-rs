"""pwdft-rs validation helpers — QE cross-checks, reference-data generators.

Submodules:
    scripts/    — numbered validation scripts (one per investigation).
    _paths      — shared filesystem layout constants.
"""

from pwdft_validation._paths import (
    CSV_REF_DIR,
    INPUTS_DIR,
    PSEUDO_DIR,
    QE_REF_DIR,
    REPO_ROOT,
    VALIDATION_DIR,
)

__all__ = [
    "CSV_REF_DIR",
    "INPUTS_DIR",
    "PSEUDO_DIR",
    "QE_REF_DIR",
    "REPO_ROOT",
    "VALIDATION_DIR",
]
