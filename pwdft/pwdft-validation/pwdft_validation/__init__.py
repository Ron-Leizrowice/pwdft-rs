"""pwdft-rs validation helpers — QE cross-checks, reference-data generators.

Submodules:
    scripts/    — numbered validation scripts (one per investigation).
    _paths      — shared filesystem layout constants.
"""

from pwdft_validation.paths import (
    CSV_REF_DIR,
    DATA_DIR,
    INPUTS_DIR,
    PROJECT_ROOT,
    PSEUDO_DIR,
    QE_REF_DIR,
)

__all__ = [
    "CSV_REF_DIR",
    "DATA_DIR",
    "INPUTS_DIR",
    "PROJECT_ROOT",
    "PSEUDO_DIR",
    "QE_REF_DIR",
]
