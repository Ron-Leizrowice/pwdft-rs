"""pwdft-rs validation helpers — QE cross-checks, reference-data generators."""

from pwdft_validation.paths import (
    CSV_REF_DIR,
    DATA_DIR,
    INPUTS_DIR,
    PROJECT_ROOT,
    PSEUDO_DIR,
    QE_REF_DIR,
)

# Alias so scripts that import REPO_ROOT still work.
REPO_ROOT = PROJECT_ROOT

__all__ = [
    "CSV_REF_DIR",
    "DATA_DIR",
    "INPUTS_DIR",
    "PROJECT_ROOT",
    "PSEUDO_DIR",
    "QE_REF_DIR",
    "REPO_ROOT",
]
