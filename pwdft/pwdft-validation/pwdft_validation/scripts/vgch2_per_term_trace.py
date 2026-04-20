#!/usr/bin/env python3
"""VGCH-2 Part A — parse per-term energies from QE 7.5 reference outputs.

Emits a CSV (stdout) with the canonical QE total-energy decomposition
printed after the ``!    total energy`` line in ``pw.x``'s stdout:

    E = E_band_sum_from_one_electron + E_hartree + E_xc + E_ewald

where (in QE's print convention, all in Ry) each contribution is
literally the number printed on the ``one-electron contribution``,
``hartree contribution``, ``xc contribution``, ``ewald contribution``
lines. For metals QE also prints a ``smearing contrib. (-TS)`` entry
and the ``!    total energy`` value is ``E − TS``. This script captures
all five for each system.

Columns:

    system, term_name, qe_value_ry, qe_value_eV

Systems covered (8): Si, C diamond, Al, Fe BCC (nspin=2 FM), Cu FCC,
GaAs, NaCl, MgO. Reads from the checked-in ``data/qe/<system>.out``
files; does NOT rerun QE.

Usage:
    uv run pwdft/pwdft-validation/pwdft_validation/scripts/vgch2_per_term_trace.py \\
        > data/csv/vgch2_per_term_trace.csv

Cross-references QE source lines that print these fields:
    qe-7.5/PW/src/electrons.f90:1612-1621 (print_energies).

Reference: VGCH-2 Part A, proposals/VGCH-2-total-energy-assembly.md.
"""

from __future__ import annotations

import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

RY_TO_EV = 13.605_693_122_994

# System -> QE output filename stem under data/qe/.
SYSTEMS: list[tuple[str, str]] = [
    ("Si", "si_scf"),
    ("C_diamond", "c_diamond_scf"),
    ("Al", "al_fcc_scf"),
    ("Fe_BCC_FM", "fe_bcc_fm_scf"),
    ("Cu_FCC", "cu_fcc_scf"),
    ("GaAs", "gaas_scf"),
    ("NaCl", "nacl_scf"),
    ("MgO", "mgo_scf"),
]

# Each regex captures a signed decimal before ``Ry`` on a dedicated line.
# QE prints to fixed precision in Ry; we don't try to parse multi-line
# sub-term breakdowns (``sum mid..``, ``local..``, ``nl-pp..``) because
# those are only printed under ``verbosity='high'`` and we treat the
# top-level terms as the canonical per-term reference.
TERM_PATTERNS: dict[str, re.Pattern[str]] = {
    "one_electron": re.compile(r"one-electron contribution\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "hartree": re.compile(r"hartree contribution\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "xc": re.compile(r"xc contribution\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "ewald": re.compile(r"ewald contribution\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "smearing_mts": re.compile(r"smearing contrib\.?\s*\(-TS\)\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "total": re.compile(r"^!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry", re.MULTILINE),
    "fermi_eV": re.compile(r"the Fermi energy is\s+(-?\d+\.\d+)\s*ev"),
}


@dataclass
class QeTerms:
    system: str
    path: Path
    values: dict[str, float]


def parse_qe_out(system: str, out_path: Path) -> QeTerms:
    """Pull the per-term breakdown from a ``pw.x`` stdout file.

    Returns the *last* match for each regex, which corresponds to the
    final SCF iteration's printout.
    """
    text = out_path.read_text()
    values: dict[str, float] = {}
    for term, pat in TERM_PATTERNS.items():
        matches = pat.findall(text)
        if not matches:
            raise RuntimeError(f"{system}: no match for {term} in {out_path}")
        # Take the last occurrence (post-convergence printout).
        values[term] = float(matches[-1])
    return QeTerms(system=system, path=out_path, values=values)


def main() -> None:
    qe_dir = Path(__file__).resolve().parents[4] / "validation" / "reference" / "qe"
    if not qe_dir.is_dir():
        print(f"error: qe_validation directory not found at {qe_dir}", file=sys.stderr)
        sys.exit(1)

    all_terms: list[QeTerms] = []
    for system, stem in SYSTEMS:
        out_path = qe_dir / f"{stem}.out"
        if not out_path.is_file():
            print(f"warning: skipping {system}: {out_path} not found", file=sys.stderr)
            continue
        all_terms.append(parse_qe_out(system, out_path))

    writer = csv.writer(sys.stdout)
    writer.writerow(["system", "term_name", "qe_value_ry", "qe_value_eV"])
    for terms in all_terms:
        for term_name, ry_value in terms.values.items():
            if term_name == "fermi_eV":
                # Fermi is already in eV in QE's stdout; Ry column is blank.
                writer.writerow([terms.system, term_name, "", f"{ry_value:.6f}"])
            else:
                ev_value = ry_value * RY_TO_EV
                writer.writerow([terms.system, term_name, f"{ry_value:.8f}", f"{ev_value:.6f}"])

    print(
        f"# parsed {len(all_terms)} systems from {qe_dir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
