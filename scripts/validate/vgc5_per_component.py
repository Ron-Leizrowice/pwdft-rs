#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""
VGC5 — Per-component energy accounting, Si diamond and Fe BCC vs QE 7.5.

Purpose
-------
VGCMP Phases 1-4 established that the pseudopotential -> Hamiltonian
assembly pipeline is bit-correct vs QE. The remaining Si 13.4 eV gap
is outside that pipeline: it must localize to one of the per-term energy
contributions (kinetic, local, non-local, Hartree, XC, Ewald) or to a
global bookkeeping term such as the V_local(G=0) background shift.

This script parses QE's standard-output energy decomposition for Si and
Fe and emits a machine-readable CSV plus a side-by-side table. The Rust
side exposes the same decomposition via `scf::EnergyComponents`; the
integration test `tests/vgc5_per_component_si.rs` pins the pwdft-rs
numbers and prints the same table for visual diffing.

QE decomposition
----------------
In QE's output (verbose printout at end of a converged run):

    one-electron contribution = E_kin + E_local + E_nonlocal + V_loc(G=0)*N
    hartree contribution      = E_H
    xc contribution           = E_xc
    ewald contribution        = E_ewald
    total energy              = one_electron + hartree + xc + ewald (- TS)

All QE values are in Ry; we convert to eV using
1 Ry = 13.605_693_122_994 eV.

pwdft-rs decomposition
----------------------
pwdft-rs computes a strict per-component breakdown via `EnergyComponents`:

    e_kinetic       = < psi | T | psi >           (computed directly)
    e_local         = int rho(r) V_local(r) dr    (G != 0 piece)
    e_local_g0_shift= V_local(G=0) * N_el         (constant background)
    e_nonlocal      = < psi | V_NL | psi >        (computed directly)
    e_hartree       = 1/2 int rho V_H dr
    e_xc            = int rho eps_xc dr           (bare, same sign as QE)
    e_ewald         = ion-ion Ewald

These sum to the total:
    E_total = e_kinetic + e_local + e_local_g0_shift + e_nonlocal
            + e_hartree + e_xc + e_ewald
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

RY_TO_EV = 13.605_693_122_994


def parse_qe_output(path: Path) -> dict[str, float]:
    """Extract QE's energy decomposition (in Ry) from a converged pw.x output."""
    text = path.read_text()

    patterns = {
        "total_energy": r"!\s*total energy\s*=\s*([-+\d.Ee]+)\s*Ry",
        "one_electron": r"one-electron contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
        "hartree":      r"hartree contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
        "xc":           r"xc contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
        "ewald":        r"ewald contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
        "smearing_ts":  r"smearing contrib\. \(-TS\)\s*=\s*([-+\d.Ee]+)\s*Ry",
        "internal_E":   r"internal energy E=F\+TS\s*=\s*([-+\d.Ee]+)\s*Ry",
        "fermi_ev":     r"the Fermi energy is\s*([-+\d.Ee]+)\s*ev",
    }

    out: dict[str, float] = {}
    for key, pat in patterns.items():
        m = re.search(pat, text)
        if m is None:
            print(f"  WARNING: {key} not found in {path.name}", file=sys.stderr)
            continue
        out[key] = float(m.group(1))
    return out


def qe_to_ev(data_ry: dict[str, float]) -> dict[str, float]:
    """Convert Ry-valued QE quantities to eV, preserving Fermi energy (already eV)."""
    return {
        k: (v * RY_TO_EV if k != "fermi_ev" else v)
        for k, v in data_ry.items()
    }


def print_table(system: str, qe_ev: dict[str, float]) -> None:
    print(f"\n=== QE reference — {system} ===")
    print(f"  total energy  = {qe_ev['total_energy']:>14.6f} eV")
    print(f"  - TS          = {qe_ev.get('smearing_ts', 0.0):>14.6f} eV")
    print(f"  internal E    = {qe_ev.get('internal_E', 0.0):>14.6f} eV")
    print(f"                  (= F + TS = F - (-TS))")
    print(f"  one-electron  = {qe_ev['one_electron']:>14.6f} eV   (= E_kin + E_loc + E_NL + V_loc(G=0)*N)")
    print(f"  hartree       = {qe_ev['hartree']:>14.6f} eV")
    print(f"  xc            = {qe_ev['xc']:>14.6f} eV")
    print(f"  ewald         = {qe_ev['ewald']:>14.6f} eV")
    # Sanity check
    s = qe_ev["one_electron"] + qe_ev["hartree"] + qe_ev["xc"] + qe_ev["ewald"]
    # QE's "total" = internal + (-TS); s should equal internal_E
    sum_err = s - qe_ev.get("internal_E", s)
    print(f"  [sum check]   sum(4 terms) = {s:.6f} eV, internal_E = "
          f"{qe_ev.get('internal_E', float('nan')):.6f} eV, Δ = {sum_err:.2e} eV")
    print(f"  fermi energy  = {qe_ev.get('fermi_ev', float('nan')):.4f} eV")


def write_csv(path: Path, system: str, qe_ev: dict[str, float]) -> None:
    rows = [
        ("total_energy", qe_ev["total_energy"]),
        ("internal_energy", qe_ev.get("internal_E", float("nan"))),
        ("smearing_ts", qe_ev.get("smearing_ts", 0.0)),
        ("one_electron", qe_ev["one_electron"]),
        ("hartree", qe_ev["hartree"]),
        ("xc", qe_ev["xc"]),
        ("ewald", qe_ev["ewald"]),
        ("fermi", qe_ev.get("fermi_ev", float("nan"))),
    ]
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["component_ev", "value_ev"])
        w.writerows(rows)
    print(f"  wrote {path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    qe_root = (here.parent.parent / "qe_validation").resolve()
    parser.add_argument("--qe-dir", type=Path, default=qe_root,
                        help="Directory containing si_scf.out and fe_bcc_fm_scf.out")
    parser.add_argument("--out-dir", type=Path, default=here,
                        help="Directory to write vgc5_qe_*.csv")
    args = parser.parse_args()

    for system, infile, outfile in [
        ("Si diamond",  "si_scf.out",        "vgc5_qe_si_components.csv"),
        ("Fe BCC (FM)", "fe_bcc_fm_scf.out", "vgc5_qe_fe_components.csv"),
    ]:
        qe_path = args.qe_dir / infile
        if not qe_path.exists():
            print(f"  SKIP {system}: {qe_path} missing", file=sys.stderr)
            continue
        ry = parse_qe_output(qe_path)
        ev = qe_to_ev(ry)
        print_table(system, ev)
        write_csv(args.out_dir / outfile, system, ev)

    return 0


if __name__ == "__main__":
    sys.exit(main())
