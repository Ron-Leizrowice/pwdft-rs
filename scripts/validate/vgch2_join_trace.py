#!/usr/bin/env python3
"""VGCH-2 Part A — join QE + pwdft per-term CSVs into a single delta table.

Reads two inputs:

- ``scripts/validate/vgch2_per_term_trace.csv``
  (produced by ``vgch2_per_term_trace.py``; columns
   ``system, term_name, qe_value_ry, qe_value_eV``)
- ``target/tmp/vgch2_per_term_trace_pwdft.csv``
  (emitted by ``tests/vgch_per_component_heavy.rs`` at Tier-2 run time;
   columns ``system, term_name, pwdft_value_eV``)

Emits a joined CSV (stdout) with columns:

    system, term_name, qe_value_eV, pwdft_value_eV, delta_meV

and a short stderr summary ranking |delta_meV| per system over the
canonical 4 QE terms (one_electron, hartree, xc, ewald) plus total.

Rows present in only one file are dropped with a stderr warning.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

QE_CANON_TERMS = ["one_electron", "hartree", "xc", "ewald", "total"]
SYSTEM_ORDER = [
    "Si",
    "C_diamond",
    "Al",
    "Fe_BCC_FM",
    "Cu_FCC",
    "GaAs",
    "NaCl",
    "MgO",
]


def load_qe(path: Path) -> dict[tuple[str, str], float]:
    """Return (system, term) -> eV."""
    data: dict[tuple[str, str], float] = {}
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            term = row["term_name"]
            if term == "smearing_mts" or term == "fermi_eV":
                continue  # not per-term of E
            ev = float(row["qe_value_eV"])
            data[(row["system"], term)] = ev
    return data


def load_pwdft(path: Path) -> dict[tuple[str, str], float]:
    data: dict[tuple[str, str], float] = {}
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            ev = float(row["pwdft_value_eV"])
            data[(row["system"], row["term_name"])] = ev
    return data


def main() -> None:
    root = Path(__file__).resolve().parent.parent.parent
    qe_path = root / "scripts" / "validate" / "vgch2_per_term_trace.csv"
    pwdft_path = root / "target" / "tmp" / "vgch2_per_term_trace_pwdft.csv"

    if not qe_path.is_file():
        print(f"error: {qe_path} missing; run vgch2_per_term_trace.py first", file=sys.stderr)
        sys.exit(1)
    if not pwdft_path.is_file():
        print(
            f"error: {pwdft_path} missing; run `cargo test --release "
            "--test vgch_per_component_heavy -- --ignored`",
            file=sys.stderr,
        )
        sys.exit(1)

    qe = load_qe(qe_path)
    pw = load_pwdft(pwdft_path)

    writer = csv.writer(sys.stdout)
    writer.writerow(
        ["system", "term_name", "qe_value_eV", "pwdft_value_eV", "delta_meV"]
    )

    deltas_by_system: dict[str, dict[str, float]] = {}
    for system in SYSTEM_ORDER:
        deltas_by_system[system] = {}
        for term in QE_CANON_TERMS:
            key = (system, term)
            if key not in qe or key not in pw:
                print(f"warning: missing key {key}", file=sys.stderr)
                continue
            delta_ev = pw[key] - qe[key]
            delta_mev = delta_ev * 1000.0
            deltas_by_system[system][term] = delta_mev
            writer.writerow(
                [
                    system,
                    term,
                    f"{qe[key]:.6f}",
                    f"{pw[key]:.6f}",
                    f"{delta_mev:+.1f}",
                ]
            )

    # stderr ranking: per system, sort |delta_meV| by term over the 4 QE terms.
    print("\n# Per-system |delta_meV| ranking (top-3 terms):", file=sys.stderr)
    print(
        f"# {'system':<12}  {'term':<14}  {'delta_meV':>12}  {'cum_frac':>8}",
        file=sys.stderr,
    )
    for system in SYSTEM_ORDER:
        deltas = deltas_by_system.get(system, {})
        # Exclude 'total' from the ranking — we want to see which component
        # carries the residual, not the residual itself.
        per_term = {t: d for t, d in deltas.items() if t != "total"}
        if not per_term:
            continue
        total_abs = sum(abs(d) for d in per_term.values())
        ranked = sorted(per_term.items(), key=lambda kv: abs(kv[1]), reverse=True)
        cum = 0.0
        for term, d in ranked[:3]:
            cum += abs(d)
            frac = cum / total_abs if total_abs > 0 else 0.0
            print(
                f"  {system:<12}  {term:<14}  {d:>+12.1f}  {frac:>8.2%}",
                file=sys.stderr,
            )
        # Print E_total residual for reference.
        e_total = deltas.get("total", 0.0)
        print(
            f"  {system:<12}  {'[total]':<14}  {e_total:>+12.1f}  {'(residual)':>8}",
            file=sys.stderr,
        )
        print("  --", file=sys.stderr)


if __name__ == "__main__":
    main()
