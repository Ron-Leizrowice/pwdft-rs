"""Per-term QE energy decomposition: extract, compare, and join traces.

Three functions:
- ``extract_components``: VGC5 per-component for Si + Fe.
- ``extract_trace``:      VGCH-2 8-system per-term trace → CSV.
- ``join_traces``:        Join QE + pwdft trace CSVs, print delta ranking.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

from pwdft_validation.units import RY_TO_EV

# ---------------------------------------------------------------------------
# VGC5 — per-component energy parsing for Si and Fe
# ---------------------------------------------------------------------------

_VGC5_PATTERNS = {
    "total_energy": r"!\s*total energy\s*=\s*([-+\d.Ee]+)\s*Ry",
    "one_electron": r"one-electron contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "hartree": r"hartree contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "xc": r"xc contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "ewald": r"ewald contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "smearing_ts": r"smearing contrib\. \(-TS\)\s*=\s*([-+\d.Ee]+)\s*Ry",
    "internal_E": r"internal energy E=F\+TS\s*=\s*([-+\d.Ee]+)\s*Ry",
    "fermi_ev": r"the Fermi energy is\s*([-+\d.Ee]+)\s*ev",
}


def _parse_vgc5(path: Path) -> dict[str, float]:
    text = path.read_text()
    out: dict[str, float] = {}
    for key, pat in _VGC5_PATTERNS.items():
        m = re.search(pat, text)
        if m is None:
            print(f"  WARNING: {key} not found in {path.name}", file=sys.stderr)
            continue
        out[key] = float(m.group(1))
    return out


def _print_vgc5(system: str, ev: dict[str, float]) -> None:
    print(f"\n=== QE reference — {system} ===")
    print(f"  total energy  = {ev['total_energy']:>14.6f} eV")
    print(f"  - TS          = {ev.get('smearing_ts', 0.0):>14.6f} eV")
    print(f"  internal E    = {ev.get('internal_E', 0.0):>14.6f} eV")
    print(f"  one-electron  = {ev['one_electron']:>14.6f} eV")
    print(f"  hartree       = {ev['hartree']:>14.6f} eV")
    print(f"  xc            = {ev['xc']:>14.6f} eV")
    print(f"  ewald         = {ev['ewald']:>14.6f} eV")
    s = ev["one_electron"] + ev["hartree"] + ev["xc"] + ev["ewald"]
    err = s - ev.get("internal_E", s)
    print(f"  [sum check]   sum={s:.6f} eV, internal_E={ev.get('internal_E', float('nan')):.6f} eV, Δ={err:.2e} eV")
    print(f"  fermi energy  = {ev.get('fermi_ev', float('nan')):.4f} eV")


def extract_components(qe_dir: Path, out_dir: Path) -> int:
    """Generate ``vgc5_qe_si_components.csv`` and ``vgc5_qe_fe_components.csv``."""
    for system, infile, outfile in [
        ("Si diamond", "si_scf.out", "vgc5_qe_si_components.csv"),
        ("Fe BCC (FM)", "fe_bcc_fm_scf.out", "vgc5_qe_fe_components.csv"),
    ]:
        qe_path = qe_dir / infile
        if not qe_path.exists():
            print(f"  SKIP {system}: {qe_path} missing", file=sys.stderr)
            continue
        ry = _parse_vgc5(qe_path)
        ev = {k: (v * RY_TO_EV if k != "fermi_ev" else v) for k, v in ry.items()}
        _print_vgc5(system, ev)
        csv_path = out_dir / outfile
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        rows = [
            ("total_energy", ev["total_energy"]),
            ("internal_energy", ev.get("internal_E", float("nan"))),
            ("smearing_ts", ev.get("smearing_ts", 0.0)),
            ("one_electron", ev["one_electron"]),
            ("hartree", ev["hartree"]),
            ("xc", ev["xc"]),
            ("ewald", ev["ewald"]),
            ("fermi", ev.get("fermi_ev", float("nan"))),
        ]
        with csv_path.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["component_ev", "value_ev"])
            w.writerows(rows)
        print(f"  wrote {csv_path}")
    return 0


# ---------------------------------------------------------------------------
# VGCH-2 Part A — 8-system per-term trace
# ---------------------------------------------------------------------------

_TRACE_SYSTEMS: list[tuple[str, str]] = [
    ("Si", "si_scf"),
    ("C_diamond", "c_diamond_scf"),
    ("Al", "al_fcc_scf"),
    ("Fe_BCC_FM", "fe_bcc_fm_scf"),
    ("Cu_FCC", "cu_fcc_scf"),
    ("GaAs", "gaas_scf"),
    ("NaCl", "nacl_scf"),
    ("MgO", "mgo_scf"),
]

_TRACE_PATTERNS: dict[str, re.Pattern[str]] = {
    "one_electron": re.compile(r"one-electron contribution\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "hartree": re.compile(r"hartree contribution\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "xc": re.compile(r"xc contribution\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "ewald": re.compile(r"ewald contribution\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "smearing_mts": re.compile(r"smearing contrib\.?\s*\(-TS\)\s*=\s*(-?\d+\.\d+)\s*Ry"),
    "total": re.compile(r"^!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry", re.MULTILINE),
    "fermi_eV": re.compile(r"the Fermi energy is\s+(-?\d+\.\d+)\s*ev"),
}


def extract_trace(qe_dir: Path, out_csv: Path) -> int:
    """Generate ``vgch2_per_term_trace.csv`` from 8 QE reference outputs."""
    all_terms = []
    for system, stem in _TRACE_SYSTEMS:
        path = qe_dir / f"{stem}.out"
        if not path.is_file():
            print(f"warning: skipping {system}: {path} not found", file=sys.stderr)
            continue
        text = path.read_text()
        values: dict[str, float] = {}
        for term, pat in _TRACE_PATTERNS.items():
            matches = pat.findall(text)
            if not matches:
                raise RuntimeError(f"{system}: no match for {term} in {path}")
            values[term] = float(matches[-1])
        all_terms.append((system, values))

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "term_name", "qe_value_ry", "qe_value_eV"])
        for system, values in all_terms:
            for term, ry in values.items():
                if term == "fermi_eV":
                    w.writerow([system, term, "", f"{ry:.6f}"])
                else:
                    w.writerow([system, term, f"{ry:.8f}", f"{ry * RY_TO_EV:.6f}"])
    print(f"# parsed {len(all_terms)} systems; wrote {out_csv}", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# VGCH-2 Part A — join QE + pwdft traces
# ---------------------------------------------------------------------------

_JOIN_CANON_TERMS = ["one_electron", "hartree", "xc", "ewald", "total"]
_JOIN_SYSTEM_ORDER = ["Si", "C_diamond", "Al", "Fe_BCC_FM", "Cu_FCC", "GaAs", "NaCl", "MgO"]


def join_traces(qe_csv: Path, pwdft_csv: Path) -> int:
    """Join QE + pwdft per-term CSVs and print delta ranking to stderr."""
    if not qe_csv.is_file():
        print(f"error: {qe_csv} missing; run 'pwdft-validate energy trace' first", file=sys.stderr)
        return 1
    if not pwdft_csv.is_file():
        print(f"error: {pwdft_csv} missing; run the Tier-2 vgch_per_component_heavy test", file=sys.stderr)
        return 1

    qe: dict[tuple[str, str], float] = {}
    with qe_csv.open() as f:
        for row in csv.DictReader(f):
            term = row["term_name"]
            if term in ("smearing_mts", "fermi_eV"):
                continue
            qe[(row["system"], term)] = float(row["qe_value_eV"])

    pw: dict[tuple[str, str], float] = {}
    with pwdft_csv.open() as f:
        for row in csv.DictReader(f):
            pw[(row["system"], row["term_name"])] = float(row["pwdft_value_eV"])

    writer = csv.writer(sys.stdout)
    writer.writerow(["system", "term_name", "qe_value_eV", "pwdft_value_eV", "delta_meV"])
    deltas_by_system: dict[str, dict[str, float]] = {}
    for system in _JOIN_SYSTEM_ORDER:
        deltas_by_system[system] = {}
        for term in _JOIN_CANON_TERMS:
            key = (system, term)
            if key not in qe or key not in pw:
                print(f"warning: missing key {key}", file=sys.stderr)
                continue
            delta_mev = (pw[key] - qe[key]) * 1000.0
            deltas_by_system[system][term] = delta_mev
            writer.writerow([system, term, f"{qe[key]:.6f}", f"{pw[key]:.6f}", f"{delta_mev:+.1f}"])

    print("\n# Per-system |delta_meV| ranking (top-3 terms):", file=sys.stderr)
    print(f"# {'system':<12}  {'term':<14}  {'delta_meV':>12}  {'cum_frac':>8}", file=sys.stderr)
    for system in _JOIN_SYSTEM_ORDER:
        per_term = {t: d for t, d in deltas_by_system.get(system, {}).items() if t != "total"}
        if not per_term:
            continue
        total_abs = sum(abs(d) for d in per_term.values())
        ranked = sorted(per_term.items(), key=lambda kv: abs(kv[1]), reverse=True)
        cum = 0.0
        for term, d in ranked[:3]:
            cum += abs(d)
            frac = cum / total_abs if total_abs > 0 else 0.0
            print(f"  {system:<12}  {term:<14}  {d:>+12.1f}  {frac:>8.2%}", file=sys.stderr)
        e_total = deltas_by_system[system].get("total", 0.0)
        print(f"  {system:<12}  {'[total]':<14}  {e_total:>+12.1f}  {'(residual)':>8}", file=sys.stderr)
        print("  --", file=sys.stderr)
    return 0
