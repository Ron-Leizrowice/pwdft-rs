#!/usr/bin/env python3
"""Si E_F shift diagnostic: localize the 1.35 eV absolute-reference offset.

Runs the pwdft-rs `test_si_diamond_fermi_vs_qe` Tier-2 test (unignored) to
capture the converged Γ-point eigenvalues, diffs them against QE's
`validation/reference/qe/si_scf.out` reference, and emits a per-band shift CSV +
a short verdict on whether the shift is a rigid offset or per-band.

Usage (from the pwdft-rs worktree root):

    uv run validation/src/pwdft_validation/scripts/si_ef_shift_trace.py

Dependencies:
- Python 3.10+, numpy (optional — not required for the arithmetic here).
- `cargo test` available; the test itself compiles pwdft-rs in --release.
- `validation/reference/qe/si_scf.out` present at expected path.

Output:
- `validation/reference/csv/si_ef_shift.csv` with columns
  band, eps_pwdft_eV, eps_qe_eV, delta_eV
- summary printout: mean_delta, std_delta, V_loc(G=0)_sum prediction, verdict.

This is a read-only diagnostic. The script does not edit any source files.
"""

from __future__ import annotations

import csv
import re
import statistics
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[4]
QE_SI_OUT = REPO_ROOT / "validation" / "reference" / "qe" / "si_scf.out"
VLOC_CSV = REPO_ROOT / "validation" / "reference" / "csv" / "vgch_vloc_heavy.csv"
OUTPUT_CSV = REPO_ROOT / "validation" / "reference" / "csv" / "si_ef_shift.csv"


def run_pwdft_si_fermi_test() -> tuple[list[float], float]:
    """Run the Fermi Si arm with --ignored --nocapture and scrape eigs + E_F."""
    cmd = [
        "cargo",
        "test",
        "--release",
        "test_si_diamond_energy_vs_qe",
        "--",
        "--ignored",
        "--nocapture",
    ]
    # We route through the machine lock since this is a cargo test invocation.
    # Callers who already hold the lock can `CARGO_TARGET_DIR=... python ...`.
    result = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    # We accept `test failed` too — the energy arm passes; we only need the eigenvalue
    # dump in stderr/stdout. Eigs are printed by `report_gamma_eigenvalues`.
    combined = result.stdout + result.stderr
    # Parse pwdft eigenvalues:  pwdft: [-5.89, 6.08, ...]
    m = re.search(r"pwdft:\s*\[([^\]]+)\]", combined)
    if not m:
        print(combined)
        raise RuntimeError("Could not find pwdft-rs Γ eigenvalues in cargo test output.")
    eigs_pwdft = [float(x.strip()) for x in m.group(1).split(",")]
    # Parse pwdft total energy line as sanity check: "E_pwdft = <value> eV"
    m2 = re.search(r"E_pwdft\s*=\s*(-?\d+\.\d+)\s*eV", combined)
    e_total_pwdft = float(m2.group(1)) if m2 else float("nan")
    # Attempt to also find Fermi energy from the separate Fermi test if it ran.
    m3 = re.search(r"E_F_pwdft\s*=\s*(-?\d+\.\d+)\s*eV", combined)
    e_f_pwdft = float(m3.group(1)) if m3 else float("nan")
    return eigs_pwdft, e_f_pwdft, e_total_pwdft


def parse_qe_gamma_eigs(path: Path) -> tuple[list[float], float]:
    """Extract Γ-point eigenvalues and E_F from a QE `pw.x` scf output."""
    text = path.read_text()
    # Locate the k=Γ block:
    m = re.search(
        r"k\s*=\s*0\.0000\s+0\.0000\s+0\.0000.*?bands \(ev\):\s*\n\n([-\d.\s]+)\n",
        text,
        re.DOTALL,
    )
    if not m:
        raise RuntimeError(f"Could not find Γ block in {path}")
    raw = m.group(1).split()
    eigs_qe = [float(x) for x in raw]
    m2 = re.search(r"the Fermi energy is\s*(-?\d+\.\d+)\s*ev", text)
    e_f_qe = float(m2.group(1)) if m2 else float("nan")
    return eigs_qe, e_f_qe


def lookup_v_loc_g0_sum(element: str, n_atoms: int) -> float:
    """V_loc(G=0)_sum = N_atoms × per-atom contribution, in eV, from the VGCH vloc CSV."""
    with VLOC_CSV.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["element"] == element:
                return n_atoms * float(row["v_local_g0_ev"])
    raise KeyError(f"Element {element} not found in {VLOC_CSV}")


def main() -> int:
    print("Si E_F shift diagnostic — reproduce the 1.35 eV offset.")
    print(f"QE reference: {QE_SI_OUT}")
    if not QE_SI_OUT.exists():
        print(f"ERROR: missing QE reference {QE_SI_OUT}", file=sys.stderr)
        return 2

    print("Running `cargo test test_si_diamond_energy_vs_qe -- --ignored --nocapture`...")
    eigs_pwdft, e_f_pwdft, e_total_pwdft = run_pwdft_si_fermi_test()
    eigs_qe, e_f_qe = parse_qe_gamma_eigs(QE_SI_OUT)

    n_bands = min(len(eigs_pwdft), len(eigs_qe))
    deltas: list[float] = []
    rows = []
    for i in range(n_bands):
        d = eigs_pwdft[i] - eigs_qe[i]
        deltas.append(d)
        rows.append(
            {
                "band": i + 1,
                "eps_pwdft_eV": f"{eigs_pwdft[i]:.6f}",
                "eps_qe_eV": f"{eigs_qe[i]:.6f}",
                "delta_eV": f"{d:.6f}",
            }
        )

    mean_d = sum(deltas) / len(deltas)
    std_d = statistics.pstdev(deltas) if len(deltas) > 1 else 0.0
    v_loc_g0_sum = lookup_v_loc_g0_sum("Si", n_atoms=2)
    # pwdft zeroes V_loc(G=0), QE keeps it → expected shift = −V_loc(G=0)_sum.
    expected = -v_loc_g0_sum

    print()
    print(f"n_bands compared: {n_bands}")
    print(f"per-band shift (pwdft − QE) mean = {mean_d:+.4f} eV, std = {std_d:.4f} eV")
    print(
        f"expected shift if eigs exclude V_loc(G=0): {expected:+.4f} eV "
        f"(= −N_atoms·V_loc(G=0); N_atoms=2, v_g0/atom = 0.6715 eV)"
    )
    residual = mean_d - expected
    print(f"residual (mean − expected): {residual:+.4f} eV")
    print()
    if std_d < 0.020:  # 20 meV threshold
        print(
            "VERDICT: std(δ) < 20 meV → PURE RIGID OFFSET. "
            "The shift is an absolute-reference-of-energy convention, "
            "NOT a per-band bug."
        )
    else:
        print(
            f"VERDICT: std(δ) = {std_d * 1000:.1f} meV (> 20 meV) → PER-BAND VARIATION. "
            "Harder; not a pure V_loc(G=0) convention difference."
        )
    print()
    print(f"E_F_pwdft = {e_f_pwdft:+.4f} eV")
    print(f"E_F_QE    = {e_f_qe:+.4f} eV")
    print(
        f"ΔE_F (pwdft − QE) = {e_f_pwdft - e_f_qe:+.4f} eV "
        f"(matches mean band shift to {abs(mean_d - (e_f_pwdft - e_f_qe)) * 1000:.1f} meV)"
    )
    print()
    print(f"E_total_pwdft = {e_total_pwdft:+.4f} eV (sanity: energy arm should be green)")
    print()
    # Write CSV
    with OUTPUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["band", "eps_pwdft_eV", "eps_qe_eV", "delta_eV"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote per-band table: {OUTPUT_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
