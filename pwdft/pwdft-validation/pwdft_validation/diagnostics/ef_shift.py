"""Si E_F shift diagnostic: localize the rigid absolute-reference offset."""

from __future__ import annotations

import csv
import math
import re
import statistics
import subprocess
import sys
from pathlib import Path

from pwdft_validation.paths import CSV_REF_DIR, PROJECT_ROOT, QE_REF_DIR


def run(repo_root: Path | None = None) -> int:
    """Run ``cargo test test_si_diamond_energy_vs_qe`` and analyse the eigenvalue shift."""
    root = repo_root or PROJECT_ROOT
    qe_si_out = QE_REF_DIR / "si_scf.out"
    vloc_csv = CSV_REF_DIR / "vgch_vloc_heavy.csv"
    output_csv = CSV_REF_DIR / "si_ef_shift.csv"

    print("Si E_F shift diagnostic — reproduce the 1.35 eV offset.")
    if not qe_si_out.exists():
        print(f"ERROR: missing QE reference {qe_si_out}", file=sys.stderr)
        return 2

    print("Running `cargo test test_si_diamond_energy_vs_qe -- --ignored --nocapture`...")
    result = subprocess.run(
        ["cargo", "test", "--release", "test_si_diamond_energy_vs_qe", "--", "--ignored", "--nocapture"],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    combined = result.stdout + result.stderr
    m = re.search(r"pwdft:\s*\[([^\]]+)\]", combined)
    if not m:
        print(combined)
        print("ERROR: Could not find pwdft-rs Γ eigenvalues in cargo test output.", file=sys.stderr)
        return 1
    eigs_pwdft = [float(x.strip()) for x in m.group(1).split(",")]
    m2 = re.search(r"E_pwdft\s*=\s*(-?\d+\.\d+)\s*eV", combined)
    e_total_pwdft = float(m2.group(1)) if m2 else float("nan")
    m3 = re.search(r"E_F_pwdft\s*=\s*(-?\d+\.\d+)\s*eV", combined)
    e_f_pwdft = float(m3.group(1)) if m3 else float("nan")

    text = qe_si_out.read_text()
    mq = re.search(
        r"k\s*=\s*0\.0000\s+0\.0000\s+0\.0000.*?bands \(ev\):\s*\n\n([-\d.\s]+)\n",
        text,
        re.DOTALL,
    )
    if not mq:
        print(f"ERROR: Could not find Γ block in {qe_si_out}", file=sys.stderr)
        return 1
    eigs_qe = [float(x) for x in mq.group(1).split()]
    mf = re.search(r"the Fermi energy is\s*(-?\d+\.\d+)\s*ev", text)
    e_f_qe = float(mf.group(1)) if mf else float("nan")

    n_bands = min(len(eigs_pwdft), len(eigs_qe))
    deltas = [eigs_pwdft[i] - eigs_qe[i] for i in range(n_bands)]
    rows = [
        {
            "band": i + 1,
            "eps_pwdft_eV": f"{eigs_pwdft[i]:.6f}",
            "eps_qe_eV": f"{eigs_qe[i]:.6f}",
            "delta_eV": f"{deltas[i]:.6f}",
        }
        for i in range(n_bands)
    ]

    mean_d = sum(deltas) / len(deltas)
    std_d = statistics.pstdev(deltas) if len(deltas) > 1 else 0.0

    # Read V_loc(G=0) prediction from the vloc CSV if available.
    v_loc_g0_sum = float("nan")
    if vloc_csv.exists():
        with vloc_csv.open() as f:
            for row in csv.DictReader(f):
                if row["element"] == "Si":
                    v_loc_g0_sum = 2 * float(row["v_local_g0_ev"])
                    break
    expected = -v_loc_g0_sum

    print(f"\nn_bands compared: {n_bands}")
    print(f"per-band shift (pwdft − QE) mean = {mean_d:+.4f} eV, std = {std_d:.4f} eV")
    if not math.isnan(expected):
        print(f"expected shift if eigs exclude V_loc(G=0): {expected:+.4f} eV")
        print(f"residual (mean − expected): {mean_d - expected:+.4f} eV")
    print()
    if std_d < 0.020:
        print("VERDICT: std(δ) < 20 meV → PURE RIGID OFFSET.")
    else:
        print(f"VERDICT: std(δ) = {std_d * 1000:.1f} meV (> 20 meV) → PER-BAND VARIATION.")
    print(f"\nE_F_pwdft = {e_f_pwdft:+.4f} eV")
    print(f"E_F_QE    = {e_f_qe:+.4f} eV")
    print(f"E_total_pwdft = {e_total_pwdft:+.4f} eV")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["band", "eps_pwdft_eV", "eps_qe_eV", "delta_eV"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote per-band table: {output_csv}")
    return 0
