#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""
PCRS — Per-Component Energy Residual Investigation.

Origin
------
VGC5 (PR #34) instrumented per-component SCF energies and observed that for
Si diamond at `conv_threshold = 1e-8`:

    |Sigma(components) - E_total| = 1.20 eV

The order-of-magnitude expectation from converged-density noise is O(1e-6 eV),
six orders of magnitude lower. This script tests whether the residual is
SCF-convergence-dependent (shrinks with tighter conv_threshold) or plateaus
(indicating a structural bookkeeping inconsistency).

Method
------
Generates YAML configs for Si SCF at a sweep of conv_threshold values
(1e-6 .. 1e-12). Runs `cargo run --release -- --input <yaml>` with
`RUST_LOG=info` for each, parses the per-component log lines emitted by
`src/scf/mod.rs:511-530`, and reports:

    conv_threshold | iters | E_total | E_sum | |E_sum - E_total|

Also computes QE's own identity closure from `validation/reference/qe/si_scf.out`
(one_electron + hartree + xc + ewald = internal_E) as the reference.

Usage
-----
    cd .claude/worktrees/<YOUR_WORKTREE>
    uv run validation/src/pwdft_validation/scripts/pcrs_residual_scan.py

Requires the machine lock for cargo commands. This script acquires it
itself via `.claude/bin/machine-lock` if present; otherwise runs raw.

See `proposals/PCRS-per-component-residual.md`.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

RY_TO_EV = 13.605_693_122_994


# ---------------------------------------------------------------------------
# YAML templates
# ---------------------------------------------------------------------------

SI_YAML_TEMPLATE = """\
# PCRS Si SCF scan; matches validation/reference/qe/si_scf.in parameters.
# a = 5.431 A, ecutwfc = 15 Ry = 204.085 eV, 4x4x4 MP, FD smearing.

system:
  lattice:
    - [0.0, 2.7155, 2.7155]
    - [2.7155, 0.0, 2.7155]
    - [2.7155, 2.7155, 0.0]
  atoms:
    - {{ symbol: Si, position: [0.00, 0.00, 0.00] }}
    - {{ symbol: Si, position: [0.25, 0.25, 0.25] }}

basis:
  ecutwfc: 204.0854
  ecutrho_ratio: 4

kpoints:
  type: monkhorst_pack
  grid: [4, 4, 4]

scf:
  max_iter: {max_iter}
  conv_threshold: {conv_threshold:.3e}
  energy_threshold: {energy_threshold:.3e}
  n_bands: 8

electrons:
  mixing_beta: 0.3
  mixing_ndim: 8
  smearing_width: 0.13605693122994   # 0.01 Ry
  mixing_mode: plain

symmetry:
  enabled: {sym_enabled}
  time_reversal: {sym_enabled}

pseudopotentials:
  Si: "{pp_path}"
"""


# Regex helpers for parsing cargo-run stderr (env_logger 'info!' format).
PC_PAT = re.compile(r"E_(band|kinetic|local|nonlocal|hartree|xc|ewald)\s*=\s*([-+\d.Ee]+)")
E_SUM_PAT = re.compile(r"E_sum\(comp\)\s*=\s*([-+\d.Ee]+)\s+\(vs E_KS\s+([-+\d.Ee]+),\s*Δ=([-+\d.Ee]+)\)")
E_LOCAL_G0_PAT = re.compile(r"E_local\(G=0\)\s*=\s*([-+\d.Ee]+)")
ITERS_PAT = re.compile(r"SCF converged after (\d+) iterations")
FINAL_DRHO_PAT = re.compile(r"SCF iter\s+\d+:.*?Δρ=([-+\d.Ee]+)\s*$", flags=re.MULTILINE)


def parse_cargo_stderr(text: str) -> dict[str, float | int | None]:
    """Extract per-component energies + identity closure from a cargo run."""
    out: dict[str, float | int | None] = {
        "e_band": None,
        "e_kinetic": None,
        "e_local": None,
        "e_local_g0_shift": None,
        "e_nonlocal": None,
        "e_hartree": None,
        "e_xc": None,
        "e_ewald": None,
        "e_sum": None,
        "e_ks": None,
        "delta": None,
        "iters": None,
        "final_drho": None,
    }
    # The `E_local` pattern matches both E_local and E_local(G=0) if we're not
    # careful — intentionally use E_local followed by a space or '=' sign.
    # Strategy: E_local(G=0) has its own pattern; collect it first and remove
    # those lines before running the generic pattern.

    g0_match = E_LOCAL_G0_PAT.search(text)
    if g0_match:
        out["e_local_g0_shift"] = float(g0_match.group(1))

    # Strip the E_local(G=0) line so it does not double-match E_local.
    text_no_g0 = E_LOCAL_G0_PAT.sub("", text)

    for m in PC_PAT.finditer(text_no_g0):
        key = f"e_{m.group(1)}"
        # Keep the last occurrence (final iteration block).
        out[key] = float(m.group(2))

    s = E_SUM_PAT.search(text)
    if s:
        out["e_sum"] = float(s.group(1))
        out["e_ks"] = float(s.group(2))
        out["delta"] = float(s.group(3))

    iters = ITERS_PAT.search(text)
    if iters:
        out["iters"] = int(iters.group(1))

    drho_matches = FINAL_DRHO_PAT.findall(text)
    if drho_matches:
        out["final_drho"] = float(drho_matches[-1])

    return out


def run_si_at_threshold(
    repo_root: Path,
    tmpdir: Path,
    conv_threshold: float,
    max_iter: int = 150,
    sym_enabled: bool = True,
) -> dict[str, float | int | None]:
    """Write a YAML, run pwdft-rs, parse output."""
    tag = "sym" if sym_enabled else "nosym"
    yaml_path = tmpdir / f"si_pcrs_{tag}_{conv_threshold:.0e}.yaml"
    # Energy threshold ~ conv_threshold * 10; loose enough to not block density.
    energy_threshold = max(conv_threshold * 10.0, 1e-10)
    pp_path = repo_root / "pseudopotentials" / "nc" / "lda" / "Si.upf"
    yaml_path.write_text(
        SI_YAML_TEMPLATE.format(
            conv_threshold=conv_threshold,
            energy_threshold=energy_threshold,
            max_iter=max_iter,
            sym_enabled="true" if sym_enabled else "false",
            pp_path=str(pp_path),
        )
    )

    env = os.environ.copy()
    env["RUST_LOG"] = "info"
    env["RUST_BACKTRACE"] = "1"

    # Invoke from repo_root so relative pseudopotential paths work.
    cmd = [
        "cargo",
        "run",
        "--release",
        "--quiet",
        "--",
        "--input",
        str(yaml_path),
    ]
    print(f"  running [{tag}] conv_threshold={conv_threshold:.0e} ... ", end="", flush=True)
    try:
        result = subprocess.run(
            cmd,
            cwd=repo_root,
            env=env,
            capture_output=True,
            text=True,
            timeout=600,
        )
    except subprocess.TimeoutExpired:
        print("TIMEOUT")
        return {"error": "timeout", "conv_threshold": conv_threshold}

    parsed = parse_cargo_stderr(result.stderr)
    parsed["conv_threshold"] = conv_threshold
    parsed["returncode"] = result.returncode
    if result.returncode != 0:
        print(f"FAILED (rc={result.returncode})")
        parsed["error"] = "non-zero returncode"
        # Save stderr for post-mortem.
        fail_log = tmpdir / f"si_pcrs_{conv_threshold:.0e}.stderr.log"
        fail_log.write_text(result.stderr)
        print(f"    saved stderr to {fail_log}")
    elif parsed.get("iters") is not None:
        iters = parsed["iters"]
        delta = parsed.get("delta")
        delta_str = f"{delta:+.3e}" if delta is not None else "?"
        print(f"OK ({iters} iters, Δ={delta_str} eV)")
    else:
        print("OK (parse failed)")
        parsed["error"] = "parse failed"
    return parsed


# ---------------------------------------------------------------------------
# QE identity closure
# ---------------------------------------------------------------------------

QE_PATS = {
    "total_energy": r"!\s*total energy\s*=\s*([-+\d.Ee]+)\s*Ry",
    "one_electron": r"one-electron contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "hartree": r"hartree contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "xc": r"xc contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "ewald": r"ewald contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "smearing_ts": r"smearing contrib\. \(-TS\)\s*=\s*([-+\d.Ee]+)\s*Ry",
    "internal_E": r"internal energy E=F\+TS\s*=\s*([-+\d.Ee]+)\s*Ry",
}


def qe_identity_closure(qe_out: Path) -> None:
    """Parse QE out; check one_el+hartree+xc+ewald == internal_E."""
    text = qe_out.read_text()
    vals_ry: dict[str, float] = {}
    for key, pat in QE_PATS.items():
        m = re.search(pat, text)
        if m is None:
            print(f"  WARNING: QE pattern {key} not found")
            continue
        vals_ry[key] = float(m.group(1))

    print(f"\n=== QE identity closure — {qe_out.name} ===")
    print(f"  one_electron   = {vals_ry['one_electron']:>16.10f} Ry  ({vals_ry['one_electron'] * RY_TO_EV:>12.6f} eV)")
    print(f"  hartree        = {vals_ry['hartree']:>16.10f} Ry")
    print(f"  xc             = {vals_ry['xc']:>16.10f} Ry")
    print(f"  ewald          = {vals_ry['ewald']:>16.10f} Ry")
    print(f"  total (F)      = {vals_ry['total_energy']:>16.10f} Ry  ({vals_ry['total_energy'] * RY_TO_EV:>12.6f} eV)")
    print(f"  -TS            = {vals_ry['smearing_ts']:>16.10f} Ry")
    print(f"  internal E     = {vals_ry['internal_E']:>16.10f} Ry  (= F - (-TS))")
    s_ry = vals_ry["one_electron"] + vals_ry["hartree"] + vals_ry["xc"] + vals_ry["ewald"]
    diff_ry = s_ry - vals_ry["internal_E"]
    diff_ev = diff_ry * RY_TO_EV
    print(f"  Sigma(4)       = {s_ry:>16.10f} Ry")
    print(f"  Sigma - intE   = {diff_ry:>+.3e} Ry  ({diff_ev:+.3e} eV)")
    if abs(diff_ev) < 1e-6:
        print("  => QE identity closes to < 1 ueV. Our ~1 eV residual is STRUCTURAL.")
    else:
        print(f"  => QE identity residual = {diff_ev:.3e} eV")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    from pwdft_validation import REPO_ROOT

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT, help="Root of the pwdft-rs repo / worktree")
    parser.add_argument(
        "--thresholds",
        nargs="+",
        type=float,
        default=[1e-6, 1e-8, 1e-9, 1e-10, 1e-11],
        help="Sweep of conv_threshold values (density RMS, e/A^3)",
    )
    parser.add_argument("--no-qe", action="store_true", help="Skip the QE identity-closure comparison")
    parser.add_argument("--keep-yaml", action="store_true", help="Keep tmp YAML configs after run")
    parser.add_argument("--max-iter", type=int, default=200, help="max_iter for pwdft-rs SCF at tightest thresholds")
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    if not (repo_root / "Cargo.toml").is_file():
        print(f"ERROR: {repo_root} is not a pwdft-rs repo", file=sys.stderr)
        return 1

    print(f"PCRS residual scan; repo = {repo_root}")
    print(f"  thresholds: {args.thresholds}")

    # QE reference identity (once).
    if not args.no_qe:
        qe_si = repo_root / "validation" / "reference" / "qe" / "si_scf.out"
        if qe_si.is_file():
            qe_identity_closure(qe_si)
        else:
            print(f"  skipping QE identity (missing {qe_si})")

    # Acquire machine lock for the whole scan if helper present.
    lock = repo_root / ".claude" / "bin" / "machine-lock"
    acquired = False
    if lock.is_file():
        rc = subprocess.call([str(lock), "acquire", "Researcher", "PCRS scan"])
        if rc == 0:
            acquired = True
            print("  machine lock acquired")
        else:
            print("  WARNING: machine lock not acquired (continuing anyway)")

    try:
        with tempfile.TemporaryDirectory(prefix="pcrs_scan_") as td:
            tmpdir = Path(td)

            results_sym: list[dict] = []
            results_nosym: list[dict] = []

            print(f"\nPCRS — Si SCF at {len(args.thresholds)} thresholds, symmetry=ON:")
            for ct in args.thresholds:
                res = run_si_at_threshold(
                    repo_root,
                    tmpdir,
                    ct,
                    max_iter=args.max_iter,
                    sym_enabled=True,
                )
                results_sym.append(res)

            print(f"\nPCRS — Si SCF at {len(args.thresholds)} thresholds, symmetry=OFF:")
            for ct in args.thresholds:
                res = run_si_at_threshold(
                    repo_root,
                    tmpdir,
                    ct,
                    max_iter=args.max_iter,
                    sym_enabled=False,
                )
                results_nosym.append(res)

            if args.keep_yaml:
                persist = repo_root / "validation" / "reference" / "csv" / "pcrs_yaml_cache"
                persist.mkdir(exist_ok=True)
                for y in tmpdir.glob("*.yaml"):
                    shutil.copy2(y, persist / y.name)
                print(f"  YAMLs saved to {persist}")

            # Summary table.
            def print_table(label: str, results: list[dict]) -> None:
                print(f"\n=== PCRS Si residual scan — {label} ===")
                print(
                    f"  {'conv_thr':>10}  {'iters':>5}  {'final_Δρ':>10}  "
                    f"{'E_KS (eV)':>14}  {'E_sum (eV)':>14}  "
                    f"{'|Δ| (eV)':>10}  notes"
                )
                print("  " + "-" * 92)
                for r in results:
                    ct = r.get("conv_threshold", float("nan"))
                    it = r.get("iters") or -1
                    drho = r.get("final_drho")
                    drho_str = f"{drho:.2e}" if drho is not None else "?"
                    e_ks = r.get("e_ks")
                    e_sum = r.get("e_sum")
                    delta = r.get("delta")
                    e_ks_s = f"{e_ks:.6f}" if e_ks is not None else "?"
                    e_sum_s = f"{e_sum:.6f}" if e_sum is not None else "?"
                    delta_s = f"{abs(delta):.2e}" if delta is not None else "?"
                    note = r.get("error", "")
                    print(
                        f"  {ct:>10.0e}  {it:>5d}  {drho_str:>10}  {e_ks_s:>14}  {e_sum_s:>14}  {delta_s:>10}  {note}"
                    )

            print_table("symmetry=ON  (48 Fd-3m ops; density symmetrized)", results_sym)
            print_table("symmetry=OFF (identity only; no density symm)", results_nosym)

            # Verdict.
            valid_sym = [r for r in results_sym if r.get("delta") is not None]
            valid_nosym = [r for r in results_nosym if r.get("delta") is not None]
            print("\n=== PCRS Verdict ===")
            if valid_sym:
                d_min = min(abs(r["delta"]) for r in valid_sym)
                d_max = max(abs(r["delta"]) for r in valid_sym)
                print(f"  symmetry=ON  residual range: [{d_min:.3e}, {d_max:.3e}] eV")
                ratio = d_min / max(d_max, 1e-30)
                if ratio > 0.99:
                    print("  -> PLATEAU (residual is CONSTANT vs conv_threshold): structural, not SCF noise.")
                else:
                    print(f"  -> residual ratio (min/max) = {ratio:.3e}")
            if valid_nosym:
                d_min = min(abs(r["delta"]) for r in valid_nosym)
                d_max = max(abs(r["delta"]) for r in valid_nosym)
                print(f"  symmetry=OFF residual range: [{d_min:.3e}, {d_max:.3e}] eV")
            if valid_sym and valid_nosym:
                d_sym = sum(abs(r["delta"]) for r in valid_sym) / len(valid_sym)
                d_nosym = sum(abs(r["delta"]) for r in valid_nosym) / len(valid_nosym)
                print(f"\n  mean residual, sym ON:  {d_sym:.3e} eV")
                print(f"  mean residual, sym OFF: {d_nosym:.3e} eV")
                if d_sym > 100 * d_nosym:
                    print(
                        f"  -> Turning symmetry OFF reduces residual by "
                        f"{d_sym / d_nosym:.1e}x; bug is in the SYMMETRY path."
                    )
                    print("     Root cause hypothesis: density symmetrization on an FFT grid")
                    print("     incompatible with the fractional translations (Si Fd-3m has τ=1/4,")
                    print("     grid 18 is not divisible by 4). `check_grid_compatibility` only")
                    print("     validates rotations, not translations.")

        return 0
    finally:
        if acquired:
            subprocess.call([str(lock), "release"])
            print("  machine lock released")


if __name__ == "__main__":
    sys.exit(main())
