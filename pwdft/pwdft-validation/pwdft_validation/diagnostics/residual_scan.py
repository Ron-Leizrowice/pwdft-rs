"""PCRS — Per-Component Energy Residual Investigation.

Sweeps SCF ``conv_threshold`` on Si and checks whether the E_sum − E_total
residual shrinks (SCF noise) or plateaus (structural bookkeeping bug).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import cyclopts

from pwdft_validation.paths import PROJECT_ROOT, PSEUDO_DIR, QE_REF_DIR
from pwdft_validation.units import RY_TO_EV

_SI_YAML = """\
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
  smearing_width: 0.13605693122994
  mixing_mode: plain

symmetry:
  enabled: {sym_enabled}
  time_reversal: {sym_enabled}

pseudopotentials:
  Si: "{pp_path}"
"""

_PC_PAT = re.compile(r"E_(band|kinetic|local|nonlocal|hartree|xc|ewald)\s*=\s*([-+\d.Ee]+)")
_E_SUM_PAT = re.compile(r"E_sum\(comp\)\s*=\s*([-+\d.Ee]+)\s+\(vs E_KS\s+([-+\d.Ee]+),\s*Δ=([-+\d.Ee]+)\)")
_E_LOCAL_G0_PAT = re.compile(r"E_local\(G=0\)\s*=\s*([-+\d.Ee]+)")
_ITERS_PAT = re.compile(r"SCF converged after (\d+) iterations")
_FINAL_DRHO_PAT = re.compile(r"SCF iter\s+\d+:.*?Δρ=([-+\d.Ee]+)\s*$", flags=re.MULTILINE)

_QE_PATS = {
    "total_energy": r"!\s*total energy\s*=\s*([-+\d.Ee]+)\s*Ry",
    "one_electron": r"one-electron contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "hartree": r"hartree contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "xc": r"xc contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "ewald": r"ewald contribution\s*=\s*([-+\d.Ee]+)\s*Ry",
    "smearing_ts": r"smearing contrib\. \(-TS\)\s*=\s*([-+\d.Ee]+)\s*Ry",
    "internal_E": r"internal energy E=F\+TS\s*=\s*([-+\d.Ee]+)\s*Ry",
}


def _parse_cargo_stderr(text: str) -> dict:
    out: dict = {
        k: None
        for k in [
            "e_band",
            "e_kinetic",
            "e_local",
            "e_local_g0_shift",
            "e_nonlocal",
            "e_hartree",
            "e_xc",
            "e_ewald",
            "e_sum",
            "e_ks",
            "delta",
            "iters",
            "final_drho",
        ]
    }
    g0 = _E_LOCAL_G0_PAT.search(text)
    if g0:
        out["e_local_g0_shift"] = float(g0.group(1))
    text_no_g0 = _E_LOCAL_G0_PAT.sub("", text)
    for m in _PC_PAT.finditer(text_no_g0):
        out[f"e_{m.group(1)}"] = float(m.group(2))
    s = _E_SUM_PAT.search(text)
    if s:
        out["e_sum"] = float(s.group(1))
        out["e_ks"] = float(s.group(2))
        out["delta"] = float(s.group(3))
    it = _ITERS_PAT.search(text)
    if it:
        out["iters"] = int(it.group(1))
    drhos = _FINAL_DRHO_PAT.findall(text)
    if drhos:
        out["final_drho"] = float(drhos[-1])
    return out


def _run_si(repo_root: Path, tmpdir: Path, conv_threshold: float, max_iter: int, sym_enabled: bool) -> dict:
    tag = "sym" if sym_enabled else "nosym"
    yaml_path = tmpdir / f"si_pcrs_{tag}_{conv_threshold:.0e}.yaml"
    pp_path = PSEUDO_DIR / "nc" / "lda" / "Si.upf"
    energy_threshold = max(conv_threshold * 10.0, 1e-10)
    yaml_path.write_text(
        _SI_YAML.format(
            conv_threshold=conv_threshold,
            energy_threshold=energy_threshold,
            max_iter=max_iter,
            sym_enabled="true" if sym_enabled else "false",
            pp_path=str(pp_path),
        )
    )
    env = {**os.environ, "RUST_LOG": "info", "RUST_BACKTRACE": "1"}
    print(f"  running [{tag}] conv_threshold={conv_threshold:.0e} ... ", end="", flush=True)
    try:
        result = subprocess.run(
            ["cargo", "run", "--release", "--quiet", "--", "--input", str(yaml_path)],
            cwd=repo_root,
            env=env,
            capture_output=True,
            text=True,
            timeout=600,
        )
    except subprocess.TimeoutExpired:
        print("TIMEOUT")
        return {"error": "timeout", "conv_threshold": conv_threshold}
    parsed = _parse_cargo_stderr(result.stderr)
    parsed["conv_threshold"] = conv_threshold
    parsed["returncode"] = result.returncode
    if result.returncode != 0:
        print(f"FAILED (rc={result.returncode})")
        parsed["error"] = "non-zero returncode"
        (tmpdir / f"si_pcrs_{conv_threshold:.0e}.stderr.log").write_text(result.stderr)
    elif parsed.get("iters") is not None:
        delta = parsed.get("delta")
        print(f"OK ({parsed['iters']} iters, Δ={delta:+.3e} eV)" if delta is not None else "OK")
    else:
        print("OK (parse failed)")
        parsed["error"] = "parse failed"
    return parsed


def _qe_identity(qe_out: Path) -> None:
    text = qe_out.read_text()
    vals: dict[str, float] = {}
    for key, pat in _QE_PATS.items():
        m = re.search(pat, text)
        if m:
            vals[key] = float(m.group(1))
    print(f"\n=== QE identity closure — {qe_out.name} ===")
    for k in ["one_electron", "hartree", "xc", "ewald", "total_energy", "smearing_ts", "internal_E"]:
        if k in vals:
            print(f"  {k:<16} = {vals[k]:>16.10f} Ry  ({vals[k] * RY_TO_EV:>12.6f} eV)")
    s = sum(vals.get(k, 0.0) for k in ["one_electron", "hartree", "xc", "ewald"])
    diff = s - vals.get("internal_E", s)
    print(f"  Sigma(4) − internal_E = {diff:+.3e} Ry  ({diff * RY_TO_EV:+.3e} eV)")
    if abs(diff * RY_TO_EV) < 1e-6:
        print("  => QE identity closes to < 1 μeV. Our residual is STRUCTURAL.")


def _print_table(label: str, results: list[dict]) -> None:
    print(f"\n=== PCRS Si residual scan — {label} ===")
    print(
        f"  {'conv_thr':>10}  {'iters':>5}  {'final_Δρ':>10}  {'E_KS (eV)':>14}  {'E_sum (eV)':>14}  {'|Δ| (eV)':>10}  notes"
    )
    print("  " + "-" * 92)
    for r in results:
        ct = r.get("conv_threshold", float("nan"))
        it = r.get("iters") or -1
        drho = r.get("final_drho")
        e_ks = r.get("e_ks")
        e_sum = r.get("e_sum")
        delta = r.get("delta")
        drho_s = f"{drho:.2e}" if drho is not None else "?"
        e_ks_s = f"{e_ks:.6f}" if e_ks is not None else "?"
        e_sum_s = f"{e_sum:.6f}" if e_sum is not None else "?"
        delta_s = f"{abs(delta):.2e}" if delta is not None else "?"
        note = r.get("error", "")
        print(f"  {ct:>10.0e}  {it:>5d}  {drho_s:>10}  {e_ks_s:>14}  {e_sum_s:>14}  {delta_s:>10}  {note}")


def run(
    repo_root: Path | None = None,
    thresholds: list[float] | None = None,
    no_qe: bool = False,
    keep_yaml: bool = False,
    max_iter: int = 200,
) -> int:
    """Run the PCRS threshold sweep."""
    root = repo_root or PROJECT_ROOT
    if not (root / "Cargo.toml").is_file():
        print(f"ERROR: {root} is not a pwdft-rs repo", file=sys.stderr)
        return 1

    thresholds = thresholds or [1e-6, 1e-8, 1e-9, 1e-10, 1e-11]
    print(f"PCRS residual scan; repo = {root}")
    print(f"  thresholds: {thresholds}")

    if not no_qe:
        qe_si = QE_REF_DIR / "si_scf.out"
        if qe_si.is_file():
            _qe_identity(qe_si)

    lock = root / ".claude" / "bin" / "machine-lock"
    acquired = False
    if lock.is_file():
        acquired = subprocess.call([str(lock), "acquire", "Researcher", "PCRS scan"]) == 0
        print(f"  machine lock {'acquired' if acquired else 'not acquired (continuing)'}")

    try:
        with tempfile.TemporaryDirectory(prefix="pcrs_scan_") as td:
            tmpdir = Path(td)
            print(f"\nPCRS — Si at {len(thresholds)} thresholds, symmetry=ON:")
            results_sym = [_run_si(root, tmpdir, ct, max_iter, True) for ct in thresholds]
            print(f"\nPCRS — Si at {len(thresholds)} thresholds, symmetry=OFF:")
            results_nosym = [_run_si(root, tmpdir, ct, max_iter, False) for ct in thresholds]
            if keep_yaml:
                persist = root / "validation" / "reference" / "csv" / "pcrs_yaml_cache"
                persist.mkdir(exist_ok=True)
                for y in tmpdir.glob("*.yaml"):
                    shutil.copy2(y, persist / y.name)
            _print_table("symmetry=ON", results_sym)
            _print_table("symmetry=OFF", results_nosym)
        return 0
    finally:
        if acquired:
            subprocess.call([str(lock), "release"])
            print("  machine lock released")


standalone_app = cyclopts.App(
    name="pwdft-residual-scan",
    help="PCRS — per-component energy residual investigation (Si SCF threshold sweep).",
)


@standalone_app.default
def _residual_scan_cmd(
    *,
    thresholds: list[float] | None = None,
    no_qe: bool = False,
    keep_yaml: bool = False,
    max_iter: int = 200,
) -> None:
    """Sweep SCF conv_threshold on Si and report per-component residuals."""
    raise SystemExit(run(thresholds=thresholds, no_qe=no_qe, keep_yaml=keep_yaml, max_iter=max_iter))


def _standalone_run() -> None:
    """Console-script entry point for ``pwdft-residual-scan``."""
    standalone_app()
