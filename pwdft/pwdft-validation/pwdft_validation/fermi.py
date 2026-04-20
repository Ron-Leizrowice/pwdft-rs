"""Fermi-level reference finder (VGCH-2 Part C).

Reproduces QE's ``efermig.f90`` bisection for Fermi-Dirac, Gaussian, and
Methfessel-Paxton order-1 smearing. Validates pwdft-rs' ``find_fermi_energy``.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from pwdft_validation.units import RY_TO_EV


@dataclass
class QeBands:
    weights: list[float]
    eigenvalues: list[list[float]]
    degauss_ry: float
    reported_ef_ev: float
    n_electrons: float


def _parse_qe_bands(path: Path, nspin: int = 1) -> QeBands:
    text = path.read_text()
    wk_pat = re.compile(r"k\(\s*(\d+)\s*\)\s*=\s*\([^)]+\),\s*wk\s*=\s*([-+]?[\d.Ee]+)")
    wk: dict[int, float] = {}
    for m in wk_pat.finditer(text):
        wk.setdefault(int(m.group(1)), float(m.group(2)))
    nk = max(wk)
    weights = [wk[i + 1] for i in range(nk)]
    if nspin == 2:
        weights = weights + weights
        nk *= 2

    end_idx = text.rfind("End of self-consistent calculation")
    if end_idx < 0:
        raise SystemExit(f"{path}: no converged-bands block found")
    block = text[end_idx:]
    eigenvalues: list[list[float]] = []
    for part in re.split(r"bands \(ev\):\s*\n", block)[1:]:
        stop = re.search(r"(?:\n\s*(?:k =|the Fermi|highest|occupation numbers|!\s*total)|\n\s*\n\s*\n)", part)
        payload = part[: stop.start()] if stop else part
        nums = re.findall(r"[-+]?\d+\.\d+", payload)
        if nums:
            eigenvalues.append([float(x) for x in nums])
    if len(eigenvalues) != nk:
        eigenvalues = eigenvalues[-nk:]
    if len(eigenvalues) != nk:
        raise SystemExit(f"{path}: parsed {len(eigenvalues)} bands blocks, expected {nk}")

    m = re.search(r"smearing, width \(Ry\)=\s*([-+]?[\d.Ee]+)", text)
    if not m:
        raise SystemExit(f"{path}: degauss not found")
    degauss_ry = float(m.group(1))

    m = re.search(r"the Fermi energy is\s*([-+]?[\d.Ee]+)\s*ev", text)
    if not m:
        raise SystemExit(f"{path}: Fermi energy not found")
    reported_ef_ev = float(m.group(1))

    m = re.search(r"number of electrons\s*=\s*([-+]?[\d.Ee]+)", text)
    if not m:
        raise SystemExit(f"{path}: number of electrons not found")

    return QeBands(
        weights=weights,
        eigenvalues=eigenvalues,
        degauss_ry=degauss_ry,
        reported_ef_ev=reported_ef_ev,
        n_electrons=float(m.group(1)),
    )


def _fd_occ(e: float, ef: float, sigma: float) -> float:
    if sigma < 1e-15:
        return 1.0 if e < ef else (0.5 if abs(e - ef) < 1e-12 else 0.0)
    x = (e - ef) / sigma
    if x > 40.0:
        return 0.0
    if x < -40.0:
        return 1.0
    return 1.0 / (1.0 + math.exp(x))


def _gauss_occ(e: float, ef: float, sigma: float) -> float:
    if sigma < 1e-15:
        return 1.0 if e < ef else 0.0
    return 0.5 * math.erfc((e - ef) / sigma)


def _mp1_occ(e: float, ef: float, sigma: float) -> float:
    if sigma < 1e-15:
        return 1.0 if e < ef else 0.0
    x = (e - ef) / sigma
    return 0.5 * math.erfc(x) - 0.5 * x * math.exp(-x * x) / math.sqrt(math.pi)


def _total_n(
    ef: float,
    bands: QeBands,
    sigma: float,
    occ_fn: Callable,
) -> float:
    return sum(w * occ_fn(e, ef, sigma) for w, evs in zip(bands.weights, bands.eigenvalues, strict=True) for e in evs)


def _find_ef_pwdft(
    bands: QeBands,
    sigma: float,
    occ_fn: Callable,
) -> float:
    e_min = min(min(e) for e in bands.eigenvalues) - 10.0 * max(sigma, 0.1)
    e_max = max(max(e) for e in bands.eigenvalues) + 10.0 * max(sigma, 0.1)
    for _ in range(200):
        e_mid = 0.5 * (e_min + e_max)
        if _total_n(e_mid, bands, sigma, occ_fn) < bands.n_electrons:
            e_min = e_mid
        else:
            e_max = e_mid
        if abs(e_max - e_min) < 1e-14:
            break
    return 0.5 * (e_min + e_max)


def _find_ef_qe(
    bands: QeBands,
    sigma: float,
    occ_fn: Callable,
) -> float:
    e_min = min(min(e) for e in bands.eigenvalues) - 10.0 * sigma
    e_max = max(max(e) for e in bands.eigenvalues) + 10.0 * sigma
    for _ in range(300):
        e_mid = 0.5 * (e_min + e_max)
        n = _total_n(e_mid, bands, sigma, occ_fn)
        if abs(n - bands.n_electrons) < 1e-10:
            return e_mid
        if n < bands.n_electrons:
            e_min = e_mid
        else:
            e_max = e_mid
    return 0.5 * (e_min + e_max)


_SYSTEMS = [
    ("Cu_FCC", "cu_fcc_scf.out", 1),
    ("Fe_BCC_FM", "fe_bcc_fm_scf.out", 2),
    ("NaCl", "nacl_scf.out", 1),
    ("C_diamond", "c_diamond_scf.out", 1),
    ("MgO", "mgo_scf.out", 1),
]

_OCC_FNS = [
    ("fermi_dirac", _fd_occ),
    ("gaussian", _gauss_occ),
    ("mp1", _mp1_occ),
]

_FINDERS = [
    ("pwdft_bracket", _find_ef_pwdft),
    ("qe_count", _find_ef_qe),
]


def compute(qe_dir: Path, out_csv: Path) -> int:
    """Generate ``vgch2c_fermi_reference.csv`` (Fermi-level bisection reference)."""
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "system",
                "smearing",
                "bisection_style",
                "ef_ev",
                "qe_reported_ef_ev",
                "ef_minus_qe_ev",
                "n_at_ef",
                "n_electrons_target",
            ]
        )
        for name, fname, nspin in _SYSTEMS:
            path = qe_dir / fname
            if not path.exists():
                print(f"# WARNING: {path} missing, skipping {name}", file=sys.stderr)
                continue
            bands = _parse_qe_bands(path, nspin=nspin)
            sigma_ev = bands.degauss_ry * RY_TO_EV
            for fn_name, occ_fn in _OCC_FNS:
                for style_name, finder in _FINDERS:
                    ef = finder(bands, sigma_ev, occ_fn)
                    n = _total_n(ef, bands, sigma_ev, occ_fn)
                    w.writerow(
                        [
                            name,
                            fn_name,
                            style_name,
                            f"{ef:.10f}",
                            f"{bands.reported_ef_ev:.6f}",
                            f"{ef - bands.reported_ef_ev:+.6f}",
                            f"{n:.10f}",
                            f"{bands.n_electrons:.6f}",
                        ]
                    )
    print(f"Wrote {out_csv}")
    return 0
