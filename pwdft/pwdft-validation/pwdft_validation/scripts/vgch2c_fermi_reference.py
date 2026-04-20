#!/usr/bin/env python3
"""VGCH-2 Part C — Python Fermi-Dirac reference Fermi-level finder.

Parses QE 7.5's final `bands (ev)` and `wk = ...` blocks from a converged
pw.x stdout, then reproduces QE's Fermi-level root-find (`efermig.f90`)
using:

 1. Fermi-Dirac (ngauss = -99): f(x) = 1/(1+exp(x)),  x = (ε−E_F)/σ.
 2. Gaussian (ngauss = 0): f(x) = 1/2 erfc(x).
 3. Methfessel-Paxton order-1 (ngauss = 1):
    f(x) = 1/2 erfc(x) − (x/(2√π)) exp(−x²).

Bisection on Σ_ik w_k · spin_factor · f_ik = N_el.

For the Cu FCC reference:
 - Reads `validation/reference/qe/cu_fcc_scf.out`.
 - nspin = 1 → spin_factor = 2; wk-sum = 2 per QE convention.
 - Cu N_el = 19.

Outputs a CSV:
    system, smearing, ef_ev, qe_reported_ef_ev, ef_minus_qe_ev, n_at_ef

This script isolates H-C1 (bisection / convergence criterion) and H-C3
(smearing-function mismatch). If pwdft-rs matches the Fermi-Dirac column
to μeV, pwdft-rs' `find_fermi_energy` is numerically correct. Differences
between the three smearing columns cap the maximum effect of H-C3.

Usage:
    uv run validation/src/pwdft_validation/scripts/vgch2c_fermi_reference.py \\
        > validation/reference/csv/vgch2c_fermi_reference.csv

Source cross-refs:
    qe-7.5/PW/src/efermig.f90 (bisection + Newton)
    qe-7.5/Modules/wgauss.f90 (occupation formulae)
    qe-7.5/PW/src/sumkg.f90  (wgauss((e-et)/degauss, ngauss))
"""

from __future__ import annotations

import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

RY_TO_EV = 13.605_693_122_994


@dataclass
class QeBands:
    # [ik] -> weight (kpt weight, summed over k gives 2 for nspin=1)
    weights: list[float]
    # [ik][ib] -> eigenvalue in eV
    eigenvalues: list[list[float]]
    # degauss in Ry (as QE prints)
    degauss_ry: float
    # QE-reported Fermi level in eV
    reported_ef_ev: float
    # N_el
    n_electrons: float


def parse_qe_output(path: Path, nspin: int = 1) -> QeBands:
    text = path.read_text()
    # K-point weights: `k(  N) = ( ... ), wk =  W` — there are typically
    # two occurrences in QE (coord-system variants). Take the first
    # sequence; skip duplicates where the same k(N) reappears.
    wk_pat = re.compile(r"k\(\s*(\d+)\s*\)\s*=\s*\([^)]+\),\s*wk\s*=\s*([-+]?[\d.Ee]+)")
    wk: dict[int, float] = {}
    for m in wk_pat.finditer(text):
        idx = int(m.group(1))
        w = float(m.group(2))
        # Only take the first occurrence per k-index (QE prints twice in
        # cart + crystal coord blocks).
        wk.setdefault(idx, w)
    # Indices are 1-based and dense 1..nk.
    nk = max(wk.keys())
    weights = [wk[i + 1] for i in range(nk)]

    # For nspin=2 QE duplicates the k-list: each k appears twice in
    # the final `bands (ev)` block (once per spin channel). The weights
    # QE reports are per-channel and do NOT include the spin-degeneracy
    # factor; sumkg integrates over both channels. Expand the weight
    # list to cover both spin blocks with identical wk.
    if nspin == 2:
        weights = weights + weights
        nk = nk * 2

    # Find the final "End of self-consistent calculation" block, then
    # read eigenvalues from the `bands (ev)` entries after it.
    end_idx = text.rfind("End of self-consistent calculation")
    if end_idx < 0:
        raise SystemExit(f"{path}: no converged-bands block found")
    block = text[end_idx:]

    # Parse per-k-point `bands (ev):` blocks. Each is followed by one or
    # more lines of floats; terminate on a blank line followed by a
    # non-float line OR another `k =` / Fermi / end-marker sentinel.
    eigenvalues: list[list[float]] = []
    # Split on 'bands (ev):' — entries alternate header / payload.
    parts = re.split(r"bands \(ev\):\s*\n", block)
    # parts[0] is pre-first; parts[1..] contain the payload + trailing.
    for part in parts[1:]:
        # Grab up to the next "k =" or Fermi or blank-line terminator.
        stop = re.search(
            r"(?:\n\s*(?:k =|the Fermi|highest|occupation numbers|!\s*total)|\n\s*\n\s*\n)",
            part,
        )
        payload = part[: stop.start()] if stop else part
        nums = re.findall(r"[-+]?\d+\.\d+", payload)
        if not nums:
            continue
        eigenvalues.append([float(x) for x in nums])

    if len(eigenvalues) != nk:
        # Sometimes QE prints extra blocks; trim to last nk.
        eigenvalues = eigenvalues[-nk:]
    if len(eigenvalues) != nk:
        raise SystemExit(f"{path}: parsed {len(eigenvalues)} bands blocks, expected {nk}")

    # Sanity: all have the same nbnd.
    nbnds = {len(e) for e in eigenvalues}
    if len(nbnds) != 1:
        raise SystemExit(f"{path}: inconsistent nbnd across k-points: {sorted(nbnds)}")

    # Extract degauss (Ry).
    m = re.search(r"smearing, width \(Ry\)=\s*([-+]?[\d.Ee]+)", text)
    if not m:
        raise SystemExit(f"{path}: degauss not found")
    degauss_ry = float(m.group(1))

    # Extract reported Fermi level (eV).
    m = re.search(r"the Fermi energy is\s*([-+]?[\d.Ee]+)\s*ev", text)
    if not m:
        raise SystemExit(f"{path}: Fermi energy not found")
    reported_ef_ev = float(m.group(1))

    # Extract N_el.
    m = re.search(r"number of electrons\s*=\s*([-+]?[\d.Ee]+)", text)
    if not m:
        raise SystemExit(f"{path}: number of electrons not found")
    n_electrons = float(m.group(1))

    return QeBands(
        weights=weights,
        eigenvalues=eigenvalues,
        degauss_ry=degauss_ry,
        reported_ef_ev=reported_ef_ev,
        n_electrons=n_electrons,
    )


def fermi_dirac_occ(e_ev: float, ef_ev: float, sigma_ev: float) -> float:
    """f(ε) = 1 / (1 + exp((ε − E_F)/σ)).

    Matches both pwdft-rs (smearing.rs:122-130) and QE
    (wgauss.f90:48-56 with x = (e − et)/degauss: wgauss = 1/(1+exp(−x))
    ⇒ 1/(1+exp((et − e)/degauss)) — equivalent when the integrand is
    Σ_b wgauss((E_F − et_b)/degauss) evaluated at the unknown E_F root).
    """
    if sigma_ev < 1e-15:
        return 1.0 if e_ev < ef_ev else (0.5 if abs(e_ev - ef_ev) < 1e-12 else 0.0)
    x = (e_ev - ef_ev) / sigma_ev
    if x > 40.0:
        return 0.0
    if x < -40.0:
        return 1.0
    return 1.0 / (1.0 + math.exp(x))


def gaussian_occ(e_ev: float, ef_ev: float, sigma_ev: float) -> float:
    if sigma_ev < 1e-15:
        return 1.0 if e_ev < ef_ev else 0.0
    x = (e_ev - ef_ev) / sigma_ev
    return 0.5 * math.erfc(x)


def mp1_occ(e_ev: float, ef_ev: float, sigma_ev: float) -> float:
    """Methfessel-Paxton order-1 occupation in [0, 1]."""
    if sigma_ev < 1e-15:
        return 1.0 if e_ev < ef_ev else 0.0
    x = (e_ev - ef_ev) / sigma_ev
    f0 = 0.5 * math.erfc(x)
    return f0 - 0.5 * x * math.exp(-x * x) / math.sqrt(math.pi)


def total_n_electrons(
    ef_ev: float,
    bands: QeBands,
    sigma_ev: float,
    occ_fn,
    spin_factor: float,
) -> float:
    """QE convention: wk already includes the spin-degeneracy factor
    (setup.f90:673, `wk *= degspin` for nspin=1). Σ wk = 2 for nspin=1.
    So `spin_factor` here is ignored against the as-printed QE weights;
    the `spin_factor` argument is retained for parity with pwdft-rs'
    `find_fermi_energy` API, which uses `spin_factor=2` on wk sums that
    equal 1.0 instead.
    """
    del spin_factor  # consumed implicitly via wk-sum = degspin.
    n = 0.0
    for w, evs in zip(bands.weights, bands.eigenvalues, strict=True):
        for e in evs:
            n += w * occ_fn(e, ef_ev, sigma_ev)
    return n


def find_fermi_pwdft(bands: QeBands, sigma_ev: float, occ_fn, spin_factor: float) -> float:
    """Bisection in pwdft-rs' style: convergence on (e_max − e_min) < 1e-14 eV."""
    e_min = min(min(e) for e in bands.eigenvalues) - 10.0 * max(sigma_ev, 0.1)
    e_max = max(max(e) for e in bands.eigenvalues) + 10.0 * max(sigma_ev, 0.1)
    target = bands.n_electrons
    for _ in range(200):
        e_mid = 0.5 * (e_min + e_max)
        n = total_n_electrons(e_mid, bands, sigma_ev, occ_fn, spin_factor)
        if n < target:
            e_min = e_mid
        else:
            e_max = e_mid
        if abs(e_max - e_min) < 1e-14:
            break
    return 0.5 * (e_min + e_max)


def find_fermi_qe(bands: QeBands, sigma_ev: float, occ_fn, spin_factor: float) -> float:
    """Bisection in QE's style: convergence on |N(ef) − N_target| < 1e-10 electrons.

    Matches `efermig.f90:47-50,289-307`. Bracket shrinks on sign of the
    electron count residual; tolerance is on the count, not the bracket.
    """
    e_min = min(min(e) for e in bands.eigenvalues) - 10.0 * sigma_ev
    e_max = max(max(e) for e in bands.eigenvalues) + 10.0 * sigma_ev
    target = bands.n_electrons
    tol_n = 1.0e-10
    for _ in range(300):
        e_mid = 0.5 * (e_min + e_max)
        n = total_n_electrons(e_mid, bands, sigma_ev, occ_fn, spin_factor)
        if abs(n - target) < tol_n:
            return e_mid
        if n < target:
            e_min = e_mid
        else:
            e_max = e_mid
    return 0.5 * (e_min + e_max)


def main() -> None:
    root = Path(__file__).resolve().parents[4]
    # Default Cu FCC path. The test could expand to more systems later.
    systems = [
        ("Cu_FCC", root / "validation" / "reference" / "qe" / "cu_fcc_scf.out", 1),
        # Tier-C triangulation: Fe, NaCl, C, MgO. Fe is nspin=2 (FM) in
        # the deck; we want to check the finder, which uses degspin=1
        # weights summing to 1.0 in QE's nspin=2 convention.
        ("Fe_BCC_FM", root / "validation" / "reference" / "qe" / "fe_bcc_fm_scf.out", 2),
        ("NaCl", root / "validation" / "reference" / "qe" / "nacl_scf.out", 1),
        ("C_diamond", root / "validation" / "reference" / "qe" / "c_diamond_scf.out", 1),
        ("MgO", root / "validation" / "reference" / "qe" / "mgo_scf.out", 1),
    ]

    writer = csv.writer(sys.stdout)
    writer.writerow(
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

    for name, path, nspin in systems:
        if not path.exists():
            print(f"# WARNING: {path} missing, skipping {name}", file=sys.stderr)
            continue
        bands = parse_qe_output(path, nspin=nspin)
        sigma_ev = bands.degauss_ry * RY_TO_EV
        spin_factor = 2.0 if nspin == 1 else 1.0

        for fn_name, occ_fn in [
            ("fermi_dirac", fermi_dirac_occ),
            ("gaussian", gaussian_occ),
            ("mp1", mp1_occ),
        ]:
            for style, finder in [
                ("pwdft_bracket", find_fermi_pwdft),
                ("qe_count", find_fermi_qe),
            ]:
                ef = finder(bands, sigma_ev, occ_fn, spin_factor)
                n = total_n_electrons(ef, bands, sigma_ev, occ_fn, spin_factor)
                writer.writerow(
                    [
                        name,
                        fn_name,
                        style,
                        f"{ef:.10f}",
                        f"{bands.reported_ef_ev:.6f}",
                        f"{ef - bands.reported_ef_ev:+.6f}",
                        f"{n:.10f}",
                        f"{bands.n_electrons:.6f}",
                    ]
                )


if __name__ == "__main__":
    main()
