#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy>=1.26",
#     "scipy>=1.11",
# ]
# ///
"""
VGCH Phase 1b Hypothesis 1 — β_l(q) projector reference for heavy-atom PPs.

For every β projector in every VGCH-scope UPF, compute the Bessel-transformed
form factor using QE's exact integration recipe
(`qe-7.5/upflib/beta_mod.f90:111-116`):

    F_l(q) = 4π · Simpson_{i}( χ_i · j_l(q·r_i) · r_i ; rab_i )

where `χ(r) = r·β(r)` is the on-disk PP_BETA array (Bohr^{-1/2}) and
`rab_i = r_i · log_step` is the log-mesh Jacobian. Output is in native QE
units: q in Bohr^{-1}, F_l(q) in Bohr^{3/2}. A second column gives the
same F_l(q) in pwdft-rs' internal Å convention (Å^{3/2}) after the
`1/√BOHR_TO_ANG` conversion applied in
`src/pseudopotential/upf/convert.rs::parse_body`.

This mirrors the already-landed `validation/src/pwdft_validation/scripts/beta_q_reference.py`
(which covers only Si) but extends to the 9 VGCH-scope elements and uses
the QE-style `Σ c_i·f_i·rab_i` quadrature — byte-identical to both QE's
`qe-7.5/upflib/simpsn.f90` and pwdft-rs' `src/numerics.rs::simpson_integrate`.

The resulting CSV at `validation/reference/csv/vgch_beta_l_heavy.csv` is then
consumed by a Rust integration test (`tests/vgch_beta_l_heavy.rs`) that
pins `NonlocalPotential`'s form-factor output against this reference to
1e-8 Bohr^{3/2} (QE's printed `tab_beta` precision).

Hypothesis-1 early verdict: if C diamond (l=0 + l=1 only, no semicore)
matches bit-perfect and its 1.45 eV E_total residual is NOT explained,
then Hypothesis 1 is not the whole story. C is the cheap case — both
coded first in the loop and reported first in stdout.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.special import spherical_jn

BOHR_TO_ANG = 0.529_177_210_903

# q-grid sampling — covers ecut=15 Ry (qmax ≈ 3.3 Bohr⁻¹) through
# ecut=50 Ry (qmax ≈ 6.3 Bohr⁻¹) plus headroom. Three values per decade.
Q_VALUES_BOHR = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

# Elements in VGCH scope: the 5 originally-ignored heavy-atom cells plus
# C (which VQEF-QC #144 showed exhibits the same opposite-sign signature).
# Si/Al are guardrails: Si passes at <40 meV, Al's gap is basis truncation.
# Both should match β_q independently and are included as cross-checks.
ELEMENTS = ["Si", "C", "Al", "Fe", "Cu", "Ga", "As", "Na", "Cl", "Mg", "O"]


# ---------------------------------------------------------------------------
# UPF v2 parser (regex-based, no external deps — matches existing scripts)
# ---------------------------------------------------------------------------


def _extract_attr(text: str, name: str) -> str:
    m = re.search(rf'{re.escape(name)}\s*=\s*"([^"]*)"', text)
    if m is None:
        raise ValueError(f"attribute {name!r} not found")
    return m.group(1).strip()


def _extract_block(text: str, tag: str) -> np.ndarray:
    m = re.search(
        rf"<{re.escape(tag)}\b[^>]*>(.*?)</{re.escape(tag)}>",
        text,
        re.DOTALL,
    )
    if m is None:
        raise ValueError(f"tag <{tag}> not found")
    body = m.group(1)
    return np.asarray(
        [float(t) for t in body.split() if t.strip()],
        dtype=np.float64,
    )


def _extract_beta_projector(text: str, idx: int, mesh: int) -> tuple[int, np.ndarray]:
    tag = f"PP_BETA.{idx}"
    m = re.search(
        rf"<{re.escape(tag)}\b([^>]*)>(.*?)</{re.escape(tag)}>",
        text,
        re.DOTALL,
    )
    if m is None:
        raise ValueError(f"tag <{tag}> not found")
    attrs = m.group(1)
    body = m.group(2)
    l_match = re.search(r'angular_momentum\s*=\s*"(-?\d+)"', attrs)
    if l_match is None:
        raise ValueError(f"{tag} has no angular_momentum attribute")
    l = int(l_match.group(1))
    values = np.asarray(
        [float(t) for t in body.split() if t.strip()],
        dtype=np.float64,
    )
    if values.size != mesh:
        raise ValueError(f"{tag}: expected {mesh} values, got {values.size}")
    return l, values


def parse_upf(path: Path) -> dict:
    text = path.read_text()
    m_hdr = re.search(r"<PP_HEADER\b([^>]*)/?>", text)
    if m_hdr is None:
        raise ValueError("PP_HEADER not found")
    hdr = m_hdr.group(1)
    z_val = float(_extract_attr(hdr, "z_valence"))
    mesh = int(_extract_attr(hdr, "mesh_size"))
    n_proj = int(_extract_attr(hdr, "number_of_proj"))

    r = _extract_block(text, "PP_R")
    rab = _extract_block(text, "PP_RAB")
    assert r.size == mesh and rab.size == mesh

    projectors = []
    for i in range(1, n_proj + 1):
        l, chi = _extract_beta_projector(text, i, mesh)
        projectors.append((l, chi))

    return {
        "z_valence": z_val,
        "mesh_size": mesh,
        "n_proj": n_proj,
        "r_bohr": r,
        "rab_bohr": rab,
        "projectors": projectors,
    }


# ---------------------------------------------------------------------------
# Simpson integrator — byte-identical to QE's `simpson(mesh, func, rab, asum)`
# and pwdft-rs's `simpson_integrate(func, rab)`.
# ---------------------------------------------------------------------------


def simpson_qe(func: np.ndarray, rab: np.ndarray) -> float:
    """QE-style Simpson: Σ c_i · func_i · rab_i with c_i ∈ {1/3, 4/3, 2/3}."""
    n = len(func)
    assert n == len(rab), "func and rab must have the same length"
    if n < 3:
        return float(np.sum(func * rab))

    # Interior: Fortran i = 2..mesh-1 → Python i = 1..n-2
    # Weight: 4 when (Fortran i) is even → (Python i) odd; 2 when Fortran odd.
    acc = 0.0
    for i in range(1, n - 1):
        weight = 4.0 if (i % 2 == 1) else 2.0
        acc += weight * func[i] * rab[i]

    if n % 2 == 1:
        return (acc + func[0] * rab[0] + func[-1] * rab[-1]) / 3.0
    # Even-mesh correction (DFTK/QE formula)
    acc += func[0] * rab[0]
    acc += -0.25 * func[n - 3] * rab[n - 3]
    acc += func[n - 2] * rab[n - 2]
    acc += 1.25 * func[n - 1] * rab[n - 1]
    return acc / 3.0


def f_l_of_q_bohr(
    r_bohr: np.ndarray,
    rab_bohr: np.ndarray,
    chi: np.ndarray,
    l: int,
    q_bohr_inv: float,
) -> float:
    """QE-convention F_l(q) in Bohr^{3/2}. `chi` = r·β(r) in Bohr^{-1/2}."""
    qr = q_bohr_inv * r_bohr
    jl = spherical_jn(l, qr)
    integrand = chi * jl * r_bohr
    integral = simpson_qe(integrand, rab_bohr)
    return 4.0 * math.pi * integral


def f_l_of_q_ang(
    r_bohr: np.ndarray,
    rab_bohr: np.ndarray,
    chi: np.ndarray,
    l: int,
    q_per_ang: float,
) -> float:
    """pwdft-rs-convention F_l(q) in Å^{3/2}.

    This is what `src/potential/nonlocal.rs::bessel_transform_projector`
    returns. The conversion from QE's Bohr units is:

        r_Å    = r_B · BOHR_TO_ANG
        rab_Å  = rab_B · BOHR_TO_ANG           (log-mesh Jacobian scales with r)
        chi_Å  = chi_B / √BOHR_TO_ANG           (Å^{-1/2} = Bohr^{-1/2} / √BOHR)
        q_Å⁻¹  = q_B⁻¹ / BOHR_TO_ANG

    Plugging in:
        F_Å(q_Å) = 4π · ∫ chi_Å · j_l(q_Å · r_Å) · r_Å · drab_Å
                 = 4π · (1/√BOHR) · BOHR² · ∫ chi_B · j_l(q_B · r_B) · r_B · drab_B
                 = BOHR^{3/2} · F_B(q_B)

    (the `j_l` argument q_Å·r_Å = q_B·r_B is dimensionless and invariant)
    """
    r_ang = r_bohr * BOHR_TO_ANG
    rab_ang = rab_bohr * BOHR_TO_ANG
    chi_ang = chi / math.sqrt(BOHR_TO_ANG)
    qr = q_per_ang * r_ang
    jl = spherical_jn(l, qr)
    integrand = chi_ang * jl * r_ang
    integral = simpson_qe(integrand, rab_ang)
    return 4.0 * math.pi * integral


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    from pwdft_validation import CSV_REF_DIR, PSEUDO_DIR

    pp_dir = PSEUDO_DIR / "nc" / "lda"
    out_csv = CSV_REF_DIR / "vgch_beta_l_heavy.csv"

    rows = []
    hdr = f"{'elem':>4} {'i':>3} {'l':>3} {'q (Bohr⁻¹)':>11} {'F_Bohr (Bohr^3/2)':>22} {'F_Å (Å^3/2)':>22}"
    print(hdr)
    print("-" * len(hdr))

    for elem in ELEMENTS:
        upf = pp_dir / f"{elem}.upf"
        if not upf.exists():
            print(f"# SKIP {elem}: {upf} missing", file=sys.stderr)
            continue
        pp = parse_upf(upf)
        print(f"# {elem}: Z_val={pp['z_valence']:.2f}, n_proj={pp['n_proj']}, ls={[l for (l, _) in pp['projectors']]}")
        for i, (l, chi) in enumerate(pp["projectors"]):
            for q_bohr in Q_VALUES_BOHR:
                f_bohr = f_l_of_q_bohr(pp["r_bohr"], pp["rab_bohr"], chi, l, q_bohr)
                # Convert q to 1/Å for the pwdft-convention integrand.
                q_ang = q_bohr / BOHR_TO_ANG
                f_ang = f_l_of_q_ang(pp["r_bohr"], pp["rab_bohr"], chi, l, q_ang)
                rows.append((elem, i, l, q_bohr, q_ang, f_bohr, f_ang))
                print(f"{elem:>4} {i:>3} {l:>3} {q_bohr:>11.6f} {f_bohr:>22.12e} {f_ang:>22.12e}")

    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "element",
                "proj_index",
                "l",
                "q_bohr_inv",
                "q_per_ang",
                "F_bohr_3halves",
                "F_ang_3halves",
            ]
        )
        for row in rows:
            w.writerow(
                [
                    row[0],
                    row[1],
                    row[2],
                    f"{row[3]:.12e}",
                    f"{row[4]:.12e}",
                    f"{row[5]:.12e}",
                    f"{row[6]:.12e}",
                ]
            )
    print(f"\nWrote {len(rows)} rows to {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
