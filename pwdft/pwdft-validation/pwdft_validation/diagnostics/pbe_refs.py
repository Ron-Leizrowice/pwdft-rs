"""Extract PBE reference data from QE output logs (GGAP Phase F)."""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

from pwdft_validation.paths import PSEUDO_DIR, QE_REF_DIR

_GAMMA_RE = re.compile(
    r"k =\s*0\.0000\s+0\.0000\s+0\.0000\s*\(\s*\d+\s*PWs\s*\)\s*bands\s*\(ev\):\s*\n\s*\n(.*?)\n\s*\n",
    re.DOTALL,
)

_SYSTEMS = [
    ("si_diamond_pbe", "si_scf_pbe.in", 5.431, 24.0, [4, 4, 4], 0.01, 1, ["Si.upf"], None),
    ("c_diamond_pbe", "c_diamond_scf_pbe.in", 3.567, 36.0, [4, 4, 4], 0.01, 1, ["C.upf"], None),
    ("al_fcc_pbe", "al_fcc_scf_pbe.in", 4.05, 24.0, [8, 8, 8], 0.02, 1, ["Al.upf"], None),
    ("fe_bcc_fm_pbe", "fe_bcc_fm_scf_pbe.in", 2.87, 60.0, [8, 8, 8], 0.02, 2, ["Fe.upf"], 0.5),
    ("gaas_zincblende_pbe", "gaas_scf_pbe.in", 5.653, 44.0, [4, 4, 4], 0.01, 1, ["Ga.upf", "As.upf"], None),
    ("cu_fcc_pbe", "cu_fcc_scf_pbe.in", 3.61, 60.0, [8, 8, 8], 0.02, 1, ["Cu.upf"], None),
    ("nacl_rocksalt_pbe", "nacl_scf_pbe.in", 5.614, 36.0, [4, 4, 4], 0.01, 1, ["Na.upf", "Cl.upf"], None),
    ("mgo_rocksalt_pbe", "mgo_scf_pbe.in", 4.212, 48.0, [4, 4, 4], 0.01, 1, ["Mg.upf", "O.upf"], None),
]


def _pp_hash(p: Path) -> str:
    with p.open("rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def _parse_out(path: Path, nspin: int = 1) -> dict:
    txt = path.read_text()
    m = re.search(r"^!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry", txt, re.MULTILINE)
    total_ry = float(m.group(1)) if m else None
    m = re.search(r"the Fermi energy is\s+(-?\d+\.\d+)\s*ev", txt)
    fermi_ev = float(m.group(1)) if m else None
    m = re.search(r"convergence has been achieved in\s+(\d+)\s+iterations", txt)
    n_iter = int(m.group(1)) if m else None
    end_scf = txt.find("End of self-consistent calculation")
    post = txt[max(end_scf, 0) :]
    matches = _GAMMA_RE.findall(post)

    def _parse_block(block: str) -> list[float]:
        return [float(t) for line in block.strip().splitlines() for t in line.split()]

    up, dn = [], []
    if nspin == 2 and len(matches) >= 2:
        up, dn = _parse_block(matches[0]), _parse_block(matches[1])
    elif matches:
        up = _parse_block(matches[0])

    mags_tot = re.findall(r"total magnetization\s*=\s*(-?\d+\.\d+)\s*Bohr", txt)
    mags_abs = re.findall(r"absolute magnetization\s*=\s*(-?\d+\.\d+)\s*Bohr", txt)
    return {
        "total_energy_ry": total_ry,
        "fermi_ev": fermi_ev,
        "n_iter": n_iter,
        "gamma_eigs_up": up,
        "gamma_eigs_dn": dn,
        "total_magnetization_mub": float(mags_tot[-1]) if mags_tot else None,
        "absolute_magnetization_mub": float(mags_abs[-1]) if mags_abs else None,
    }


def _fmt_eigs(evs: list[float], n_max: int = 8) -> str:
    return "[" + ", ".join(f"{v:.4f}" for v in evs[:n_max]) + "]"


def extract(qe_dir: Path | None = None, pseudo_dir: Path | None = None) -> int:
    """Print TOML snippet for all 8 PBE reference systems to stdout."""
    qe_dir = qe_dir or QE_REF_DIR
    pbe_pseudo_dir = (pseudo_dir or PSEUDO_DIR) / "nc" / "pbe"

    lines = [
        "# --- PBE reference data (GGAP Phase F-pre) ---",
        "#",
        "# Generated with QE 7.5 (qe-7.5/build/bin/pw.x), 16 MPI ranks,",
        "# PseudoDojo ONCV NC/PBE v0.4 .standard pseudopotentials.",
        "",
    ]
    summary_rows = []

    for section, in_file, a_ang, ecut_ry, k_grid, degauss_ry, nspin, pp_files, start_mag in _SYSTEMS:
        out_path = qe_dir / in_file.replace(".in", ".out")
        parsed = _parse_out(out_path, nspin=nspin)
        pp_hashes = {}
        for p in pp_files:
            pp_path = pbe_pseudo_dir / p
            if pp_path.exists():
                pp_hashes[p] = _pp_hash(pp_path)[:16]

        lines += [
            f"[{section}]",
            f'input_file        = "{in_file}"',
            f"lattice_param_ang = {a_ang}",
            f"ecutwfc_ry        = {ecut_ry}",
            f"k_grid            = {k_grid}",
            f"degauss_ry        = {degauss_ry}",
            f"nspin             = {nspin}",
        ]
        if start_mag is not None:
            lines.append(f"starting_magnetization = {start_mag}")
        lines += [
            f"total_energy_ry   = {parsed['total_energy_ry']:.8f}",
            f"fermi_energy_ev   = {parsed['fermi_ev']:.4f}",
            f"n_iterations      = {parsed['n_iter']}",
        ]
        if nspin == 2:
            lines += [
                f"total_magnetization_mub = {parsed['total_magnetization_mub']:.2f}",
                f"absolute_magnetization_mub = {parsed['absolute_magnetization_mub']:.2f}",
                f"gamma_eigenvalues_up_ev = {_fmt_eigs(parsed['gamma_eigs_up'])}",
                f"gamma_eigenvalues_dn_ev = {_fmt_eigs(parsed['gamma_eigs_dn'])}",
            ]
        else:
            lines.append(f"gamma_eigenvalues_ev = {_fmt_eigs(parsed['gamma_eigs_up'])}")
        for p, h in pp_hashes.items():
            lines.append(f'# pp_sha256_{p.replace(".", "_")} = "{h}..."')
        lines.append("")

        e_ev = parsed["total_energy_ry"] * 13.605693122994 if parsed["total_energy_ry"] else float("nan")
        summary_rows.append(
            (
                section,
                ecut_ry,
                k_grid,
                parsed["total_energy_ry"],
                e_ev,
                parsed["fermi_ev"],
                parsed["n_iter"],
                pp_files,
                parsed.get("absolute_magnetization_mub"),
            )
        )

    print("\n".join(lines))
    print("=== SUMMARY ===", file=sys.stderr)
    print(
        f"{'system':<22} {'ecut':>6} {'kgrid':<12} {'E_tot_Ry':>14} {'E_tot_eV':>14} {'E_F_eV':>8} {'it':>3} PP",
        file=sys.stderr,
    )
    for section, ecut_ry, k_grid, e_ry, e_ev, ef, ni, pp_files, am in summary_rows:
        mag_str = f" |M|={am:.2f} μB" if am is not None else ""
        print(
            f"{section:<22} {ecut_ry:>6.1f} {k_grid!s:<12} {e_ry:>14.6f} {e_ev:>14.4f} {ef:>8.4f} {ni:>3} {','.join(pp_files)}{mag_str}",
            file=sys.stderr,
        )
    return 0
