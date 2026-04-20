"""Extract PBE reference data from QE output logs for GGAP Phase F.

Parses each *_pbe.out file for:
  - total energy (Ry)
  - Fermi energy (eV)
  - convergence iteration count
  - Gamma-point eigenvalues (eV; first 8 bands)
  - (spin-polarized only) total/absolute magnetization (Bohr mag/cell)

Emits a TOML snippet suitable for appending to reference_data.toml.
"""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

from pwdft_validation import PSEUDO_DIR, QE_REF_DIR

GAMMA_RE = re.compile(
    r"k =\s*0\.0000\s+0\.0000\s+0\.0000\s*\(\s*\d+\s*PWs\s*\)\s*bands\s*\(ev\):\s*\n\s*\n(.*?)\n\s*\n",
    re.DOTALL,
)


def pp_hash(pp_file: Path) -> str:
    with pp_file.open("rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def parse_out(path: Path, nspin: int = 1):
    """Return dict of parsed fields from a QE pw.x SCF output."""
    txt = path.read_text()

    # Total energy line: "!    total energy              =     -16.91056535 Ry"
    m = re.search(r"^!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry", txt, re.M)
    total_energy_ry = float(m.group(1)) if m else None

    # Fermi energy: "the Fermi energy is     6.5649 ev"
    m = re.search(r"the Fermi energy is\s+(-?\d+\.\d+)\s*ev", txt)
    fermi_ev = float(m.group(1)) if m else None

    # Convergence iterations: "convergence has been achieved in   7 iterations"
    m = re.search(r"convergence has been achieved in\s+(\d+)\s+iterations", txt)
    n_iter = int(m.group(1)) if m else None

    # Gamma eigenvalues. Find the first k-block at k=0,0,0 AFTER "End of self-consistent"
    # For nspin=2 there are two: spin up and spin down (two sections labeled "---" etc).
    end_scf = txt.find("End of self-consistent calculation")
    if end_scf < 0:
        end_scf = 0
    post_scf = txt[end_scf:]

    # Match "k = 0.0000 0.0000 0.0000 (   NNN PWs)   bands (ev):" followed by eigenvalue lines
    # until the next blank line.
    gamma_eigs_up = []
    gamma_eigs_dn = []

    matches = GAMMA_RE.findall(post_scf)

    def parse_block(block: str):
        vals = []
        for line in block.strip().splitlines():
            for tok in line.split():
                vals.append(float(tok))
        return vals

    if nspin == 2:
        if len(matches) >= 2:
            gamma_eigs_up = parse_block(matches[0])
            gamma_eigs_dn = parse_block(matches[1])
    else:
        if matches:
            gamma_eigs_up = parse_block(matches[0])

    # Magnetization (nspin=2 only)
    tot_mag = None
    abs_mag = None
    # "total magnetization       =     2.34 Bohr mag/cell"
    # Use the LAST occurrence (converged).
    mags_tot = re.findall(r"total magnetization\s*=\s*(-?\d+\.\d+)\s*Bohr", txt)
    mags_abs = re.findall(r"absolute magnetization\s*=\s*(-?\d+\.\d+)\s*Bohr", txt)
    if mags_tot:
        tot_mag = float(mags_tot[-1])
    if mags_abs:
        abs_mag = float(mags_abs[-1])

    return {
        "total_energy_ry": total_energy_ry,
        "fermi_ev": fermi_ev,
        "n_iter": n_iter,
        "gamma_eigs_up": gamma_eigs_up,
        "gamma_eigs_dn": gamma_eigs_dn,
        "total_magnetization_mub": tot_mag,
        "absolute_magnetization_mub": abs_mag,
    }


# Schema: (toml_section, input_file, lattice_ang, ecutwfc_ry, k_grid, degauss_ry, nspin, pp_files, starting_mag)
SYSTEMS = [
    ("si_diamond_pbe", "si_scf_pbe.in", 5.431, 24.0, [4, 4, 4], 0.01, 1, ["Si.upf"], None),
    ("c_diamond_pbe", "c_diamond_scf_pbe.in", 3.567, 36.0, [4, 4, 4], 0.01, 1, ["C.upf"], None),
    ("al_fcc_pbe", "al_fcc_scf_pbe.in", 4.05, 24.0, [8, 8, 8], 0.02, 1, ["Al.upf"], None),
    ("fe_bcc_fm_pbe", "fe_bcc_fm_scf_pbe.in", 2.87, 60.0, [8, 8, 8], 0.02, 2, ["Fe.upf"], 0.5),
    ("gaas_zincblende_pbe", "gaas_scf_pbe.in", 5.653, 44.0, [4, 4, 4], 0.01, 1, ["Ga.upf", "As.upf"], None),
    ("cu_fcc_pbe", "cu_fcc_scf_pbe.in", 3.61, 60.0, [8, 8, 8], 0.02, 1, ["Cu.upf"], None),
    ("nacl_rocksalt_pbe", "nacl_scf_pbe.in", 5.614, 36.0, [4, 4, 4], 0.01, 1, ["Na.upf", "Cl.upf"], None),
    ("mgo_rocksalt_pbe", "mgo_scf_pbe.in", 4.212, 48.0, [4, 4, 4], 0.01, 1, ["Mg.upf", "O.upf"], None),
]


def fmt_eigs(evs, n_max=8):
    return "[" + ", ".join(f"{v:.4f}" for v in evs[:n_max]) + "]"


def main():
    lines = []
    summary_rows = []
    pseudo_dir = PSEUDO_DIR / "nc" / "pbe"

    lines.append("# --- PBE reference data (GGAP Phase F-pre) ---")
    lines.append("#")
    lines.append("# Generated 2026-04-19 with QE 7.5 (qe-7.5/build/bin/pw.x), 16 MPI ranks,")
    lines.append("# PseudoDojo ONCV NC/PBE v0.4 .standard pseudopotentials.")
    lines.append("# Each *_pbe section mirrors its LDA sibling's cell/k-grid/smearing;")
    lines.append("# ecut is lifted to the PseudoDojo PBE .standard recommendation per")
    lines.append("# species-max (Ha → Ry = ×2) so the reference is basis-converged.")
    lines.append("")

    for section, in_file, a_ang, ecut_ry, k_grid, degauss_ry, nspin, pp_files, start_mag in SYSTEMS:
        out_path = QE_REF_DIR / in_file.replace(".in", ".out")
        parsed = parse_out(out_path, nspin=nspin)

        # PP hashes
        pp_hashes = {p: pp_hash(pseudo_dir / p)[:16] for p in pp_files}

        lines.append(f"[{section}]")
        lines.append(f'input_file        = "{in_file}"')
        lines.append(f"lattice_param_ang = {a_ang}")
        lines.append(f"ecutwfc_ry        = {ecut_ry}")
        lines.append(f"k_grid            = {k_grid}")
        lines.append(f"degauss_ry        = {degauss_ry}")
        lines.append(f"nspin             = {nspin}")
        if start_mag is not None:
            lines.append(f"starting_magnetization = {start_mag}")
        lines.append(f"total_energy_ry   = {parsed['total_energy_ry']:.8f}")
        lines.append(f"fermi_energy_ev   = {parsed['fermi_ev']:.4f}")
        lines.append(f"n_iterations      = {parsed['n_iter']}")
        if nspin == 2:
            lines.append(f"total_magnetization_mub = {parsed['total_magnetization_mub']:.2f}")
            lines.append(f"absolute_magnetization_mub = {parsed['absolute_magnetization_mub']:.2f}")
            lines.append(f"gamma_eigenvalues_up_ev = {fmt_eigs(parsed['gamma_eigs_up'])}")
            lines.append(f"gamma_eigenvalues_dn_ev = {fmt_eigs(parsed['gamma_eigs_dn'])}")
        else:
            lines.append(f"gamma_eigenvalues_ev = {fmt_eigs(parsed['gamma_eigs_up'])}")

        # PP hash comments
        for pp, h in pp_hashes.items():
            lines.append(f'# pp_sha256_{pp.replace(".", "_")} = "{h}..."')
        lines.append("")

        e_ev = parsed["total_energy_ry"] * 13.605693122994
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
                parsed.get("total_magnetization_mub"),
                parsed.get("absolute_magnetization_mub"),
            )
        )

    print("\n".join(lines))
    print()
    print("=== SUMMARY ===", file=sys.stderr)
    print(
        f"{'system':<22} {'ecut':>6} {'kgrid':<12} {'E_tot_Ry':>14} {'E_tot_eV':>14} {'E_F_eV':>8} {'it':>3} PP",
        file=sys.stderr,
    )
    for row in summary_rows:
        section, ecut_ry, k_grid, e_ry, e_ev, ef, ni, pp_files, _tm, am = row
        mag_str = f" |M|={am:.2f} μB" if am is not None else ""
        print(
            f"{section:<22} {ecut_ry:>6.1f} {k_grid!s:<12} {e_ry:>14.6f} {e_ev:>14.4f} {ef:>8.4f} {ni:>3} {','.join(pp_files)}{mag_str}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
