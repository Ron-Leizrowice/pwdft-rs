"""SAD initial-density reference generator (VGCH Phase 1c Hypothesis 2).

Reproduces QE's ``atomic_rho.f90`` pipeline for 7 VGCH systems and emits
shell-averaged ρ(r) CSVs for comparison with pwdft-rs.
"""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pwdft_validation.integrate import simpson_qe
from pwdft_validation.units import BOHR_TO_ANG
from pwdft_validation.upf import UpfData, parse_upf


@dataclass
class AtomSpec:
    element: str
    crystal_frac: tuple[float, float, float]


@dataclass
class SystemSpec:
    name: str
    lattice_type: str
    celldm1_bohr: float
    atoms: list[AtomSpec]
    grid_dims: tuple[int, int, int]
    ecutwfc_ry: float
    ecutrho_ratio: int


SYSTEMS: list[SystemSpec] = [
    SystemSpec(
        "c_diamond",
        "fcc",
        6.7409,
        [AtomSpec("C", (0.0, 0.0, 0.0)), AtomSpec("C", (0.25, 0.25, 0.25))],
        (32, 32, 32),
        30.0,
        4,
    ),
    SystemSpec("al_fcc", "fcc", 7.6527, [AtomSpec("Al", (0.0, 0.0, 0.0))], (24, 24, 24), 24.0, 4),
    SystemSpec("fe_bcc", "bcc", 5.4235, [AtomSpec("Fe", (0.0, 0.0, 0.0))], (24, 24, 24), 15.0, 4),
    SystemSpec("cu_fcc", "fcc", 6.8219, [AtomSpec("Cu", (0.0, 0.0, 0.0))], (24, 24, 24), 25.0, 4),
    SystemSpec(
        "gaas",
        "fcc",
        10.6829,
        [AtomSpec("Ga", (0.0, 0.0, 0.0)), AtomSpec("As", (0.25, 0.25, 0.25))],
        (32, 32, 32),
        20.0,
        4,
    ),
    SystemSpec(
        "nacl",
        "rocksalt_fcc",
        10.6078,
        [AtomSpec("Na", (0.0, 0.0, 0.0)), AtomSpec("Cl", (0.5, 0.5, 0.5))],
        (32, 32, 32),
        25.0,
        4,
    ),
    SystemSpec(
        "mgo",
        "rocksalt_fcc",
        7.9586,
        [AtomSpec("Mg", (0.0, 0.0, 0.0)), AtomSpec("O", (0.5, 0.5, 0.5))],
        (24, 24, 24),
        30.0,
        4,
    ),
]


def _build_lattice(spec: SystemSpec) -> tuple[np.ndarray, list[tuple[str, np.ndarray]]]:
    a = spec.celldm1_bohr * BOHR_TO_ANG
    if spec.lattice_type in ("fcc", "rocksalt_fcc"):
        A = (a / 2.0) * np.array([[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]])
    elif spec.lattice_type == "bcc":
        A = (a / 2.0) * np.array([[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]])
    else:
        raise ValueError(f"unknown lattice {spec.lattice_type!r}")
    atoms = [(aspec.element, np.asarray(aspec.crystal_frac) @ A) for aspec in spec.atoms]
    return A, atoms


def _reciprocal(A: np.ndarray) -> np.ndarray:
    a1, a2, a3 = A[0], A[1], A[2]
    omega = np.dot(a1, np.cross(a2, a3))
    b1 = 2 * math.pi * np.cross(a2, a3) / omega
    b2 = 2 * math.pi * np.cross(a3, a1) / omega
    b3 = 2 * math.pi * np.cross(a1, a2) / omega
    return np.stack([b1, b2, b3])


def _miller_grid(dims: tuple[int, int, int]) -> np.ndarray:
    nx, ny, nz = dims
    idx = np.arange(nx * ny * nz)
    i1, i2, i3 = idx // (ny * nz), (idx // nz) % ny, idx % nz
    n1 = np.where(i1 > nx // 2, i1 - nx, i1)
    n2 = np.where(i2 > ny // 2, i2 - ny, i2)
    n3 = np.where(i3 > nz // 2, i3 - nz, i3)
    return np.stack([n1, n2, n3], axis=-1)


def _rho_at_of_g(pp: UpfData, g_norms: np.ndarray, omega_ang3: float) -> np.ndarray:
    assert pp.rho_at_ebohr is not None
    r = pp.r_bohr * BOHR_TO_ANG
    rab = pp.rab_bohr * BOHR_TO_ANG
    rho = pp.rho_at_ebohr / BOHR_TO_ANG
    out = np.empty_like(g_norms)
    for k, g in enumerate(g_norms):
        if g < 1e-12:
            integrand = rho
        else:
            gr = g * r
            j0 = np.where(gr < 1e-10, 1.0 - gr * gr / 6.0, np.sin(gr) / np.where(gr == 0.0, 1.0, gr))
            integrand = rho * j0
        out[k] = simpson_qe(integrand, rab)
    return out / omega_ang3


def _build_rho_init(
    spec: SystemSpec, pp_cache: dict[str, UpfData], renormalize: bool
) -> tuple[np.ndarray, np.ndarray, list, float]:
    A, atoms = _build_lattice(spec)
    omega = float(np.dot(A[0], np.cross(A[1], A[2])))
    B = _reciprocal(A)
    dims = spec.grid_dims
    nx, ny, nz = dims
    n_grid = nx * ny * nz
    n_miller = _miller_grid(dims)
    G = n_miller.astype(np.float64) @ B
    G_norms = np.linalg.norm(G, axis=1)
    G_uniq, inv_idx = np.unique(G_norms.round(decimals=10), return_inverse=True)
    rho_g = np.zeros(n_grid, dtype=np.complex128)
    species: dict[str, list[np.ndarray]] = {}
    for elem, cart in atoms:
        species.setdefault(elem, []).append(cart)
    for elem, taus in species.items():
        rho_at_uniq = _rho_at_of_g(pp_cache[elem], G_uniq, omega)
        rho_at_full = rho_at_uniq[inv_idx]
        taus_arr = np.stack(taus)
        phase = -G @ taus_arr.T
        strf = np.sum(np.exp(1j * phase), axis=1)
        rho_g += strf * rho_at_full
    if renormalize:
        n_el = sum(pp_cache[a.element].z_valence for a in spec.atoms)
        charge = float(omega * rho_g[0].real)
        if abs(charge) > 1e-8:
            rho_g *= n_el / charge
    rho_r = np.fft.ifftn(rho_g.reshape(nx, ny, nz)).real * n_grid
    return rho_r, A, atoms, omega


def _shell_average(
    rho_3d: np.ndarray, dims: tuple, A: np.ndarray, tau: np.ndarray, r_edges: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    nx, ny, nz = dims
    ix, iy, iz = np.meshgrid(np.arange(nx) / nx, np.arange(ny) / ny, np.arange(nz) / nz, indexing="ij")
    frac = np.stack([ix, iy, iz], axis=-1).astype(np.float64)
    B = _reciprocal(A)
    tau_frac = (tau @ B.T) / (2 * math.pi)
    df = frac - tau_frac
    df = df - np.floor(df + 0.5)
    dist = np.linalg.norm(df @ A, axis=-1)
    bin_idx = np.digitize(dist.ravel(), r_edges) - 1
    n_bins = len(r_edges) - 1
    rho_flat = rho_3d.ravel()
    valid = (bin_idx >= 0) & (bin_idx < n_bins)
    counts = np.bincount(bin_idx[valid], minlength=n_bins)
    rho_sum = np.bincount(bin_idx[valid], weights=rho_flat[valid], minlength=n_bins)
    rho_avg = np.where(counts > 0, rho_sum / counts, 0.0)
    return rho_avg, counts


def generate(pseudo_dir: Path, out_csv: Path) -> int:
    """Generate ``vgch_sad_heavy.csv`` and ``vgch_sad_heavy_samples.csv``."""
    r_edges = np.linspace(0.0, 2.0, 51)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    rows: list[list] = []

    print(f"{'system':>12} {'atom':>4} {'dims':>12} {'Ω (Å³)':>10} {'∫ρdr':>8} {'N_el':>6} {'rho(r≈0)':>12}")
    print("-" * 80)

    for spec in SYSTEMS:
        unique_elems = sorted({a.element for a in spec.atoms})
        pp_cache: dict[str, UpfData] = {}
        for elem in unique_elems:
            upf = pseudo_dir / f"{elem}.upf"
            if not upf.exists():
                print(f"# SKIP {spec.name}: {upf} missing", file=sys.stderr)
                pp_cache = {}
                break
            pp_cache[elem] = parse_upf(upf)
        if not pp_cache:
            continue

        rho_3d, A, atoms, omega = _build_rho_init(spec, pp_cache, renormalize=True)
        dvol = omega / rho_3d.size
        n_el = sum(pp_cache[elem].z_valence for elem in (a.element for a in spec.atoms))

        for atom_idx, (elem, cart) in enumerate(atoms):
            rho_avg, counts = _shell_average(rho_3d, spec.grid_dims, A, cart, r_edges)
            label = f"{elem}{atom_idx}"
            if atom_idx == 0:
                integrated = float(rho_3d.sum() * dvol)
                dims = spec.grid_dims
                print(
                    f"{spec.name:>12} {label:>4} {dims[0]}×{dims[1]}×{dims[2]:>3} "
                    f"{omega:>10.3f} {integrated:>8.3f} {n_el:>6.2f} {rho_avg[0]:>12.4e}"
                )
            for i in range(len(r_centers)):
                if counts[i] == 0:
                    continue
                rows.append([spec.name, label, f"{r_centers[i]:.8f}", f"{rho_avg[i]:.10e}", str(int(counts[i]))])

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "atom_label", "r_bin_center_ang", "rho_avg_e_per_ang3", "n_bin_points"])
        w.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {out_csv}", file=sys.stderr)

    # Raw-grid samples (200 per system).
    sample_rows: list[list] = []
    for spec in SYSTEMS:
        unique_elems = sorted({a.element for a in spec.atoms})
        pp_cache2: dict[str, UpfData] = {}
        for elem in unique_elems:
            upf = pseudo_dir / f"{elem}.upf"
            if not upf.exists():
                pp_cache2 = {}
                break
            pp_cache2[elem] = parse_upf(upf)
        if not pp_cache2:
            continue
        rho_3d, _, _, _ = _build_rho_init(spec, pp_cache2, renormalize=True)
        nx, ny, nz = spec.grid_dims
        step = max(1, (nx * ny * nz) // 200)
        for flat_idx in range(0, nx * ny * nz, step):
            i1 = flat_idx // (ny * nz)
            i2 = (flat_idx // nz) % ny
            i3 = flat_idx % nz
            sample_rows.append([spec.name, str(i1), str(i2), str(i3), f"{rho_3d[i1, i2, i3]:.10e}"])

    sample_csv = out_csv.with_name("vgch_sad_heavy_samples.csv")
    with sample_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "i1", "i2", "i3", "rho_e_per_ang3"])
        w.writerows(sample_rows)
    print(f"Wrote {len(sample_rows)} rows to {sample_csv}", file=sys.stderr)
    return 0
