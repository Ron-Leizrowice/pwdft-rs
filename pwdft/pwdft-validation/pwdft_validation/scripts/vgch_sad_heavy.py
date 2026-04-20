#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy>=1.26",
#     "scipy>=1.11",
# ]
# ///
"""
VGCH Phase 1c Hypothesis 2 — SAD (Superposition of Atomic Densities)
initial-density reference.

For each of the 7 VGCH-class systems (C, Al, Fe, Cu, GaAs, NaCl, MgO),
compute the initial charge density ρ(r) on a fixed FFT grid by faithfully
reproducing the QE `atomic_rho.f90` pipeline:

  1. For each species type, Bessel-transform the UPF PP_RHOATOM
     radial atomic density via QE's Simpson rule
     (`qe-7.5/upflib/rhoat_mod.f90:101-111`):

         ρ_at(G)  =  (1/Ω) · Simpson_{i}( rho_at_i · sin(q r_i)/(q r_i) ; rab_i )

     where rho_at stored in UPF is 4π r² ρ(r) already (so no r²
     weighting here — QE confirms this at lines 103-106).

  2. Sum over all atoms of the species with the structure factor
     `strf(G,nt) = Σ_{I ∈ species} exp(-iG·τ_I)`:

         ρ(G) = Σ_nt strf(G,nt) · ρ_at_nt(G)

  3. Inverse FFT to real space:

         ρ(r) = Σ_G ρ(G) exp(+iG·r)

     (pwdft-rs' unnormalized inverse convention; matches QE's
     `rho_g2r` up to the same Fourier convention.)

  4. Compute shell-average ρ(r) vs r around each atom by binning
     grid points by distance-to-atom (minimum image convention).
     Write `data/csv/vgch_sad_heavy.csv` with rows
     `(system, atom_label, r_bin_center_ang, rho_avg_e_per_ang3,
     n_bin_points)`.

The Rust test `tests/vgch_sad_heavy.rs` calls
`pwdft_rs::scf::initial_density::build_sad_density_for_diagnostic`
on the *same* FFT grid dims and the same crystal, computes the same
shell-average, and compares bin-by-bin.

Key design choices:
  - **Grid dims fixed per system** (hardcoded below) to isolate
    Bessel-transform + structure-factor + FFT-convention agreement
    from grid-selection disagreement.
  - **No clamp, no renormalization in Python.** pwdft-rs' current code
    clamps `ρ<0 → 0` in real space and renormalizes to N_el (see
    `src/scf/initial_density.rs:109-123`). QE skips both steps (see
    `qe-7.5/PW/src/atomic_rho.f90:186-188` for QE's explicit
    "negative charge will re-appear" comment). If the diagnostic
    reveals large Δρ(r) near atom cores, this asymmetry is the
    suspect. A second Python variant applies QE's G=0 renormalization
    (`rho%of_g = rho%of_g * nelec / charge` on line 223 of
    `potinit.f90`) and is reported alongside as an optional toggle.
  - **C diamond is the critical light-atom control.** If C's SAD
    agrees bit-perfect but C still has a 1.45 eV E_total residual,
    H2 is cleared — the bug is either in the mixer basin (H3) or in
    the energy-assembly code (V_loc(G=0) compensation).

The 7 VGCH systems mirror `data/qe/*.in`. GaAs, NaCl, MgO are
two-species; others are one-species.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

BOHR_TO_ANG = 0.529_177_210_903


# ---------------------------------------------------------------------------
# System registry — mirrors `data/qe/*.in` and
# `tests/vgch_per_component_heavy.rs`. Grid dims are chosen small enough
# to make shell-averaging well-sampled near the atom cores (nx·ny·nz ≈
# a few thousand points per Å³ cell) while large enough to resolve the
# PP_RHOATOM Bessel transform at production ecut.
#
# `lattice_type` is:
#   - "fcc"    : FCC conventional side a; primitive basis
#                (0,a/2,a/2), (a/2,0,a/2), (a/2,a/2,0).
#   - "bcc"    : BCC conventional side a; primitive basis
#                (-a/2,a/2,a/2), (a/2,-a/2,a/2), (a/2,a/2,-a/2).
#   - "rocksalt_fcc" : FCC with two-atom basis (species A at
#                0,0,0; species B at 0.5,0.5,0.5 in crystal coords).
# `ecutwfc` and `ecutrho_ratio` are not used by Python (we pass grid
# dims explicitly) but are echoed in the CSV header for audit.
# ---------------------------------------------------------------------------


@dataclass
class AtomSpec:
    element: str  # symbol matching `pseudopotentials/nc/lda/<Symbol>.upf`
    crystal_frac: tuple[float, float, float]


@dataclass
class SystemSpec:
    name: str
    lattice_type: str
    celldm1_bohr: float  # QE celldm(1) = side a in Bohr for ibrav=2 (FCC) or ibrav=3 (BCC)
    atoms: list[AtomSpec]
    grid_dims: tuple[int, int, int]
    ecutwfc_ry: float  # diagnostic; not used in Python
    ecutrho_ratio: int  # diagnostic; not used in Python


SYSTEMS: list[SystemSpec] = [
    # C diamond: 2-atom FCC, light atom, no semicore — critical control
    SystemSpec(
        name="c_diamond",
        lattice_type="fcc",
        celldm1_bohr=6.7409,  # from data/qe/c_diamond_scf.in
        atoms=[
            AtomSpec("C", (0.00, 0.00, 0.00)),
            AtomSpec("C", (0.25, 0.25, 0.25)),
        ],
        # 32³ chosen deliberately over 24³: the two C atoms sit at grid
        # points (0,0,0) and (8,8,8) on a 32³ grid, which cleanly
        # preserves the diamond glide-translation τ=(1/4,1/4,1/4) as a
        # integer grid translation. At 24³ the atom would land at
        # (6,6,6) — still integer, so the nint-glide issue is not the
        # cause of the asymmetry observed in the diagnostic; it's
        # inversion-symmetry breaking between the bin-around-A and
        # bin-around-B discrete shell-average buckets.
        grid_dims=(32, 32, 32),
        ecutwfc_ry=30.0,
        ecutrho_ratio=4,
    ),
    # Al FCC: 1-atom FCC, light-medium, no semicore
    SystemSpec(
        name="al_fcc",
        lattice_type="fcc",
        celldm1_bohr=7.6527,
        atoms=[AtomSpec("Al", (0.00, 0.00, 0.00))],
        grid_dims=(24, 24, 24),
        ecutwfc_ry=24.0,
        ecutrho_ratio=4,
    ),
    # Fe BCC: 1-atom BCC, 3s/3p/3d semicore
    SystemSpec(
        name="fe_bcc",
        lattice_type="bcc",
        celldm1_bohr=5.4235,
        atoms=[AtomSpec("Fe", (0.00, 0.00, 0.00))],
        grid_dims=(24, 24, 24),
        ecutwfc_ry=15.0,
        ecutrho_ratio=4,
    ),
    # Cu FCC: 1-atom FCC, 3s/3p/3d semicore
    SystemSpec(
        name="cu_fcc",
        lattice_type="fcc",
        celldm1_bohr=6.8219,
        atoms=[AtomSpec("Cu", (0.00, 0.00, 0.00))],
        grid_dims=(24, 24, 24),
        ecutwfc_ry=25.0,
        ecutrho_ratio=4,
    ),
    # GaAs zinc-blende: 2-atom FCC, Ga 3d semicore, As 3d semicore
    SystemSpec(
        name="gaas",
        lattice_type="fcc",
        celldm1_bohr=10.6829,
        atoms=[
            AtomSpec("Ga", (0.00, 0.00, 0.00)),
            AtomSpec("As", (0.25, 0.25, 0.25)),
        ],
        grid_dims=(32, 32, 32),
        ecutwfc_ry=20.0,
        ecutrho_ratio=4,
    ),
    # NaCl rocksalt: 2-atom FCC, Na 2s/2p semicore (ONCVPSP z=9)
    SystemSpec(
        name="nacl",
        lattice_type="rocksalt_fcc",
        celldm1_bohr=10.6078,
        atoms=[
            AtomSpec("Na", (0.00, 0.00, 0.00)),
            AtomSpec("Cl", (0.50, 0.50, 0.50)),
        ],
        grid_dims=(32, 32, 32),
        ecutwfc_ry=25.0,
        ecutrho_ratio=4,
    ),
    # MgO rocksalt: 2-atom FCC, Mg 2s/2p semicore
    SystemSpec(
        name="mgo",
        lattice_type="rocksalt_fcc",
        celldm1_bohr=7.9586,
        atoms=[
            AtomSpec("Mg", (0.00, 0.00, 0.00)),
            AtomSpec("O", (0.50, 0.50, 0.50)),
        ],
        grid_dims=(24, 24, 24),
        ecutwfc_ry=30.0,
        ecutrho_ratio=4,
    ),
]


# ---------------------------------------------------------------------------
# UPF v2 parser — regex, consistent with `vgch_beta_l_heavy.py`.
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


@dataclass
class UpfAtomicRho:
    element: str
    z_valence: float
    mesh: int
    r_bohr: np.ndarray  # Bohr
    rab_bohr: np.ndarray  # Bohr (log-mesh Jacobian)
    rho_at_ebohr: np.ndarray  # e/Bohr (4π r² ρ(r), stored convention)


def parse_upf_rhoatom(upf_path: Path) -> UpfAtomicRho:
    text = upf_path.read_text()
    m_hdr = re.search(r"<PP_HEADER\b([^>]*)/?>", text)
    if m_hdr is None:
        raise ValueError(f"{upf_path}: PP_HEADER not found")
    hdr = m_hdr.group(1)
    element = _extract_attr(hdr, "element")
    z_val = float(_extract_attr(hdr, "z_valence"))
    mesh = int(_extract_attr(hdr, "mesh_size"))
    r = _extract_block(text, "PP_R")
    rab = _extract_block(text, "PP_RAB")
    rho_at = _extract_block(text, "PP_RHOATOM")
    assert r.size == mesh and rab.size == mesh and rho_at.size == mesh
    return UpfAtomicRho(
        element=element,
        z_valence=z_val,
        mesh=mesh,
        r_bohr=r,
        rab_bohr=rab,
        rho_at_ebohr=rho_at,
    )


# ---------------------------------------------------------------------------
# Simpson integrator — byte-identical to QE `simpson` (even-mesh correction
# included) and pwdft-rs `simpson_integrate`.
# ---------------------------------------------------------------------------


def simpson_qe(func: np.ndarray, rab: np.ndarray) -> float:
    n = len(func)
    assert n == len(rab)
    if n < 3:
        return float(np.sum(func * rab))
    acc = 0.0
    for i in range(1, n - 1):
        weight = 4.0 if (i % 2 == 1) else 2.0
        acc += weight * func[i] * rab[i]
    if n % 2 == 1:
        return (acc + func[0] * rab[0] + func[-1] * rab[-1]) / 3.0
    acc += func[0] * rab[0]
    acc += -0.25 * func[n - 3] * rab[n - 3]
    acc += func[n - 2] * rab[n - 2]
    acc += 1.25 * func[n - 1] * rab[n - 1]
    return acc / 3.0


# ---------------------------------------------------------------------------
# Lattice builders — primitive cells matching pwdft-rs conventions.
# Return (a_ang, lattice_rows_ang, atom_cartesians_ang).
# ---------------------------------------------------------------------------


def build_lattice_and_atoms(
    sys_spec: SystemSpec,
) -> tuple[float, np.ndarray, list[tuple[str, np.ndarray]]]:
    a_ang = sys_spec.celldm1_bohr * BOHR_TO_ANG
    if sys_spec.lattice_type in ("fcc", "rocksalt_fcc"):
        A = (a_ang / 2.0) * np.array(
            [
                [0.0, 1.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 0.0],
            ]
        )
    elif sys_spec.lattice_type == "bcc":
        A = (a_ang / 2.0) * np.array(
            [
                [-1.0, 1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, -1.0],
            ]
        )
    else:
        raise ValueError(f"unknown lattice type {sys_spec.lattice_type!r}")

    atoms = []
    for aspec in sys_spec.atoms:
        frac = np.asarray(aspec.crystal_frac, dtype=np.float64)
        # row-vector convention: cart = frac · A  (matches Crystal.cart_position)
        cart = frac @ A
        atoms.append((aspec.element, cart))
    return a_ang, A, atoms


def reciprocal_lattice(A: np.ndarray) -> np.ndarray:
    """Compute reciprocal lattice B from direct A (rows a,b,c) such that
    a_i · b_j = 2π δ_ij. Returns B as rows."""
    a1, a2, a3 = A[0], A[1], A[2]
    omega = np.dot(a1, np.cross(a2, a3))
    b1 = 2.0 * math.pi * np.cross(a2, a3) / omega
    b2 = 2.0 * math.pi * np.cross(a3, a1) / omega
    b3 = 2.0 * math.pi * np.cross(a1, a2) / omega
    return np.stack([b1, b2, b3], axis=0)


# ---------------------------------------------------------------------------
# FFT grid indexing — matches `src/scf/grid.rs::g_vector_at_dims`:
#   Flat idx = i1·(ny·nz) + i2·nz + i3   (row-major, C order)
#   Miller n_i = i_i               if i_i <= n_i/2
#             = i_i - n_i          if i_i  > n_i/2
# ---------------------------------------------------------------------------


def miller_grid(dims: tuple[int, int, int]) -> np.ndarray:
    """Return (N, 3) array of Miller indices for each flat grid point,
    in the pwdft-rs row-major flat order."""
    nx, ny, nz = dims
    idx = np.arange(nx * ny * nz)
    i1 = idx // (ny * nz)
    i2 = (idx // nz) % ny
    i3 = idx % nz
    n1 = np.where(i1 > nx // 2, i1 - nx, i1)
    n2 = np.where(i2 > ny // 2, i2 - ny, i2)
    n3 = np.where(i3 > nz // 2, i3 - nz, i3)
    return np.stack([n1, n2, n3], axis=-1)


# ---------------------------------------------------------------------------
# Bessel transform of one species' PP_RHOATOM at each unique |G|.
# Internal units: convert r, rab, rho_at to Å before integrating
# (r_Å = r_B · BOHR, rab_Å = rab_B · BOHR, rho_at_Å = rho_at_B / BOHR).
# ---------------------------------------------------------------------------


def rho_at_of_g(
    pp: UpfAtomicRho,
    g_norms_per_ang: np.ndarray,
    omega_ang3: float,
) -> np.ndarray:
    """Return ρ_at(G) / Ω in e/Å³ for each |G| in `g_norms_per_ang`.
    ρ_at(G) = Simpson( [4πr²ρ(r)] · sin(Gr)/(Gr) ; rab ) in Å-units."""
    r_ang = pp.r_bohr * BOHR_TO_ANG
    rab_ang = pp.rab_bohr * BOHR_TO_ANG
    rho_ang = pp.rho_at_ebohr / BOHR_TO_ANG  # e/Å  (stored as 4π r² ρ)

    out = np.empty_like(g_norms_per_ang, dtype=np.float64)
    for k, g in enumerate(g_norms_per_ang):
        if g < 1e-12:
            # j0 = 1 everywhere
            integrand = rho_ang
        else:
            gr = g * r_ang
            # Taylor fallback for gr → 0 avoids divide-by-zero at r[0]
            j0 = np.where(gr < 1e-10, 1.0 - gr * gr / 6.0, np.sin(gr) / np.where(gr == 0.0, 1.0, gr))
            integrand = rho_ang * j0
        out[k] = simpson_qe(integrand, rab_ang)
    return out / omega_ang3


# ---------------------------------------------------------------------------
# SAD ρ(G) assembly and IFFT.
# ---------------------------------------------------------------------------


def build_rho_init_real(
    sys_spec: SystemSpec,
    pp_cache: dict[str, UpfAtomicRho],
    renormalize: bool,
) -> tuple[np.ndarray, np.ndarray, list[tuple[str, np.ndarray]], float]:
    """Compute ρ_init(r) on the fixed FFT grid via the QE-convention
    atomic_rho recipe. Returns (rho_r_3d, A_rows_ang, atoms, omega_ang3).

    `renormalize=True` applies QE's G=0 renormalization:
        charge = Ω · ρ(G=0)
        ρ(G) ← ρ(G) · (N_el / charge)
    which matches `qe-7.5/PW/src/potinit.f90:218-223`.
    """
    _a_ang, A_rows, atoms = build_lattice_and_atoms(sys_spec)
    omega_ang3 = np.dot(A_rows[0], np.cross(A_rows[1], A_rows[2]))
    B_rows = reciprocal_lattice(A_rows)

    dims = sys_spec.grid_dims
    nx, ny, nz = dims
    n_grid = nx * ny * nz

    # G vectors at each flat index (in 1/Å).
    n = miller_grid(dims)  # (N, 3)
    # G = n1·b1 + n2·b2 + n3·b3 (row vectors)
    G = n.astype(np.float64) @ B_rows  # (N, 3) in 1/Å
    G_norms = np.linalg.norm(G, axis=1)

    # Unique |G| → interpolation lookup reduces Simpson work by 4-8×
    # per grid (depends on cell symmetry).
    G_norms_unique, inverse_idx = np.unique(G_norms.round(decimals=10), return_inverse=True)

    rho_g = np.zeros(n_grid, dtype=np.complex128)

    # Group atoms by species; sum structure factor within each group.
    species_groups: dict[str, list[np.ndarray]] = {}
    for elem, cart in atoms:
        species_groups.setdefault(elem, []).append(cart)

    for elem, taus in species_groups.items():
        pp = pp_cache[elem]
        # Species-level Bessel transform at unique G.
        rho_at_unique = rho_at_of_g(pp, G_norms_unique, omega_ang3)  # (nU,)
        rho_at_full = rho_at_unique[inverse_idx]  # (N,)
        # Structure factor for this species: Σ_I exp(-iG·τ_I).
        # G @ tau is dot product of each G-row with each tau-row.
        taus_arr = np.stack(taus, axis=0)  # (n_I, 3)
        phase = -G @ taus_arr.T  # (N, n_I)
        strf = np.sum(np.exp(1j * phase), axis=1)  # (N,)
        rho_g += strf * rho_at_full

    # Optional G=0 renormalization — QE-equivalent.
    n_electrons = sum(pp_cache[aspec.element].z_valence for aspec in sys_spec.atoms)
    if renormalize:
        charge = float(omega_ang3 * rho_g[0].real)  # G=0 → flat idx 0 (Miller (0,0,0))
        if abs(charge) > 1e-8:
            rho_g *= n_electrons / charge

    # Reshape rho_g from pwdft-rs flat index order to (nx, ny, nz) for
    # numpy.fft.ifftn, which uses the same row-major C order.
    rho_g_3d = rho_g.reshape(nx, ny, nz)

    # Inverse FFT. pwdft-rs uses "unnormalized inverse": ρ(r) = Σ_G ρ(G) e^{+iGr}.
    # numpy.fft.ifftn = (1/N) Σ_G ρ(G) e^{+2πi n·k/N}, so multiply by N:
    rho_r_3d = np.fft.ifftn(rho_g_3d) * n_grid

    rho_r = rho_r_3d.real  # imag part is FFT roundoff

    return rho_r, A_rows, atoms, float(omega_ang3)


# ---------------------------------------------------------------------------
# Shell-average ρ(r) around each atom.
#
# For each atom, for each FFT grid point, compute the minimum-image
# distance to the atom in cartesian coords and bin.
#
# `r_edges_ang` default covers 0 → 2.0 Å in 50 bins of 0.04 Å width
# — dense enough to resolve the PP's r_grid structure near the core,
# spanning the atomic-density support for all elements in the set.
# ---------------------------------------------------------------------------


def shell_average(
    rho_r_3d: np.ndarray,
    dims: tuple[int, int, int],
    A_rows: np.ndarray,
    tau_cart: np.ndarray,
    r_edges_ang: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (rho_avg, n_points_per_bin).

    For each grid point at fractional coord f_grid = (i/nx, j/ny, k/nz),
    compute cart r_grid = f_grid · A. Displacement δ = r_grid - τ, wrapped
    to [-1/2, 1/2) in fractional coords (minimum image). Distance |δ|
    determines which bin the grid point contributes to."""
    nx, ny, nz = dims
    # Cartesian grid coords
    ix, iy, iz = np.meshgrid(
        np.arange(nx, dtype=np.float64) / nx,
        np.arange(ny, dtype=np.float64) / ny,
        np.arange(nz, dtype=np.float64) / nz,
        indexing="ij",
    )
    frac = np.stack([ix, iy, iz], axis=-1)  # (nx, ny, nz, 3)

    # Atom fractional coord
    B_rows = reciprocal_lattice(A_rows)
    tau_frac = (tau_cart @ B_rows.T) / (2.0 * math.pi)

    # Displacement in fractional coords, wrapped to [-0.5, +0.5).
    # Use floor-based wrap (not np.round, which rounds half-to-even) so
    # this matches the Rust-side shell_average exactly at d = ±0.5.
    # On 2-atom FCC cells with τ=(1/4,1/4,1/4), grid points where
    # d = ±0.5 exist and pick up a signed cartesian displacement that
    # differs between "round-half-to-even" and "round-half-away-from
    # zero"; the disagreement shows up as a ~20% shell-average asymmetry.
    df = frac - tau_frac
    df = df - np.floor(df + 0.5)

    # Cartesian displacement
    dc = df @ A_rows
    dist = np.linalg.norm(dc, axis=-1)  # (nx, ny, nz)

    # Bin
    bin_idx = np.digitize(dist.ravel(), r_edges_ang) - 1
    n_bins = len(r_edges_ang) - 1
    rho_flat = rho_r_3d.ravel()

    rho_avg = np.zeros(n_bins, dtype=np.float64)
    counts = np.zeros(n_bins, dtype=np.int64)
    # Use bincount for efficiency
    valid = (bin_idx >= 0) & (bin_idx < n_bins)
    counts = np.bincount(bin_idx[valid], minlength=n_bins)
    rho_sum = np.bincount(bin_idx[valid], weights=rho_flat[valid], minlength=n_bins)
    nz_mask = counts > 0
    rho_avg[nz_mask] = rho_sum[nz_mask] / counts[nz_mask]

    return rho_avg, counts


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    from pwdft_validation import CSV_REF_DIR, PSEUDO_DIR

    pp_dir = PSEUDO_DIR / "nc" / "lda"
    out_csv = CSV_REF_DIR / "vgch_sad_heavy.csv"

    # Shell bins: 50 bins from 0 → 2.0 Å, 0.04 Å wide
    r_edges = np.linspace(0.0, 2.0, 51)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])

    rows: list[list[str]] = []

    print(f"{'system':>12} {'atom':>4} {'dims':>12} {'Ω (Å³)':>10} {'∫ρdr':>8} {'N_el':>6} {'rho(r≈0)':>12}")
    print("-" * 80)

    for sys_spec in SYSTEMS:
        # Cache pseudopotentials for each unique species in the system
        unique_elements = sorted({a.element for a in sys_spec.atoms})
        pp_cache: dict[str, UpfAtomicRho] = {}
        for elem in unique_elements:
            upf = pp_dir / f"{elem}.upf"
            if not upf.exists():
                print(f"# SKIP {sys_spec.name}: {upf} missing", file=sys.stderr)
                pp_cache = {}
                break
            pp_cache[elem] = parse_upf_rhoatom(upf)
        if not pp_cache:
            continue

        # Build ρ_init(r) — match pwdft-rs by renormalizing (clamp not
        # applied here; the test compares against pwdft-rs' post-clamp
        # ρ so the "no-clamp" Python truth is the correct baseline for
        # the diagnostic).
        rho_r_3d, A_rows, atoms, omega_ang3 = build_rho_init_real(sys_spec, pp_cache, renormalize=True)

        # Integrated charge (sanity check)
        dvol = omega_ang3 / rho_r_3d.size
        integrated = float(rho_r_3d.sum() * dvol)
        n_el = sum(pp_cache[elem].z_valence for elem in (a.element for a in sys_spec.atoms))

        # Shell-average around each atom
        for atom_idx, (elem, cart) in enumerate(atoms):
            rho_avg, counts = shell_average(rho_r_3d, sys_spec.grid_dims, A_rows, cart, r_edges)
            atom_label = f"{elem}{atom_idx}"

            # First bin center (near the atom core) for quick eyeball
            if atom_idx == 0:
                first_rho = rho_avg[0]
                print(
                    f"{sys_spec.name:>12} {atom_label:>4} "
                    f"{sys_spec.grid_dims[0]}×{sys_spec.grid_dims[1]}×{sys_spec.grid_dims[2]:>3} "
                    f"{omega_ang3:>10.3f} {integrated:>8.3f} {n_el:>6.2f} {first_rho:>12.4e}"
                )

            for i in range(len(r_centers)):
                if counts[i] == 0:
                    continue  # no sample points in this bin
                rows.append(
                    [
                        sys_spec.name,
                        atom_label,
                        f"{r_centers[i]:.8f}",
                        f"{rho_avg[i]:.10e}",
                        str(int(counts[i])),
                    ]
                )

    # Write CSV
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "atom_label", "r_bin_center_ang", "rho_avg_e_per_ang3", "n_bin_points"])
        w.writerows(rows)
    print(f"\nwrote {len(rows)} rows to {out_csv}", file=sys.stderr)

    # Also dump a per-system raw-grid sample so the Rust test can do
    # pointwise diffs, not just shell averages. Limit to 200 points per
    # system to keep the file small.
    sample_rows = []
    for sys_spec in SYSTEMS:
        unique_elements = sorted({a.element for a in sys_spec.atoms})
        pp_cache2: dict[str, UpfAtomicRho] = {}
        for elem in unique_elements:
            upf = pp_dir / f"{elem}.upf"
            if not upf.exists():
                pp_cache2 = {}
                break
            pp_cache2[elem] = parse_upf_rhoatom(upf)
        if not pp_cache2:
            continue
        rho_r_3d, _A, _atoms, _omega = build_rho_init_real(sys_spec, pp_cache2, renormalize=True)
        nx, ny, nz = sys_spec.grid_dims
        # Deterministic 200-point sampling: every (nx·ny·nz // 200)-th index
        step = max(1, (nx * ny * nz) // 200)
        for idx in range(0, nx * ny * nz, step):
            i1 = idx // (ny * nz)
            i2 = (idx // nz) % ny
            i3 = idx % nz
            sample_rows.append(
                [
                    sys_spec.name,
                    str(i1),
                    str(i2),
                    str(i3),
                    f"{rho_r_3d[i1, i2, i3]:.10e}",
                ]
            )
    sample_csv = out_csv.with_name("vgch_sad_heavy_samples.csv")
    with sample_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "i1", "i2", "i3", "rho_e_per_ang3"])
        w.writerows(sample_rows)
    print(f"wrote {len(sample_rows)} rows to {sample_csv}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
