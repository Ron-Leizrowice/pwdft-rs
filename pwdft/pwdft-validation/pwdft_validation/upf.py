"""Unified UPF v2 parser.

Single regex-based parser covering all blocks needed by validation scripts.
No external UPF library required; handles ONCVPSP v2 format as produced by
QE 7.5's pseudopotential generation tools.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Low-level XML helpers
# ---------------------------------------------------------------------------


def _extract_attr(text: str, name: str) -> str:
    m = re.search(rf'{re.escape(name)}\s*=\s*"([^"]*)"', text)
    if m is None:
        raise ValueError(f"attribute {name!r} not found")
    return m.group(1).strip()


def _extract_block(text: str, tag: str) -> np.ndarray:
    m = re.search(rf"<{re.escape(tag)}\b[^>]*>(.*?)</{re.escape(tag)}>", text, re.DOTALL)
    if m is None:
        raise ValueError(f"tag <{tag}> not found")
    return np.asarray([float(t) for t in m.group(1).split() if t.strip()], dtype=np.float64)


def _extract_beta(text: str, idx: int, mesh: int) -> tuple[int, np.ndarray]:
    tag = f"PP_BETA.{idx}"
    m = re.search(rf"<{re.escape(tag)}\b([^>]*)>(.*?)</{re.escape(tag)}>", text, re.DOTALL)
    if m is None:
        raise ValueError(f"tag <{tag}> not found")
    attrs, body = m.group(1), m.group(2)
    lm = re.search(r'angular_momentum\s*=\s*"(-?\d+)"', attrs)
    if lm is None:
        raise ValueError(f"{tag} has no angular_momentum attribute")
    l = int(lm.group(1))
    values = np.asarray([float(t) for t in body.split() if t.strip()], dtype=np.float64)
    if values.size != mesh:
        raise ValueError(f"{tag}: expected {mesh} values, got {values.size}")
    return l, values


# ---------------------------------------------------------------------------
# UpfData: parsed UPF v2 file
# ---------------------------------------------------------------------------


@dataclass
class UpfData:
    """All blocks from a UPF v2 pseudopotential file.

    Fields are None when the block is absent (e.g. PP_NLCC on a PP without
    core correction, PP_LOCAL on a PAW PP, etc.).
    """

    element: str
    z_valence: float
    mesh_size: int
    has_nlcc: bool

    r_bohr: np.ndarray
    rab_bohr: np.ndarray

    v_local_ry: np.ndarray | None = field(default=None)
    projectors: list[tuple[int, np.ndarray]] = field(default_factory=list)
    dij_ry: np.ndarray | None = field(default=None)
    rho_at_ebohr: np.ndarray | None = field(default=None)
    rho_core_ebohr3: np.ndarray | None = field(default=None)

    @property
    def n_proj(self) -> int:
        return len(self.projectors)


def parse_upf(path: Path) -> UpfData:
    """Parse all available blocks from *path* (UPF v2 XML format)."""
    text = path.read_text()

    m_hdr = re.search(r"<PP_HEADER\b([^>]*)/?>", text)
    if m_hdr is None:
        raise ValueError(f"{path}: PP_HEADER not found")
    hdr = m_hdr.group(1)

    element = _extract_attr(hdr, "element") if "element" in hdr else path.stem
    z_val = float(_extract_attr(hdr, "z_valence"))
    mesh = int(_extract_attr(hdr, "mesh_size"))
    n_proj_hdr = int(_extract_attr(hdr, "number_of_proj")) if "number_of_proj" in hdr else 0
    has_nlcc = _extract_attr(hdr, "core_correction").lower() == "t" if "core_correction" in hdr else False

    r = _extract_block(text, "PP_R")
    rab = _extract_block(text, "PP_RAB")
    for name, arr in (("PP_R", r), ("PP_RAB", rab)):
        if arr.size != mesh:
            raise ValueError(f"{path}: {name}: expected {mesh} values, got {arr.size}")

    # PP_LOCAL
    v_local: np.ndarray | None = None
    try:
        v_local = _extract_block(text, "PP_LOCAL")
        if v_local.size != mesh:
            raise ValueError(f"{path}: PP_LOCAL: expected {mesh} values, got {v_local.size}")
    except ValueError:
        pass

    # PP_BETA projectors
    projectors: list[tuple[int, np.ndarray]] = []
    for i in range(1, n_proj_hdr + 1):
        l, chi = _extract_beta(text, i, mesh)
        projectors.append((l, chi))

    # PP_DIJ
    dij: np.ndarray | None = None
    if n_proj_hdr > 0:
        try:
            dij_flat = _extract_block(text, "PP_DIJ")
            if dij_flat.size == n_proj_hdr * n_proj_hdr:
                dij = dij_flat.reshape(n_proj_hdr, n_proj_hdr)
        except ValueError:
            pass

    # PP_RHOATOM
    rho_at: np.ndarray | None = None
    try:
        rho_at = _extract_block(text, "PP_RHOATOM")
        if rho_at.size != mesh:
            rho_at = None
    except ValueError:
        pass

    # PP_NLCC
    rho_core: np.ndarray | None = None
    if has_nlcc:
        try:
            rho_core = _extract_block(text, "PP_NLCC")
            if rho_core.size != mesh:
                rho_core = None
        except ValueError:
            pass

    return UpfData(
        element=element,
        z_valence=z_val,
        mesh_size=mesh,
        has_nlcc=has_nlcc,
        r_bohr=r,
        rab_bohr=rab,
        v_local_ry=v_local,
        projectors=projectors,
        dij_ry=dij,
        rho_at_ebohr=rho_at,
        rho_core_ebohr3=rho_core,
    )
