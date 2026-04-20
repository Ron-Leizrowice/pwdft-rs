"""Parse QE's converged charge density (VGCH-2 Part B).

Reads Fortran sequential-access ``charge-density.dat`` and emits a portable
flat binary bundle (``VGCH2BIN`` format) for Rust-side consumption.
"""

from __future__ import annotations

import struct
import sys
from pathlib import Path
from typing import BinaryIO

import cyclopts
import numpy as np


def _read_record(fp: BinaryIO) -> bytes:
    head = fp.read(4)
    if len(head) != 4:
        raise EOFError("unexpected EOF reading record head")
    n = struct.unpack("<i", head)[0]
    body = fp.read(n)
    if len(body) != n:
        raise EOFError(f"unexpected EOF reading record body ({len(body)}/{n})")
    tail = fp.read(4)
    if len(tail) != 4:
        raise EOFError("unexpected EOF reading record tail")
    n2 = struct.unpack("<i", tail)[0]
    if n != n2:
        raise ValueError(f"Fortran record length markers disagree: head={n}, tail={n2}")
    return body


def parse_charge_density(path: Path) -> dict:
    """Parse QE ``charge-density.dat`` into numpy arrays.

    Returns: gamma_only, nspin, ngm, b1/b2/b3, mill, rho_g.
    """
    with path.open("rb") as fp:
        rec1 = _read_record(fp)
        if len(rec1) != 12:
            raise ValueError(f"rec 1 size {len(rec1)} != 12")
        gamma_only_raw, ngm_g, nspin = struct.unpack("<iii", rec1)
        gamma_only = bool(gamma_only_raw)
        if gamma_only:
            raise NotImplementedError(
                "gamma_only=.TRUE. densities are not supported — re-run QE with a k-grid to write the full G-sphere"
            )

        rec2 = _read_record(fp)
        if len(rec2) != 72:
            raise ValueError(f"rec 2 size {len(rec2)} != 72")
        bg = np.frombuffer(rec2, dtype="<f8").reshape(3, 3)
        b1, b2, b3 = bg[0].copy(), bg[1].copy(), bg[2].copy()

        rec3 = _read_record(fp)
        if len(rec3) != 3 * ngm_g * 4:
            raise ValueError(f"rec 3 size {len(rec3)} != {3 * ngm_g * 4}")
        mill = np.frombuffer(rec3, dtype="<i4").reshape(ngm_g, 3, order="C").copy()

        rho_g_list = []
        for ispin in range(nspin):
            recn = _read_record(fp)
            if len(recn) != ngm_g * 16:
                raise ValueError(f"rho rec {ispin + 1} size {len(recn)} != {ngm_g * 16}")
            rho_g_list.append(np.frombuffer(recn, dtype="<c16").copy())
        rho_g = np.stack(rho_g_list, axis=-1) if nspin > 1 else rho_g_list[0]

    return {
        "gamma_only": gamma_only,
        "nspin": int(nspin),
        "ngm": int(ngm_g),
        "b1": b1,
        "b2": b2,
        "b3": b3,
        "mill": mill,
        "rho_g": rho_g,
    }


def _sanity_print(d: dict) -> None:
    print(f"gamma_only = {d['gamma_only']}")
    print(f"nspin      = {d['nspin']}")
    print(f"ngm        = {d['ngm']}")
    print(f"b1 (1/Bohr)= {d['b1']}")
    print(f"b2 (1/Bohr)= {d['b2']}")
    print(f"b3 (1/Bohr)= {d['b3']}")
    print(
        f"mill range per-axis: ({d['mill'][:, 0].min():+d}..{d['mill'][:, 0].max():+d}) "
        f"({d['mill'][:, 1].min():+d}..{d['mill'][:, 1].max():+d}) "
        f"({d['mill'][:, 2].min():+d}..{d['mill'][:, 2].max():+d})"
    )
    g0_mask = np.all(d["mill"] == 0, axis=1)
    g0_idx = np.where(g0_mask)[0]
    if len(g0_idx) == 1:
        rho0 = d["rho_g"][g0_idx[0]] if d["nspin"] == 1 else d["rho_g"][g0_idx[0], 0]
        print(f"rho(G=0)    = {rho0.real:.8e} + {rho0.imag:.2e}j  (e/Bohr^3)")
    mag = np.abs(d["rho_g"] if d["nspin"] == 1 else d["rho_g"][:, 0])
    print(f"|rho_g| min/median/max = {mag.min():.3e} / {np.median(mag):.3e} / {mag.max():.3e}")


def write_vgch2bin(d: dict, out: Path) -> None:
    """Write parsed density to ``VGCH2BIN`` flat binary bundle."""
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as fp:
        fp.write(b"VGCH2BIN")
        fp.write(struct.pack("<I", 1))
        fp.write(struct.pack("<B", 1 if d["gamma_only"] else 0))
        fp.write(struct.pack("<i", int(d["nspin"])))
        fp.write(struct.pack("<i", int(d["ngm"])))
        fp.write(d["b1"].astype("<f8", copy=False).tobytes())
        fp.write(d["b2"].astype("<f8", copy=False).tobytes())
        fp.write(d["b3"].astype("<f8", copy=False).tobytes())
        fp.write(np.ascontiguousarray(d["mill"], dtype="<i4").tobytes())
        fp.write(np.ascontiguousarray(d["rho_g"], dtype="<c16").tobytes())


def parse(rho_path: Path, out_path: Path, verbose: bool = False) -> int:
    """Parse ``charge-density.dat`` at *rho_path* and write VGCH2BIN to *out_path*."""
    if not rho_path.is_file():
        print(f"error: {rho_path} not found", file=sys.stderr)
        return 1
    d = parse_charge_density(rho_path)
    if verbose:
        _sanity_print(d)
    write_vgch2bin(d, out_path)
    print(f"wrote {out_path} (ngm={d['ngm']}, nspin={d['nspin']})")
    return 0


standalone_app = cyclopts.App(
    name="pwdft-density-parse",
    help="Parse QE charge-density.dat → VGCH2BIN flat binary bundle.",
)


@standalone_app.default
def _density_parse_cmd(rho: Path, out: Path, *, verbose: bool = False) -> None:
    """Parse charge-density.dat at *rho* and write VGCH2BIN to *out*."""
    raise SystemExit(parse(rho, out, verbose=verbose))


def _standalone_parse() -> None:
    """Console-script entry point for ``pwdft-density-parse``."""
    standalone_app()
