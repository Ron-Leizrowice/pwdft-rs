#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy>=1.26",
# ]
# ///
"""
VGCH-2 Part B — parse QE's converged charge density.

QE 7.5 writes the self-consistent charge density to
`<prefix>.save/charge-density.dat` (binary Fortran sequential-access
file) or `charge-density.hdf5` when QE is built with HDF5. This
script parses the `.dat` variant (the default when QE was built
without HDF5, as at the time of this change) and emits a portable
`.npz` archive suitable for the Rust transplant harness to consume.

Record layout (from `qe-7.5/Modules/io_base.f90:482-590`):

    rec 1 : gamma_only (logical), ngm_g (int32), nspin (int32)
    rec 2 : b1(3), b2(3), b3(3) — reciprocal vectors in a.u. (1/Bohr)
    rec 3 : mill_g(3, ngm_g)    — Miller indices, int32
    rec 4 : rho_g(ngm_g)        — complex128, spin 1
    [rec 5: rho_g(ngm_g)        — complex128, spin 2, if nspin == 2]

Each Fortran sequential record is preceded and followed by a 4-byte
length marker (standard ifort/gfortran convention).

Output is a flat little-endian binary bundle (`.bin`) with a fixed
header + byte-packed payloads. Rust-side parsing avoids a ZIP dep.

Layout:
    magic (8 bytes, ASCII): "VGCH2BIN"
    version: u32 = 1
    gamma_only: u8 (0 or 1)
    nspin: i32
    ngm: i32
    b1: 3 × f64
    b2: 3 × f64
    b3: 3 × f64
    mill: ngm × 3 × i32  (C-order: mill[ig, 0..3])
    rho_g: ngm × nspin × (f64, f64)  complex128, e/Bohr³
           (convention: ρ(G) = (1/Ω) ∫ ρ(r) exp(-i G·r) d³r,
            QE's normalization, matches pwdft-rs' (1/N) forward-FFT)

Usage:
    uv run validation/src/pwdft_validation/scripts/vgch2_parse_qe_density.py \\
        --rho path/to/prefix.save/charge-density.dat \\
        --out path/to/cu_rho_qe.bin

Notes:
 - Does NOT convert to e/Å³; downstream consumers apply the
   BOHR_TO_ANG⁻³ factor themselves so the boundary stays explicit.
 - Does NOT map onto the pwdft-rs FFT grid; the consumer indexes
   pwdft-rs' grid via Miller indices directly (mill[i,0],
   mill[i,1], mill[i,2]).
 - Assumes gamma_only == False (full G-sphere written). VGCH-2
   heavy-atom systems all use Γ-centered k-grids with multiple k-points;
   gamma_only would be rare here. The parser errors out if
   gamma_only == True so the caller knows to re-run QE without it.
"""

from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path

import numpy as np


def _read_record(fp) -> bytes:
    """Read one Fortran sequential-access record (4-byte head + body + 4-byte tail)."""
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
    """Parse QE's `charge-density.dat` and return a dict of numpy arrays.

    Returns keys: gamma_only (bool), nspin (int), ngm (int),
    b1/b2/b3 ((3,) float64, 1/Bohr), mill ((ngm,3) int32),
    rho_g ((ngm,) or (ngm,nspin) complex128, e/Bohr³).
    """
    with open(path, "rb") as fp:
        # rec 1: gamma_only (logical*4), ngm_g (int32), nspin (int32)
        rec1 = _read_record(fp)
        if len(rec1) != 12:
            raise ValueError(f"rec 1 size {len(rec1)} != 12 (expected logical*4 + 2 int32)")
        gamma_only_raw, ngm_g, nspin = struct.unpack("<iii", rec1)
        # Fortran LOGICAL: 0 = .FALSE., any other = .TRUE. (gfortran: 1 = .TRUE.).
        gamma_only = bool(gamma_only_raw)
        if gamma_only:
            raise NotImplementedError(
                "gamma_only=.TRUE. densities are not supported by this parser — "
                "re-run QE with a k-grid (even 1x1x1 unshifted Γ) to write "
                "the full G-sphere"
            )

        # rec 2: b1(3), b2(3), b3(3) — 9 doubles
        rec2 = _read_record(fp)
        if len(rec2) != 72:
            raise ValueError(f"rec 2 size {len(rec2)} != 72 (9 doubles)")
        bg = np.frombuffer(rec2, dtype="<f8").reshape(3, 3)  # rows are b1,b2,b3
        b1, b2, b3 = bg[0], bg[1], bg[2]

        # rec 3: mill_g(3, ngm_g) — Fortran column-major storage. In memory
        # the layout is [m1_g1, m2_g1, m3_g1, m1_g2, m2_g2, m3_g2, …] because
        # the first Fortran dimension (the coordinate axis, length 3) varies
        # fastest. Reading this as a C-order (ngm_g, 3) array gives
        # `mill[ig, :] = [m1, m2, m3]` for each g-vector directly. NOTE that
        # `order="F"` here would give an *incorrect* mapping (swaps mill[ig,:]
        # columns across g-vectors); verified against the G=0 sanity check
        # `rho(G=0) ≈ N_el/Ω = 19/79.38 = 0.2394 e/Bohr³` for Cu FCC, which
        # only lands at `mill == [0, 0, 0]` under the C-order interpretation.
        rec3 = _read_record(fp)
        if len(rec3) != 3 * ngm_g * 4:
            raise ValueError(f"rec 3 size {len(rec3)} != 3*ngm_g*4 = {3 * ngm_g * 4}")
        mill = np.frombuffer(rec3, dtype="<i4").reshape(ngm_g, 3, order="C").copy()

        # rec 4 [+ rec 5]: rho_g(ngm_g) as complex*16 (float64 real + float64 imag)
        rho_g_list = []
        for ispin in range(nspin):
            recn = _read_record(fp)
            if len(recn) != ngm_g * 16:
                raise ValueError(f"rho rec {ispin + 1} size {len(recn)} != ngm_g*16 = {ngm_g * 16}")
            arr = np.frombuffer(recn, dtype="<c16").copy()
            rho_g_list.append(arr)
        rho_g = np.stack(rho_g_list, axis=-1) if nspin > 1 else rho_g_list[0]

    return {
        "gamma_only": gamma_only,
        "nspin": int(nspin),
        "ngm": int(ngm_g),
        "b1": b1.copy(),
        "b2": b2.copy(),
        "b3": b3.copy(),
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
    # Find G=0 entry and print rho(G=0) — should equal N_el / Ω in e/Bohr³.
    g0_mask = np.all(d["mill"] == 0, axis=1)
    g0_idx = np.where(g0_mask)[0]
    if len(g0_idx) == 1:
        rho0 = d["rho_g"][g0_idx[0]] if d["nspin"] == 1 else d["rho_g"][g0_idx[0], 0]
        print(f"rho(G=0)    = {rho0.real:.8e} + {rho0.imag:.2e}j  (e/Bohr^3)")
    else:
        print(f"rho(G=0) not found uniquely: {len(g0_idx)} matches")
    # rho(G) integrated charge (sum over shell 0 components times Ω would recover
    # total charge, but we don't have Ω here — only print the magnitude distribution).
    mag = np.abs(d["rho_g"] if d["nspin"] == 1 else d["rho_g"][:, 0])
    print(f"|rho_g| min/median/max = {mag.min():.3e} / {np.median(mag):.3e} / {mag.max():.3e}")


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description="Parse QE charge-density.dat")
    p.add_argument("--rho", required=True, type=Path, help="path to charge-density.dat")
    p.add_argument(
        "--out",
        required=True,
        type=Path,
        help="path to output .bin (VGCH2BIN flat binary bundle)",
    )
    p.add_argument("--verbose", action="store_true", help="print parsed summary to stdout")
    args = p.parse_args(argv)

    if not args.rho.is_file():
        print(f"error: {args.rho} not found", file=sys.stderr)
        return 1

    d = parse_charge_density(args.rho)
    if args.verbose:
        _sanity_print(d)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    # Flat binary bundle, little-endian, self-describing header.
    # See module docstring § Layout.
    with open(args.out, "wb") as fp:
        fp.write(b"VGCH2BIN")  # magic
        fp.write(struct.pack("<I", 1))  # version
        fp.write(struct.pack("<B", 1 if d["gamma_only"] else 0))
        fp.write(struct.pack("<i", int(d["nspin"])))
        fp.write(struct.pack("<i", int(d["ngm"])))
        fp.write(d["b1"].astype("<f8", copy=False).tobytes())
        fp.write(d["b2"].astype("<f8", copy=False).tobytes())
        fp.write(d["b3"].astype("<f8", copy=False).tobytes())
        # mill is (ngm, 3) int32 C-order
        fp.write(np.ascontiguousarray(d["mill"], dtype="<i4").tobytes())
        # rho_g is (ngm,) or (ngm, nspin) complex128 — unit nspin writes as
        # a single block of ngm complex128.
        rho_flat = np.ascontiguousarray(d["rho_g"], dtype="<c16")
        fp.write(rho_flat.tobytes())
    print(f"wrote {args.out} (ngm={d['ngm']}, nspin={d['nspin']})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
