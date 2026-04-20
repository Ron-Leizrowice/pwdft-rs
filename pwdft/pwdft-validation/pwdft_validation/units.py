"""Physical constants and unit-conversion factors.

Single source of truth for all scripts. Native units follow QE convention
(Ry, Bohr) at the UPF boundary; pwdft-rs internal units are eV and Å.
"""

from __future__ import annotations

BOHR_TO_ANG: float = 0.529_177_210_903
BOHR3_TO_ANG3: float = BOHR_TO_ANG**3
RY_TO_EV: float = 13.605_693_122_994
# e² in Rydberg atomic units: e²/r gives Ry when r in Bohr (e²=2 Ry·Bohr).
E2_RY_BOHR: float = 2.0
