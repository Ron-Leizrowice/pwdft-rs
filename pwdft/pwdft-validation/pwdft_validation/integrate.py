"""Quadrature helpers matching QE and pwdft-rs conventions."""

from __future__ import annotations

import numpy as np


def simpson_qe(func: np.ndarray, rab: np.ndarray) -> float:
    """Simpson integrator byte-identical to QE ``simpson(mesh, func, rab, asum)``.

    Weights c_i ∈ {1/3, 4/3, 2/3} with the even-mesh correction from
    ``qe-7.5/upflib/simpsn.f90``. Matches pwdft-rs ``numerics::simpson_integrate``.
    """
    n = len(func)
    assert n == len(rab), "func and rab must have the same length"
    if n < 3:
        return float(np.sum(func * rab))

    acc = 0.0
    for i in range(1, n - 1):
        weight = 4.0 if (i % 2 == 1) else 2.0
        acc += weight * func[i] * rab[i]

    if n % 2 == 1:
        return (acc + func[0] * rab[0] + func[-1] * rab[-1]) / 3.0

    # Even-mesh correction (DFTK/QE formula).
    acc += func[0] * rab[0]
    acc += -0.25 * func[n - 3] * rab[n - 3]
    acc += func[n - 2] * rab[n - 2]
    acc += 1.25 * func[n - 1] * rab[n - 1]
    return acc / 3.0
