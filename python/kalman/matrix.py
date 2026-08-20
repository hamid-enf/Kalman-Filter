"""Minimal numerical helpers shared by the KF/EKF/UKF ports.

The C library is dependency-free; the Python port uses numpy for the matrix
arithmetic (the idiomatic, high-performance choice on a host). The *numerical
strategy* is identical to the C core:

  * no explicit matrix inversions where a linear solve will do
    (``np.linalg.solve`` / ``np.linalg.cholesky``),
  * Joseph-form covariance update,
  * covariance symmetrisation after every step.
"""

from __future__ import annotations

import numpy as np


def symmetrize(a: np.ndarray) -> np.ndarray:
    """Return (a + a.T) / 2 in place, suppressing accumulated rounding error."""
    a[:] = 0.5 * (a + a.T)
    return a


def is_finite(a: np.ndarray) -> bool:
    """True if every element of ``a`` is finite."""
    return bool(np.all(np.isfinite(a)))
