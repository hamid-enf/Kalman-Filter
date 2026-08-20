"""Kalman filtering library — Python (numpy) port of the C core.

Provides the same three estimators as the C library, with identical formulas:

    from kalman import KalmanFilter, ExtendedKalmanFilter, UnscentedKalmanFilter

plus the KF extensions (innovation gating, adaptive R, RTS smoother) and the
``rts_smooth_step`` / ``rts_smoother`` helpers.

The Python port uses numpy for matrix arithmetic; the numerical strategy
(Cholesky solves instead of inversions where sensible, Joseph-form update,
symmetrisation) is kept identical to the dependency-free C core.
"""

from .kf import KalmanFilter, rts_smooth_step, rts_smoother
from .ekf import ExtendedKalmanFilter
from .ukf import UnscentedKalmanFilter

__all__ = [
    "KalmanFilter",
    "ExtendedKalmanFilter",
    "UnscentedKalmanFilter",
    "rts_smooth_step",
    "rts_smoother",
]

__version__ = "1.0.0"
