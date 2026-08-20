"""Extended Kalman filter (numpy port of the C ``kf_ekf_*`` API).

Nonlinear state-transition and measurement models, linearised with the user's
Jacobians:

    x(k+1) = f(x, u, dt) + w,    z = h(x) + v
    F = df/dx,  H = dh/dx   (evaluated at the current estimate)
"""

from __future__ import annotations

import numpy as np

from .matrix import symmetrize


class ExtendedKalmanFilter:
    """Extended Kalman filter.

    Parameters
    ----------
    n, m : int
        State / measurement dimension.
    f : callable(x, u, dt) -> x_out
        State-transition function.
    F_jac : callable(x, u, dt) -> F (n, n)
        Jacobian of f w.r.t. the state.
    h : callable(x) -> z_out
        Measurement function.
    H_jac : callable(x) -> H (m, n)
        Jacobian of h w.r.t. the state.
    """

    def __init__(self, n, m, f=None, F_jac=None, h=None, H_jac=None):
        self.n = int(n)
        self.m = int(m)
        self.x = np.zeros(self.n)
        self.P = np.eye(self.n)
        self.Q = np.zeros((self.n, self.n))
        self.R = np.zeros((self.m, self.m))

        self.f = f
        self.F_jac = F_jac
        self.h = h
        self.H_jac = H_jac

        self._y = np.zeros(self.m)

    def set_models(self, f, F_jac, h, H_jac) -> None:
        """Register the four nonlinear model callbacks."""
        self.f, self.F_jac, self.h, self.H_jac = f, F_jac, h, H_jac

    # --- state / covariance / noise (same semantics as KalmanFilter) ---
    def set_state(self, x) -> None:
        self.x[:] = np.asarray(x, dtype=float)

    def get_state(self) -> np.ndarray:
        return self.x.copy()

    def set_covariance(self, P) -> None:
        self.P[:] = np.asarray(P, dtype=float)

    def set_covariance_diagonal(self, diag) -> None:
        self.P[:] = np.diag(np.asarray(diag, dtype=float))

    def set_covariance_scalar(self, p) -> None:
        self.P[:] = float(p) * np.eye(self.n)

    def get_covariance(self) -> np.ndarray:
        return self.P.copy()

    def set_process_noise(self, Q) -> None:
        self.Q[:] = np.asarray(Q, dtype=float)

    def set_process_noise_scalar(self, q) -> None:
        self.Q[:] = float(q) * np.eye(self.n)

    def set_measurement_noise(self, R) -> None:
        self.R[:] = np.asarray(R, dtype=float)

    def set_measurement_noise_scalar(self, r) -> None:
        self.R[:] = float(r) * np.eye(self.m)

    # --- filtering ---
    def predict(self, u=None, dt=0.0) -> None:
        if self.f is None or self.F_jac is None:
            raise RuntimeError("models not set; call set_models() first")
        F = np.asarray(self.F_jac(self.x, u, dt), dtype=float)
        if F.shape != (self.n, self.n):
            raise ValueError("F_jac must return (%d,%d)" % (self.n, self.n))
        x_new = np.asarray(self.f(self.x, u, dt), dtype=float)
        self.x[:] = x_new
        self.P[:] = symmetrize(F @ self.P @ F.T + self.Q)

    def update(self, z) -> None:
        if self.h is None or self.H_jac is None:
            raise RuntimeError("models not set; call set_models() first")
        z = np.asarray(z, dtype=float)
        H = np.asarray(self.H_jac(self.x), dtype=float)
        if H.shape != (self.m, self.n):
            raise ValueError("H_jac must return (%d,%d)" % (self.m, self.n))
        z_hat = np.asarray(self.h(self.x), dtype=float)
        y = z - z_hat
        self._y[:] = y

        P, R = self.P, self.R
        S = H @ P @ H.T + R
        K = np.linalg.solve(S, (P @ H.T).T).T
        self.x = self.x + K @ y
        IKH = np.eye(self.n) - K @ H
        self.P[:] = symmetrize(IKH @ P @ IKH.T + K @ R @ K.T)
