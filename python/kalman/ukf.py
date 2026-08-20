"""Unscented Kalman filter (numpy port of the C ``kf_ukf_*`` API).

Nonlinear models without Jacobians: the distribution is sampled at 2n+1 sigma
points (from the Cholesky factor of P), propagated through f/h, and the mean
and covariance reconstructed from the results.

Sigma-point parameters: alpha (spread), beta (distribution prior), kappa
(secondary scaling). Defaults alpha=1, beta=2, kappa=0 are robust
(all-positive weights).
"""

from __future__ import annotations

import numpy as np

from .matrix import symmetrize


class UnscentedKalmanFilter:
    """Unscented Kalman filter.

    Parameters
    ----------
    n, m : int
        State / measurement dimension.
    f : callable(x, u, dt) -> x_out
        State-transition function.
    h : callable(x) -> z_out
        Measurement function.
    """

    def __init__(self, n, m, f=None, h=None):
        self.n = int(n)
        self.m = int(m)
        self.x = np.zeros(self.n)
        self.P = np.eye(self.n)
        self.Q = np.zeros((self.n, self.n))
        self.R = np.zeros((self.m, self.m))

        self.f = f
        self.h = h

        self.alpha = 1.0
        self.beta = 2.0
        self.kappa = 0.0
        # lambda = alpha^2 (n + kappa) - n  = 0 for the defaults (matches C)
        self.lambda_ = self.alpha * self.alpha * (self.n + self.kappa) - self.n
        self._recompute_weights()

    def set_models(self, f, h) -> None:
        self.f, self.h = f, h

    def set_parameters(self, alpha, beta, kappa) -> None:
        if alpha <= 0 or beta < 0:
            raise ValueError("alpha must be > 0 and beta >= 0")
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.kappa = float(kappa)
        self.lambda_ = alpha**2 * (self.n + kappa) - self.n
        self._recompute_weights()

    def _recompute_weights(self) -> None:
        n, lam = self.n, self.lambda_
        denom = n + lam
        nsig = 2 * n + 1
        self.Wm = np.full(nsig, 0.5 / denom)
        self.Wc = np.full(nsig, 0.5 / denom)
        self.Wm[0] = lam / denom
        self.Wc[0] = lam / denom + (1.0 - self.alpha**2 + self.beta)

    # --- state / covariance / noise ---
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

    # --- sigma points ---
    def _sigma_points(self):
        n = self.n
        L = np.linalg.cholesky(self.P)      # lower; P = L L^T
        c = np.sqrt(n + self.lambda_)
        sig = np.empty((2 * n + 1, n))
        sig[0] = self.x
        for i in range(n):
            sig[i + 1] = self.x + c * L[:, i]
            sig[i + 1 + n] = self.x - c * L[:, i]
        return sig

    # --- filtering ---
    def predict(self, u=None, dt=0.0) -> None:
        if self.f is None:
            raise RuntimeError("model f not set; call set_models() first")
        sig = self._sigma_points()
        X = np.array([self.f(si, u, dt) for si in sig])
        xm = self.Wm @ X
        d = X - xm
        self.P[:] = symmetrize(d.T @ (self.Wc[:, None] * d) + self.Q)
        self.x[:] = xm

    def update(self, z) -> None:
        if self.h is None:
            raise RuntimeError("model h not set; call set_models() first")
        z = np.asarray(z, dtype=float)
        sig = self._sigma_points()
        Z = np.array([self.h(si) for si in sig])
        zm = self.Wm @ Z

        dz = Z - zm
        S = dz.T @ (self.Wc[:, None] * dz) + self.R

        dx = sig - self.x
        Pxz = dx.T @ (self.Wc[:, None] * dz)     # n x m

        K = Pxz @ np.linalg.inv(S)               # n x m
        self.x = self.x + K @ (z - zm)
        self.P[:] = symmetrize(self.P - K @ Pxz.T)
