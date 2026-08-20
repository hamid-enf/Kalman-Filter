"""Linear Kalman filter (numpy port of the C ``kf_kf_*`` API).

Mirrors the C core exactly (same formulas, same numerically-stable choices):

    predict : x = F x + u ;  P = F P F^T + Q
    update  : y = z - H x ;  S = H P H^T + R ;  K = (P H^T) S^-1 (solved)
              x = x + K y ;  P = (I-KH) P (I-KH)^T + K R K^T  (Joseph form)

Also provides the same optional extensions as the C library: innovation
gating (NIS outlier rejection), adaptive measurement noise (residual-based
covariance matching), and the Rauch-Tung-Striebel smoother.
"""

from __future__ import annotations

import numpy as np

from .matrix import symmetrize, is_finite


class KalmanFilter:
    """Linear Kalman filter.

    Parameters
    ----------
    n : int
        State dimension.
    m : int
        Measurement dimension.
    """

    def __init__(self, n: int, m: int):
        if n <= 0 or m <= 0:
            raise ValueError("dimensions must be positive")
        self.n = int(n)
        self.m = int(m)

        self.x = np.zeros(self.n)                       # state
        self.P = np.eye(self.n)                         # covariance
        self.Q = np.zeros((self.n, self.n))             # process noise
        self.F = np.eye(self.n)                         # transition
        self.R = np.zeros((self.m, self.m))             # measurement noise
        self.H = np.zeros((self.m, self.n))             # measurement

        # gating / adaptive-R state
        self.gate_threshold = 0.0
        self.last_nis = 0.0
        self._y = np.zeros(self.m)                      # innovation
        self._resid = np.zeros(self.m)                  # residual r = z - H x+

    # ------------------------------------------------------------------
    # Convenience constructors (mirror kf_kf_init_1d / _constant_velocity)
    # ------------------------------------------------------------------
    @classmethod
    def constant_signal(cls, q: float, r: float) -> "KalmanFilter":
        """1-D constant-signal filter: F=[1], H=[1], P=1."""
        f = cls(1, 1)
        f.F[0, 0] = 1.0
        f.H[0, 0] = 1.0
        f.Q[0, 0] = float(q)
        f.R[0, 0] = float(r)
        f.P[0, 0] = 1.0
        return f

    @classmethod
    def constant_velocity(cls, dt: float, q_accel: float, r: float,
                          p0_pos: float, p0_vel: float) -> "KalmanFilter":
        """Constant-velocity model x = [pos, vel], measuring position only.

        F = [[1, dt],[0, 1]]; Q is the white-acceleration discrete covariance
        [[dt^4/4, dt^3/2],[dt^3/2, dt^2]] * q_accel.
        """
        f = cls(2, 1)
        dt = float(dt)
        f.F = np.array([[1.0, dt], [0.0, 1.0]])
        f.H = np.array([[1.0, 0.0]])
        f.Q = q_accel * np.array([[dt**4 / 4, dt**3 / 2],
                                  [dt**3 / 2, dt**2]])
        f.R[0, 0] = float(r)
        f.P = np.diag([float(p0_pos), float(p0_vel)])
        return f

    # ------------------------------------------------------------------
    # State / covariance / noise access
    # ------------------------------------------------------------------
    def set_state(self, x) -> None:
        x = np.asarray(x, dtype=float)
        if x.shape != (self.n,):
            raise ValueError("state must have shape (%d,)" % self.n)
        self.x[:] = x

    def get_state(self) -> np.ndarray:
        return self.x.copy()

    def set_covariance(self, P) -> None:
        P = np.asarray(P, dtype=float)
        if P.shape != (self.n, self.n):
            raise ValueError("covariance must be %dx%d" % (self.n, self.n))
        if not is_finite(P):
            raise ValueError("covariance contains non-finite values")
        self.P[:] = P

    def set_covariance_diagonal(self, diag) -> None:
        diag = np.asarray(diag, dtype=float)
        if diag.shape != (self.n,):
            raise ValueError("diagonal must have length %d" % self.n)
        self.P[:] = np.diag(diag)

    def set_covariance_scalar(self, p: float) -> None:
        self.P[:] = float(p) * np.eye(self.n)

    def get_covariance(self) -> np.ndarray:
        return self.P.copy()

    def set_process_noise(self, Q) -> None:
        Q = np.asarray(Q, dtype=float)
        if Q.shape != (self.n, self.n):
            raise ValueError("Q must be %dx%d" % (self.n, self.n))
        self.Q[:] = Q

    def set_process_noise_diagonal(self, diag) -> None:
        diag = np.asarray(diag, dtype=float)
        self.Q[:] = np.diag(diag)

    def set_process_noise_scalar(self, q: float) -> None:
        self.Q[:] = float(q) * np.eye(self.n)

    def get_process_noise(self) -> np.ndarray:
        return self.Q.copy()

    def set_measurement_noise(self, R) -> None:
        R = np.asarray(R, dtype=float)
        if R.shape != (self.m, self.m):
            raise ValueError("R must be %dx%d" % (self.m, self.m))
        self.R[:] = R

    def set_measurement_noise_diagonal(self, diag) -> None:
        diag = np.asarray(diag, dtype=float)
        self.R[:] = np.diag(diag)

    def set_measurement_noise_scalar(self, r: float) -> None:
        self.R[:] = float(r) * np.eye(self.m)

    def get_measurement_noise(self) -> np.ndarray:
        return self.R.copy()

    def set_transition_matrix(self, F) -> None:
        F = np.asarray(F, dtype=float)
        if F.shape != (self.n, self.n):
            raise ValueError("F must be %dx%d" % (self.n, self.n))
        self.F[:] = F

    def set_measurement_matrix(self, H) -> None:
        H = np.asarray(H, dtype=float)
        if H.shape != (self.m, self.n):
            raise ValueError("H must be %dx%d" % (self.m, self.n))
        self.H[:] = H

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------
    def predict(self, u=None, dt: float = 0.0) -> None:
        """Time update: x = F x + u, P = F P F^T + Q.

        ``u`` is an optional control input (n,) already in state space; pass
        None for no input. ``dt`` is unused by the linear KF (it is encoded in
        F/Q) and kept only for API symmetry with the EKF/UKF.
        """
        F, P, Q = self.F, self.P, self.Q
        self.x = F @ self.x
        if u is not None:
            u = np.asarray(u, dtype=float)
            if u.shape != (self.n,):
                raise ValueError("control input must have shape (%d,)" % self.n)
            self.x = self.x + u
        P[:] = symmetrize(F @ P @ F.T + Q)

    def _update(self, z: np.ndarray, gate: bool):
        H, P, R = self.H, self.P, self.R
        z = np.asarray(z, dtype=float)
        if z.shape != (self.m,):
            raise ValueError("measurement must have shape (%d,)" % self.m)

        y = z - H @ self.x                       # innovation
        self._y[:] = y
        S = H @ P @ H.T + R                      # innovation covariance

        if gate and self.gate_threshold > 0.0:
            nis = float(y @ np.linalg.solve(S, y))
            self.last_nis = nis
            if nis > self.gate_threshold:
                return False                      # rejected: state unchanged

        # K = (P H^T) S^-1  (solve, don't invert)
        K = np.linalg.solve(S, (P @ H.T).T).T
        self.x = self.x + K @ y

        # residual r = z - H x+  (used by adaptive R)
        self._resid[:] = z - H @ self.x

        IKH = np.eye(self.n) - K @ H
        P[:] = symmetrize(IKH @ P @ IKH.T + K @ R @ K.T)   # Joseph form
        return True

    def update(self, z) -> None:
        """Measurement update. Raises ``numpy.linalg.LinAlgError`` on a
        singular innovation covariance."""
        self._update(z, gate=False)

    # ------------------------------------------------------------------
    # Extensions: gating, adaptive R, smoother
    # ------------------------------------------------------------------
    def set_gate_threshold(self, chi2: float) -> None:
        """Set the NIS (chi-square) gate; 0 disables gating."""
        if chi2 < 0:
            raise ValueError("gate threshold must be non-negative")
        self.gate_threshold = float(chi2)

    @property
    def nis(self) -> float:
        """Normalized innovation squared of the most recent update."""
        return self.last_nis

    def update_gated(self, z) -> bool:
        """Gated update. Returns True if accepted, False if rejected as an
        outlier (state left unchanged)."""
        return bool(self._update(z, gate=True))

    def adapt_r(self, gamma: float, r_min: float) -> None:
        """Adapt R online (residual-based covariance matching):

            R <- gamma R + (1 - gamma) (r r^T + H P H^T)

        Call after an *accepted* update (never after a gated rejection).
        ``gamma`` in (0,1) is the forgetting factor; ``r_min`` floors the
        diagonal to keep R positive definite.
        """
        if not (0.0 < gamma < 1.0):
            raise ValueError("gamma must be in (0, 1)")
        H, P = self.H, self.P
        s = H @ P @ H.T + np.outer(self._resid, self._resid)
        self.R[:] = gamma * self.R + (1.0 - gamma) * s
        np.fill_diagonal(self.R, np.maximum(np.diag(self.R), r_min))


def rts_smooth_step(n, x_filt, P_filt, F, x_pred, P_pred,
                    x_smooth_next, P_smooth_next):
    """One backward step of the RTS fixed-interval smoother.

    Returns ``(x_smooth, P_smooth)`` for step k given the forward-pass data at
    k (filtered) and the already-smoothed estimate at k+1.
    """
    A = F @ P_filt                    # = F P_k  (P_k symmetric)
    C = np.linalg.solve(P_pred, A).T  # C = P_k F^T P_pred^-1
    d = x_smooth_next - x_pred
    x_smooth = x_filt + C @ d
    D = P_smooth_next - P_pred
    P_smooth = symmetrize(P_filt + C @ D @ C.T)
    return x_smooth, P_smooth


def rts_smoother(x_filt, P_filt, F, x_pred, P_pred):
    """Run a full RTS backward pass over a stored forward trajectory.

    ``x_filt``/``P_filt`` are (N, n) / (N, n, n) filtered estimates,
    ``x_pred``/``P_pred`` the predicted ones, ``F`` the (n, n) transition.
    Returns the smoothed ``(x_smooth, P_smooth)`` arrays. Seed the pass with
    the final filtered estimate (done internally).
    """
    N = x_filt.shape[0]
    xs = x_filt.copy()
    Ps = P_filt.copy()
    for k in range(N - 2, -1, -1):
        xs[k], Ps[k] = rts_smooth_step(
            x_filt.shape[1], x_filt[k], P_filt[k], F,
            x_pred[k + 1], P_pred[k + 1], xs[k + 1], Ps[k + 1])
    return xs, Ps
