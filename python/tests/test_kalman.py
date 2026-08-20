"""Numerical tests for the Python port.

Mirrors the key checks of the C suite: convergence, UKF==KF equivalence on a
linear system, RTS smoother vs. an exact scalar reference, gating, and EKF
nonlinear measurement. Run with:  python -m pytest tests/  (or python tests/test_kalman.py)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from kalman import (KalmanFilter, ExtendedKalmanFilter, UnscentedKalmanFilter,
                    rts_smooth_step)


# deterministic PRNG (xorshift32) so the tests are reproducible
_rng = 0x12345678


def _rand():
    global _rng
    x = _rng
    x ^= (x << 13) & 0xFFFFFFFF
    x ^= x >> 17
    x ^= (x << 5) & 0xFFFFFFFF
    _rng = x & 0xFFFFFFFF
    return x


def _noise():
    return ((int(_rand() % 2000001) - 1000000) / 1000000.0)


def _assert_close(a, b, tol, msg=""):
    assert abs(a - b) <= tol, "%s: %r vs %r (tol %r)" % (msg, a, b, tol)


def test_1d_convergence():
    f = KalmanFilter.constant_signal(q=1e-4, r=1.0)
    for _ in range(500):
        f.predict()
        z = 10.0 + _noise()
        f.update(np.array([z]))
    _assert_close(f.get_state()[0], 10.0, 0.3, "1D convergence")


def test_constant_velocity():
    f = KalmanFilter.constant_velocity(dt=1.0, q_accel=0.05, r=1.0,
                                       p0_pos=100.0, p0_vel=100.0)
    pos = 0.0
    for _ in range(500):
        pos += 2.0
        f.predict()
        z = pos + 0.5 * _noise()
        f.update(np.array([z]))
    _assert_close(f.get_state()[0], pos, 5.0, "CV position")
    _assert_close(f.get_state()[1], 2.0, 0.2, "CV velocity")


def test_ukf_matches_kf_linear():
    """The unscented transform is exact for linear systems, so a UKF with a
    linear model must match the linear KF to high accuracy."""
    F = np.array([[1.0, 1.0], [0.0, 1.0]])
    H = np.array([[1.0, 0.0]])
    Q = np.array([[0.0, 0.0], [0.0, 0.05]])
    P0 = np.diag([10.0, 10.0])

    kf = KalmanFilter(2, 1)
    kf.set_transition_matrix(F)
    kf.set_measurement_matrix(H)
    kf.set_process_noise(Q)
    kf.set_measurement_noise_scalar(1.0)
    kf.set_covariance(P0)

    ukf = UnscentedKalmanFilter(2, 1,
                                f=lambda x, u, dt: F @ x,
                                h=lambda x: H @ x)
    ukf.set_process_noise(Q)
    ukf.set_measurement_noise_scalar(1.0)
    ukf.set_covariance(P0)

    for i in range(200):
        z = 2.0 * (i + 1) + 0.5 * _noise()
        kf.predict()
        kf.update(np.array([z]))
        ukf.predict()
        ukf.update(np.array([z]))

    _assert_close(kf.get_state()[0], ukf.get_state()[0], 0.05, "UKF==KF pos")
    _assert_close(kf.get_state()[1], ukf.get_state()[1], 0.02, "UKF==KF vel")


def test_gating():
    f = KalmanFilter.constant_signal(q=1e-4, r=1.0)
    for _ in range(100):
        f.predict()
        f.update(np.array([10.0 + _noise()]))

    f.set_gate_threshold(6.63)   # ~99% chi-square, 1 dof

    f.predict()
    accepted = f.update_gated(np.array([10.0 + _noise()]))
    assert accepted and f.nis < 6.63, "normal reading should be accepted"

    before = f.get_state()[0]
    f.predict()
    accepted = f.update_gated(np.array([1000.0]))
    assert not accepted, "spike should be rejected"
    _assert_close(f.get_state()[0], before, 1e-4, "state unchanged after gate")


def test_adaptive_r():
    f = KalmanFilter.constant_signal(q=1e-4, r=0.05)
    for _ in range(3000):
        f.predict()
        z = 10.0 + math.sqrt(2.0) * _noise()   # variance 2 (uniform var 1/3 *2? )
        f.update(np.array([z]))
        f.adapt_r(gamma=0.995, r_min=0.001)
    # converge toward the true sensor variance (within a loose band)
    r = f.get_measurement_noise()[0, 0]
    assert 0.3 < r < 4.0, "adaptive R drifted: %r" % r


def test_rts_smoother_matches_reference():
    """1-D RTS smoother vs. the exact scalar recursion C = P_k / P_{k+1|k}."""
    N = 60
    f = KalmanFilter.constant_signal(q=0.01, r=1.0)
    xf = np.zeros(N + 1); Pf = np.zeros(N + 1)
    xp = np.zeros(N + 1); Pp = np.zeros(N + 1)
    xf[0], Pf[0] = 0.0, 1.0
    for k in range(N):
        f.predict()
        xp[k + 1] = f.get_state()[0]
        Pp[k + 1] = f.get_covariance()[0, 0]
        f.update(np.array([_noise()]))
        xf[k + 1] = f.get_state()[0]
        Pf[k + 1] = f.get_covariance()[0, 0]

    xs = np.zeros(N + 1); Ps = np.zeros(N + 1)
    xs[N], Ps[N] = xf[N], Pf[N]
    for k in range(N - 1, -1, -1):
        x, P = rts_smooth_step(1, np.array([xf[k]]), np.array([[Pf[k]]]),
                               np.array([[1.0]]),
                               np.array([xp[k + 1]]), np.array([[Pp[k + 1]]]),
                               np.array([xs[k + 1]]), np.array([[Ps[k + 1]]]))
        xs[k], Ps[k] = x[0], P[0, 0]
    # scalar reference
    for k in range(N - 1, -1, -1):
        C = Pf[k] / Pp[k + 1]
        xr = xf[k] + C * (xs[k + 1] - xp[k + 1])
        Pr = Pf[k] + C * C * (Ps[k + 1] - Pp[k + 1])
        _assert_close(xs[k], xr, 1e-4, "RTS x")
        _assert_close(Ps[k], Pr, 1e-4, "RTS P")


def test_ekf_nonlinear_measurement():
    ekf = ExtendedKalmanFilter(1, 1,
                               f=lambda x, u, dt: x,
                               F_jac=lambda x, u, dt: np.array([[1.0]]),
                               h=lambda x: np.array([x[0] ** 2]),
                               H_jac=lambda x: np.array([[2.0 * x[0]]]))
    ekf.set_state(np.array([2.0]))
    ekf.set_process_noise_scalar(1e-4)
    ekf.set_measurement_noise_scalar(1.0)
    ekf.set_covariance_scalar(1.0)
    for _ in range(300):
        ekf.predict()
        ekf.update(np.array([9.0 + 2.0 * _noise()]))
    _assert_close(ekf.get_state()[0], 3.0, 0.15, "EKF h=x^2")


def test_ukf_nonlinear_measurement():
    ukf = UnscentedKalmanFilter(1, 1,
                                f=lambda x, u, dt: x,
                                h=lambda x: np.array([x[0] ** 2]))
    ukf.set_state(np.array([2.0]))
    ukf.set_process_noise_scalar(1e-4)
    ukf.set_measurement_noise_scalar(1.0)
    ukf.set_covariance_scalar(1.0)
    for _ in range(300):
        ukf.predict()
        ukf.update(np.array([9.0 + 2.0 * _noise()]))
    _assert_close(ukf.get_state()[0], 3.0, 0.15, "UKF h=x^2")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print("  PASS  %s" % t.__name__)
        except AssertionError as e:
            failed += 1
            print("  FAIL  %s: %s" % (t.__name__, e))
    print("\n%s (%d failed)" % ("ALL TESTS PASSED" if failed == 0 else "FAILED", failed))
    return failed


if __name__ == "__main__":
    sys.exit(main())
