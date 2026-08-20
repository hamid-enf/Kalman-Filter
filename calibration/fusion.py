"""Cross-language calibration scenario — Python implementation.

See README.md for the spec. Identical xorshift32 + noise to fusion.c so the
printed numbers match the C reference exactly.
"""

import numpy as np

from kalman import KalmanFilter

TRUE_VEL = 2.0
DT = 0.01
GNSS_EVERY = 100

_rng = 0xC0FFEE


def xorshift32():
    global _rng
    x = _rng
    x ^= (x << 13) & 0xFFFFFFFF
    x ^= x >> 17
    x ^= (x << 5) & 0xFFFFFFFF
    _rng = x & 0xFFFFFFFF
    return _rng


def uniform():
    return (int(xorshift32() % 2000001) - 1000000) / 1000000.0


def gauss():
    return (uniform() + uniform() + uniform()) / 3.0


def main():
    F = np.array([[1.0, DT], [0.0, 1.0]])
    Hvel = np.array([[0.0, 1.0]])
    Hpos = np.array([[1.0, 0.0]])
    Q = np.array([[0.0, 0.0], [0.0, 0.05]])
    P0 = np.diag([10.0, 10.0])

    kf = KalmanFilter(2, 1)
    kf.set_transition_matrix(F)
    kf.set_process_noise(Q)
    kf.set_covariance(P0)

    true_pos = 0.0
    fast = 0

    print("  t(s) | fused pos | fused vel")
    for step in range(400):
        odom = TRUE_VEL + 0.2 * gauss()
        true_pos += TRUE_VEL * DT

        kf.predict()

        kf.set_measurement_matrix(Hvel)
        kf.set_measurement_noise_scalar(0.04)
        kf.update(np.array([odom]))

        fast += 1
        if fast >= GNSS_EVERY:
            gnss = true_pos + 2.0 * gauss()
            fast = 0
            kf.set_measurement_matrix(Hpos)
            kf.set_measurement_noise_scalar(4.0)
            kf.update(np.array([gnss]))
            x = kf.get_state()
            print(" %5.2f | %9.5f | %9.5f" % (step * DT, x[0], x[1]))


if __name__ == "__main__":
    main()
