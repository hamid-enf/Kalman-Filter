"""Example 4 (nonlinear) — EKF: estimate an angle from sin(angle) readings.

Mirrors the C example 09. The measurement is nonlinear (h(x) = sin(x)), so the
linear KF does not apply; the EKF linearises with the Jacobian H = cos(x).
"""

import math

import numpy as np

from kalman import ExtendedKalmanFilter

TRUE = 0.8
rng = 9


def rand(state):
    x = state
    x ^= (x << 13) & 0xFFFFFFFF
    x ^= x >> 17
    x ^= (x << 5) & 0xFFFFFFFF
    return x & 0xFFFFFFFF


def gauss(state):
    n = lambda: ((int(rand(state) % 2000001) - 1000000) / 1000000.0)
    return (n() + n() + n()) / 3.0


def main():
    ekf = ExtendedKalmanFilter(
        1, 1,
        f=lambda x, u, dt: x,
        F_jac=lambda x, u, dt: np.array([[1.0]]),
        h=lambda x: np.array([math.sin(x[0])]),
        H_jac=lambda x: np.array([[math.cos(x[0])]]))
    ekf.set_state(np.array([0.3]))
    ekf.set_process_noise_scalar(1e-4)
    ekf.set_measurement_noise_scalar(0.01)
    ekf.set_covariance_scalar(0.5)

    print("step | measurement | estimate")
    for i in range(40):
        ekf.predict()
        z = math.sin(TRUE) + 0.1 * gauss(rng)
        ekf.update(np.array([z]))
        if i % 5 == 0:
            print("%4d | %11.3f | %8.3f" % (i, z, ekf.get_state()[0]))
    print("\nEstimated angle: %.3f rad (true %.3f)" % (ekf.get_state()[0], TRUE))


if __name__ == "__main__":
    main()
