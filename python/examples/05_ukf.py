"""Example 5 (nonlinear) — UKF: same sin(angle) problem, no Jacobians.

Mirrors the C example 10. Only the raw nonlinear functions f and h are needed;
the UKF samples sigma points instead of linearising.
"""

import math

import numpy as np

from kalman import UnscentedKalmanFilter

TRUE = 0.8
rng = 10


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
    ukf = UnscentedKalmanFilter(1, 1,
                                f=lambda x, u, dt: x,
                                h=lambda x: np.array([math.sin(x[0])]))
    ukf.set_state(np.array([0.3]))
    ukf.set_process_noise_scalar(1e-4)
    ukf.set_measurement_noise_scalar(0.01)
    ukf.set_covariance_scalar(0.5)

    print("step | measurement | estimate")
    for i in range(40):
        ukf.predict()
        z = math.sin(TRUE) + 0.1 * gauss(rng)
        ukf.update(np.array([z]))
        if i % 5 == 0:
            print("%4d | %11.3f | %8.3f" % (i, z, ukf.get_state()[0]))
    print("\nEstimated angle: %.3f rad (true %.3f)" % (ukf.get_state()[0], TRUE))


if __name__ == "__main__":
    main()
