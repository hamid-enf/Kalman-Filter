"""Example 3 (advanced) — fuse odometry (velocity) with GNSS (position).

Mirrors examples/stm32_h743/03_sensor_fusion.c. Two sensors at different rates
and qualities, fused via sequential updates with different H and R.
"""

import numpy as np

from kalman import KalmanFilter

TRUE_VEL = 2.0
FAST_DT = 0.01
GNSS_EVERY = 100
rng = 0xAB0BA9C


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
    F = np.array([[1.0, FAST_DT], [0.0, 1.0]])
    H_vel = np.array([[0.0, 1.0]])
    H_pos = np.array([[1.0, 0.0]])
    Q = np.array([[0.0, 0.0], [0.0, 0.05]])
    P0 = np.diag([10.0, 10.0])

    kf = KalmanFilter(2, 1)
    kf.set_transition_matrix(F)
    kf.set_process_noise(Q)
    kf.set_covariance(P0)

    pos = 0.0
    fast_count = 0
    print(" t(s) | odom vel | gnss pos | fused pos | fused vel")
    for step in range(400):
        odom = TRUE_VEL + 0.2 * gauss(rng)
        pos += TRUE_VEL * FAST_DT

        kf.predict()

        # fast update: correct velocity with odometry
        kf.set_measurement_matrix(H_vel)
        kf.set_measurement_noise_scalar(0.04)
        kf.update(np.array([odom]))

        fast_count += 1
        if fast_count >= GNSS_EVERY:
            gnss = pos + 2.0 * gauss(rng)
            fast_count = 0
            kf.set_measurement_matrix(H_pos)
            kf.set_measurement_noise_scalar(4.0)
            kf.update(np.array([gnss]))

            t = step * FAST_DT
            print("%5.2f | %8.2f | %8.2f | %9.2f | %9.2f"
                  % (t, odom, gnss, kf.get_state()[0], kf.get_state()[1]))
    print("\nFinal: fused pos %.2f (true %.2f), vel %.2f (true %.1f)"
          % (kf.get_state()[0], pos, kf.get_state()[1], TRUE_VEL))


if __name__ == "__main__":
    main()
