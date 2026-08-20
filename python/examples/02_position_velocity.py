"""Example 2 (medium) — position + velocity from position-only measurements.

Mirrors examples/stm32_h743/02_position_velocity.c. A cart moves at 1.5 m/s;
we measure only position (sigma 0.3 m). The filter infers the velocity.
"""

import numpy as np

from kalman import KalmanFilter

TRUE_VEL = 1.5
DT = 0.01
rng = 0x5EED5EED


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
    kf = KalmanFilter.constant_velocity(dt=DT, q_accel=0.1, r=0.09,
                                        p0_pos=1.0, p0_vel=1.0)
    pos = 0.0
    print(" t(s) | true pos | measured | filt pos | filt vel")
    for step in range(100):
        pos += TRUE_VEL * DT
        z = pos + 0.3 * gauss(rng)
        kf.predict()
        kf.update(np.array([z]))
        if step % 10 == 0:
            t = step * DT
            print("%5.2f | %8.2f | %8.2f | %8.2f | %8.2f"
                  % (t, pos, z, kf.get_state()[0], kf.get_state()[1]))
    print("\nFinal: pos %.2f (true %.2f), vel %.2f (true %.1f)"
          % (kf.get_state()[0], pos, kf.get_state()[1], TRUE_VEL))


if __name__ == "__main__":
    main()
