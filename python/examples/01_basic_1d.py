"""Example 1 (simple) — filter a noisy 1-D signal.

Mirrors examples/stm32_h743/01_simple_temperature.c. A constant signal of 25.0
is read by a noisy sensor (sigma 0.5) that also spikes +20 every 15th sample.
The filter cleans the noise and the innovation gate rejects the spikes.
"""

import numpy as np

from kalman import KalmanFilter

TRUE = 25.0
rng = 0x1A2B3C4D


def rand(state):
    x = state
    x ^= (x << 13) & 0xFFFFFFFF
    x ^= x >> 17
    x ^= (x << 5) & 0xFFFFFFFF
    return x & 0xFFFFFFFF


def noise(state):
    return ((int(rand(state) % 2000001) - 1000000) / 1000000.0)


def gauss(state):
    return (noise(state) + noise(state) + noise(state)) / 3.0


def read_temperature(state, step):
    z = TRUE + 0.5 * gauss(state)
    if step % 15 == 5:
        z += 20.0
    return z


def main():
    kf = KalmanFilter.constant_signal(q=1e-3, r=0.25)

    # warm-up WITHOUT gating (the estimate starts at 0, far from 25)
    for step in range(20):
        z = read_temperature(rng, step)
        kf.predict()
        kf.update(np.array([z]))

    kf.set_gate_threshold(6.63)   # ~99% chi-square, 1 dof

    print("step |  raw  | filtered | note")
    for step in range(20, 45):
        z = read_temperature(rng, step)
        kf.predict()
        accepted = kf.update_gated(np.array([z]))
        note = "" if accepted else " <-- spike rejected"
        print("%4d | %5.2f | %7.2f |%s" % (step, z, kf.get_state()[0], note))
    print("\nFinal estimate: %.2f (true %.1f)" % (kf.get_state()[0], TRUE))


if __name__ == "__main__":
    main()
