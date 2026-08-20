# Kalman filtering — MATLAB port

A MATLAB port of the [C library](../README.md), with identical formulas and
numerical strategy (Joseph-form update, Cholesky-based solves, symmetrisation,
innovation gating, adaptive R, RTS smoother).

> **Honesty note:** this port was written to mirror the tested C and Python
> implementations exactly, but the MATLAB files were **reviewed by inspection,
> not executed** (no MATLAB/Octave runtime was available in the build
> environment). Run `tests/test_kalman.m` in MATLAB/Octave to verify on your
> side.

## Usage

Add the package folder to the path and use the `kalman.` classes:

```matlab
addpath('matlab');           % makes +kalman/ visible

kf = kalman.KF(1, 1);
kf.set_constant_signal(1e-3, 0.25);   % F=[1], H=[1], P=1

for z = noisy_readings
    kf.predict();
    kf.update(z);
    x = kf.get_state();
end
```

## API (mirrors the C `kf_*` API)

| C | MATLAB |
|---|--------|
| `kf_kf_init(n, m)` | `kalman.KF(n, m)` |
| `kf_kf_init_1d / _constant_velocity` | `set_constant_signal / set_constant_velocity` |
| `kf_kf_set_state / get_state` | `set_state / get_state` |
| `kf_kf_set_covariance* / get_covariance` | `set_covariance* / get_covariance` |
| `kf_kf_set_process_noise* / set_measurement_noise*` | `set_process_noise* / set_measurement_noise*` |
| `kf_kf_set_transition_matrix / set_measurement_matrix` | `set_transition_matrix / set_measurement_matrix` |
| `kf_kf_predict / update` | `predict / update` |
| `kf_kf_set_gate_threshold / nis / update_gated` | `set_gate_threshold / nis / update_gated` |
| `kf_kf_adapt_r` | `adapt_r` |
| `kf_kf_smooth_step` | `kalman.KF.rts_smooth_step(...)` (static) |
| `kf_ekf_*` | `kalman.EKF` (function handles f, F_jac, h, H_jac) |
| `kf_ukf_*` | `kalman.UKF` (function handles f, h) |

## Differences from the C core (intentional)

- **Errors are MATLAB errors** (`error(...)`), not status codes; `update_gated`
  returns a logical (true = accepted, false = rejected).
- **MATLAB `\` and `chol`** are used for the solves / Cholesky factorisation —
  these are the idiomatic, numerically-stable MATLAB equivalents of the C
  `kf_matrix_cholesky` + `kf_matrix_cholesky_solve`.

## Examples

| File | Content |
|------|---------|
| `example01_basic_1d.m` | 1-D temperature filter with spike rejection (simple) |
| `example02_position_velocity.m` | constant velocity (medium) |
| `example03_sensor_fusion.m` | odometry + GNSS fusion (advanced) |
| `example04_ekf.m` | EKF: angle from sin(angle) |
| `example05_ukf.m` | UKF: same problem, no Jacobians |

## Tests

`tests/test_kalman.m` mirrors the C/Python suite's key numerical checks
(convergence, UKF≡KF equivalence, RTS smoother vs an exact scalar reference,
gating, adaptive-R, EKF/UKF nonlinear). Run it with `addpath('..'); test_kalman()`.
