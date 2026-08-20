# Kalman filtering — Python (numpy) port

A Python port of the [C library](../README.md), with **identical formulas and
numerical strategy** (Joseph-form covariance update, Cholesky-based solves
instead of naive inversions, covariance symmetrisation, innovation gating,
adaptive R, RTS smoother).

Use this for prototyping and validating your filter on the PC before porting
the same model to the STM32 — the semantics match 1:1.

## Install / run

```sh
cd python
pip install -e .          # install the `kalman` package
python tests/test_kalman.py
PYTHONPATH=. python examples/01_basic_1d.py
```

Only dependency: **numpy** (see `requirements.txt`).

## Quick start (1-D signal)

```python
import numpy as np
from kalman import KalmanFilter

kf = KalmanFilter.constant_signal(q=1e-3, r=0.25)   # F=[1], H=[1], P=1

for z in noisy_readings:                            # z: scalar
    kf.predict()
    kf.update(np.array([z]))
    x = kf.get_state()[0]                           # filtered estimate
```

## API (mirrors the C `kf_*` API)

| C | Python |
|---|--------|
| `kf_kf_init(n, m)` | `KalmanFilter(n, m)` |
| `kf_kf_init_1d / _constant_velocity` | `KalmanFilter.constant_signal / .constant_velocity` |
| `kf_kf_set_state / get_state` | `set_state / get_state` |
| `kf_kf_set_covariance* / get_covariance` | `set_covariance* / get_covariance` |
| `kf_kf_set_process_noise* / set_measurement_noise*` | `set_process_noise* / set_measurement_noise*` |
| `kf_kf_set_transition_matrix / set_measurement_matrix` | `set_transition_matrix / set_measurement_matrix` |
| `kf_kf_predict / update` | `predict / update` |
| `kf_kf_set_gate_threshold / nis / update_gated` | `set_gate_threshold / nis / update_gated` |
| `kf_kf_adapt_r` | `adapt_r` |
| `kf_kf_smooth_step` | `rts_smooth_step` (plus `rts_smoother` for a full pass) |
| `kf_ekf_*` | `ExtendedKalmanFilter` (callables f, F_jac, h, H_jac) |
| `kf_ukf_*` | `UnscentedKalmanFilter` (callables f, h) |

## Differences from the C core (intentional)

- **Errors are exceptions**, not status codes: bad dimensions raise
  `ValueError`; a singular innovation covariance raises
  `numpy.linalg.LinAlgError`. `update_gated` returns a bool (True = accepted,
  False = rejected outlier) instead of a `KF_WARN_GATED` status.
- **numpy is used** for matrix arithmetic (the idiomatic host choice); the C
  core remains dependency-free. The algorithms are the same.

## Examples

| File | Content |
|------|---------|
| `01_basic_1d.py` | 1-D temperature filter with spike rejection (simple) |
| `02_position_velocity.py` | constant-velocity; velocity inferred from position (medium) |
| `03_sensor_fusion.py` | odometry (velocity) + GNSS (position) fusion (advanced) |
| `04_ekf.py` | EKF: estimate an angle from sin(angle) |
| `05_ukf.py` | UKF: same problem, no Jacobians |

## Tests

`tests/test_kalman.py` mirrors the C suite's key numerical checks
(convergence, UKF≡KF equivalence on a linear system, RTS smoother vs an exact
scalar reference, gating, adaptive-R convergence, EKF/UKF nonlinear
measurement). All pass.
