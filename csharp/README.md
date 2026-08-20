# Kalman filtering — C# (.NET) port

A C# port of the [C library](../README.md), with identical formulas and
numerical strategy (Joseph-form update, Cholesky-based solves instead of naive
inversions, symmetrisation, innovation gating, adaptive R, RTS smoother).

> **Honesty note:** this port was written to mirror the tested C and Python
> implementations exactly, but the C# files were **reviewed by inspection, not
> compiled** (no .NET SDK was available in the build environment). Run
> `dotnet run --project Tests` to verify on your side.

## Build & run

```sh
cd csharp
dotnet run --project Examples        # run the examples
dotnet run --project Tests           # run the numerical tests
```

Requires .NET SDK 6.0+.

## Structure

| Project | Contents |
|---------|----------|
| `Kalman/` | the library: `LinAlg`, `KalmanFilter`, `ExtendedKalmanFilter`, `UnscentedKalmanFilter` |
| `Examples/` | 3 examples (1-D, constant-velocity, sensor fusion) |
| `Tests/` | numerical tests mirroring the C/Python suite |

## API (mirrors the C `kf_*` API)

| C | C# |
|---|----|
| `kf_kf_init(n, m)` | `new KalmanFilter(n, m)` |
| `kf_kf_init_1d / _constant_velocity` | `KalmanFilter.ConstantSignal / .ConstantVelocity` |
| `kf_kf_set_state / get_state` | `SetState / GetState` |
| `kf_kf_set_covariance* / get_covariance` | `SetCovariance* / GetCovariance` |
| `kf_kf_set_process_noise* / set_measurement_noise*` | `SetProcessNoise* / SetMeasurementNoise*` |
| `kf_kf_set_transition_matrix / set_measurement_matrix` | `SetTransitionMatrix / SetMeasurementMatrix` |
| `kf_kf_predict / update` | `Predict / Update` |
| `kf_kf_set_gate_threshold / nis / update_gated` | `SetGateThreshold / Nis() / UpdateGated` |
| `kf_kf_adapt_r` | `AdaptR` |
| `kf_kf_smooth_step` | `KalmanFilter.RtsSmoothStep(...)` (static) |
| `kf_ekf_*` | `ExtendedKalmanFilter` (delegates f, F_jac, h, H_jac) |
| `kf_ukf_*` | `UnscentedKalmanFilter` (delegates f, h) |

## Differences from the C core (intentional)

- **Errors are .NET exceptions** (`ArgumentException`, `ArgumentOutOfRangeException`,
  `InvalidOperationException`), not status codes; `UpdateGated` returns `bool`
  (true = accepted, false = rejected outlier).
- **A small internal `LinAlg` class** provides the matrix arithmetic (the
  standard .NET library has no dense-matrix type); the algorithms are the same
  Cholesky/solve operations as the C `kf_matrix_*` engine.
