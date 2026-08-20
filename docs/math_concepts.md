# Mathematical concepts (plain-language)

You do not need deep mathematics to use this library well, but understanding a
few ideas will make tuning far easier.

## What the filter does

A Kalman filter fuses two things that are both imperfect:

1. A **model** of how the system evolves (it "predicts" the next state), and
2. **Measurements** from sensors (which correct that prediction).

It keeps a **Gaussian belief** about the true state: a best estimate `x` (the
mean) and an uncertainty `P` (the covariance). Each step it *predicts* (the
uncertainty grows) and then *updates* (the uncertainty shrinks as a
measurement arrives). The result is the best possible estimate for linear
systems with Gaussian noise — and a very good estimate in practice otherwise.

## The core quantities

| Symbol | Meaning | In this library |
|--------|---------|-----------------|
| `x` | The state vector — everything you want to estimate | `kf_*_get_state()` |
| `P` | Covariance — *how sure* you are about `x` | `kf_*_get_covariance()` |
| `F` | State-transition matrix — the model of how `x` evolves (`x ← F x`) | `kf_kf_set_transition_matrix()` |
| `H` | Measurement matrix — how the state maps to what the sensor reads (`z ≈ H x`) | `kf_kf_set_measurement_matrix()` |
| `Q` | Process-noise covariance — how much you trust the *model* | `kf_*_set_process_noise*()` |
| `R` | Measurement-noise covariance — how much you trust the *sensor* | `kf_*_set_measurement_noise*()` |
| `K` | Kalman gain — how strongly a measurement corrects the state | `kf_*_get_gain()` |
| `y` | Innovation — the difference between the measurement and its prediction (`z − Hx`) | `kf_*_get_innovation()` |

### State `x`

The state is the list of quantities you want to know. Examples:

- `x = [temperature]`
- `x = [position, velocity]`
- `x = [position, velocity, acceleration]`

You choose the state; the filter estimates it.

### Prediction (the `predict` step)

`x ← F x` moves the state forward in time according to your model, and
`P ← F P Fᵀ + Q` grows the uncertainty (because the model is not perfect — `Q`
quantifies how imperfect). For example, with a constant-velocity model,
`F = [[1, dt],[0, 1]]` adds `velocity × dt` to the position each step.

### Measurement (the `update` step)

A sensor reading `z` arrives. `H x` is what the model *predicts* the sensor
should read; `y = z − H x` is the surprise (the innovation). The Kalman gain
`K` decides how much of that surprise to trust, and the filter does
`x ← x + K y`, then shrinks `P`. When the sensor is noisy (large `R`), `K` is
small and the filter trusts its model; when the sensor is precise (small `R`),
`K` is large and the filter trusts the measurement.

### The Kalman gain `K`

`K` is the filter's automatic weighting between model and measurement. You
never set it by hand — it is computed from `P`, `H`, and `R`. If you want the
filter to trust measurements more, lower `R` (or raise `Q`); to trust the
model more, raise `R` (or lower `Q`).

## Why the EKF exists

The linear KF assumes `F` and `H` are linear. Many real models are not — e.g.
measuring `sin(angle)`, or `distance = sqrt(x² + y²)`. The **E**xtended KF
linearises the model around the current estimate using the Jacobian (the
derivative matrix) of the nonlinear functions, then applies the linear KF. It
works well when the nonlinearity is mild over one step; it can diverge for
strongly nonlinear or poorly initialised problems.

## Why the UKF exists

The **U**nscented KF avoids linearisation entirely. Instead of a Jacobian, it
picks a small deterministic set of *sigma points* that represent the
distribution, pushes those points through the nonlinear function, and
reconstructs the mean/covariance from the results. This captures the mean and
covariance accurately to 2nd order (3rd for Gaussians) without derivatives,
and it needs no Jacobians — only the raw `f` and `h`. The cost is evaluating
`f`/`h` at `2n+1` points each step.

### Choosing KF vs EKF vs UKF

| Situation | Recommendation |
|-----------|----------------|
| Linear model | KF (simplest, cheapest, optimal) |
| Mild nonlinearity, derivatives available | EKF |
| Strong nonlinearity, or no easy Jacobians | UKF |

## Further reading

The [Tuning guide](tuning_guide.md) explains how to choose `Q`, `R`, `P` — the
three numbers that most affect behaviour.
