# Examples

Ten host examples progressively demonstrate the library. Each prints results;
build and run them with:

```sh
make run-examples
```

| # | Example | State | Measurement | Demonstrates |
|---|---------|-------|-------------|--------------|
| 01 | `01_basic_1d` | scalar value | scalar | The simplest possible KF |
| 02 | `02_position` | position | position | Random-walk tracking |
| 03 | `03_position_velocity` | [pos, vel] | position | Inferring velocity from position |
| 04 | `04_position_velocity_acceleration` | [pos, vel, acc] | position | Third-order model |
| 05 | `05_multiple_measurements` | [x, y] | [x, y] | Vector measurement update |
| 06 | `06_two_sensors` | scalar | two sensors | Sequential sensor fusion |
| 07 | `07_predict_update_rates` | [pos, vel] | position (slow) | Predict ≠ update rate |
| 08 | `08_variable_dt` | [pos, vel] | position | Irregular sampling |
| 09 | `09_ekf_basic` | angle | sin(angle) | EKF with a nonlinear measurement |
| 10 | `10_ukf_basic` | angle | sin(angle) | Same, but UKF (no Jacobians) |

## 01 — 1D signal filtering
A constant signal observed through noise. `F = [1]`, `H = [1]`, scalar `Q`, `R`,
`P₀`. Shows the absolute minimum API surface.

## 02 — Position estimation
A position that random-walks. Same 1-D setup, but the "truth" moves, so `Q`
must be large enough to follow it — the classic smoothing-vs-lag trade-off.

## 03 — Position + velocity
`x = [p, v]`, `F = [[1, dt],[0, 1]]`, `H = [1, 0]`. Only position is measured,
yet the filter estimates velocity too, and can predict position between
measurements.

## 04 — Position + velocity + acceleration
`x = [p, v, a]`, `F` has `dt²/2` and `dt` terms. Position-only measurements
recover all three quantities. Shows how the constant-acceleration model
generalises example 03.

## 05 — Multiple measurements
A 2-D position updated with a 2-element measurement vector in one call. `R` is
diagonal with *different* noise per axis — demonstrating per-axis trust.

## 06 — Two sensors, different noise
Two sensors measuring the same quantity with different `R` values, applied as
two sequential `update` calls. Mathematically equivalent to a combined update
for uncorrelated sensors, but simpler and rate-flexible.

## 07 — Different predict/update rates
A fast loop predicts 10× per measurement (the IMU/GNSS pattern). Shows that
`predict` and `update` are fully decoupled.

## 08 — Variable dt
`F` is rebuilt with the current `dt` each step. The filter tracks a trajectory
despite irregular sampling.

## 09 — EKF, nonlinear measurement
Estimate an angle from `sin(angle)` measurements. Registers the four EKF
callbacks (`f`, `F`, `h`, `H`) and converges to the true angle.

## 10 — UKF, nonlinear measurement
The same problem solved with the UKF: only `f` and `h` are needed (no
Jacobians). Also shows `kf_ukf_set_parameters`.

The STM32-specific skeletons live in `examples/stm32/` (see
[STM32 integration](stm32_integration.md)).
