# Q / R / P tuning guide

Behaviour of the filter is dominated by three covariances. Getting them
approximately right matters far more than any other detail.

## The three knobs

| Knob | Question it answers | Too small | Too large |
|------|--------------------|-----------|-----------|
| `R` (measurement noise) | "How much do I trust the sensor?" | Filter overreacts to every noisy sample (jittery estimate) | Filter ignores the sensor and lags badly |
| `Q` (process noise) | "How much do I trust the model?" | Estimate lags real changes (too smooth) | Estimate is noisy, tracks model noise |
| `P₀` (initial covariance) | "How uncertain is my starting guess?" | Slow initial convergence (if the guess is actually bad) | Large initial transients, then fast convergence |

The Kalman gain is the *ratio*: `K ≈ P Hᵀ / (H P Hᵀ + R)`. It is the relative
trust in model vs. measurement that matters, not the absolute scale of `Q` and
`R` individually in most cases.

## Practical recipe

1. **Set `R` from sensor reality.** Measure your sensor's actual variance:
   record N samples of a fixed reference, compute their variance, and set `R`
   to that. A diagonal `R` (independent axes) is almost always fine.
2. **Set `P₀` honestly.** If you genuinely don't know the initial state, use a
   large diagonal (e.g. 10–100× the expected squared range). If you know it,
   use a small value. `P₀` mainly affects startup; it is forgotten quickly.
3. **Tune `Q` last.** Start small and increase it until the filter tracks real
   changes without getting visibly noisy. If you know the process noise
   physically (e.g. the variance of the acceleration driving a
   constant-velocity model), use it directly.

## Model-specific `Q`

The process noise should be placed on the quantities you expect to be driven by
unmodelled disturbances:

- **Constant-signal model** (`F = I`): `Q` is a small scalar; it lets the
  filter follow slow drift.
- **Constant-velocity model** (`x = [p, v]`): put `Q` on the *velocity*, not
  the position — unknown accelerations enter through velocity. With discrete
  time and a white acceleration of variance `σ²`, a common choice is
  `Q = σ² [[dt⁴/4, dt³/2],[dt³/2, dt²]]`.
- **Constant-acceleration model** (`x = [p, v, a]`): put `Q` on the
  *acceleration* (jerk).

## Rules of thumb

- **Oscillating / noisy but fast** estimate → `R` too small or `Q` too large.
- **Lagging / smooth but slow** estimate → `R` too large or `Q` too small.
- **Never converges** from a bad start → `P₀` too small.
- **Explodes or goes NaN** → check that `Q`/`R` are positive (semi)definite and
  that `R` is not zero while the initial covariance is also zero (singular
  innovation); see [Troubleshooting](troubleshooting.md).

## Scalar vs. diagonal vs. full

For convenience the library accepts noise in three forms:

```c
kf_kf_set_process_noise_scalar(&kf, 0.01f);              /* Q = 0.01 * I      */
kf_kf_set_process_noise_diagonal(&kf, diag);             /* independent axes  */
kf_kf_set_process_noise(&kf, Q);                         /* full matrix       */
```

Start with scalar or diagonal; use the full matrix only when you have a good
reason (e.g. cross-coupled process noise). The same applies to `R`.

## A worked example

Filtering a temperature that drifts slowly, sampled by a sensor with a standard
deviation of 0.5 °C:

```c
kf_kf_set_measurement_noise_scalar(&kf, 0.5f * 0.5f);   /* R = variance     */
kf_kf_set_process_noise_scalar(&kf, 1e-3f);             /* Q: allow drift   */
kf_kf_set_covariance_scalar(&kf, 10.0f);                /* P0: unknown start*/
```

If the estimate lags real changes, raise `Q` to `1e-2`. If it looks noisy,
lower `Q` back. This two-knob loop is the whole art of Kalman tuning in most
embedded projects.
