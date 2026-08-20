# Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| Estimate is very noisy / oscillates | `R` too small, or `Q` too large | Raise `R`, lower `Q` (see [Tuning](tuning_guide.md)) |
| Estimate lags real changes | `R` too large, or `Q` too small | Lower `R`, raise `Q` |
| Filter never converges from a bad start | `P₀` too small | Increase the initial covariance |
| Estimate starts fine then diverges / goes NaN | Covariance blow-up, singular `S`, or a bad model | Enable diagnostics; check `Q`,`R` positive-definite; check `R` not zero with zero `P₀` |
| `kf_kf_update` returns `KF_ERROR_NOT_POSITIVE_DEFINITE` | Singular innovation covariance `S = HPHᵀ + R` | Ensure `R` is positive (semi)definite and not zero while `P` is zero; check `H`/dimensions |
| `kf_kf_init` returns `KF_ERROR_INVALID_DIMENSION` | `n`/`m` is 0 or exceeds `KF_MAX_*_DIM` | Fix the call or raise the max in `kalman_config.h` |
| `KF_ERROR_NOT_INITIALIZED` | Used the instance before `init`, or after zeroing it | Call `kf_*_init()` first |
| `KF_ERROR_FEATURE_DISABLED` / symbols missing | The filter was compiled out | Enable `KF_ENABLE_KF/EKF/UKF` in `kalman_config.h` (or via `-D`) |
| `KF_ERROR_NULL_POINTER` | Passed `NULL` (e.g. a `NULL` measurement or control pointer, or a `NULL` instance) | Check the call |
| EKF/UKF return `KF_ERROR_INVALID_PARAMETER` | Model callbacks not registered | Call `kf_ekf_set_models()` / `kf_ukf_set_models()` |
| EKF diverges | Linearisation error (too nonlinear / bad init) | Initialise near the true state, use smaller `dt`, or switch to UKF |
| UKF gives wild results | Sigma-point parameters inappropriate | Use defaults (`alpha=1, beta=2, kappa=0`); very small `alpha` produces a large negative `Wm[0]` |
| Build fails with "C11 required" | Old compiler | Use a C11 compiler (GCC ≥ 5, ARM GCC ≥ 5, Clang) |
| Unexpected "unknown error" from `kf_status_str` | Custom status value | Only values returned by the library are meaningful |

## A debugging checklist for divergence

1. Build with `-DKF_ENABLE_DIAGNOSTICS=1` and call `kf_kf_check(&kf)` (or the
   EKF/UKF equivalent) each iteration; log any non-`KF_OK` result.
2. Verify the model matrices: dimensions (`F` is `n×n`, `H` is `m×n`),
   row-major layout, and that `F`/`H` are the matrices you intend.
3. Verify units are consistent (e.g. `dt` in seconds everywhere).
4. Verify `Q` and `R` are positive semi-definite (diagonal entries ≥ 0).
5. Log the innovation `y` and the gain `K` (advanced API): a monotonically
   growing innovation or a gain that collapses to zero are tell-tale signs.

## Sanity-checking the model

For a linear KF, a good quick check is to run the filter with `Q = 0` and a
large `R` on a *constant* input and confirm the estimate converges to a smooth
constant. Then re-enable `Q` gradually. If the model matrices are wrong, the
filter will visibly misbehave even on clean data — fix the model before
tuning noise.
