# Testing

The test suite is host-runnable and deterministic (a seeded PRNG generates
reproducible noise). Run it with:

```sh
make test
```

A CMake build registers the suite with CTest (`ctest`).

## Structure

| File | Covers |
|------|--------|
| `tests/test_framework.{h,c}` | Tiny assertion harness (no dependencies) |
| `tests/test_matrix.c` | Matrix engine correctness against hand-computed values |
| `tests/test_kf.c` | Linear KF behaviour and error handling |
| `tests/test_ekf.c` | EKF nonlinear models and error handling |
| `tests/test_ukf.c` | UKF nonlinear models, KF equivalence, error handling |
| `tests/test_extensions.c` | Gating/NIS, adaptive R, RTS smoother |
| `tests/test_reference.c` | Exact comparison against naive double-precision reference implementations |
| `tests/test_boundary.c` | Control-input path, maximum dimensions, UKF small-alpha |
| `tests/test_stress.c` | Property-based: P symmetric/positive-definite/finite across random configs |
| `tests/main.c` | Runner |

## What is tested

**Matrix engine**
- Construction, element access, dimension/NULL error paths.
- Identity, add/sub/scale, `add_scaled`.
- Matrix multiply plus the two transpose-helper products, verified against
  hand-computed results; dimension-mismatch errors.
- Transpose (rectangular) and symmetrisation.
- Cholesky factorisation of a known matrix; singular and indefinite matrices
  return `KF_ERROR_NOT_POSITIVE_DEFINITE`.
- Cholesky solve of a known linear system.
- Vector ops and matrix-vector multiply.

**Linear KF**
- 1-D convergence to a constant under noise.
- Position+velocity tracking (velocity inferred from position only).
- Multi-measurement (2×2) update.
- Sequential sensor fusion (two noise levels).
- Zero-noise convergence.
- Variable `dt`.
- Error handling: NULL pointers, out-of-range dimensions, uninitialised use,
  invalid element index.
- NaN/Inf inputs do not crash (deterministic).
- Reset behaviour (state re-zeroed, `P` back to identity).
- Singular innovation covariance returns a clean error.

**EKF**
- Nonlinear measurement `h(x)=x²` converges to the true value.
- Nonlinear transition `f(x)=sin(x)` tracked.
- Missing model callbacks and invalid dimensions are rejected.

**UKF**
- Nonlinear measurement and transition (same models as EKF).
- **Equivalence with the linear KF** on a linear system — since the unscented
  transform is exact for linear systems, the UKF must match the KF; this is a
  strong check of the sigma points, weights, and covariance reconstruction.
- Error handling incl. invalid sigma-point parameters.

**Extensions**
- Gating: NIS thresholding rejects spikes and leaves the state unchanged.
- Adaptive R: converges toward the true (simulated) sensor-noise variance.
- RTS smoother: matches an exact scalar reference to float precision, and
  `P_smooth <= P_filt` holds throughout.

**Reference comparisons**
- The linear KF (3-state / 2-measurement) is run in lock-step against a naive
  double-precision implementation and must match to float precision.
- The UKF's sigma-point machinery is verified to reproduce the covariance
  `A P A^T` exactly for a linear map.

**Boundary / untested paths**
- The control-input `u` term is verified to accumulate exactly (KF) and to be
  forwarded to the EKF transition callback.
- The filters are exercised at `n = KF_MAX_STATE_DIM`,
  `m = KF_MAX_MEASUREMENT_DIM` to flush out any off-by-one in the internal
  scratch/indexing (these run under AddressSanitizer).
- The UKF is run with a small `alpha` (large negative `Wm[0]`).

**Property-based stress**
- Across 40 random KF and 30 random UKF configurations, after every step the
  covariance `P` is checked to be symmetric, positive-definite (Cholesky
  succeeds), and finite.

## Dynamic analysis

The suite is also run under AddressSanitizer + UndefinedBehaviorSanitizer
(including `bounds` and `pointer-overflow`) at `-O1`, which detects buffer
overruns, leaks, and undefined behaviour. See `tools/` and the CI entry points
for how this is invoked.

## Coverage philosophy

The tests exercise **numerical behaviour**, not just compilation: convergence
tolerances, tracking of known trajectories, and equivalence against a trusted
reference (the linear KF). Stochastic tests use a fixed seed so results are
bit-reproducible.

## Extending the suite

Add a `void test_xxx(void)` function, call it from `tests/main.c`, and use the
`CHECK` / `CHECK_NEAR` macros. Keep tests deterministic (seeded PRNG) and avoid
timing-dependent assertions.
