# Architecture

## Layers

The library is layered so that the algorithm, the linear algebra, the
configuration, and the user API stay independent:

```
┌────────────────────────────────────────────────────────────┐
│  User application                                          │
│  (high-level: kf_kf_* / kf_ekf_* / kf_ukf_*  … or advanced)│
├────────────────────────────────────────────────────────────┤
│  Filter layers          kalman_kf / kalman_ekf / kalman_ukf │
├────────────────────────────────────────────────────────────┤
│  Matrix engine          kalman_matrix (views, no ownership) │
├────────────────────────────────────────────────────────────┤
│  Types & config         kalman_types / kalman_config        │
└────────────────────────────────────────────────────────────┘
          (no dependency on any HAL / RTOS / CMSIS / libc beyond sqrt)
```

- **Configuration** (`kalman_config.h`) is processed first and defines the
  scalar type, the maximum dimensions, and every feature switch.
- **Types** (`kalman_types.h`) defines `kf_status_t` and the scalar type.
- **Matrix engine** (`kalman_matrix.{h,c}`) is the only place that touches
  arithmetic on matrices. It is a *view-based* design: a `kf_matrix_t` holds
  only `{rows, cols, data*}` and never owns memory.
- **Filters** (`kalman_kf/ekf/ukf.{h,c}`) contain the estimation algorithms and
  own their state. They call the matrix engine for linear algebra.
- **Examples / STM32 / tests / benchmarks** are entirely outside the core.

## Memory model

- Every filter instance is a plain `struct` whose size is fixed at compile time
  (`sizeof(kf_kf_t)`, `sizeof(kf_ekf_t)`, `sizeof(kf_ukf_t)`). There is **no
  dynamic allocation** anywhere.
- Instances are provided by the caller: as a `static`/global (zero-initialised,
  then `init`ed), on the stack, or in a user-managed pool.
- Within an instance, the state `x`, covariance `P`, noise `Q`/`R`, model
  `F`/`H`, and a private `scratch` buffer are all embedded. The scratch buffer
  is reused by `predict` and `update` and is sized by a compile-time macro that
  is verified against the layout with a `#error` guard.

The RAM cost is therefore fully predictable from `n`, `m`, the enabled filters,
and the scalar type — see [Performance](performance.md).

## Ownership & reentrancy

- No hidden global mutable state. All functions are reentrant and touch only
  the instance passed in.
- A single instance is **not** safe to use concurrently from multiple contexts
  (ISR + task + main loop) without external synchronisation; each call is a
  read-modify-write of the whole instance. Distinct instances are independent
  and may be used from distinct contexts freely.
- The library never creates threads, locks, or interrupts.

## Numerical strategy

1. **No explicit matrix inversions.** The Kalman gain `K = (P Hᵀ) S⁻¹` is
   computed by solving `S Kᵀ = (P Hᵀ)ᵀ` via a Cholesky factorisation
   (`kf_matrix_cholesky` + `kf_matrix_cholesky_solve`). Cholesky is the most
   stable practical method for symmetric positive-definite systems.
2. **Joseph-form covariance update** (default): `P = (I−KH)P(I−KH)ᵀ + KRKᵀ`.
   Unlike the simplified `P = (I−KH)P`, this form is symmetric by construction
   and keeps `P` positive semi-definite even with rounding.
3. **Symmetrisation.** Covariances are symmetrised (`(P+Pᵀ)/2`) after each
   update to suppress the tiny asymmetry that floating-point rounding
   introduces.
4. **Sigma points from Cholesky.** The UKF derives its sigma points from
   `sqrt(P)` via `P = L Lᵀ`, avoiding an eigendecomposition.
5. **Explicit singular/indefinite detection.** A non-positive Cholesky pivot
   returns `KF_ERROR_NOT_POSITIVE_DEFINITE` instead of producing NaNs.

## Error handling

Every fallible function returns `kf_status_t`. Errors are deterministic and
specific (`NULL` pointer, bad dimension, singular matrix, non-PD matrix,
non-finite, …). See the [API reference](api_reference.md) for the full enum.
