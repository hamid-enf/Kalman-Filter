# API reference

Naming is consistent across the three filters: the prefix identifies the
filter (`kf_kf_`, `kf_ekf_`, `kf_ukf_`) and the suffix is the operation. Only
the linear KF is listed fully here; the EKF and UKF expose the same
state/covariance/noise setters with their own prefix, plus the model-specific
calls shown at the end.

All functions that can fail return a `kf_status_t` (`KF_OK` = 0, negative =
error). Check it — the library does not silently ignore errors.

## Common types

| Type | Description |
|------|-------------|
| `kf_real_t` | The scalar type (`float` by default). |
| `kf_status_t` | Status code enum (see below). |
| `kf_matrix_t` | A matrix *view*: `{ rows, cols, data }` (row-major, caller-owned storage). |

### Status codes

`KF_OK`, `KF_ERROR_NULL_POINTER`, `KF_ERROR_INVALID_DIMENSION`,
`KF_ERROR_INVALID_PARAMETER`, `KF_ERROR_NOT_INITIALIZED`,
`KF_ERROR_FEATURE_DISABLED`, `KF_ERROR_SINGULAR_MATRIX`,
`KF_ERROR_NOT_POSITIVE_DEFINITE`, `KF_ERROR_NUMERICAL`,
`KF_ERROR_NON_FINITE`, `KF_ERROR_DIVERGED`.

`const char *kf_status_str(kf_status_t s)` returns a static, human-readable
description (safe to use for logging; do not free).

## Linear Kalman filter (`kf_kf_*`)

### Lifecycle

```c
kf_status_t kf_kf_init(kf_kf_t *kf, uint16_t n, uint16_t m);
kf_status_t kf_kf_reset(kf_kf_t *kf);
```

`kf_kf_init` sets `x = 0`, `F = I`, `H = 0`, `P = I`, `Q = 0`, `R = 0` and
marks the instance initialised. `kf_kf_reset` restores those defaults while
keeping the dimensions.

### State access

```c
kf_status_t kf_kf_set_state(kf_kf_t *kf, const kf_real_t *x);   /* n values  */
const kf_real_t *kf_kf_get_state(const kf_kf_t *kf);            /* returns x */
kf_status_t kf_kf_set_state_element(kf_kf_t *kf, uint16_t i, kf_real_t v);
```

### Covariance P

```c
kf_status_t kf_kf_set_covariance(kf_kf_t *kf, const kf_real_t *P);           /* n*n */
kf_status_t kf_kf_set_covariance_diagonal(kf_kf_t *kf, const kf_real_t *d);  /* n    */
kf_status_t kf_kf_set_covariance_scalar(kf_kf_t *kf, kf_real_t p);           /* p*I  */
const kf_real_t *kf_kf_get_covariance(const kf_kf_t *kf);                    /* n*n  */
```

### Process noise Q / measurement noise R

```c
kf_status_t kf_kf_set_process_noise(kf_kf_t *kf, const kf_real_t *Q);           /* n*n */
kf_status_t kf_kf_set_process_noise_diagonal(kf_kf_t *kf, const kf_real_t *d);  /* n   */
kf_status_t kf_kf_set_process_noise_scalar(kf_kf_t *kf, kf_real_t q);           /* q*I */

kf_status_t kf_kf_set_measurement_noise(kf_kf_t *kf, const kf_real_t *R);        /* m*m */
kf_status_t kf_kf_set_measurement_noise_diagonal(kf_kf_t *kf, const kf_real_t *d);
kf_status_t kf_kf_set_measurement_noise_scalar(kf_kf_t *kf, kf_real_t r);

const kf_real_t *kf_kf_get_process_noise(const kf_kf_t *kf);     /* n*n, read Q */
const kf_real_t *kf_kf_get_measurement_noise(const kf_kf_t *kf); /* m*m, read R */
```

### Model matrices

```c
kf_status_t kf_kf_set_transition_matrix(kf_kf_t *kf, const kf_real_t *F);  /* n*n */
kf_status_t kf_kf_set_measurement_matrix(kf_kf_t *kf, const kf_real_t *H); /* m*n */
```

### Real-time path

```c
kf_status_t kf_kf_predict(kf_kf_t *kf, const kf_real_t *u, kf_real_t dt);
kf_status_t kf_kf_update(kf_kf_t *kf, const kf_real_t *z);              /* m values */
```

- `u` is an optional n-vector control input, already expressed in state space
  (i.e. pre-multiplied by any `B` and by `dt`). Pass `NULL` for none.
- `dt` is accepted for API symmetry with the EKF/UKF. For a linear KF the step
  size is encoded in `F` (and `Q`); for variable `dt`, re-set `F`/`Q` before
  `predict` (see example 08).
- `z` is the measurement vector; it may be `NULL` only if you want prediction
  without correction (though calling only `predict` is the normal way to do
  that).

### Convenience constructors

```c
kf_status_t kf_kf_init_1d(kf_kf_t *kf, kf_real_t q, kf_real_t r);
kf_status_t kf_kf_init_constant(kf_kf_t *kf, uint16_t n,
                                kf_real_t q, kf_real_t r, kf_real_t p0);
kf_status_t kf_kf_init_constant_velocity(kf_kf_t *kf, kf_real_t dt,
                                         kf_real_t q_accel, kf_real_t r,
                                         kf_real_t p0_pos, kf_real_t p0_vel);
kf_status_t kf_kf_init_constant_acceleration(kf_kf_t *kf, kf_real_t dt,
                                             kf_real_t q_jerk, kf_real_t r,
                                             kf_real_t p0);
```

These configure the whole filter (F, H, Q, R, P) in one call for the most
common models: 1-D / N-D constant signal, constant velocity, and constant
acceleration. They build the model matrices for you, so beginners never have to
write matrix algebra.

### Extensions (optional feature sets)

```c
/* Innovation gating / NIS (KF_ENABLE_GATING) */
kf_status_t kf_kf_set_gate_threshold(kf_kf_t *kf, kf_real_t chi2);
kf_real_t   kf_kf_nis(const kf_kf_t *kf);
kf_status_t kf_kf_update_gated(kf_kf_t *kf, const kf_real_t *z);

/* Adaptive measurement noise (KF_ENABLE_ADAPTIVE_R) */
kf_status_t kf_kf_adapt_r(kf_kf_t *kf, kf_real_t gamma, kf_real_t r_min);

/* RTS fixed-interval smoother (KF_ENABLE_SMOOTHER) */
kf_status_t kf_kf_smooth_step(kf_kf_t *kf,
                              const kf_real_t *x_filt, const kf_real_t *P_filt,
                              const kf_real_t *F,
                              const kf_real_t *x_pred, const kf_real_t *P_pred,
                              const kf_real_t *x_smooth_next,
                              const kf_real_t *P_smooth_next,
                              kf_real_t *x_smooth, kf_real_t *P_smooth);
```

- `kf_kf_update_gated` returns `KF_WARN_GATED` (positive) when a measurement is
  rejected as an outlier; the state is left unchanged in that case.
- `kf_kf_nis` returns the normalized innovation squared of the last update.
- `kf_kf_adapt_r` updates `R` in place from the last innovation;
  `kf_kf_get_measurement_noise()` returns the adapted `R`.
- `kf_kf_smooth_step` performs one backward smoothing step; sweep it from the
  end of a stored forward trajectory toward the start (see example 13).

### Advanced access (`KF_ENABLE_ADVANCED_API`)

```c
kf_status_t kf_kf_get_innovation(const kf_kf_t *kf, kf_real_t *y);  /* m, y = z - Hx */
kf_status_t kf_kf_get_gain(const kf_kf_t *kf, kf_real_t *K);        /* n*m, Kalman gain */
```

These return the innovation and gain from the *most recent* update.

### Diagnostics (`KF_ENABLE_DIAGNOSTICS`)

```c
kf_status_t kf_kf_check(const kf_kf_t *kf);  /* finiteness + P symmetry */
```

## Extended Kalman filter (`kf_ekf_*`)

Same lifecycle/state/covariance/noise setters as the KF (`kf_ekf_set_state`,
`kf_ekf_set_covariance*`, `kf_ekf_set_process_noise*`,
`kf_ekf_set_measurement_noise*`), plus:

```c
typedef void (*kf_ekf_f_fn)(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                            kf_real_t *x_out, void *ctx);
typedef void (*kf_ekf_h_fn)(const kf_real_t *x, kf_real_t *z_out, void *ctx);
typedef void (*kf_ekf_jacobian_f_fn)(const kf_real_t *x, const kf_real_t *u,
                                     kf_real_t dt, kf_matrix_t *F, void *ctx);
typedef void (*kf_ekf_jacobian_h_fn)(const kf_real_t *x, kf_matrix_t *H, void *ctx);

kf_status_t kf_ekf_set_models(kf_ekf_t *ekf,
                              kf_ekf_f_fn f, kf_ekf_jacobian_f_fn F,
                              kf_ekf_h_fn h, kf_ekf_jacobian_h_fn H);
void        kf_ekf_set_context_f(kf_ekf_t *ekf, void *ctx);
void        kf_ekf_set_context_h(kf_ekf_t *ekf, void *ctx);

kf_status_t kf_ekf_predict(kf_ekf_t *ekf, const kf_real_t *u, kf_real_t dt);
kf_status_t kf_ekf_update(kf_ekf_t *ekf, const kf_real_t *z);
```

The Jacobian callbacks receive a pre-sized `kf_matrix_t` bound to internal
storage; fill it with `kf_matrix_set(&F, i, j, value)`. Do not resize it or
change its `data` pointer.

## Unscented Kalman filter (`kf_ukf_*`)

Same lifecycle/state/covariance/noise setters as the KF, plus:

```c
typedef void (*kf_ukf_f_fn)(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                            kf_real_t *x_out, void *ctx);
typedef void (*kf_ukf_h_fn)(const kf_real_t *x, kf_real_t *z_out, void *ctx);

kf_status_t kf_ukf_set_models(kf_ukf_t *ukf, kf_ukf_f_fn f, kf_ukf_h_fn h);
void        kf_ukf_set_context_f(kf_ukf_t *ukf, void *ctx);
void        kf_ukf_set_context_h(kf_ukf_t *ukf, void *ctx);

kf_status_t kf_ukf_set_parameters(kf_ukf_t *ukf, kf_real_t alpha,
                                  kf_real_t beta, kf_real_t kappa);

kf_status_t kf_ukf_predict(kf_ukf_t *ukf, const kf_real_t *u, kf_real_t dt);
kf_status_t kf_ukf_update(kf_ukf_t *ukf, const kf_real_t *z);
```

`alpha` (spread), `beta` (distribution prior), `kappa` (secondary scaling)
control the sigma points; `kf_ukf_set_parameters` recomputes the weights.
Defaults: `alpha = 1`, `beta = 2`, `kappa = 0`.

## Matrix engine (`kf_matrix_*`, advanced)

A minimal, dependency-free matrix toolkit used internally. Useful when you want
to build custom models or inspect internals directly.

```c
kf_status_t kf_matrix_init(kf_matrix_t *m, uint16_t rows, uint16_t cols, kf_real_t *data);
kf_real_t   kf_matrix_get(const kf_matrix_t *m, uint16_t r, uint16_t c);
void        kf_matrix_set(kf_matrix_t *m, uint16_t r, uint16_t c, kf_real_t v);

kf_status_t kf_matrix_zero(kf_matrix_t *m);
kf_status_t kf_matrix_identity(kf_matrix_t *m);
kf_status_t kf_matrix_copy(kf_matrix_t *dst, const kf_matrix_t *src);
kf_status_t kf_matrix_scale(kf_matrix_t *m, kf_real_t s);
kf_status_t kf_matrix_add(kf_matrix_t *out, const kf_matrix_t *a, const kf_matrix_t *b);
kf_status_t kf_matrix_add_scaled(kf_matrix_t *out, const kf_matrix_t *a, const kf_matrix_t *b, kf_real_t s);
kf_status_t kf_matrix_sub(kf_matrix_t *out, const kf_matrix_t *a, const kf_matrix_t *b);
kf_status_t kf_matrix_symmetrize(kf_matrix_t *m);          /* m = (m + m^T)/2 */

kf_status_t kf_matrix_mul(kf_matrix_t *out, const kf_matrix_t *a, const kf_matrix_t *b);
kf_status_t kf_matrix_mul_transpose_b(kf_matrix_t *out, const kf_matrix_t *a, const kf_matrix_t *b); /* a * b^T */
kf_status_t kf_matrix_mul_transpose_a(kf_matrix_t *out, const kf_matrix_t *a, const kf_matrix_t *b); /* a^T * b */
kf_status_t kf_matrix_transpose(kf_matrix_t *out, const kf_matrix_t *a);

kf_status_t kf_matrix_cholesky(kf_matrix_t *m);            /* in-place, A = L L^T */
kf_status_t kf_matrix_cholesky_solve(const kf_matrix_t *L, kf_matrix_t *B); /* L L^T X = B */

/* vectors (plain kf_real_t arrays) */
void kf_vec_copy/kf_vec_zero/kf_vec_scale/kf_vec_add/kf_vec_sub/kf_vec_axpy(...);
kf_real_t kf_vec_dot(const kf_real_t *a, const kf_real_t *b, uint16_t n);
kf_status_t kf_mat_vec_mul(kf_real_t *out, const kf_matrix_t *A, const kf_real_t *x);
```

**Aliasing rules**: element-wise operations may alias the output with an input;
the `*_mul`/`*_transpose_*` functions must not alias the output with either
input (use a scratch buffer).

**Note on inversion**: the engine deliberately provides `cholesky`/`solve`
rather than a general `invert`, because the filters never need an explicit
inverse — they solve `S X = B`. If you need an inverse yourself, call
`kf_matrix_cholesky` on a symmetric positive-definite matrix and then
`kf_matrix_cholesky_solve` with an identity right-hand side.
