# Quick start

This page gets you filtering in a few minutes. It assumes you have the library
sources on your include/source path (see [Installation](installation.md)).

## 1. Linear Kalman filter (the common case)

The mental model is simple:

1. Initialise with a **one-call constructor** for the model you want (1-D,
   constant velocity, constant acceleration, …).
2. Loop: `kf_kf_predict()` (advance time) then `kf_kf_update()` (feed a
   measurement), and read the result with `kf_kf_get_state()`.

```c
#include "kalman.h"

static kf_kf_t kf;                      /* static => no heap, known size */

void init(void) {
    /* position+velocity, dt=1; builds F, H, Q, R, P for you */
    kf_kf_init_constant_velocity(&kf, /*dt=*/1.0f, /*q_accel=*/0.1f,
                                 /*r=*/1.0f, /*p0_pos=*/10.0f, /*p0_vel=*/10.0f);
}

void step(kf_real_t measurement, kf_real_t dt) {
    kf_kf_predict(&kf, NULL, dt);
    kf_kf_update(&kf, &measurement);

    const kf_real_t *x = kf_kf_get_state(&kf);
    /* x[0] = filtered position, x[1] = filtered velocity */
}
```

Available constructors:

| Constructor | Model |
|-------------|-------|
| `kf_kf_init_1d(q, r)` | 1-D constant signal |
| `kf_kf_init_constant(n, q, r, p0)` | N-D constant signal |
| `kf_kf_init_constant_velocity(dt, q, r, p0_pos, p0_vel)` | position + velocity |
| `kf_kf_init_constant_acceleration(dt, q, r, p0)` | position + velocity + acceleration |

For a custom model, use the general `kf_kf_init()` plus the granular setters
(`kf_kf_set_transition_matrix`, `kf_kf_set_measurement_matrix`,
`kf_kf_set_process_noise`, …).

## 2. Extended Kalman filter

The EKF needs four callbacks: the transition function `f`, its Jacobian `F`,
the measurement function `h`, and its Jacobian `H`.

```c
#include "kalman.h"
#include <math.h>

static kf_ekf_t ekf;

static void f_fn(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                 kf_real_t *xo, void *ctx) { xo[0] = x[0] + dt * x[1]; }

static void F_fn(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                 kf_matrix_t *F, void *ctx) {
    kf_matrix_set(F, 0, 0, 1); kf_matrix_set(F, 0, 1, dt);
    kf_matrix_set(F, 1, 0, 0); kf_matrix_set(F, 1, 1, 1);
}

static void h_fn(const kf_real_t *x, kf_real_t *z, void *ctx) {
    z[0] = x[0] * x[0];                /* e.g. power from an amplitude */
}

static void H_fn(const kf_real_t *x, kf_matrix_t *H, void *ctx) {
    kf_matrix_set(H, 0, 0, 2 * x[0]);  /* dh/dx */
}

void init(void) {
    kf_ekf_init(&ekf, 2, 1);
    kf_ekf_set_models(&ekf, f_fn, F_fn, h_fn, H_fn);
    kf_ekf_set_process_noise_scalar(&ekf, 1e-4f);
    kf_ekf_set_measurement_noise_scalar(&ekf, 1.0f);
    kf_ekf_set_covariance_scalar(&ekf, 1.0f);
}
```

Then `kf_ekf_predict(&ekf, u, dt)` / `kf_ekf_update(&ekf, &z)` in the loop.

## 3. Unscented Kalman filter

The UKF needs only `f` and `h` — no Jacobians:

```c
static kf_ukf_t ukf;

static void f_fn(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                 kf_real_t *xo, void *ctx) { xo[0] = x[0] + dt * x[1]; }

static void h_fn(const kf_real_t *x, kf_real_t *z, void *ctx) { z[0] = x[0]; }

void init(void) {
    kf_ukf_init(&ukf, 2, 1);
    kf_ukf_set_models(&ukf, f_fn, h_fn);
    /* optional: kf_ukf_set_parameters(&ukf, 1.0f, 2.0f, 0.0f); */
    kf_ukf_set_process_noise_scalar(&ukf, 1e-4f);
    kf_ukf_set_measurement_noise_scalar(&ukf, 1.0f);
    kf_ukf_set_covariance_scalar(&ukf, 1.0f);
}
```

## Next steps

- See [Math concepts](math_concepts.md) to understand what the matrices *mean*.
- See [Tuning guide](tuning_guide.md) to choose `Q`, `R`, `P` sensibly.
- See [Examples](examples.md) for ten progressively more advanced examples.
