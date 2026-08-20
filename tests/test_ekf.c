/**
 * @file    test_ekf.c
 * @brief   Extended Kalman filter unit tests.
 *
 * Covers: nonlinear measurement (h = x^2), nonlinear state transition
 * (f = sin), Jacobian correctness, and error handling.
 */

#include "test_framework.h"
#include "kalman.h"

#include <stddef.h>
#include <math.h>

#if KF_ENABLE_EKF

static uint32_t rng_state = 0xA5A5A5A5u;

static uint32_t rng_next(void)
{
    uint32_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_state = x;
    return x;
}

static kf_real_t noise_unit(void)
{
    return (kf_real_t)((int32_t)(rng_next() % 2000001u) - 1000000) / 1000000.0f;
}

/* ---- Model A: identity transition, quadratic measurement h(x) = x^2 ---- */
static void f_ident(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                    kf_real_t *x_out, void *ctx)
{
    (void)u; (void)dt; (void)ctx;
    x_out[0] = x[0];
}

static void F_ident(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                    kf_matrix_t *F, void *ctx)
{
    (void)x; (void)u; (void)dt; (void)ctx;
    kf_matrix_set(F, 0, 0, 1.0f);
}

static void h_square(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = x[0] * x[0];
}

static void H_square(const kf_real_t *x, kf_matrix_t *H, void *ctx)
{
    (void)ctx;
    kf_matrix_set(H, 0, 0, 2.0f * x[0]);
}

/* ---- Model B: nonlinear transition f(x) = sin(x), direct measurement ---- */
static void f_sin(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                  kf_real_t *x_out, void *ctx)
{
    (void)u; (void)dt; (void)ctx;
    x_out[0] = (kf_real_t)sin((double)x[0]);
}

static void F_sin(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                  kf_matrix_t *F, void *ctx)
{
    (void)u; (void)dt; (void)ctx;
    kf_matrix_set(F, 0, 0, (kf_real_t)cos((double)x[0]));
}

static void h_ident(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = x[0];
}

static void H_ident(const kf_real_t *x, kf_matrix_t *H, void *ctx)
{
    (void)x; (void)ctx;
    kf_matrix_set(H, 0, 0, 1.0f);
}

static void test_ekf_nonlinear_measurement(void)
{
    kf_ekf_t ekf;
    kf_real_t x0[1] = {2.0f};   /* start below the true value 3 */
    kf_real_t z;
    int i;

    t_begin("EKF: nonlinear measurement h(x)=x^2");
    CHECK(kf_ekf_init(&ekf, 1, 1) == KF_OK);
    CHECK(kf_ekf_set_models(&ekf, f_ident, F_ident, h_square, H_square) == KF_OK);
    kf_ekf_set_process_noise_scalar(&ekf, 1e-4f);
    kf_ekf_set_measurement_noise_scalar(&ekf, 1.0f);
    kf_ekf_set_covariance_scalar(&ekf, 1.0f);
    kf_ekf_set_state(&ekf, x0);

    for (i = 0; i < 300; i++) {
        CHECK(kf_ekf_predict(&ekf, NULL, 1.0f) == KF_OK);
        z = 9.0f + 2.0f * noise_unit();   /* true x = 3, z = 9 + noise */
        CHECK(kf_ekf_update(&ekf, &z) == KF_OK);
    }
    CHECK_NEAR(kf_ekf_get_state(&ekf)[0], 3.0, 0.15);
    t_end();
}

static void test_ekf_nonlinear_transition(void)
{
    kf_ekf_t ekf;
    kf_real_t x0[1] = {1.0f};
    kf_real_t z, truth = 1.0f;
    int i;

    t_begin("EKF: nonlinear transition f(x)=sin(x)");
    CHECK(kf_ekf_init(&ekf, 1, 1) == KF_OK);
    CHECK(kf_ekf_set_models(&ekf, f_sin, F_sin, h_ident, H_ident) == KF_OK);
    kf_ekf_set_process_noise_scalar(&ekf, 1e-6f);
    kf_ekf_set_measurement_noise_scalar(&ekf, 0.01f);
    kf_ekf_set_covariance_scalar(&ekf, 0.1f);
    kf_ekf_set_state(&ekf, x0);

    for (i = 0; i < 100; i++) {
        kf_ekf_predict(&ekf, NULL, 1.0f);
        truth = (kf_real_t)sin((double)truth);
        z = truth + 0.1f * noise_unit();
        kf_ekf_update(&ekf, &z);
    }
    CHECK_NEAR(kf_ekf_get_state(&ekf)[0], truth, 0.1);
    t_end();
}

static void test_ekf_error_handling(void)
{
    kf_ekf_t ekf;
    kf_real_t z = 1.0f;

    t_begin("EKF: error handling");
    CHECK(kf_ekf_init(&ekf, 1, 1) == KF_OK);

    /* No models registered -> predict/update must fail. */
    CHECK(kf_ekf_predict(&ekf, NULL, 1.0f) == KF_ERROR_INVALID_PARAMETER);
    CHECK(kf_ekf_update(&ekf, &z) == KF_ERROR_INVALID_PARAMETER);

#if KF_ENABLE_RUNTIME_CHECKS
    CHECK(kf_ekf_init(NULL, 1, 1) == KF_ERROR_NULL_POINTER);
#endif
    CHECK(kf_ekf_init(&ekf, KF_MAX_STATE_DIM + 1, 1) == KF_ERROR_INVALID_DIMENSION);
    t_end();
}

void test_ekf(void)
{
    test_ekf_nonlinear_measurement();
    test_ekf_nonlinear_transition();
    test_ekf_error_handling();
}

#else /* KF_ENABLE_EKF */

void test_ekf(void)
{
}

#endif /* KF_ENABLE_EKF */
