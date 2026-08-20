/**
 * @file    test_ukf.c
 * @brief   Unscented Kalman filter unit tests.
 *
 * Covers: nonlinear measurement, nonlinear transition, equivalence with the
 * linear KF on a linear problem (the unscented transform is exact for linear
 * systems, so this is a strong correctness check), and error handling.
 */

#include "test_framework.h"
#include "kalman.h"

#include <stddef.h>
#include <math.h>

#if KF_ENABLE_UKF

static uint32_t rng_state = 0x5A5A5A5Au;

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

static void h_square(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = x[0] * x[0];
}

/* ---- Model B: nonlinear transition f(x) = sin(x), direct measurement ---- */
static void f_sin(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                  kf_real_t *x_out, void *ctx)
{
    (void)u; (void)dt; (void)ctx;
    x_out[0] = (kf_real_t)sin((double)x[0]);
}

static void h_ident(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = x[0];
}

/* ---- Model C: linear constant-velocity model, F=[[1,1],[0,1]], H=[1,0] ---- */
static void f_linear(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                     kf_real_t *x_out, void *ctx)
{
    (void)u; (void)dt; (void)ctx;
    x_out[0] = x[0] + x[1];
    x_out[1] = x[1];
}

static void h_linear(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = x[0];
}

static void test_ukf_nonlinear_measurement(void)
{
    kf_ukf_t ukf;
    kf_real_t x0[1] = {2.0f};
    kf_real_t z;
    int i;

    t_begin("UKF: nonlinear measurement h(x)=x^2");
    CHECK(kf_ukf_init(&ukf, 1, 1) == KF_OK);
    CHECK(kf_ukf_set_models(&ukf, f_ident, h_square) == KF_OK);
    kf_ukf_set_process_noise_scalar(&ukf, 1e-4f);
    kf_ukf_set_measurement_noise_scalar(&ukf, 1.0f);
    kf_ukf_set_covariance_scalar(&ukf, 1.0f);
    kf_ukf_set_state(&ukf, x0);

    for (i = 0; i < 300; i++) {
        CHECK(kf_ukf_predict(&ukf, NULL, 1.0f) == KF_OK);
        z = 9.0f + 2.0f * noise_unit();
        CHECK(kf_ukf_update(&ukf, &z) == KF_OK);
    }
    CHECK_NEAR(kf_ukf_get_state(&ukf)[0], 3.0, 0.15);
    t_end();
}

static void test_ukf_nonlinear_transition(void)
{
    kf_ukf_t ukf;
    kf_real_t x0[1] = {1.0f};
    kf_real_t z, truth = 1.0f;
    int i;

    t_begin("UKF: nonlinear transition f(x)=sin(x)");
    CHECK(kf_ukf_init(&ukf, 1, 1) == KF_OK);
    CHECK(kf_ukf_set_models(&ukf, f_sin, h_ident) == KF_OK);
    kf_ukf_set_process_noise_scalar(&ukf, 1e-6f);
    kf_ukf_set_measurement_noise_scalar(&ukf, 0.01f);
    kf_ukf_set_covariance_scalar(&ukf, 0.1f);
    kf_ukf_set_state(&ukf, x0);

    for (i = 0; i < 100; i++) {
        kf_ukf_predict(&ukf, NULL, 1.0f);
        truth = (kf_real_t)sin((double)truth);
        z = truth + 0.1f * noise_unit();
        kf_ukf_update(&ukf, &z);
    }
    CHECK_NEAR(kf_ukf_get_state(&ukf)[0], truth, 0.1);
    t_end();
}

static void test_ukf_matches_kf_linear(void)
{
    /* The unscented transform reproduces the linear KF exactly, so a UKF with
       linear models must agree with the linear KF to high accuracy. */
    kf_ukf_t ukf;
    kf_kf_t kf;
    kf_real_t F[4] = {1, 1, 0, 1};   /* constant velocity, dt = 1 */
    kf_real_t H[2] = {1, 0};
    kf_real_t Q[4] = {0, 0, 0, 0.05f};
    kf_real_t P0[4] = {10, 0, 0, 10};
    kf_real_t x0[2] = {0, 0};
    int i;

    t_begin("UKF: equivalence with linear KF on a linear system");

    CHECK(kf_kf_init(&kf, 2, 1) == KF_OK);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise(&kf, Q);
    kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
    kf_kf_set_covariance(&kf, P0);
    kf_kf_set_state(&kf, x0);

    CHECK(kf_ukf_init(&ukf, 2, 1) == KF_OK);
    CHECK(kf_ukf_set_models(&ukf, f_linear, h_linear) == KF_OK);
    kf_ukf_set_process_noise(&ukf, Q);
    kf_ukf_set_measurement_noise_scalar(&ukf, 1.0f);
    kf_ukf_set_covariance(&ukf, P0);
    kf_ukf_set_state(&ukf, x0);

    for (i = 0; i < 200; i++) {
        kf_real_t z = 2.0f * (kf_real_t)(i + 1) + 0.5f * noise_unit();

        kf_kf_predict(&kf, NULL, 1.0f);
        kf_kf_update(&kf, &z);

        kf_ukf_predict(&ukf, NULL, 1.0f);
        kf_ukf_update(&ukf, &z);
    }

    CHECK_NEAR(kf_kf_get_state(&kf)[0], kf_ukf_get_state(&ukf)[0], 0.05);
    CHECK_NEAR(kf_kf_get_state(&kf)[1], kf_ukf_get_state(&ukf)[1], 0.02);
    t_end();
}

static void test_ukf_error_handling(void)
{
    kf_ukf_t ukf;
    kf_real_t z = 1.0f;

    t_begin("UKF: error handling");
    CHECK(kf_ukf_init(&ukf, 1, 1) == KF_OK);

    /* No models registered -> predict/update fail. */
    CHECK(kf_ukf_predict(&ukf, NULL, 1.0f) == KF_ERROR_INVALID_PARAMETER);
    CHECK(kf_ukf_update(&ukf, &z) == KF_ERROR_INVALID_PARAMETER);

#if KF_ENABLE_RUNTIME_CHECKS
    CHECK(kf_ukf_init(NULL, 1, 1) == KF_ERROR_NULL_POINTER);
#endif
    CHECK(kf_ukf_init(&ukf, KF_MAX_STATE_DIM + 1, 1) == KF_ERROR_INVALID_DIMENSION);

    /* Invalid sigma-point parameters. */
    CHECK(kf_ukf_init(&ukf, 1, 1) == KF_OK);
    CHECK(kf_ukf_set_parameters(&ukf, 0.0f, 2.0f, 0.0f) == KF_ERROR_INVALID_PARAMETER);
    CHECK(kf_ukf_set_parameters(&ukf, -1.0f, 2.0f, 0.0f) == KF_ERROR_INVALID_PARAMETER);
    CHECK(kf_ukf_set_parameters(&ukf, 1.0f, -1.0f, 0.0f) == KF_ERROR_INVALID_PARAMETER);
    t_end();
}

void test_ukf(void)
{
    test_ukf_nonlinear_measurement();
    test_ukf_nonlinear_transition();
    test_ukf_matches_kf_linear();
    test_ukf_error_handling();
}

#else /* KF_ENABLE_UKF */

void test_ukf(void)
{
}

#endif /* KF_ENABLE_UKF */
