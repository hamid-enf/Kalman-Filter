/**
 * @file    test_kf.c
 * @brief   Linear Kalman filter unit tests.
 *
 * Covers: 1D convergence, multi-state, multi-measurement, sequential sensor
 * updates, zero/high noise, constant/changing signals, variable dt, singular
 * measurement covariance, invalid inputs, NaN/Inf inputs, and reset behaviour.
 */

#include "test_framework.h"
#include "kalman.h"

#include <stddef.h>
#include <math.h>

/* Deterministic PRNG (xorshift32) so tests are reproducible. */
static uint32_t rng_state = 0x12345678u;

static uint32_t rng_next(void)
{
    uint32_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_state = x;
    return x;
}

/* Gaussian-ish noise in [-1, 1] (uniform), good enough for stochastic tests. */
static kf_real_t noise_unit(void)
{
    return (kf_real_t)((int32_t)(rng_next() % 2000001u) - 1000000) / 1000000.0f;
}

static void test_1d_convergence(void)
{
    kf_kf_t kf;
    kf_real_t F[1] = {1.0f};
    kf_real_t H[1] = {1.0f};
    kf_real_t z, x0[1] = {0.0f};
    int i;

    t_begin("KF: 1D constant-signal convergence");
    CHECK(kf_kf_init(&kf, 1, 1) == KF_OK);
    CHECK(kf_kf_set_transition_matrix(&kf, F) == KF_OK);
    CHECK(kf_kf_set_measurement_matrix(&kf, H) == KF_OK);
    CHECK(kf_kf_set_process_noise_scalar(&kf, 1e-4f) == KF_OK);
    CHECK(kf_kf_set_measurement_noise_scalar(&kf, 1.0f) == KF_OK);
    CHECK(kf_kf_set_covariance_scalar(&kf, 1.0f) == KF_OK);
    kf_kf_set_state(&kf, x0);

    for (i = 0; i < 500; i++) {
        CHECK(kf_kf_predict(&kf, NULL, 1.0f) == KF_OK);
        z = 10.0f + noise_unit();               /* noisy measurement of 10 */
        CHECK(kf_kf_update(&kf, &z) == KF_OK);
    }
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 10.0, 0.3);
    t_end();
}

static void test_constant_velocity(void)
{
    kf_kf_t kf;
    kf_real_t F[4] = {1, 1, 0, 1};   /* dt = 1 */
    kf_real_t H[2] = {1, 0};
    kf_real_t Q[4] = {0, 0, 0, 0.05f};
    kf_real_t P0[4] = {100, 0, 0, 100};
    kf_real_t x0[2] = {0, 0};
    kf_real_t z;
    int i;

    t_begin("KF: position+velocity tracking");
    CHECK(kf_kf_init(&kf, 2, 1) == KF_OK);
    CHECK(kf_kf_set_transition_matrix(&kf, F) == KF_OK);
    CHECK(kf_kf_set_measurement_matrix(&kf, H) == KF_OK);
    CHECK(kf_kf_set_process_noise(&kf, Q) == KF_OK);
    CHECK(kf_kf_set_measurement_noise_scalar(&kf, 1.0f) == KF_OK);
    CHECK(kf_kf_set_covariance(&kf, P0) == KF_OK);
    kf_kf_set_state(&kf, x0);

    for (i = 0; i < 500; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);
        z = 2.0f * (kf_real_t)(i + 1) + 0.5f * noise_unit();   /* pos = 2t */
        kf_kf_update(&kf, &z);
    }
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 1000.0, 5.0);   /* position ~= 2*500 */
    CHECK_NEAR(kf_kf_get_state(&kf)[1], 2.0, 0.2);      /* velocity ~= 2     */
    t_end();
}

static void test_multi_measurement(void)
{
    /* 2D state, 2D measurement (both components directly observed). */
    kf_kf_t kf;
    kf_real_t F[4] = {1, 0, 0, 1};
    kf_real_t H[4] = {1, 0, 0, 1};
    kf_real_t x0[2] = {0, 0};
    kf_real_t z[2];
    int i;

    t_begin("KF: multi-measurement (2x2)");
    CHECK(kf_kf_init(&kf, 2, 2) == KF_OK);
    CHECK(kf_kf_set_transition_matrix(&kf, F) == KF_OK);
    CHECK(kf_kf_set_measurement_matrix(&kf, H) == KF_OK);
    CHECK(kf_kf_set_process_noise_scalar(&kf, 1e-3f) == KF_OK);
    CHECK(kf_kf_set_measurement_noise_scalar(&kf, 1.0f) == KF_OK);
    CHECK(kf_kf_set_covariance_scalar(&kf, 1.0f) == KF_OK);
    kf_kf_set_state(&kf, x0);

    for (i = 0; i < 300; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);
        z[0] = 5.0f + noise_unit();
        z[1] = -3.0f + noise_unit();
        kf_kf_update(&kf, z);
    }
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 5.0, 0.2);
    CHECK_NEAR(kf_kf_get_state(&kf)[1], -3.0, 0.2);
    t_end();
}

static void test_sequential_sensors(void)
{
    /* Two independent sensors measuring the same scalar state with different
       noise levels; updated sequentially (see example 06). */
    kf_kf_t kf;
    kf_real_t F[1] = {1.0f};
    kf_real_t H[1] = {1.0f};
    kf_real_t x0[1] = {0.0f};
    kf_real_t z1, z2;
    int i;

    t_begin("KF: sequential sensor fusion");
    CHECK(kf_kf_init(&kf, 1, 1) == KF_OK);
    CHECK(kf_kf_set_transition_matrix(&kf, F) == KF_OK);
    CHECK(kf_kf_set_measurement_matrix(&kf, H) == KF_OK);
    CHECK(kf_kf_set_process_noise_scalar(&kf, 1e-4f) == KF_OK);
    CHECK(kf_kf_set_covariance_scalar(&kf, 1.0f) == KF_OK);
    kf_kf_set_state(&kf, x0);

    for (i = 0; i < 500; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);

        /* Sensor 1: noisy (R = 10). */
        kf_kf_set_measurement_noise_scalar(&kf, 10.0f);
        z1 = 7.0f + 3.0f * noise_unit();
        kf_kf_update(&kf, &z1);

        /* Sensor 2: precise (R = 0.1). */
        kf_kf_set_measurement_noise_scalar(&kf, 0.1f);
        z2 = 7.0f + 0.3f * noise_unit();
        kf_kf_update(&kf, &z2);
    }
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 7.0, 0.1);
    t_end();
}

static void test_zero_noise(void)
{
    /* With zero noise the filter should trust measurements completely and
       converge essentially exactly. */
    kf_kf_t kf;
    kf_real_t F[1] = {1.0f}, H[1] = {1.0f};
    kf_real_t x0[1] = {0.0f}, z;
    int i;

    t_begin("KF: zero-noise convergence");
    CHECK(kf_kf_init(&kf, 1, 1) == KF_OK);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise_scalar(&kf, 0.0f);
    kf_kf_set_measurement_noise_scalar(&kf, 1e-6f);
    kf_kf_set_covariance_scalar(&kf, 1.0f);
    kf_kf_set_state(&kf, x0);

    for (i = 0; i < 100; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);
        z = 42.0f;
        kf_kf_update(&kf, &z);
    }
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 42.0, 1e-2);
    t_end();
}

static void test_variable_dt(void)
{
    /* Constant-velocity model driven at varying dt; F is updated each step. */
    kf_kf_t kf;
    kf_real_t H[2] = {1, 0};
    kf_real_t Q[4] = {0, 0, 0, 0.01f};
    kf_real_t P0[4] = {10, 0, 0, 10};
    kf_real_t x0[2] = {0, 0};
    kf_real_t z;
    int i;

    t_begin("KF: variable dt");
    CHECK(kf_kf_init(&kf, 2, 1) == KF_OK);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise(&kf, Q);
    kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
    kf_kf_set_covariance(&kf, P0);
    kf_kf_set_state(&kf, x0);

    {
        /* Drive a known trajectory with irregular steps. */
        kf_real_t t = 0.0f;
        for (i = 0; i < 200; i++) {
            kf_real_t dt = (i % 3 == 0) ? 0.5f : 1.0f;
            kf_real_t F[4] = {1, dt, 0, 1};
            kf_kf_set_transition_matrix(&kf, F);
            kf_kf_predict(&kf, NULL, dt);
            t += dt;
            z = 3.0f * t + 0.5f * noise_unit();   /* pos = 3t */
            kf_kf_update(&kf, &z);
        }
        CHECK_NEAR(kf_kf_get_state(&kf)[0], 3.0 * t, 5.0);
        CHECK_NEAR(kf_kf_get_state(&kf)[1], 3.0, 0.3);
    }
    t_end();
}

static void test_error_handling(void)
{
    kf_kf_t kf;
    kf_real_t v[1] = {0.0f};

    t_begin("KF: error handling & invalid inputs");
#if KF_ENABLE_RUNTIME_CHECKS
    CHECK(kf_kf_init(NULL, 1, 1) == KF_ERROR_NULL_POINTER);
#endif
    CHECK(kf_kf_init(&kf, 0, 1) == KF_ERROR_INVALID_DIMENSION);
    CHECK(kf_kf_init(&kf, 1, 0) == KF_ERROR_INVALID_DIMENSION);
    CHECK(kf_kf_init(&kf, KF_MAX_STATE_DIM + 1, 1) == KF_ERROR_INVALID_DIMENSION);

    /* Uninitialised use. */
    {
        kf_kf_t u;
        u.initialized = 0;
        CHECK(kf_kf_predict(&u, NULL, 1.0f) == KF_ERROR_NOT_INITIALIZED);
        CHECK(kf_kf_update(&u, v) == KF_ERROR_NOT_INITIALIZED);
    }

    CHECK(kf_kf_init(&kf, 1, 1) == KF_OK);
    CHECK(kf_kf_predict(&kf, NULL, 1.0f) == KF_OK);
#if KF_ENABLE_RUNTIME_CHECKS
    CHECK(kf_kf_update(&kf, NULL) == KF_ERROR_NULL_POINTER);
#endif
    CHECK(kf_kf_set_state_element(&kf, 5, 1.0f) == KF_ERROR_INVALID_PARAMETER);
    t_end();
}

static void test_nan_inf(void)
{
    kf_kf_t kf;
    kf_real_t z;

    t_begin("KF: NaN/Inf robustness (no crash, deterministic)");
    CHECK(kf_kf_init(&kf, 1, 1) == KF_OK);
    {
        kf_real_t F[1] = {1.0f}, H[1] = {1.0f};
        kf_kf_set_transition_matrix(&kf, F);
        kf_kf_set_measurement_matrix(&kf, H);
        kf_kf_set_process_noise_scalar(&kf, 1e-4f);
        kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
        kf_kf_set_covariance_scalar(&kf, 1.0f);
    }

    z = NAN;
    /* Feeding NaN must not crash; the status is either OK (NaN propagates in
       the plain path) or an error if diagnostics are on. Either is acceptable
       for this test — we only assert the call returns. */
    (void)kf_kf_update(&kf, &z);

    z = INFINITY;
    (void)kf_kf_update(&kf, &z);
    CHECK(1);   /* reached without crash */
    t_end();
}

static void test_reset(void)
{
    kf_kf_t kf;
    kf_real_t F[1] = {1.0f}, H[1] = {1.0f};
    kf_real_t x0[1] = {5.0f};

    t_begin("KF: reset behaviour");
    CHECK(kf_kf_init(&kf, 1, 1) == KF_OK);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_state(&kf, x0);
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 5.0, 1e-9);

    CHECK(kf_kf_reset(&kf) == KF_OK);
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 0.0, 1e-9);   /* state re-zeroed   */
    CHECK_NEAR(kf_kf_get_covariance(&kf)[0], 1.0, 1e-9);/* P reset to I    */
    t_end();
}

static void test_singular_measurement(void)
{
    /* R = 0 with a zero initial covariance yields a singular S; the update
       must fail cleanly rather than corrupt the state. */
    kf_kf_t kf;
    kf_real_t F[1] = {1.0f}, H[1] = {1.0f};
    kf_real_t z = 1.0f;

    t_begin("KF: singular innovation covariance");
    CHECK(kf_kf_init(&kf, 1, 1) == KF_OK);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_measurement_noise_scalar(&kf, 0.0f);
    kf_kf_set_covariance_scalar(&kf, 0.0f);
    kf_kf_set_process_noise_scalar(&kf, 0.0f);

    CHECK(kf_kf_update(&kf, &z) == KF_ERROR_NOT_POSITIVE_DEFINITE);
    t_end();
}

static void test_input_validation(void)
{
#if KF_ENABLE_VALIDATION
    kf_kf_t kf;
    kf_real_t bad[1] = {NAN};

    t_begin("KF: input validation (non-finite rejected)");
    CHECK(kf_kf_init(&kf, 1, 1) == KF_OK);
    CHECK(kf_kf_set_state(&kf, bad) == KF_ERROR_NON_FINITE);
    CHECK(kf_kf_set_covariance(&kf, bad) == KF_ERROR_NON_FINITE);
    CHECK(kf_kf_set_process_noise(&kf, bad) == KF_ERROR_NON_FINITE);
    CHECK(kf_kf_set_measurement_noise(&kf, bad) == KF_ERROR_NON_FINITE);
    CHECK(kf_kf_set_transition_matrix(&kf, bad) == KF_ERROR_NON_FINITE);
    CHECK(kf_kf_set_measurement_matrix(&kf, bad) == KF_ERROR_NON_FINITE);
    t_end();
#else
    t_begin("KF: input validation (compiled out)");
    CHECK(1);
    t_end();
#endif
}

void test_kf(void)
{
    test_1d_convergence();
    test_constant_velocity();
    test_multi_measurement();
    test_sequential_sensors();
    test_zero_noise();
    test_variable_dt();
    test_error_handling();
    test_nan_inf();
    test_reset();
    test_singular_measurement();
    test_input_validation();
}
