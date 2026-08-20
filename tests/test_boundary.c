/**
 * @file    test_boundary.c
 * @brief   Boundary / untested-path coverage: control-input (u) forwarding,
 *          maximum dimensions, and UKF small-alpha sigma points.
 *
 * These run under AddressSanitizer in CI, so any off-by-one in the internal
 * scratch/indexing at the maximum configured dimensions is caught as a memory
 * error rather than silently corrupting adjacent state.
 */

#include "test_framework.h"
#include "kalman.h"

#include <stddef.h>
#include <math.h>

static uint32_t rng_state = 0xB0B0B0B0u;

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

/* --------------------------------------------------------------------------
 * Control input u: x <- F x + u must accumulate u exactly for F = I.
 * ------------------------------------------------------------------------ */
static void test_control_input(void)
{
    kf_kf_t kf;
    kf_real_t F[4] = {1, 0, 0, 1};   /* F = I */
    kf_real_t H[2] = {1, 0};
    kf_real_t u[2] = {0.5f, -0.25f};
    kf_real_t x0[2] = {1.0f, 2.0f};
    int i;

    t_begin("boundary: control input (u) accumulates in predict");

    CHECK(kf_kf_init(&kf, 2, 1) == KF_OK);
    CHECK(kf_kf_set_transition_matrix(&kf, F) == KF_OK);
    CHECK(kf_kf_set_measurement_matrix(&kf, H) == KF_OK);
    CHECK(kf_kf_set_process_noise_scalar(&kf, 0.0f) == KF_OK);
    CHECK(kf_kf_set_measurement_noise_scalar(&kf, 1.0f) == KF_OK);
    CHECK(kf_kf_set_state(&kf, x0) == KF_OK);

    for (i = 0; i < 10; i++) {
        CHECK(kf_kf_predict(&kf, u, 1.0f) == KF_OK);
    }
    /* x = x0 + 10 * u  exactly. */
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 1.0f + 10 * 0.5f, 1e-5);
    CHECK_NEAR(kf_kf_get_state(&kf)[1], 2.0f + 10 * (-0.25f), 1e-5);

    /* NULL control input is a no-op (x unchanged by input term). */
    CHECK(kf_kf_predict(&kf, NULL, 1.0f) == KF_OK);
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 1.0f + 10 * 0.5f, 1e-5);
    t_end();
}

#if KF_ENABLE_EKF
/* EKF: verify the control input is forwarded to the transition callback. */
static void ekf_f_accel(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                        kf_real_t *xo, void *ctx)
{
    (void)ctx;
    xo[0] = x[0] + x[1] * dt;
    xo[1] = x[1] + ((u != NULL) ? u[1] * dt : 0.0f);   /* a = u[1] */
}

static void ekf_F_accel(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                        kf_matrix_t *F, void *ctx)
{
    (void)x; (void)u; (void)ctx;
    kf_matrix_set(F, 0, 0, 1); kf_matrix_set(F, 0, 1, dt);
    kf_matrix_set(F, 1, 0, 0); kf_matrix_set(F, 1, 1, 1);
}

static void ekf_h_pos(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = x[0];
}

static void ekf_H_pos(const kf_real_t *x, kf_matrix_t *H, void *ctx)
{
    (void)x; (void)ctx;
    kf_matrix_set(H, 0, 0, 1); kf_matrix_set(H, 0, 1, 0);
}

static void test_ekf_control_input(void)
{
    kf_ekf_t ekf;
    kf_real_t x0[2] = {0, 0};
    kf_real_t u[2] = {0, 1.0f};      /* constant acceleration = 1 */
    int i;

    t_begin("boundary: EKF forwards control input u");
    CHECK(kf_ekf_init(&ekf, 2, 1) == KF_OK);
    CHECK(kf_ekf_set_models(&ekf, ekf_f_accel, ekf_F_accel, ekf_h_pos, ekf_H_pos) == KF_OK);
    kf_ekf_set_process_noise_scalar(&ekf, 1e-6f);
    kf_ekf_set_measurement_noise_scalar(&ekf, 0.01f);
    kf_ekf_set_state(&ekf, x0);

    for (i = 0; i < 20; i++) {
        CHECK(kf_ekf_predict(&ekf, u, 1.0f) == KF_OK);
    }
    /* v = 20 * 1.0 = 20; p = 0+1+...+19 = 190 (velocity starts at 0). */
    CHECK_NEAR(kf_ekf_get_state(&ekf)[1], 20.0, 1e-3);
    CHECK_NEAR(kf_ekf_get_state(&ekf)[0], 190.0, 0.1);
    t_end();
}
#endif /* KF_ENABLE_EKF */

/* --------------------------------------------------------------------------
 * Maximum dimensions: run a KF at n=KF_MAX_STATE_DIM, m=KF_MAX_MEASUREMENT_DIM
 * and a UKF at n=KF_MAX_STATE_DIM to exercise the full scratch/index range.
 * ------------------------------------------------------------------------ */
static void test_max_dimensions(void)
{
    kf_kf_t kf;
    kf_real_t F[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
    kf_real_t H[KF_MAX_MEASUREMENT_DIM * KF_MAX_STATE_DIM];
    kf_real_t z[KF_MAX_MEASUREMENT_DIM];
    uint16_t i;
    int k;

    t_begin("boundary: KF at maximum dimensions");

    CHECK(kf_kf_init(&kf, KF_MAX_STATE_DIM, KF_MAX_MEASUREMENT_DIM) == KF_OK);

    for (i = 0; i < (uint16_t)(KF_MAX_STATE_DIM * KF_MAX_STATE_DIM); i++) {
        F[i] = 0.0f;
    }
    for (i = 0; i < KF_MAX_STATE_DIM; i++) {
        F[(size_t)i * KF_MAX_STATE_DIM + i] = 1.0f;
    }
    for (i = 0; i < (uint16_t)(KF_MAX_MEASUREMENT_DIM * KF_MAX_STATE_DIM); i++) {
        H[i] = 0.0f;
    }
    for (i = 0; i < KF_MAX_MEASUREMENT_DIM; i++) {
        H[(size_t)i * KF_MAX_STATE_DIM + i] = 1.0f;   /* measure first m states */
    }

    CHECK(kf_kf_set_transition_matrix(&kf, F) == KF_OK);
    CHECK(kf_kf_set_measurement_matrix(&kf, H) == KF_OK);
    CHECK(kf_kf_set_process_noise_scalar(&kf, 1e-3f) == KF_OK);
    CHECK(kf_kf_set_measurement_noise_scalar(&kf, 1.0f) == KF_OK);

    for (k = 0; k < 100; k++) {
        for (i = 0; i < KF_MAX_MEASUREMENT_DIM; i++) {
            z[i] = 1.0f + noise_unit();
        }
        CHECK(kf_kf_predict(&kf, NULL, 1.0f) == KF_OK);
        CHECK(kf_kf_update(&kf, z) == KF_OK);
    }
    /* First KF_MAX_MEASUREMENT_DIM states converge toward 1.0. */
    for (i = 0; i < KF_MAX_MEASUREMENT_DIM; i++) {
        CHECK_NEAR(kf_kf_get_state(&kf)[i], 1.0, 0.3);
    }
    t_end();
}

#if KF_ENABLE_UKF
static void ukf_f_ident(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                        kf_real_t *xo, void *ctx)
{
    uint16_t i;
    (void)u; (void)dt; (void)ctx;
    for (i = 0; i < KF_MAX_STATE_DIM; i++) xo[i] = x[i];
}

static void ukf_h_first(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = x[0];
}

static void test_ukf_max_dimensions(void)
{
    kf_ukf_t ukf;
    kf_real_t z[1];
    int k;

    t_begin("boundary: UKF at maximum state dimension");
    CHECK(kf_ukf_init(&ukf, KF_MAX_STATE_DIM, 1) == KF_OK);
    CHECK(kf_ukf_set_models(&ukf, ukf_f_ident, ukf_h_first) == KF_OK);
    kf_ukf_set_process_noise_scalar(&ukf, 1e-4f);
    kf_ukf_set_measurement_noise_scalar(&ukf, 1.0f);

    for (k = 0; k < 50; k++) {
        z[0] = 2.0f + noise_unit();
        CHECK(kf_ukf_predict(&ukf, NULL, 1.0f) == KF_OK);
        CHECK(kf_ukf_update(&ukf, &z[0]) == KF_OK);
    }
    CHECK_NEAR(kf_ukf_get_state(&ukf)[0], 2.0, 0.4);
    t_end();
}

static void test_ukf_small_alpha(void)
{
    kf_ukf_t ukf;
    kf_real_t x0[1] = {0.0f};
    kf_real_t z;
    int i;

    t_begin("boundary: UKF with small alpha (tight sigma points)");
    CHECK(kf_ukf_init(&ukf, 1, 1) == KF_OK);
    CHECK(kf_ukf_set_models(&ukf, ukf_f_ident, ukf_h_first) == KF_OK);
    /* Small alpha -> large negative Wm[0]; must still run without NaN. */
    CHECK(kf_ukf_set_parameters(&ukf, 1e-2f, 2.0f, 0.0f) == KF_OK);
    kf_ukf_set_process_noise_scalar(&ukf, 1e-4f);
    kf_ukf_set_measurement_noise_scalar(&ukf, 1.0f);
    kf_ukf_set_covariance_scalar(&ukf, 1.0f);
    kf_ukf_set_state(&ukf, x0);

    for (i = 0; i < 200; i++) {
        z = 5.0f + noise_unit();
        CHECK(kf_ukf_predict(&ukf, NULL, 1.0f) == KF_OK);
        CHECK(kf_ukf_update(&ukf, &z) == KF_OK);
    }
    CHECK_NEAR(kf_ukf_get_state(&ukf)[0], 5.0, 0.4);
    t_end();
}
#endif /* KF_ENABLE_UKF */

void test_boundary(void)
{
    test_control_input();
#if KF_ENABLE_EKF
    test_ekf_control_input();
#endif
    test_max_dimensions();
#if KF_ENABLE_UKF
    test_ukf_max_dimensions();
    test_ukf_small_alpha();
#endif
}
