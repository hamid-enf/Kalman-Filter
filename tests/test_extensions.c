/**
 * @file    test_extensions.c
 * @brief   Unit tests for the KF extensions: innovation gating (NIS), adaptive
 *          measurement noise (R), and the RTS smoother.
 */

#include "test_framework.h"
#include "kalman.h"

#include <stddef.h>
#include <math.h>

#if KF_ENABLE_GATING || KF_ENABLE_ADAPTIVE_R || KF_ENABLE_SMOOTHER

static uint32_t rng_state = 0xDEADBEEFu;

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

#endif

/* --------------------------------------------------------------------------
 * Innovation gating
 * ------------------------------------------------------------------------ */

#if KF_ENABLE_GATING
static void test_gating(void)
{
    kf_kf_t kf;
    kf_real_t z;
    int i;

    t_begin("KF ext: innovation gating / NIS");
    CHECK(kf_kf_init_1d(&kf, 1e-4f, 1.0f) == KF_OK);

    /* Converge on a constant WITHOUT gating first. */
    for (i = 0; i < 100; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);
        z = 10.0f + noise_unit();
        CHECK(kf_kf_update(&kf, &z) == KF_OK);
    }
    CHECK_NEAR(kf_kf_get_state(&kf)[0], 10.0, 0.3);

    /* Enable gating only now. */
    CHECK(kf_kf_set_gate_threshold(&kf, 6.63f) == KF_OK);  /* ~99% chi-square */

    /* Normal measurement -> small NIS. */
    kf_kf_predict(&kf, NULL, 1.0f);
    z = 10.0f + noise_unit();
    CHECK(kf_kf_update_gated(&kf, &z) == KF_OK);
    CHECK(kf_kf_nis(&kf) < 6.63f);

    /* Huge outlier -> large NIS -> rejected, state unchanged. */
    {
        kf_real_t before = kf_kf_get_state(&kf)[0];
        kf_kf_predict(&kf, NULL, 1.0f);
        z = 1000.0f;
        CHECK(kf_kf_update_gated(&kf, &z) == KF_WARN_GATED);
        CHECK(kf_kf_nis(&kf) > 6.63f);
        CHECK_NEAR(kf_kf_get_state(&kf)[0], before, 1e-4);  /* unchanged */
    }
    t_end();
}
#endif

/* --------------------------------------------------------------------------
 * Adaptive R
 * ------------------------------------------------------------------------ */

#if KF_ENABLE_ADAPTIVE_R
static void test_adaptive_r(void)
{
    kf_kf_t kf;
    kf_real_t z;

    t_begin("KF ext: adaptive measurement noise");
    /* R starts small (0.01); the sensor is actually noisy (variance ~1). */
    CHECK(kf_kf_init_1d(&kf, 1e-4f, 0.01f) == KF_OK);

    /* Large innovation should drive R upward. */
    kf_kf_predict(&kf, NULL, 1.0f);
    z = 10.0f;                       /* far from x=0 -> big innovation */
    CHECK(kf_kf_update(&kf, &z) == KF_OK);
    CHECK(kf_kf_adapt_r(&kf, 0.9f, 0.01f) == KF_OK);

    CHECK(kf_kf_get_measurement_noise(&kf)[0] > 0.01f);  /* R grew       */
    CHECK(kf_kf_get_measurement_noise(&kf)[0] < 200.0f); /* but bounded  */

    /* r_min floor is enforced. */
    {
        kf_real_t R[1] = {0.001f};
        (void)kf_kf_set_measurement_noise(&kf, R);
        (void)kf_kf_adapt_r(&kf, 0.5f, 0.05f);
        CHECK(kf_kf_get_measurement_noise(&kf)[0] >= 0.05f - 1e-6f);
    }

    /* Invalid forgetting factor rejected. */
    CHECK(kf_kf_adapt_r(&kf, 1.5f, 0.01f) == KF_ERROR_INVALID_PARAMETER);
    CHECK(kf_kf_adapt_r(&kf, 0.0f, 0.01f) == KF_ERROR_INVALID_PARAMETER);
    t_end();
}
#endif

/* --------------------------------------------------------------------------
 * RTS smoother
 * ------------------------------------------------------------------------ */

#if KF_ENABLE_SMOOTHER
#define SMOOTH_N 60

static void test_smoother(void)
{
    kf_kf_t kf;
    kf_real_t x_filt[SMOOTH_N + 1];
    kf_real_t P_filt[SMOOTH_N + 1];
    kf_real_t x_pred[SMOOTH_N + 1];
    kf_real_t P_pred[SMOOTH_N + 1];
    kf_real_t x_smooth[SMOOTH_N + 1];
    kf_real_t P_smooth[SMOOTH_N + 1];
    kf_real_t F[1] = {1.0f};
    kf_real_t truth = 0.0f;
    kf_real_t filt_err = 0.0f;
    kf_real_t smooth_err = 0.0f;
    int k;

    t_begin("KF ext: RTS smoother");

    CHECK(kf_kf_init_1d(&kf, 0.01f, 1.0f) == KF_OK);
    kf_kf_set_state_element(&kf, 0u, 0.0f);

    /* Forward pass: random walk, store the trajectory. */
    x_filt[0] = 0.0f;
    P_filt[0] = 1.0f;
    for (k = 0; k < SMOOTH_N; k++) {
        kf_real_t z;
        truth += 0.5f * noise_unit();      /* true signal random-walks */
        kf_kf_predict(&kf, NULL, 1.0f);
        x_pred[k + 1] = kf_kf_get_state(&kf)[0];
        P_pred[k + 1] = kf_kf_get_covariance(&kf)[0];
        z = truth + noise_unit();
        kf_kf_update(&kf, &z);
        x_filt[k + 1] = kf_kf_get_state(&kf)[0];
        P_filt[k + 1] = kf_kf_get_covariance(&kf)[0];
        filt_err += (x_filt[k + 1] - truth) * (x_filt[k + 1] - truth);
    }

    /* Backward pass: seed with the final filtered estimate. */
    x_smooth[SMOOTH_N] = x_filt[SMOOTH_N];
    P_smooth[SMOOTH_N] = P_filt[SMOOTH_N];
    for (k = SMOOTH_N - 1; k >= 0; k--) {
        kf_status_t st = kf_kf_smooth_step(&kf,
                                           &x_filt[k], &P_filt[k], F,
                                           &x_pred[k + 1], &P_pred[k + 1],
                                           &x_smooth[k + 1], &P_smooth[k + 1],
                                           &x_smooth[k], &P_smooth[k]);
        CHECK(st == KF_OK);
    }
    smooth_err = 0.0f;
    for (k = 0; k <= SMOOTH_N; k++) {
        smooth_err += (x_smooth[k] - truth) * (x_smooth[k] - truth);
    }

    /* Smoothed covariance never exceeds filtered (smoothing reduces variance). */
    for (k = 0; k <= SMOOTH_N; k++) {
        CHECK(P_smooth[k] <= P_filt[k] + 1e-6f);
    }

    /* Smoothed trajectory has lower total squared error than filtered. */
    CHECK(smooth_err < filt_err);

    /* Error-handling: NULL inputs. */
#if KF_ENABLE_RUNTIME_CHECKS
    CHECK(kf_kf_smooth_step(&kf, NULL, &P_filt[0], F, &x_pred[1], &P_pred[1],
                            &x_smooth[1], &P_smooth[1],
                            &x_smooth[0], &P_smooth[0]) == KF_ERROR_NULL_POINTER);
#endif
    t_end();
}
#endif

void test_extensions(void)
{
#if KF_ENABLE_GATING
    test_gating();
#endif
#if KF_ENABLE_ADAPTIVE_R
    test_adaptive_r();
#endif
#if KF_ENABLE_SMOOTHER
    test_smoother();
#endif
#if !KF_ENABLE_GATING && !KF_ENABLE_ADAPTIVE_R && !KF_ENABLE_SMOOTHER
    t_begin("KF extensions: (all compiled out)");
    CHECK(1);
    t_end();
#endif
}
