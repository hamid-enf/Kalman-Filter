/**
 * @file    test_stress.c
 * @brief   Property-based stress tests: run the filters across many random
 *          configurations and verify the mathematical invariants that must
 *          hold after every step:
 *            - the covariance P stays symmetric,
 *            - P stays positive-definite (a Cholesky factorisation succeeds),
 *            - all state/covariance entries stay finite.
 *
 * These catch numerical instabilities that hand-picked examples miss.
 */

#include "test_framework.h"
#include "kalman.h"

#include <stddef.h>
#include <math.h>

static uint32_t rng_state = 0x13579BDFu;

static uint32_t rng_next(void)
{
    uint32_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_state = x;
    return x;
}

/* Uniform in [-1, 1]. */
static kf_real_t randu(void)
{
    return (kf_real_t)((int32_t)(rng_next() % 2000001u) - 1000000) / 1000000.0f;
}

/* --------------------------------------------------------------------------
 * Matrix-engine fuzz: random dimensions/data must never crash or produce
 * non-finite output (the functions return a status instead).
 * ------------------------------------------------------------------------ */
static void test_matrix_fuzz(void)
{
    kf_matrix_t A, B, C;
    kf_real_t Ab[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
    kf_real_t Bb[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
    kf_real_t Cb[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
    int trial;

    t_begin("stress: matrix engine fuzz (no crash / no non-finite output)");

    for (trial = 0; trial < 200; trial++) {
        uint16_t p = (uint16_t)(1u + (rng_next() % KF_MAX_STATE_DIM));
        uint16_t q = (uint16_t)(1u + (rng_next() % KF_MAX_STATE_DIM));
        uint16_t r = (uint16_t)(1u + (rng_next() % KF_MAX_STATE_DIM));
        uint16_t i;

        for (i = 0; i < (uint16_t)(p * q); i++) Ab[i] = randu();
        for (i = 0; i < (uint16_t)(q * r); i++) Bb[i] = randu();

        A.rows = p; A.cols = q; A.data = Ab;
        B.rows = q; B.cols = r; B.data = Bb;
        C.rows = p; C.cols = r; C.data = Cb;

        if (kf_matrix_mul(&C, &A, &B) == KF_OK) {
            CHECK(kf_matrix_is_finite(&C));
        }

        /* Symmetric random matrix: cholesky must return OK or NOT_POSITIVE_DEFINITE. */
        {
            kf_matrix_t S;
            uint16_t j;
            for (i = 0; i < p; i++) {
                for (j = 0; j < p; j++) {
                    Ab[(size_t)i * p + j] = (i == j) ? (kf_real_t)(1.0 + randu()) : randu();
                    Ab[(size_t)j * p + i] = Ab[(size_t)i * p + j];
                }
            }
            S.rows = p; S.cols = p; S.data = Ab;
            {
                kf_status_t st = kf_matrix_cholesky(&S);
                CHECK(st == KF_OK || st == KF_ERROR_NOT_POSITIVE_DEFINITE);
            }
        }
    }
    t_end();
}

/* Build a random symmetric positive-definite n x n matrix A = M M^T + eps I. */
static void rand_spd(kf_real_t *A, uint16_t n, kf_real_t scale, kf_real_t eps)
{
    kf_real_t M[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
    uint16_t i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            M[(size_t)i * n + j] = scale * randu();
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            kf_real_t s = 0.0f;
            for (k = 0; k < n; k++) {
                s += M[(size_t)i * n + k] * M[(size_t)j * n + k];
            }
            A[(size_t)i * n + j] = s + ((i == j) ? eps : 0.0f);
        }
    }
}

/* Verify P is symmetric and positive-definite (Cholesky succeeds). */
static int cov_ok(const kf_real_t *P, uint16_t n)
{
    kf_matrix_t L;
    kf_real_t buf[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
    uint16_t i;

    if (!kf_mat_is_symmetric(P, n)) {
        return 0;
    }
    L.rows = n;
    L.cols = n;
    L.data = buf;
    for (i = 0; i < (uint16_t)(n * n); i++) {
        buf[i] = P[i];
    }
    return (kf_matrix_cholesky(&L) == KF_OK) ? 1 : 0;
}

static void test_kf_stress(void)
{
    kf_kf_t kf;
    int trial;

    t_begin("stress: KF invariants across random configurations");

    for (trial = 0; trial < 40; trial++) {
        uint16_t n = (uint16_t)(1u + (rng_next() % KF_MAX_STATE_DIM));
        uint16_t m = (uint16_t)(1u + (rng_next() % KF_MAX_MEASUREMENT_DIM));
        kf_real_t F[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
        kf_real_t H[KF_MAX_MEASUREMENT_DIM * KF_MAX_STATE_DIM];
        kf_real_t Q[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
        kf_real_t R[KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM];
        kf_real_t P0[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
        kf_real_t z[KF_MAX_MEASUREMENT_DIM];
        uint16_t i, j;
        int step;

        CHECK(kf_kf_init(&kf, n, m) == KF_OK);

        /* Near-identity F with small off-diagonal coupling. */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                F[(size_t)i * n + j] = (i == j) ? 1.0f : 0.1f * randu();
            }
        }
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                H[(size_t)i * n + j] = randu();
            }
        }
        rand_spd(Q, n, 0.1f, 1e-3f);
        rand_spd(R, m, 0.3f, 1e-2f);
        rand_spd(P0, n, 1.0f, 0.1f);

        CHECK(kf_kf_set_transition_matrix(&kf, F) == KF_OK);
        CHECK(kf_kf_set_measurement_matrix(&kf, H) == KF_OK);
        CHECK(kf_kf_set_process_noise(&kf, Q) == KF_OK);
        CHECK(kf_kf_set_measurement_noise(&kf, R) == KF_OK);
        CHECK(kf_kf_set_covariance(&kf, P0) == KF_OK);

        for (step = 0; step < 30; step++) {
            for (i = 0; i < m; i++) z[i] = randu();
            CHECK(kf_kf_predict(&kf, NULL, 0.1f) == KF_OK);
            CHECK(kf_kf_update(&kf, z) == KF_OK);
            CHECK(cov_ok(kf_kf_get_covariance(&kf), n));
            CHECK(kf_vec_is_finite(kf_kf_get_state(&kf), n));
        }
    }
    t_end();
}

#if KF_ENABLE_UKF
static void ukf_f_linear(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                         kf_real_t *xo, void *ctx)
{
    /* F = I for the stress test (keeps the system stable). */
    uint16_t i;
    (void)u; (void)dt; (void)ctx;
    for (i = 0; i < KF_MAX_STATE_DIM; i++) xo[i] = x[i];
}

static void ukf_h_first(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = x[0];
}

static void test_ukf_stress(void)
{
    kf_ukf_t ukf;
    int trial;

    t_begin("stress: UKF invariants across random configurations");

    for (trial = 0; trial < 30; trial++) {
        uint16_t n = (uint16_t)(1u + (rng_next() % KF_MAX_STATE_DIM));
        kf_real_t P0[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
        kf_real_t z;
        int step;

        CHECK(kf_ukf_init(&ukf, n, 1) == KF_OK);
        CHECK(kf_ukf_set_models(&ukf, ukf_f_linear, ukf_h_first) == KF_OK);
        rand_spd(P0, n, 1.0f, 0.1f);
        kf_ukf_set_covariance(&ukf, P0);
        kf_ukf_set_process_noise_scalar(&ukf, 1e-3f);
        kf_ukf_set_measurement_noise_scalar(&ukf, 1.0f);

        for (step = 0; step < 30; step++) {
            z = randu();
            CHECK(kf_ukf_predict(&ukf, NULL, 0.1f) == KF_OK);
            CHECK(kf_ukf_update(&ukf, &z) == KF_OK);
            CHECK(cov_ok(kf_ukf_get_covariance(&ukf), n));
            CHECK(kf_vec_is_finite(kf_ukf_get_state(&ukf), n));
        }
    }
    t_end();
}
#endif /* KF_ENABLE_UKF */

void test_stress(void)
{
    test_matrix_fuzz();
    test_kf_stress();
#if KF_ENABLE_UKF
    test_ukf_stress();
#endif
}
