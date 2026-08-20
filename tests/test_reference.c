/**
 * @file    test_reference.c
 * @brief   Reference-comparison tests: the library is checked against naive
 *          double-precision implementations of the same recursions, to catch
 *          any subtle algorithmic/numerical bug that convergence tests might
 *          miss (e.g. a wrong Joseph form, gain solve, or UKF reconstruction).
 */

#include "test_framework.h"
#include "kalman.h"

#include <stddef.h>
#include <string.h>

/* --------------------------------------------------------------------------
 * Naive double-precision KF reference (predict + Joseph update).
 * ------------------------------------------------------------------------ */
static void ref_kf_step(int n, int m,
                        double *x, double *P,
                        const double *F, const double *Q,
                        const double *H, const double *R,
                        const double *z)
{
    double xp[8];
    double P1[64], P2[64], PHt[64], S[64], Sinv[64], K[64];
    double KH[64], T[64], Pp[64], KR[64], KRK[64], y[8];
    double Saug[8][16];
    int i, j, k;

    for (i = 0; i < n; i++) { double s = 0; for (k = 0; k < n; k++) s += F[i*n+k]*x[k]; xp[i] = s; }
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
        double s = 0; for (k = 0; k < n; k++) s += F[i*n+k]*P[k*n+j]; P1[i*n+j] = s;
    }
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
        double s = 0; for (k = 0; k < n; k++) s += P1[i*n+k]*F[j*n+k]; P2[i*n+j] = s + Q[i*n+j];
    }
    for (i = 0; i < m; i++) { double s = 0; for (k = 0; k < n; k++) s += H[i*n+k]*xp[k]; y[i] = z[i] - s; }
    for (i = 0; i < n; i++) for (j = 0; j < m; j++) {
        double s = 0; for (k = 0; k < n; k++) s += P2[i*n+k]*H[j*n+k]; PHt[i*m+j] = s;
    }
    for (i = 0; i < m; i++) for (j = 0; j < m; j++) {
        double s = 0; for (k = 0; k < n; k++) s += H[i*n+k]*PHt[k*m+j]; S[i*m+j] = s + R[i*m+j];
    }
    /* Gauss-Jordan inverse of S (m x m). */
    memset(Saug, 0, sizeof(Saug));
    for (i = 0; i < m; i++) for (j = 0; j < m; j++) Saug[i][j] = S[i*m+j];
    for (i = 0; i < m; i++) Saug[i][i+m] = 1.0;
    for (i = 0; i < m; i++) {
        double piv = Saug[i][i];
        for (j = 0; j < 2*m; j++) Saug[i][j] /= piv;
        for (k = 0; k < m; k++) if (k != i) {
            double f = Saug[k][i];
            for (j = 0; j < 2*m; j++) Saug[k][j] -= f * Saug[i][j];
        }
    }
    for (i = 0; i < m; i++) for (j = 0; j < m; j++) Sinv[i*m+j] = Saug[i][j+m];
    for (i = 0; i < n; i++) for (j = 0; j < m; j++) {
        double s = 0; for (k = 0; k < m; k++) s += PHt[i*m+k]*Sinv[k*m+j]; K[i*m+j] = s;
    }
    for (i = 0; i < n; i++) { double s = 0; for (k = 0; k < m; k++) s += K[i*m+k]*y[k]; x[i] = xp[i] + s; }
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
        double s = 0; for (k = 0; k < m; k++) s += K[i*m+k]*H[k*n+j];
        KH[i*n+j] = (i == j ? 1.0 : 0.0) - s;
    }
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
        double s = 0; for (k = 0; k < n; k++) s += KH[i*n+k]*P2[k*n+j]; T[i*n+j] = s;
    }
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
        double s = 0; for (k = 0; k < n; k++) s += T[i*n+k]*KH[j*n+k]; Pp[i*n+j] = s;
    }
    for (i = 0; i < n; i++) for (j = 0; j < m; j++) {
        double s = 0; for (k = 0; k < m; k++) s += K[i*m+k]*R[k*m+j]; KR[i*m+j] = s;
    }
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
        double s = 0; for (k = 0; k < m; k++) s += KR[i*m+k]*K[j*m+k]; KRK[i*n+j] = s;
    }
    for (i = 0; i < n*n; i++) P[i] = Pp[i] + KRK[i];
}

/* --------------------------------------------------------------------------
 * Test: linear KF (n=3, m=2) must match the double-precision reference.
 * ------------------------------------------------------------------------ */
static void test_kf_reference(void)
{
    kf_kf_t kf;
    double x[3] = {0.5, -0.2, 1.3}, P[9] = {2, 0.1, 0, 0.1, 3, 0, 0, 0, 1};
    double F[9] = {1, 0.1, 0, 0, 1, 0.1, 0, 0, 1};
    double Q[9] = {0.01, 0, 0, 0, 0.02, 0, 0, 0, 0.03};
    double H[6] = {1, 0, 0, 0, 1, 0};
    double R[4] = {0.5, 0, 0, 0.8};
    kf_real_t xf[3], Pf[9], Ff[9], Qf[9], Hf[6], Rf[4];
    int i, k;

    t_begin("reference: KF (3x2) vs double-precision implementation");

    for (i = 0; i < 3; i++) xf[i] = (kf_real_t)x[i];
    for (i = 0; i < 9; i++) { Pf[i] = (kf_real_t)P[i]; Ff[i] = (kf_real_t)F[i]; Qf[i] = (kf_real_t)Q[i]; }
    for (i = 0; i < 6; i++) Hf[i] = (kf_real_t)H[i];
    for (i = 0; i < 4; i++) Rf[i] = (kf_real_t)R[i];

    CHECK(kf_kf_init(&kf, 3, 2) == KF_OK);
    CHECK(kf_kf_set_state(&kf, xf) == KF_OK);
    CHECK(kf_kf_set_covariance(&kf, Pf) == KF_OK);
    CHECK(kf_kf_set_transition_matrix(&kf, Ff) == KF_OK);
    CHECK(kf_kf_set_process_noise(&kf, Qf) == KF_OK);
    CHECK(kf_kf_set_measurement_matrix(&kf, Hf) == KF_OK);
    CHECK(kf_kf_set_measurement_noise(&kf, Rf) == KF_OK);

    for (k = 0; k < 50; k++) {
        double z[2] = {1.0 + 0.01 * k, -0.5 + 0.005 * k};
        kf_real_t zf[2] = {(kf_real_t)z[0], (kf_real_t)z[1]};
        CHECK(kf_kf_predict(&kf, NULL, 0.1) == KF_OK);
        CHECK(kf_kf_update(&kf, zf) == KF_OK);
        ref_kf_step(3, 2, x, P, F, Q, H, R, z);
    }

    {
        const kf_real_t *xs = kf_kf_get_state(&kf);
        const kf_real_t *Ps = kf_kf_get_covariance(&kf);
        for (i = 0; i < 3; i++) CHECK_NEAR(xs[i], x[i], 2e-3);
        for (i = 0; i < 9; i++) CHECK_NEAR(Ps[i], P[i], 2e-2);
    }
    t_end();
}

/* --------------------------------------------------------------------------
 * Test: UKF reconstructs P exactly for a linear map (unscented transform is
 * exact for linear systems, so a single predict with f(x)=A x, Q=0 must give
 * P = A P A^T to machine precision).
 * ------------------------------------------------------------------------ */
#if KF_ENABLE_UKF
static void f_scale2(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                     kf_real_t *xo, void *ctx)
{
    (void)u; (void)dt; (void)ctx;
    xo[0] = 2.0f * x[0] + 0.5f * x[1];
    xo[1] = -1.0f * x[1];
}

static void h_skip(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = x[0];
}

static void test_ukf_sigma_reconstruction(void)
{
    kf_ukf_t ukf;
    kf_real_t x0[2] = {1.0f, -2.0f};
    kf_real_t P0[4] = {3.0f, 0.4f, 0.4f, 2.0f};
    /* A = [[2, 0.5],[0, -1]]; expected P' = A P A^T. */
    kf_real_t Pexp[4];
    kf_real_t A[4] = {2.0f, 0.5f, 0.0f, -1.0f};

    t_begin("reference: UKF sigma points reconstruct covariance exactly");

    /* Pexp = A P A^T */
    {
        int i, j, kk;
        for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) {
            double s = 0;
            for (kk = 0; kk < 2; kk++) {
                double t = 0; int l;
                for (l = 0; l < 2; l++) t += (double)A[i*2+l] * (double)P0[l*2+kk];
                s += t * (double)A[j*2+kk];
            }
            Pexp[i*2+j] = (kf_real_t)s;
        }
    }

    CHECK(kf_ukf_init(&ukf, 2, 1) == KF_OK);
    CHECK(kf_ukf_set_models(&ukf, f_scale2, h_skip) == KF_OK);
    CHECK(kf_ukf_set_state(&ukf, x0) == KF_OK);
    CHECK(kf_ukf_set_covariance(&ukf, P0) == KF_OK);
    CHECK(kf_ukf_set_process_noise_scalar(&ukf, 0.0f) == KF_OK);
    CHECK(kf_ukf_set_measurement_noise_scalar(&ukf, 1.0f) == KF_OK);

    /* One predict (no update): P must become A P A^T exactly. */
    CHECK(kf_ukf_predict(&ukf, NULL, 1.0f) == KF_OK);
    {
        const kf_real_t *Ps = kf_ukf_get_covariance(&ukf);
        int i;
        for (i = 0; i < 4; i++) CHECK_NEAR(Ps[i], Pexp[i], 1e-3);
    }
    /* State mean must be A x exactly. */
    CHECK_NEAR(kf_ukf_get_state(&ukf)[0], 2.0f * 1.0f + 0.5f * (-2.0f), 1e-4);
    CHECK_NEAR(kf_ukf_get_state(&ukf)[1], -1.0f * (-2.0f), 1e-4);
    t_end();
}
#endif /* KF_ENABLE_UKF */

void test_reference(void)
{
    test_kf_reference();
#if KF_ENABLE_UKF
    test_ukf_sigma_reconstruction();
#endif
}
