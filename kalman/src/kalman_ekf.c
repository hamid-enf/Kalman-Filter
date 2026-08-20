/**
 * @file    kalman_ekf.c
 * @brief   Extended Kalman filter implementation.
 */

#include "kalman_ekf.h"

#if KF_ENABLE_EKF
#if KF_ENABLE_VALIDATION
#define KF_REQUIRE_FINITE(arr, len) do { \
        kf_status_t _kf_vs = kf_vec_validate_finite((arr), (uint16_t)(len)); \
        if (_kf_vs != KF_OK) { return _kf_vs; } \
    } while (0)
#else
#define KF_REQUIRE_FINITE(arr, len) ((void)0)
#endif


#define KF_NN  (KF_MAX_STATE_DIM * KF_MAX_STATE_DIM)
#define KF_NM  (KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM)
#define KF_MM  (KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM)

#define OFF_N1    (0u)
#define OFF_N2    (OFF_N1 + KF_MAX_STATE_DIM)
#define OFF_M1    (OFF_N2 + KF_MAX_STATE_DIM)
#define OFF_M2    (OFF_M1 + KF_MAX_MEASUREMENT_DIM)
#define OFF_NN1   (OFF_M2 + KF_MAX_MEASUREMENT_DIM)
#define OFF_NN2   (OFF_NN1 + KF_NN)
#define OFF_NM1   (OFF_NN2 + KF_NN)
#define OFF_NM2   (OFF_NM1 + KF_NM)
#define OFF_MN    (OFF_NM2 + KF_NM)
#define OFF_MM    (OFF_MN + KF_NM)
#define OFF_END   (OFF_MM + KF_MM)

#if OFF_END != KF_EKF_SCRATCH_FLOATS
#error "kalman_ekf.c: scratch layout does not match KF_EKF_SCRATCH_FLOATS"
#endif

static void kf_view(kf_matrix_t *m, kf_real_t *data, uint16_t rows, uint16_t cols)
{
    m->rows = rows;
    m->cols = cols;
    m->data = data;
}

/* --------------------------------------------------------------------------
 * Initialisation / reset
 * ------------------------------------------------------------------------ */

kf_status_t kf_ekf_init(kf_ekf_t *ekf, uint16_t n, uint16_t m)
{
    KF_NULL_CHECK(ekf);

    if (n == 0u || n > KF_MAX_STATE_DIM) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    if (m == 0u || m > KF_MAX_MEASUREMENT_DIM) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    ekf->n = n;
    ekf->m = m;
    ekf->f = NULL;
    ekf->h = NULL;
    ekf->F_jac = NULL;
    ekf->H_jac = NULL;
    ekf->ctx_f = NULL;
    ekf->ctx_h = NULL;
    ekf->initialized = 1u;
    return kf_ekf_reset(ekf);
}

kf_status_t kf_ekf_reset(kf_ekf_t *ekf)
{
    kf_matrix_t M;
    uint16_t n, m;

    KF_NULL_CHECK(ekf);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    n = ekf->n;
    m = ekf->m;

    kf_vec_zero(ekf->x, n);
    kf_view(&M, ekf->P, n, n);  kf_matrix_identity(&M);
    kf_view(&M, ekf->Q, n, n);  kf_matrix_zero(&M);
    kf_view(&M, ekf->F, n, n);  kf_matrix_identity(&M);
    kf_view(&M, ekf->R, m, m);  kf_matrix_zero(&M);
    kf_view(&M, ekf->H, m, n);  kf_matrix_zero(&M);

#if KF_ENABLE_ADVANCED_API
    {
        uint16_t i;
        kf_vec_zero(ekf->y, m);
        for (i = 0u; i < (uint16_t)(n * m); i++) {
            ekf->K[i] = (kf_real_t)0;
        }
    }
#endif

    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Model registration
 * ------------------------------------------------------------------------ */

kf_status_t kf_ekf_set_models(kf_ekf_t *ekf,
                              kf_ekf_f_fn f, kf_ekf_jacobian_f_fn F_jac,
                              kf_ekf_h_fn h, kf_ekf_jacobian_h_fn H_jac)
{
    KF_NULL_CHECK(ekf);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    ekf->f = f;
    ekf->F_jac = F_jac;
    ekf->h = h;
    ekf->H_jac = H_jac;
    return KF_OK;
}

void kf_ekf_set_context_f(kf_ekf_t *ekf, void *ctx)
{
    if (ekf != NULL) {
        ekf->ctx_f = ctx;
    }
}

void kf_ekf_set_context_h(kf_ekf_t *ekf, void *ctx)
{
    if (ekf != NULL) {
        ekf->ctx_h = ctx;
    }
}

/* --------------------------------------------------------------------------
 * State / covariance / noise setters
 * ------------------------------------------------------------------------ */

kf_status_t kf_ekf_set_state(kf_ekf_t *ekf, const kf_real_t *x)
{
    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(x);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(x, ekf->n);
    kf_vec_copy(ekf->x, x, ekf->n);
    return KF_OK;
}

const kf_real_t *kf_ekf_get_state(const kf_ekf_t *ekf)
{
    return (ekf == NULL) ? NULL : ekf->x;
}

kf_status_t kf_ekf_set_covariance(kf_ekf_t *ekf, const kf_real_t *P)
{
    uint16_t i;
    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(P);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(P, ekf->n * ekf->n);
    for (i = 0u; i < (uint16_t)(ekf->n * ekf->n); i++) {
        ekf->P[i] = P[i];
    }
    return KF_OK;
}

kf_status_t kf_ekf_set_covariance_diagonal(kf_ekf_t *ekf, const kf_real_t *diag)
{
    kf_matrix_t M;
    uint16_t i;
    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(diag);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ekf->P, ekf->n, ekf->n);
    kf_matrix_zero(&M);
    for (i = 0u; i < ekf->n; i++) {
        kf_matrix_set(&M, i, i, diag[i]);
    }
    return KF_OK;
}

kf_status_t kf_ekf_set_covariance_scalar(kf_ekf_t *ekf, kf_real_t p)
{
    kf_matrix_t M;
    KF_NULL_CHECK(ekf);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ekf->P, ekf->n, ekf->n);
    kf_matrix_identity(&M);
    kf_matrix_scale(&M, p);
    return KF_OK;
}

const kf_real_t *kf_ekf_get_covariance(const kf_ekf_t *ekf)
{
    return (ekf == NULL) ? NULL : ekf->P;
}

kf_status_t kf_ekf_set_process_noise(kf_ekf_t *ekf, const kf_real_t *Q)
{
    uint16_t i;
    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(Q);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(Q, ekf->n * ekf->n);
    for (i = 0u; i < (uint16_t)(ekf->n * ekf->n); i++) {
        ekf->Q[i] = Q[i];
    }
    return KF_OK;
}

kf_status_t kf_ekf_set_process_noise_diagonal(kf_ekf_t *ekf, const kf_real_t *diag)
{
    kf_matrix_t M;
    uint16_t i;
    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(diag);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ekf->Q, ekf->n, ekf->n);
    kf_matrix_zero(&M);
    for (i = 0u; i < ekf->n; i++) {
        kf_matrix_set(&M, i, i, diag[i]);
    }
    return KF_OK;
}

kf_status_t kf_ekf_set_process_noise_scalar(kf_ekf_t *ekf, kf_real_t q)
{
    kf_matrix_t M;
    KF_NULL_CHECK(ekf);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ekf->Q, ekf->n, ekf->n);
    kf_matrix_identity(&M);
    kf_matrix_scale(&M, q);
    return KF_OK;
}

kf_status_t kf_ekf_set_measurement_noise(kf_ekf_t *ekf, const kf_real_t *R)
{
    uint16_t i;
    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(R);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(R, ekf->m * ekf->m);
    for (i = 0u; i < (uint16_t)(ekf->m * ekf->m); i++) {
        ekf->R[i] = R[i];
    }
    return KF_OK;
}

kf_status_t kf_ekf_set_measurement_noise_diagonal(kf_ekf_t *ekf, const kf_real_t *diag)
{
    kf_matrix_t M;
    uint16_t i;
    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(diag);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ekf->R, ekf->m, ekf->m);
    kf_matrix_zero(&M);
    for (i = 0u; i < ekf->m; i++) {
        kf_matrix_set(&M, i, i, diag[i]);
    }
    return KF_OK;
}

kf_status_t kf_ekf_set_measurement_noise_scalar(kf_ekf_t *ekf, kf_real_t r)
{
    kf_matrix_t M;
    KF_NULL_CHECK(ekf);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ekf->R, ekf->m, ekf->m);
    kf_matrix_identity(&M);
    kf_matrix_scale(&M, r);
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Predict
 * ------------------------------------------------------------------------ */

kf_status_t kf_ekf_predict(kf_ekf_t *ekf, const kf_real_t *u, kf_real_t dt)
{
    kf_matrix_t Fm, Pm, Qm, tmp;
    kf_real_t *s;
    uint16_t n;

    KF_NULL_CHECK(ekf);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    if (ekf->f == NULL || ekf->F_jac == NULL) {
        return KF_ERROR_INVALID_PARAMETER;
    }
    n = ekf->n;
    s = ekf->scratch;

    /* Linearise: F = df/dx evaluated at the current state. */
    kf_view(&Fm, ekf->F, n, n);
    ekf->F_jac(ekf->x, u, dt, &Fm, ekf->ctx_f);

    /* Propagate the state: x = f(x, u, dt). */
    ekf->f(ekf->x, u, dt, &s[OFF_N1], ekf->ctx_f);

    /* Covariance: P = F P F^T + Q. */
    kf_view(&Pm, ekf->P, n, n);
    kf_view(&Qm, ekf->Q, n, n);
    kf_view(&tmp, &s[OFF_NN1], n, n);
    kf_matrix_mul(&tmp, &Fm, &Pm);
    kf_matrix_mul_transpose_b(&Pm, &tmp, &Fm);
    kf_matrix_add(&Pm, &Pm, &Qm);
    kf_matrix_symmetrize(&Pm);

    /* Commit the propagated state. */
    kf_vec_copy(ekf->x, &s[OFF_N1], n);

    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Update
 * ------------------------------------------------------------------------ */

kf_status_t kf_ekf_update(kf_ekf_t *ekf, const kf_real_t *z)
{
    kf_matrix_t Hm, Pm, Rm, Bm, Sm, Km, Tm, T2, KR, RHS;
    kf_real_t *s;
    kf_status_t st;
    uint16_t n, m, r, c;

    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(z);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    if (ekf->h == NULL || ekf->H_jac == NULL) {
        return KF_ERROR_INVALID_PARAMETER;
    }
    n = ekf->n;
    m = ekf->m;
    s = ekf->scratch;

    /* Linearise: H = dh/dx evaluated at the current state. */
    kf_view(&Hm, ekf->H, m, n);
    ekf->H_jac(ekf->x, &Hm, ekf->ctx_h);

    /* Predicted measurement z_hat = h(x), then innovation y = z - z_hat. */
    ekf->h(ekf->x, &s[OFF_M1], ekf->ctx_h);
    kf_vec_sub(&s[OFF_M1], z, &s[OFF_M1], m);
#if KF_ENABLE_ADVANCED_API
    kf_vec_copy(ekf->y, &s[OFF_M1], m);
#endif

    kf_view(&Pm, ekf->P, n, n);
    kf_view(&Rm, ekf->R, m, m);

    kf_view(&Bm, &s[OFF_NM1], n, m);
    kf_matrix_mul_transpose_b(&Bm, &Pm, &Hm);   /* B = P H^T        */

    kf_view(&Sm, &s[OFF_MM], m, m);
    kf_matrix_mul(&Sm, &Hm, &Bm);               /* S = H B          */
    kf_matrix_add(&Sm, &Sm, &Rm);               /* S += R           */

    st = kf_matrix_cholesky(&Sm);
    if (st != KF_OK) {
        return st;
    }
    kf_view(&RHS, &s[OFF_MN], m, n);
    kf_matrix_transpose(&RHS, &Bm);
    st = kf_matrix_cholesky_solve(&Sm, &RHS);   /* RHS = K^T        */
    if (st != KF_OK) {
        return st;
    }
    kf_view(&Km, &s[OFF_NM1], n, m);
    kf_matrix_transpose(&Km, &RHS);             /* K = (K^T)^T      */

#if KF_ENABLE_ADVANCED_API
    {
        uint16_t i;
        for (i = 0u; i < (uint16_t)(n * m); i++) {
            ekf->K[i] = s[OFF_NM1 + i];
        }
    }
#endif

    /* x = x + K y */
    kf_mat_vec_mul(&s[OFF_N2], &Km, &s[OFF_M1]);
    kf_vec_add(ekf->x, ekf->x, &s[OFF_N2], n);

    /* Covariance update */
    kf_view(&Tm, &s[OFF_NN1], n, n);
    kf_matrix_mul(&Tm, &Km, &Hm);
    for (r = 0u; r < n; r++) {
        for (c = 0u; c < n; c++) {
            kf_real_t ident = (r == c) ? (kf_real_t)1 : (kf_real_t)0;
            kf_matrix_set(&Tm, r, c, ident - kf_matrix_get(&Tm, r, c));
        }
    }

#if KF_USE_JOSEPH_FORM
    kf_view(&T2, &s[OFF_NN2], n, n);
    kf_matrix_mul(&T2, &Tm, &Pm);
    kf_matrix_mul_transpose_b(&Pm, &T2, &Tm);
    kf_view(&KR, &s[OFF_NM2], n, m);
    kf_matrix_mul(&KR, &Km, &Rm);
    kf_matrix_mul_transpose_b(&T2, &KR, &Km);
    kf_matrix_add(&Pm, &Pm, &T2);
#else
    kf_view(&T2, &s[OFF_NN2], n, n);
    kf_matrix_mul(&T2, &Tm, &Pm);
    kf_matrix_copy(&Pm, &T2);
#endif
    kf_matrix_symmetrize(&Pm);

    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Advanced / diagnostics
 * ------------------------------------------------------------------------ */

#if KF_ENABLE_ADVANCED_API
kf_status_t kf_ekf_get_innovation(const kf_ekf_t *ekf, kf_real_t *y_out)
{
    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(y_out);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_vec_copy(y_out, ekf->y, ekf->m);
    return KF_OK;
}

kf_status_t kf_ekf_get_gain(const kf_ekf_t *ekf, kf_real_t *K_out)
{
    uint16_t i;
    KF_NULL_CHECK(ekf);
    KF_NULL_CHECK(K_out);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    for (i = 0u; i < (uint16_t)(ekf->n * ekf->m); i++) {
        K_out[i] = ekf->K[i];
    }
    return KF_OK;
}
#endif

#if KF_ENABLE_DIAGNOSTICS
kf_status_t kf_ekf_check(const kf_ekf_t *ekf)
{
    kf_matrix_t M;
    uint16_t i;

    KF_NULL_CHECK(ekf);

    if (ekf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    for (i = 0u; i < ekf->n; i++) {
        if (!kf_isfinite(ekf->x[i])) {
            return KF_ERROR_NON_FINITE;
        }
    }
    kf_view(&M, (kf_real_t *)ekf->P, ekf->n, ekf->n);
    if (!kf_matrix_is_finite(&M)) {
        return KF_ERROR_NON_FINITE;
    }
    if (!kf_matrix_is_symmetric(&M)) {
        return KF_ERROR_NUMERICAL;
    }
    return KF_OK;
}
#endif

#endif /* KF_ENABLE_EKF */
