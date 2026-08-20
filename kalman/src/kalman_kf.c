/**
 * @file    kalman_kf.c
 * @brief   Linear Kalman filter implementation.
 */

#include "kalman_kf.h"

#if KF_ENABLE_KF
#if KF_ENABLE_VALIDATION
#define KF_REQUIRE_FINITE(arr, len) do { \
        kf_status_t _kf_vs = kf_vec_validate_finite((arr), (uint16_t)(len)); \
        if (_kf_vs != KF_OK) { return _kf_vs; } \
    } while (0)
#else
#define KF_REQUIRE_FINITE(arr, len) ((void)0)
#endif


/* --------------------------------------------------------------------------
 * Scratch buffer layout (offsets in kf_real_t elements, sized for the maximum
 * configured dimensions so all offsets are compile-time constants).
 * ------------------------------------------------------------------------ */

#define KF_NN  (KF_MAX_STATE_DIM * KF_MAX_STATE_DIM)
#define KF_NM  (KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM)
#define KF_MM  (KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM)

#define OFF_VEC_N  (0u)
#define OFF_VEC_M  (OFF_VEC_N + KF_MAX_STATE_DIM)
#define OFF_NN1    (OFF_VEC_M + KF_MAX_MEASUREMENT_DIM)
#define OFF_NN2    (OFF_NN1 + KF_NN)
#define OFF_NM1    (OFF_NN2 + KF_NN)
#define OFF_NM2    (OFF_NM1 + KF_NM)
#define OFF_MN     (OFF_NM2 + KF_NM)
#define OFF_MM     (OFF_MN + KF_NM)
#define OFF_END    (OFF_MM + KF_MM)

#if OFF_END != KF_KF_SCRATCH_FLOATS
#error "kalman_kf.c: scratch layout does not match KF_KF_SCRATCH_FLOATS"
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

kf_status_t kf_kf_init(kf_kf_t *kf, uint16_t n, uint16_t m)
{
    KF_NULL_CHECK(kf);

    if (n == 0u || n > KF_MAX_STATE_DIM) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    if (m == 0u || m > KF_MAX_MEASUREMENT_DIM) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    kf->n = n;
    kf->m = m;
    kf->initialized = 1u;
    return kf_kf_reset(kf);
}

kf_status_t kf_kf_reset(kf_kf_t *kf)
{
    kf_matrix_t M;
    uint16_t n;
    uint16_t m;

    KF_NULL_CHECK(kf);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    n = kf->n;
    m = kf->m;

    kf_vec_zero(kf->x, n);

    kf_view(&M, kf->P, n, n);
    (void)kf_matrix_identity(&M);
    kf_view(&M, kf->Q, n, n);
    (void)kf_matrix_zero(&M);
    kf_view(&M, kf->F, n, n);
    (void)kf_matrix_identity(&M);
    kf_view(&M, kf->R, m, m);
    (void)kf_matrix_zero(&M);
    kf_view(&M, kf->H, m, n);
    (void)kf_matrix_zero(&M);

#if KF_ENABLE_ADVANCED_API
    kf_vec_zero(kf->y, m);
    {
        uint16_t i;
        for (i = 0u; i < (uint16_t)(n * m); i++) {
            kf->K[i] = (kf_real_t)0;
        }
    }
#endif

    return KF_OK;
}

/* --------------------------------------------------------------------------
 * State access
 * ------------------------------------------------------------------------ */

kf_status_t kf_kf_set_state(kf_kf_t *kf, const kf_real_t *x)
{
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(x);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(x, kf->n);
    kf_vec_copy(kf->x, x, kf->n);
    return KF_OK;
}

const kf_real_t *kf_kf_get_state(const kf_kf_t *kf)
{
    if (kf == NULL) {
        return NULL;
    }
    return kf->x;
}

kf_status_t kf_kf_set_state_element(kf_kf_t *kf, uint16_t index, kf_real_t value)
{
    KF_NULL_CHECK(kf);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    if (index >= kf->n) {
        return KF_ERROR_INVALID_PARAMETER;
    }
    kf->x[index] = value;
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Covariance access
 * ------------------------------------------------------------------------ */

kf_status_t kf_kf_set_covariance(kf_kf_t *kf, const kf_real_t *P)
{
    uint16_t i;
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(P);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(P, kf->n * kf->n);
    for (i = 0u; i < (uint16_t)(kf->n * kf->n); i++) {
        kf->P[i] = P[i];
    }
    return KF_OK;
}

kf_status_t kf_kf_set_covariance_diagonal(kf_kf_t *kf, const kf_real_t *diag)
{
    kf_matrix_t M;
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(diag);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, kf->P, kf->n, kf->n);
    (void)kf_matrix_zero(&M);
    {
        uint16_t i;
        for (i = 0u; i < kf->n; i++) {
            kf_matrix_set(&M, i, i, diag[i]);
        }
    }
    return KF_OK;
}

kf_status_t kf_kf_set_covariance_scalar(kf_kf_t *kf, kf_real_t p)
{
    kf_matrix_t M;
    KF_NULL_CHECK(kf);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, kf->P, kf->n, kf->n);
    (void)kf_matrix_identity(&M);
    (void)kf_matrix_scale(&M, p);
    return KF_OK;
}

const kf_real_t *kf_kf_get_covariance(const kf_kf_t *kf)
{
    if (kf == NULL) {
        return NULL;
    }
    return kf->P;
}

/* --------------------------------------------------------------------------
 * Noise configuration
 * ------------------------------------------------------------------------ */

kf_status_t kf_kf_set_process_noise(kf_kf_t *kf, const kf_real_t *Q)
{
    uint16_t i;
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(Q);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(Q, kf->n * kf->n);
    for (i = 0u; i < (uint16_t)(kf->n * kf->n); i++) {
        kf->Q[i] = Q[i];
    }
    return KF_OK;
}

kf_status_t kf_kf_set_process_noise_diagonal(kf_kf_t *kf, const kf_real_t *diag)
{
    kf_matrix_t M;
    uint16_t i;
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(diag);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, kf->Q, kf->n, kf->n);
    (void)kf_matrix_zero(&M);
    for (i = 0u; i < kf->n; i++) {
        kf_matrix_set(&M, i, i, diag[i]);
    }
    return KF_OK;
}

kf_status_t kf_kf_set_process_noise_scalar(kf_kf_t *kf, kf_real_t q)
{
    kf_matrix_t M;
    KF_NULL_CHECK(kf);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, kf->Q, kf->n, kf->n);
    (void)kf_matrix_identity(&M);
    (void)kf_matrix_scale(&M, q);
    return KF_OK;
}

kf_status_t kf_kf_set_measurement_noise(kf_kf_t *kf, const kf_real_t *R)
{
    uint16_t i;
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(R);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(R, kf->m * kf->m);
    for (i = 0u; i < (uint16_t)(kf->m * kf->m); i++) {
        kf->R[i] = R[i];
    }
    return KF_OK;
}

kf_status_t kf_kf_set_measurement_noise_diagonal(kf_kf_t *kf, const kf_real_t *diag)
{
    kf_matrix_t M;
    uint16_t i;
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(diag);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, kf->R, kf->m, kf->m);
    (void)kf_matrix_zero(&M);
    for (i = 0u; i < kf->m; i++) {
        kf_matrix_set(&M, i, i, diag[i]);
    }
    return KF_OK;
}

kf_status_t kf_kf_set_measurement_noise_scalar(kf_kf_t *kf, kf_real_t r)
{
    kf_matrix_t M;
    KF_NULL_CHECK(kf);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, kf->R, kf->m, kf->m);
    (void)kf_matrix_identity(&M);
    (void)kf_matrix_scale(&M, r);
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Model configuration
 * ------------------------------------------------------------------------ */

kf_status_t kf_kf_set_transition_matrix(kf_kf_t *kf, const kf_real_t *F)
{
    uint16_t i;
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(F);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(F, kf->n * kf->n);
    for (i = 0u; i < (uint16_t)(kf->n * kf->n); i++) {
        kf->F[i] = F[i];
    }
    return KF_OK;
}

kf_status_t kf_kf_set_measurement_matrix(kf_kf_t *kf, const kf_real_t *H)
{
    uint16_t i;
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(H);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(H, kf->m * kf->n);
    for (i = 0u; i < (uint16_t)(kf->m * kf->n); i++) {
        kf->H[i] = H[i];
    }
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Predict
 * ------------------------------------------------------------------------ */

kf_status_t kf_kf_predict(kf_kf_t *kf, const kf_real_t *u, kf_real_t dt)
{
    kf_matrix_t Fm, Pm, Qm, tmp;
    kf_real_t *s;
    uint16_t n;

    (void)dt; /* encoded in F/Q for the linear case */

    KF_NULL_CHECK(kf);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    n = kf->n;
    s = kf->scratch;

    kf_view(&Fm, kf->F, n, n);
    kf_view(&Pm, kf->P, n, n);
    kf_view(&Qm, kf->Q, n, n);

    /* x = F x + u */
    (void)kf_mat_vec_mul(&s[OFF_VEC_N], &Fm, kf->x);
    if (u != NULL) {
        kf_vec_add(&s[OFF_VEC_N], &s[OFF_VEC_N], u, n);
    }
    kf_vec_copy(kf->x, &s[OFF_VEC_N], n);

    /* P = F P F^T + Q */
    kf_view(&tmp, &s[OFF_NN1], n, n);
    (void)kf_matrix_mul(&tmp, &Fm, &Pm);            /* tmp = F P           */
    (void)kf_matrix_mul_transpose_b(&Pm, &tmp, &Fm);/* P = (F P) F^T       */
    (void)kf_matrix_add(&Pm, &Pm, &Qm);             /* P += Q              */
    (void)kf_matrix_symmetrize(&Pm);

    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Update
 * ------------------------------------------------------------------------ */

kf_status_t kf_kf_update(kf_kf_t *kf, const kf_real_t *z)
{
    kf_matrix_t Hm, Pm, Rm, Bm, Sm, Km, Tm, T2, KR, RHS;
    kf_real_t *s;
    kf_status_t st;
    uint16_t n;
    uint16_t m;
    uint16_t r;
    uint16_t c;

    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(z);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    n = kf->n;
    m = kf->m;
    s = kf->scratch;

    kf_view(&Hm, kf->H, m, n);
    kf_view(&Pm, kf->P, n, n);
    kf_view(&Rm, kf->R, m, m);

    /* Innovation y = z - H x */
    (void)kf_mat_vec_mul(&s[OFF_VEC_M], &Hm, kf->x);   /* Hx            */
    kf_vec_sub(&s[OFF_VEC_M], z, &s[OFF_VEC_M], m); /* y = z - Hx */
#if KF_ENABLE_ADVANCED_API
    kf_vec_copy(kf->y, &s[OFF_VEC_M], m);
#endif

    /* B = P H^T  (n x m) */
    kf_view(&Bm, &s[OFF_NM1], n, m);
    (void)kf_matrix_mul_transpose_b(&Bm, &Pm, &Hm);

    /* S = H B + R  (m x m) */
    kf_view(&Sm, &s[OFF_MM], m, m);
    (void)kf_matrix_mul(&Sm, &Hm, &Bm);
    (void)kf_matrix_add(&Sm, &Sm, &Rm);

    /* Solve for the gain without an explicit inverse:
       S = L L^T, then S K^T = B^T  ->  K = (B^T S^-1)^T. */
    st = kf_matrix_cholesky(&Sm);
    if (st != KF_OK) {
        return st;
    }
    kf_view(&RHS, &s[OFF_MN], m, n);
    (void)kf_matrix_transpose(&RHS, &Bm);       /* RHS = B^T          */
    st = kf_matrix_cholesky_solve(&Sm, &RHS); /* RHS = K^T      */
    if (st != KF_OK) {
        return st;
    }
    kf_view(&Km, &s[OFF_NM1], n, m);
    (void)kf_matrix_transpose(&Km, &RHS);       /* K = (K^T)^T        */

#if KF_ENABLE_ADVANCED_API
    {
        uint16_t i;
        for (i = 0u; i < (uint16_t)(n * m); i++) {
            kf->K[i] = s[OFF_NM1 + i];
        }
    }
#endif

    /* x = x + K y */
    (void)kf_mat_vec_mul(&s[OFF_VEC_N], &Km, &s[OFF_VEC_M]);
    kf_vec_add(kf->x, kf->x, &s[OFF_VEC_N], n);

    /* Covariance update */
    kf_view(&Tm, &s[OFF_NN1], n, n);
    (void)kf_matrix_mul(&Tm, &Km, &Hm);          /* Tm = K H           */
    for (r = 0u; r < n; r++) {
        for (c = 0u; c < n; c++) {
            kf_real_t ident = (r == c) ? (kf_real_t)1 : (kf_real_t)0;
            kf_matrix_set(&Tm, r, c, ident - kf_matrix_get(&Tm, r, c));
        }
    }                                        /* Tm = I - K H      */

#if KF_USE_JOSEPH_FORM
    /* P = (I-KH) P (I-KH)^T + K R K^T   (Joseph form) */
    kf_view(&T2, &s[OFF_NN2], n, n);
    (void)kf_matrix_mul(&T2, &Tm, &Pm);            /* T2 = (I-KH) P    */
    (void)kf_matrix_mul_transpose_b(&Pm, &T2, &Tm);/* P = T2 (I-KH)^T  */
    kf_view(&KR, &s[OFF_NM2], n, m);
    (void)kf_matrix_mul(&KR, &Km, &Rm);            /* KR = K R         */
    (void)kf_matrix_mul_transpose_b(&T2, &KR, &Km);/* T2 = K R K^T     */
    (void)kf_matrix_add(&Pm, &Pm, &T2);            /* P += K R K^T     */
#else
    /* P = (I-KH) P   (simple form) */
    kf_view(&T2, &s[OFF_NN2], n, n);
    (void)kf_matrix_mul(&T2, &Tm, &Pm);
    (void)kf_matrix_copy(&Pm, &T2);
#endif
    (void)kf_matrix_symmetrize(&Pm);

    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Advanced / diagnostics
 * ------------------------------------------------------------------------ */

#if KF_ENABLE_ADVANCED_API
kf_status_t kf_kf_get_innovation(const kf_kf_t *kf, kf_real_t *y_out)
{
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(y_out);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_vec_copy(y_out, kf->y, kf->m);
    return KF_OK;
}

kf_status_t kf_kf_get_gain(const kf_kf_t *kf, kf_real_t *K_out)
{
    uint16_t i;
    KF_NULL_CHECK(kf);
    KF_NULL_CHECK(K_out);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    for (i = 0u; i < (uint16_t)(kf->n * kf->m); i++) {
        K_out[i] = kf->K[i];
    }
    return KF_OK;
}
#endif /* KF_ENABLE_ADVANCED_API */

#if KF_ENABLE_DIAGNOSTICS
kf_status_t kf_kf_check(const kf_kf_t *kf)
{
    uint16_t i;

    KF_NULL_CHECK(kf);

    if (kf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    if (!kf_vec_is_finite(kf->x, kf->n)) {
        return KF_ERROR_NON_FINITE;
    }
    for (i = 0u; i < (uint16_t)(kf->n * kf->n); i++) {
        if (!kf_isfinite(kf->P[i])) {
            return KF_ERROR_NON_FINITE;
        }
    }
    if (!kf_mat_is_symmetric(kf->P, kf->n)) {
        return KF_ERROR_NUMERICAL;
    }
    return KF_OK;
}
#endif /* KF_ENABLE_DIAGNOSTICS */

#endif /* KF_ENABLE_KF */
