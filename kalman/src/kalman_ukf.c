/**
 * @file    kalman_ukf.c
 * @brief   Unscented Kalman filter implementation.
 */

#include "kalman_ukf.h"

#if KF_ENABLE_UKF
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

#define OFF_PCOPY (0u)
#define OFF_NN1   (OFF_PCOPY + KF_NN)
#define OFF_NM1   (OFF_NN1 + KF_NN)
#define OFF_NM2   (OFF_NM1 + KF_NM)
#define OFF_MN1   (OFF_NM2 + KF_NM)
#define OFF_MM1   (OFF_MN1 + KF_NM)
#define OFF_VN1   (OFF_MM1 + KF_MM)
#define OFF_VN2   (OFF_VN1 + KF_MAX_STATE_DIM)
#define OFF_VM1   (OFF_VN2 + KF_MAX_STATE_DIM)
#define OFF_VM2   (OFF_VM1 + KF_MAX_MEASUREMENT_DIM)
#define OFF_END   (OFF_VM2 + KF_MAX_MEASUREMENT_DIM)

#if OFF_END != KF_UKF_SCRATCH_FLOATS
#error "kalman_ukf.c: scratch layout does not match KF_UKF_SCRATCH_FLOATS"
#endif

static void kf_view(kf_matrix_t *m, kf_real_t *data, uint16_t rows, uint16_t cols)
{
    m->rows = rows;
    m->cols = cols;
    m->data = data;
}

/* M += w * outer(a, b)  (M is rows x cols, row-major). */
static void kf_outer_add(kf_real_t *M, uint16_t rows, uint16_t cols,
                         kf_real_t w, const kf_real_t *a, const kf_real_t *b)
{
    uint16_t r;
    for (r = 0u; r < rows; r++) {
        kf_real_t ar = w * a[r];
        uint16_t c;
        for (c = 0u; c < cols; c++) {
            M[(size_t)r * cols + c] += ar * b[c];
        }
    }
}

/* --------------------------------------------------------------------------
 * Sigma-point generation
 * ------------------------------------------------------------------------ */

/* Builds the 2n+1 sigma points from (x, P) using the Cholesky factor of P.
 * Uses the OFF_PCOPY scratch region for the factorisation. */
static kf_status_t ukf_sigma_points(kf_ukf_t *ukf)
{
    kf_matrix_t L;
    kf_real_t c;
    uint16_t n = ukf->n;
    uint16_t i;
    uint16_t j;
    kf_real_t scale = kf_sqrt((kf_real_t)n + ukf->lambda);

    kf_view(&L, &ukf->scratch[OFF_PCOPY], n, n);
    for (i = 0u; i < (uint16_t)(n * n); i++) {
        L.data[i] = ukf->P[i];
    }
    {
        kf_status_t st = kf_matrix_cholesky(&L);
        if (st != KF_OK) {
            return st;
        }
    }

    /* sigma_0 = x */
    for (j = 0u; j < n; j++) {
        ukf->sigma[j] = ukf->x[j];
    }
    /* sigma_i = x + c * L[:, i-1]; sigma_{i+n} = x - c * L[:, i-1] */
    c = scale;
    for (i = 1u; i <= n; i++) {
        uint16_t col = (uint16_t)(i - 1u);
        size_t row_pos = (size_t)i * (size_t)n;
        size_t row_neg = ((size_t)i + (size_t)n) * (size_t)n;
        for (j = 0u; j < n; j++) {
            kf_real_t lv = (j >= col) ? L.data[(size_t)j * (size_t)n + (size_t)col] : (kf_real_t)0;
            ukf->sigma[row_pos + (size_t)j] = ukf->x[j] + c * lv;
            ukf->sigma[row_neg + (size_t)j] = ukf->x[j] - c * lv;
        }
    }
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Initialisation / reset
 * ------------------------------------------------------------------------ */

static void ukf_recompute_weights(kf_ukf_t *ukf)
{
    uint16_t i;
    uint16_t n = ukf->n;
    kf_real_t nlambda = (kf_real_t)n + ukf->lambda;
    uint16_t nsig = (uint16_t)(2u * n + 1u);

    ukf->Wm[0] = ukf->lambda / nlambda;
    ukf->Wc[0] = ukf->lambda / nlambda +
                 ((kf_real_t)1 - ukf->alpha * ukf->alpha + ukf->beta);
    for (i = 1u; i < nsig; i++) {
        ukf->Wm[i] = (kf_real_t)0.5 / nlambda;
        ukf->Wc[i] = (kf_real_t)0.5 / nlambda;
    }
}

kf_status_t kf_ukf_init(kf_ukf_t *ukf, uint16_t n, uint16_t m)
{
    KF_NULL_CHECK(ukf);

    if (n == 0u || n > KF_MAX_STATE_DIM) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    if (m == 0u || m > KF_MAX_MEASUREMENT_DIM) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    ukf->n = n;
    ukf->m = m;
    ukf->f = NULL;
    ukf->h = NULL;
    ukf->ctx_f = NULL;
    ukf->ctx_h = NULL;
    ukf->initialized = 1u;

    ukf->alpha = (kf_real_t)1.0;
    ukf->beta  = (kf_real_t)2.0;
    ukf->kappa = (kf_real_t)0.0;
    ukf->lambda = ukf->alpha * ukf->alpha * ((kf_real_t)n + ukf->kappa) - (kf_real_t)n;
    ukf_recompute_weights(ukf);

    return kf_ukf_reset(ukf);
}

kf_status_t kf_ukf_reset(kf_ukf_t *ukf)
{
    kf_matrix_t M;
    uint16_t n;
    uint16_t m;
    uint16_t i;

    KF_NULL_CHECK(ukf);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    n = ukf->n;
    m = ukf->m;

    kf_vec_zero(ukf->x, n);
    kf_view(&M, ukf->P, n, n);
    (void)kf_matrix_identity(&M);
    kf_view(&M, ukf->Q, n, n);
    (void)kf_matrix_zero(&M);
    kf_view(&M, ukf->R, m, m);
    (void)kf_matrix_zero(&M);

    for (i = 0u; i < (uint16_t)((2u * n + 1u) * n); i++) {
        ukf->sigma[i] = (kf_real_t)0;
        ukf->sigma_prop[i] = (kf_real_t)0;
    }
    for (i = 0u; i < (uint16_t)((2u * n + 1u) * m); i++) {
        ukf->sigma_z[i] = (kf_real_t)0;
    }

#if KF_ENABLE_ADVANCED_API
    kf_vec_zero(ukf->y, m);
    for (i = 0u; i < (uint16_t)(n * m); i++) {
        ukf->K[i] = (kf_real_t)0;
    }
#endif

    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Model / parameters
 * ------------------------------------------------------------------------ */

kf_status_t kf_ukf_set_models(kf_ukf_t *ukf, kf_ukf_f_fn f, kf_ukf_h_fn h)
{
    KF_NULL_CHECK(ukf);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    ukf->f = f;
    ukf->h = h;
    return KF_OK;
}

void kf_ukf_set_context_f(kf_ukf_t *ukf, void *ctx)
{
    if (ukf != NULL) {
        ukf->ctx_f = ctx;
    }
}

void kf_ukf_set_context_h(kf_ukf_t *ukf, void *ctx)
{
    if (ukf != NULL) {
        ukf->ctx_h = ctx;
    }
}

kf_status_t kf_ukf_set_parameters(kf_ukf_t *ukf, kf_real_t alpha,
                                  kf_real_t beta, kf_real_t kappa)
{
    KF_NULL_CHECK(ukf);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    if (alpha <= (kf_real_t)0 || beta < (kf_real_t)0) {
        return KF_ERROR_INVALID_PARAMETER;
    }
    ukf->alpha = alpha;
    ukf->beta  = beta;
    ukf->kappa = kappa;
    ukf->lambda = alpha * alpha * ((kf_real_t)ukf->n + kappa) - (kf_real_t)ukf->n;
    ukf_recompute_weights(ukf);
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * State / covariance / noise setters
 * ------------------------------------------------------------------------ */

kf_status_t kf_ukf_set_state(kf_ukf_t *ukf, const kf_real_t *x)
{
    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(x);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(x, ukf->n);
    kf_vec_copy(ukf->x, x, ukf->n);
    return KF_OK;
}

const kf_real_t *kf_ukf_get_state(const kf_ukf_t *ukf)
{
    return (ukf == NULL) ? NULL : ukf->x;
}

kf_status_t kf_ukf_set_covariance(kf_ukf_t *ukf, const kf_real_t *P)
{
    uint16_t i;
    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(P);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(P, ukf->n * ukf->n);
    for (i = 0u; i < (uint16_t)(ukf->n * ukf->n); i++) {
        ukf->P[i] = P[i];
    }
    return KF_OK;
}

kf_status_t kf_ukf_set_covariance_diagonal(kf_ukf_t *ukf, const kf_real_t *diag)
{
    kf_matrix_t M;
    uint16_t i;
    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(diag);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ukf->P, ukf->n, ukf->n);
    (void)kf_matrix_zero(&M);
    for (i = 0u; i < ukf->n; i++) {
        kf_matrix_set(&M, i, i, diag[i]);
    }
    return KF_OK;
}

kf_status_t kf_ukf_set_covariance_scalar(kf_ukf_t *ukf, kf_real_t p)
{
    kf_matrix_t M;
    KF_NULL_CHECK(ukf);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ukf->P, ukf->n, ukf->n);
    (void)kf_matrix_identity(&M);
    (void)kf_matrix_scale(&M, p);
    return KF_OK;
}

const kf_real_t *kf_ukf_get_covariance(const kf_ukf_t *ukf)
{
    return (ukf == NULL) ? NULL : ukf->P;
}

kf_status_t kf_ukf_set_process_noise(kf_ukf_t *ukf, const kf_real_t *Q)
{
    uint16_t i;
    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(Q);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(Q, ukf->n * ukf->n);
    for (i = 0u; i < (uint16_t)(ukf->n * ukf->n); i++) {
        ukf->Q[i] = Q[i];
    }
    return KF_OK;
}

kf_status_t kf_ukf_set_process_noise_diagonal(kf_ukf_t *ukf, const kf_real_t *diag)
{
    kf_matrix_t M;
    uint16_t i;
    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(diag);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ukf->Q, ukf->n, ukf->n);
    (void)kf_matrix_zero(&M);
    for (i = 0u; i < ukf->n; i++) {
        kf_matrix_set(&M, i, i, diag[i]);
    }
    return KF_OK;
}

kf_status_t kf_ukf_set_process_noise_scalar(kf_ukf_t *ukf, kf_real_t q)
{
    kf_matrix_t M;
    KF_NULL_CHECK(ukf);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ukf->Q, ukf->n, ukf->n);
    (void)kf_matrix_identity(&M);
    (void)kf_matrix_scale(&M, q);
    return KF_OK;
}

kf_status_t kf_ukf_set_measurement_noise(kf_ukf_t *ukf, const kf_real_t *R)
{
    uint16_t i;
    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(R);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    KF_REQUIRE_FINITE(R, ukf->m * ukf->m);
    for (i = 0u; i < (uint16_t)(ukf->m * ukf->m); i++) {
        ukf->R[i] = R[i];
    }
    return KF_OK;
}

kf_status_t kf_ukf_set_measurement_noise_diagonal(kf_ukf_t *ukf, const kf_real_t *diag)
{
    kf_matrix_t M;
    uint16_t i;
    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(diag);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ukf->R, ukf->m, ukf->m);
    (void)kf_matrix_zero(&M);
    for (i = 0u; i < ukf->m; i++) {
        kf_matrix_set(&M, i, i, diag[i]);
    }
    return KF_OK;
}

kf_status_t kf_ukf_set_measurement_noise_scalar(kf_ukf_t *ukf, kf_real_t r)
{
    kf_matrix_t M;
    KF_NULL_CHECK(ukf);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_view(&M, ukf->R, ukf->m, ukf->m);
    (void)kf_matrix_identity(&M);
    (void)kf_matrix_scale(&M, r);
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Predict
 * ------------------------------------------------------------------------ */

kf_status_t kf_ukf_predict(kf_ukf_t *ukf, const kf_real_t *u, kf_real_t dt)
{
    kf_status_t st;
    kf_real_t *x_mean;
    kf_real_t *diff;
    uint16_t n = ukf->n;
    uint16_t nsig = (uint16_t)(2u * n + 1u);
    uint16_t i;

    KF_NULL_CHECK(ukf);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    if (ukf->f == NULL) {
        return KF_ERROR_INVALID_PARAMETER;
    }

    st = ukf_sigma_points(ukf);
    if (st != KF_OK) {
        return st;
    }

    /* Propagate each sigma point through f. */
    for (i = 0u; i < nsig; i++) {
        ukf->f(&ukf->sigma[(size_t)i * n], u, dt,
               &ukf->sigma_prop[(size_t)i * n], ukf->ctx_f);
    }

    /* Weighted mean. */
    x_mean = &ukf->scratch[OFF_VN1];
    kf_vec_zero(x_mean, n);
    for (i = 0u; i < nsig; i++) {
        kf_vec_axpy(x_mean, ukf->Wm[i], &ukf->sigma_prop[(size_t)i * n], n);
    }

    /* Covariance: P = Q + sum Wc[i] (X_i - x)(X_i - x)^T. */
    {
        uint16_t j;
        for (j = 0u; j < (uint16_t)(n * n); j++) {
            ukf->P[j] = ukf->Q[j];
        }
    }
    diff = &ukf->scratch[OFF_VN2];
    for (i = 0u; i < nsig; i++) {
        kf_vec_sub(diff, &ukf->sigma_prop[(size_t)i * n], x_mean, n);
        kf_outer_add(ukf->P, n, n, ukf->Wc[i], diff, diff);
    }
    {
        kf_matrix_t Pm;
        kf_view(&Pm, ukf->P, n, n);
        (void)kf_matrix_symmetrize(&Pm);
    }

    kf_vec_copy(ukf->x, x_mean, n);
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Update
 * ------------------------------------------------------------------------ */

kf_status_t kf_ukf_update(kf_ukf_t *ukf, const kf_real_t *z)
{
    kf_status_t st;
    kf_matrix_t Sm, Km, Pxzm, RHSm, Tm;
    kf_real_t *z_mean;
    kf_real_t *diff_z;
    kf_real_t *diff_x;
    uint16_t n = ukf->n;
    uint16_t m = ukf->m;
    uint16_t nsig = (uint16_t)(2u * n + 1u);
    uint16_t i;

    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(z);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    if (ukf->h == NULL) {
        return KF_ERROR_INVALID_PARAMETER;
    }

    st = ukf_sigma_points(ukf);
    if (st != KF_OK) {
        return st;
    }

    /* Propagate sigma points through h. */
    for (i = 0u; i < nsig; i++) {
        ukf->h(&ukf->sigma[(size_t)i * n], &ukf->sigma_z[(size_t)i * m], ukf->ctx_h);
    }

    /* Predicted measurement mean. */
    z_mean = &ukf->scratch[OFF_VM1];
    kf_vec_zero(z_mean, m);
    for (i = 0u; i < nsig; i++) {
        kf_vec_axpy(z_mean, ukf->Wm[i], &ukf->sigma_z[(size_t)i * m], m);
    }

    /* Innovation covariance S = R + sum Wc[i] (Z_i - z)(Z_i - z)^T. */
    kf_view(&Sm, &ukf->scratch[OFF_MM1], m, m);
    {
        uint16_t j;
        for (j = 0u; j < (uint16_t)(m * m); j++) {
            Sm.data[j] = ukf->R[j];
        }
    }
    diff_z = &ukf->scratch[OFF_VM2];
    for (i = 0u; i < nsig; i++) {
        kf_vec_sub(diff_z, &ukf->sigma_z[(size_t)i * m], z_mean, m);
        kf_outer_add(Sm.data, m, m, ukf->Wc[i], diff_z, diff_z);
    }

    /* Cross-covariance Pxz = sum Wc[i] (X_i - x)(Z_i - z)^T  (n x m). */
    kf_view(&Pxzm, &ukf->scratch[OFF_NM1], n, m);
    (void)kf_matrix_zero(&Pxzm);
    diff_x = &ukf->scratch[OFF_VN1];
    for (i = 0u; i < nsig; i++) {
        kf_vec_sub(diff_x, &ukf->sigma[(size_t)i * n], ukf->x, n);
        kf_vec_sub(diff_z, &ukf->sigma_z[(size_t)i * m], z_mean, m);
        kf_outer_add(Pxzm.data, n, m, ukf->Wc[i], diff_x, diff_z);
    }

    /* Gain: K = Pxz S^-1  ->  solve S K^T = Pxz^T. */
    st = kf_matrix_cholesky(&Sm);
    if (st != KF_OK) {
        return st;
    }
    kf_view(&RHSm, &ukf->scratch[OFF_MN1], m, n);
    (void)kf_matrix_transpose(&RHSm, &Pxzm);
    st = kf_matrix_cholesky_solve(&Sm, &RHSm);
    if (st != KF_OK) {
        return st;
    }
    kf_view(&Km, &ukf->scratch[OFF_NM2], n, m);
    (void)kf_matrix_transpose(&Km, &RHSm);

#if KF_ENABLE_ADVANCED_API
    {
        uint16_t j;
        for (j = 0u; j < (uint16_t)(n * m); j++) {
            ukf->K[j] = Km.data[j];
        }
        kf_vec_sub(ukf->y, z, z_mean, m);
    }
#endif

    /* Innovation y = z - z_mean, then x = x + K y. */
    kf_vec_sub(&ukf->scratch[OFF_VM2], z, z_mean, m);
    (void)kf_mat_vec_mul(&ukf->scratch[OFF_VN2], &Km, &ukf->scratch[OFF_VM2]);
    kf_vec_add(ukf->x, ukf->x, &ukf->scratch[OFF_VN2], n);

    /* Covariance: P = P - K S K^T = P - K Pxz^T. */
    kf_view(&Tm, &ukf->scratch[OFF_NN1], n, n);
    (void)kf_matrix_mul_transpose_b(&Tm, &Km, &Pxzm);   /* Tm = K Pxz^T     */
    {
        kf_matrix_t Pm;
        kf_view(&Pm, ukf->P, n, n);
        (void)kf_matrix_sub(&Pm, &Pm, &Tm);
        (void)kf_matrix_symmetrize(&Pm);
    }

    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Advanced / diagnostics
 * ------------------------------------------------------------------------ */

#if KF_ENABLE_ADVANCED_API
kf_status_t kf_ukf_get_innovation(const kf_ukf_t *ukf, kf_real_t *y_out)
{
    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(y_out);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    kf_vec_copy(y_out, ukf->y, ukf->m);
    return KF_OK;
}

kf_status_t kf_ukf_get_gain(const kf_ukf_t *ukf, kf_real_t *K_out)
{
    uint16_t i;
    KF_NULL_CHECK(ukf);
    KF_NULL_CHECK(K_out);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    for (i = 0u; i < (uint16_t)(ukf->n * ukf->m); i++) {
        K_out[i] = ukf->K[i];
    }
    return KF_OK;
}
#endif

#if KF_ENABLE_DIAGNOSTICS
kf_status_t kf_ukf_check(const kf_ukf_t *ukf)
{
    uint16_t i;

    KF_NULL_CHECK(ukf);

    if (ukf->initialized == 0u) {
        return KF_ERROR_NOT_INITIALIZED;
    }
    if (!kf_vec_is_finite(ukf->x, ukf->n)) {
        return KF_ERROR_NON_FINITE;
    }
    for (i = 0u; i < (uint16_t)(ukf->n * ukf->n); i++) {
        if (!kf_isfinite(ukf->P[i])) {
            return KF_ERROR_NON_FINITE;
        }
    }
    if (!kf_mat_is_symmetric(ukf->P, ukf->n)) {
        return KF_ERROR_NUMERICAL;
    }
    return KF_OK;
}
#endif

#endif /* KF_ENABLE_UKF */
