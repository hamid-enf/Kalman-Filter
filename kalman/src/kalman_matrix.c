/**
 * @file    kalman_matrix.c
 * @brief   Implementation of the minimal matrix engine.
 *
 * All routines are hand-written, allocation-free, and use explicit integer
 * types. Loops are written to be friendly to the compiler's auto-vectorisation
 * on Cortex-M4/M7/M33 (contiguous row-major access where possible).
 */

#include "kalman_matrix.h"

#include <stddef.h>

/* --------------------------------------------------------------------------
 * Construction / access
 * ------------------------------------------------------------------------ */

kf_status_t kf_matrix_init(kf_matrix_t *m, uint16_t rows, uint16_t cols,
                           kf_real_t *data)
{
    KF_NULL_CHECK(m);
    KF_NULL_CHECK(data);

    if (rows == 0u || cols == 0u) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    m->rows = rows;
    m->cols = cols;
    m->data = data;
    return KF_OK;
}

kf_real_t kf_matrix_get(const kf_matrix_t *m, uint16_t r, uint16_t c)
{
    return m->data[(size_t)r * m->cols + c];
}

void kf_matrix_set(kf_matrix_t *m, uint16_t r, uint16_t c, kf_real_t v)
{
    m->data[(size_t)r * m->cols + c] = v;
}

/* --------------------------------------------------------------------------
 * Basic element-wise operations
 * ------------------------------------------------------------------------ */

kf_status_t kf_matrix_zero(kf_matrix_t *m)
{
    uint16_t r;
    KF_NULL_CHECK(m);
    KF_NULL_CHECK(m->data);

    for (r = 0u; r < m->rows; r++) {
        uint16_t c;
        for (c = 0u; c < m->cols; c++) {
            kf_matrix_set(m, r, c, (kf_real_t)0);
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_identity(kf_matrix_t *m)
{
    uint16_t r, c;
    KF_NULL_CHECK(m);
    KF_NULL_CHECK(m->data);

    if (m->rows != m->cols) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    for (r = 0u; r < m->rows; r++) {
        for (c = 0u; c < m->cols; c++) {
            kf_matrix_set(m, r, c, (r == c) ? (kf_real_t)1 : (kf_real_t)0);
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_copy(kf_matrix_t *dst, const kf_matrix_t *src)
{
    uint16_t r;
    KF_NULL_CHECK(dst);
    KF_NULL_CHECK(src);
    KF_NULL_CHECK(dst->data);
    KF_NULL_CHECK(src->data);

    if (dst->rows != src->rows || dst->cols != src->cols) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    for (r = 0u; r < src->rows; r++) {
        uint16_t c;
        for (c = 0u; c < src->cols; c++) {
            kf_matrix_set(dst, r, c, kf_matrix_get(src, r, c));
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_scale(kf_matrix_t *m, kf_real_t s)
{
    uint16_t r;
    KF_NULL_CHECK(m);
    KF_NULL_CHECK(m->data);

    for (r = 0u; r < m->rows; r++) {
        uint16_t c;
        for (c = 0u; c < m->cols; c++) {
            kf_matrix_set(m, r, c, kf_matrix_get(m, r, c) * s);
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_add(kf_matrix_t *out, const kf_matrix_t *a,
                          const kf_matrix_t *b)
{
    uint16_t r;
    KF_NULL_CHECK(out);
    KF_NULL_CHECK(a);
    KF_NULL_CHECK(b);
    KF_NULL_CHECK(out->data);
    KF_NULL_CHECK(a->data);
    KF_NULL_CHECK(b->data);

    if (out->rows != a->rows || out->cols != a->cols ||
        a->rows != b->rows || a->cols != b->cols) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    for (r = 0u; r < a->rows; r++) {
        uint16_t c;
        for (c = 0u; c < a->cols; c++) {
            kf_matrix_set(out, r, c, kf_matrix_get(a, r, c) + kf_matrix_get(b, r, c));
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_add_scaled(kf_matrix_t *out, const kf_matrix_t *a,
                                 const kf_matrix_t *b, kf_real_t s)
{
    uint16_t r;
    KF_NULL_CHECK(out);
    KF_NULL_CHECK(a);
    KF_NULL_CHECK(b);
    KF_NULL_CHECK(out->data);
    KF_NULL_CHECK(a->data);
    KF_NULL_CHECK(b->data);

    if (out->rows != a->rows || out->cols != a->cols ||
        a->rows != b->rows || a->cols != b->cols) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    for (r = 0u; r < a->rows; r++) {
        uint16_t c;
        for (c = 0u; c < a->cols; c++) {
            kf_matrix_set(out, r, c,
                          kf_matrix_get(a, r, c) + s * kf_matrix_get(b, r, c));
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_sub(kf_matrix_t *out, const kf_matrix_t *a,
                          const kf_matrix_t *b)
{
    uint16_t r;
    KF_NULL_CHECK(out);
    KF_NULL_CHECK(a);
    KF_NULL_CHECK(b);
    KF_NULL_CHECK(out->data);
    KF_NULL_CHECK(a->data);
    KF_NULL_CHECK(b->data);

    if (out->rows != a->rows || out->cols != a->cols ||
        a->rows != b->rows || a->cols != b->cols) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    for (r = 0u; r < a->rows; r++) {
        uint16_t c;
        for (c = 0u; c < a->cols; c++) {
            kf_matrix_set(out, r, c, kf_matrix_get(a, r, c) - kf_matrix_get(b, r, c));
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_symmetrize(kf_matrix_t *m)
{
    uint16_t r;
    KF_NULL_CHECK(m);
    KF_NULL_CHECK(m->data);

    if (m->rows != m->cols) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    for (r = 0u; r < m->rows; r++) {
        uint16_t c;
        for (c = r + 1u; c < m->cols; c++) {
            kf_real_t avg = (kf_matrix_get(m, r, c) + kf_matrix_get(m, c, r)) * (kf_real_t)0.5;
            kf_matrix_set(m, r, c, avg);
            kf_matrix_set(m, c, r, avg);
        }
    }
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Multiplication
 * ------------------------------------------------------------------------ */

kf_status_t kf_matrix_mul(kf_matrix_t *out, const kf_matrix_t *a,
                          const kf_matrix_t *b)
{
    uint16_t i;
    uint16_t p, q, r;

    KF_NULL_CHECK(out);
    KF_NULL_CHECK(a);
    KF_NULL_CHECK(b);
    KF_NULL_CHECK(out->data);
    KF_NULL_CHECK(a->data);
    KF_NULL_CHECK(b->data);

    if (a->cols != b->rows) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    p = a->rows;
    q = a->cols;
    r = b->cols;
    if (out->rows != p || out->cols != r) {
        return KF_ERROR_INVALID_DIMENSION;
    }

    for (i = 0u; i < p; i++) {
        uint16_t j;
        for (j = 0u; j < r; j++) {
            kf_real_t sum = (kf_real_t)0;
            uint16_t k;
            for (k = 0u; k < q; k++) {
                sum += kf_matrix_get(a, i, k) * kf_matrix_get(b, k, j);
            }
            kf_matrix_set(out, i, j, sum);
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_mul_transpose_b(kf_matrix_t *out, const kf_matrix_t *a,
                                      const kf_matrix_t *b)
{
    uint16_t i;
    uint16_t p, q, r;

    KF_NULL_CHECK(out);
    KF_NULL_CHECK(a);
    KF_NULL_CHECK(b);
    KF_NULL_CHECK(out->data);
    KF_NULL_CHECK(a->data);
    KF_NULL_CHECK(b->data);

    if (a->cols != b->cols) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    p = a->rows;     /* out rows              */
    q = a->cols;     /* shared inner dim      */
    r = b->rows;     /* b^T has r columns     */
    if (out->rows != p || out->cols != r) {
        return KF_ERROR_INVALID_DIMENSION;
    }

    for (i = 0u; i < p; i++) {
        uint16_t j;
        for (j = 0u; j < r; j++) {
            kf_real_t sum = (kf_real_t)0;
            uint16_t k;
            for (k = 0u; k < q; k++) {
                sum += kf_matrix_get(a, i, k) * kf_matrix_get(b, j, k);
            }
            kf_matrix_set(out, i, j, sum);
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_mul_transpose_a(kf_matrix_t *out, const kf_matrix_t *a,
                                      const kf_matrix_t *b)
{
    uint16_t i;
    uint16_t p, q, r;

    KF_NULL_CHECK(out);
    KF_NULL_CHECK(a);
    KF_NULL_CHECK(b);
    KF_NULL_CHECK(out->data);
    KF_NULL_CHECK(a->data);
    KF_NULL_CHECK(b->data);

    if (a->rows != b->rows) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    q = a->rows;     /* shared inner dim      */
    p = a->cols;     /* a^T has p rows        */
    r = b->cols;
    if (out->rows != p || out->cols != r) {
        return KF_ERROR_INVALID_DIMENSION;
    }

    for (i = 0u; i < p; i++) {
        uint16_t j;
        for (j = 0u; j < r; j++) {
            kf_real_t sum = (kf_real_t)0;
            uint16_t k;
            for (k = 0u; k < q; k++) {
                sum += kf_matrix_get(a, k, i) * kf_matrix_get(b, k, j);
            }
            kf_matrix_set(out, i, j, sum);
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_transpose(kf_matrix_t *out, const kf_matrix_t *a)
{
    KF_NULL_CHECK(out);
    KF_NULL_CHECK(a);
    KF_NULL_CHECK(out->data);
    KF_NULL_CHECK(a->data);

    if (out->rows != a->cols || out->cols != a->rows) {
        return KF_ERROR_INVALID_DIMENSION;
    }
    if (out == a) {
        /* In-place transpose: only valid for square matrices. */
        uint16_t r;
        if (out->rows != out->cols) {
            return KF_ERROR_INVALID_DIMENSION;
        }
        for (r = 0u; r < out->rows; r++) {
            uint16_t c;
            for (c = r + 1u; c < out->cols; c++) {
                kf_real_t t = kf_matrix_get(out, r, c);
                kf_matrix_set(out, r, c, kf_matrix_get(out, c, r));
                kf_matrix_set(out, c, r, t);
            }
        }
    } else {
        uint16_t r;
        for (r = 0u; r < a->rows; r++) {
            uint16_t c;
            for (c = 0u; c < a->cols; c++) {
                kf_matrix_set(out, c, r, kf_matrix_get(a, r, c));
            }
        }
    }
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Cholesky decomposition & solve
 * ------------------------------------------------------------------------ */

kf_status_t kf_matrix_cholesky(kf_matrix_t *m)
{
    uint16_t i;

    KF_NULL_CHECK(m);
    KF_NULL_CHECK(m->data);

    if (m->rows != m->cols) {
        return KF_ERROR_INVALID_DIMENSION;
    }

    for (i = 0u; i < m->rows; i++) {
        uint16_t j;
        for (j = 0u; j <= i; j++) {
            kf_real_t sum = kf_matrix_get(m, i, j);
            uint16_t k;
            for (k = 0u; k < j; k++) {
                sum -= kf_matrix_get(m, i, k) * kf_matrix_get(m, j, k);
            }
            if (i > j) {
                kf_real_t diag = kf_matrix_get(m, j, j);
                if (diag <= KF_MIN_PIVOT || !kf_isfinite(diag)) {
                    return KF_ERROR_NOT_POSITIVE_DEFINITE;
                }
                kf_matrix_set(m, i, j, sum / diag);
            } else {
                if (sum <= KF_MIN_PIVOT || !kf_isfinite(sum)) {
                    return KF_ERROR_NOT_POSITIVE_DEFINITE;
                }
                kf_matrix_set(m, i, j, kf_sqrt(sum));
            }
        }
    }
    return KF_OK;
}

kf_status_t kf_matrix_cholesky_solve(const kf_matrix_t *L, kf_matrix_t *B)
{
    uint16_t col;

    KF_NULL_CHECK(L);
    KF_NULL_CHECK(B);
    KF_NULL_CHECK(L->data);
    KF_NULL_CHECK(B->data);

    if (L->rows != L->cols || B->rows != L->cols) {
        return KF_ERROR_INVALID_DIMENSION;
    }

    for (col = 0u; col < B->cols; col++) {
        uint16_t i;

        /* Forward substitution: L y = b  (in place). */
        for (i = 0u; i < L->rows; i++) {
            kf_real_t sum = kf_matrix_get(B, i, col);
            uint16_t k;
            for (k = 0u; k < i; k++) {
                sum -= kf_matrix_get(L, i, k) * kf_matrix_get(B, k, col);
            }
            kf_matrix_set(B, i, col, sum / kf_matrix_get(L, i, i));
        }

        /* Backward substitution: L^T x = y  (in place). */
        for (i = L->rows; i-- > 0u; ) {
            kf_real_t sum = kf_matrix_get(B, i, col);
            uint16_t k;
            for (k = i + 1u; k < L->rows; k++) {
                sum -= kf_matrix_get(L, k, i) * kf_matrix_get(B, k, col);
            }
            kf_matrix_set(B, i, col, sum / kf_matrix_get(L, i, i));
        }
    }
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Vector operations
 * ------------------------------------------------------------------------ */

void kf_vec_copy(kf_real_t *dst, const kf_real_t *src, uint16_t n)
{
    uint16_t i;
    for (i = 0u; i < n; i++) {
        dst[i] = src[i];
    }
}

void kf_vec_zero(kf_real_t *v, uint16_t n)
{
    uint16_t i;
    for (i = 0u; i < n; i++) {
        v[i] = (kf_real_t)0;
    }
}

void kf_vec_scale(kf_real_t *v, kf_real_t s, uint16_t n)
{
    uint16_t i;
    for (i = 0u; i < n; i++) {
        v[i] *= s;
    }
}

void kf_vec_add(kf_real_t *out, const kf_real_t *a, const kf_real_t *b, uint16_t n)
{
    uint16_t i;
    for (i = 0u; i < n; i++) {
        out[i] = a[i] + b[i];
    }
}

void kf_vec_sub(kf_real_t *out, const kf_real_t *a, const kf_real_t *b, uint16_t n)
{
    uint16_t i;
    for (i = 0u; i < n; i++) {
        out[i] = a[i] - b[i];
    }
}

void kf_vec_axpy(kf_real_t *y, kf_real_t a, const kf_real_t *x, uint16_t n)
{
    uint16_t i;
    for (i = 0u; i < n; i++) {
        y[i] += a * x[i];
    }
}

kf_real_t kf_vec_dot(const kf_real_t *a, const kf_real_t *b, uint16_t n)
{
    kf_real_t sum = (kf_real_t)0;
    uint16_t i;
    for (i = 0u; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

kf_status_t kf_mat_vec_mul(kf_real_t *out, const kf_matrix_t *A, const kf_real_t *x)
{
    uint16_t i;

    KF_NULL_CHECK(out);
    KF_NULL_CHECK(A);
    KF_NULL_CHECK(x);
    KF_NULL_CHECK(A->data);

    for (i = 0u; i < A->rows; i++) {
        kf_real_t sum = (kf_real_t)0;
        uint16_t k;
        for (k = 0u; k < A->cols; k++) {
            sum += kf_matrix_get(A, i, k) * x[k];
        }
        out[i] = sum;
    }
    return KF_OK;
}

kf_status_t kf_vec_validate_finite(const kf_real_t *v, uint16_t n)
{
    uint16_t i;
    KF_NULL_CHECK(v);

    for (i = 0u; i < n; i++) {
        if (!kf_isfinite(v[i])) {
            return KF_ERROR_NON_FINITE;
        }
    }
    return KF_OK;
}

/* --------------------------------------------------------------------------
 * Diagnostics
 * ------------------------------------------------------------------------ */

int kf_matrix_is_finite(const kf_matrix_t *m)
{
    uint16_t r;
    if (m == NULL || m->data == NULL) {
        return 0;
    }
    for (r = 0u; r < m->rows; r++) {
        uint16_t c;
        for (c = 0u; c < m->cols; c++) {
            if (!kf_isfinite(kf_matrix_get(m, r, c))) {
                return 0;
            }
        }
    }
    return 1;
}

int kf_matrix_is_symmetric(const kf_matrix_t *m)
{
    uint16_t r;
    if (m == NULL || m->data == NULL) {
        return 0;
    }
    if (m->rows != m->cols) {
        return 0;
    }
    for (r = 0u; r < m->rows; r++) {
        uint16_t c;
        for (c = r + 1u; c < m->cols; c++) {
            kf_real_t diff = kf_matrix_get(m, r, c) - kf_matrix_get(m, c, r);
            if (kf_fabs(diff) > (kf_real_t)1e-4) {
                return 0;
            }
        }
    }
    return 1;
}
