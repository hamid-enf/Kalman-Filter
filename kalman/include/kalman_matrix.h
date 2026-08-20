/**
 * @file    kalman_matrix.h
 * @brief   Minimal, dependency-free matrix engine used internally by the
 *          KF/EKF/UKF implementations.
 *
 * Design notes
 * ------------
 *  - Matrices are stored **row-major**.
 *  - A `kf_matrix_t` is a *view*: it only holds dimensions and a pointer to
 *    caller-supplied storage. The library itself never allocates memory, so
 *    views let the filters share a single scratch buffer and avoid copying.
 *  - Only the operations actually required by the filters are provided. This
 *    is intentionally NOT a general-purpose linear algebra library.
 *  - The filters prefer numerically-stable linear-system solves over explicit
 *    matrix inversion; see kf_matrix_cholesky()/kf_matrix_cholesky_solve().
 *
 * Aliasing rules
 * --------------
 *  - Element-wise operations (add, sub, scale, zero, identity, copy,
 *    symmetrize, transpose when dst != src) may alias inputs with the output.
 *  - Matrix multiplication (kf_matrix_mul and the *_transpose_* helpers) MUST
 *    NOT alias the output with either input. Use a scratch buffer.
 *
 * All functions are reentrant and deterministic. None of them allocate
 * memory, print anything, or depend on any RTOS / HAL / CMSIS facility.
 */

#ifndef KALMAN_MATRIX_H
#define KALMAN_MATRIX_H

#include "kalman_config.h"
#include "kalman_types.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------------
 * Matrix view type
 * ------------------------------------------------------------------------ */

typedef struct {
    uint16_t    rows;   /**< Number of rows                                    */
    uint16_t    cols;   /**< Number of columns                                 */
    kf_real_t  *data;   /**< Row-major storage, length rows*cols (caller-owned)*/
} kf_matrix_t;

/* --------------------------------------------------------------------------
 * Numeric helpers (inline, no allocation)
 * ------------------------------------------------------------------------ */

/** @brief Square root of a kf_real_t (uses the correct precision). */
static inline kf_real_t kf_sqrt(kf_real_t x)
{
#if defined(KF_ENABLE_FLOAT64) && (KF_ENABLE_FLOAT64 == 1)
    return (kf_real_t)sqrt((double)x);
#else
    return (kf_real_t)sqrtf((float)x);
#endif
}

/** @brief Absolute value of a kf_real_t. */
static inline kf_real_t kf_fabs(kf_real_t x)
{
    return (x < (kf_real_t)0) ? -x : x;
}

/** @brief Returns 1 if x is finite (not NaN, not +/-Inf), 0 otherwise. */
static inline int kf_isfinite(kf_real_t x)
{
    return isfinite((double)x) ? 1 : 0;
}

/* --------------------------------------------------------------------------
 * Construction / access
 * ------------------------------------------------------------------------ */

/** @brief Bind a view to caller-provided storage. */
kf_status_t kf_matrix_init(kf_matrix_t *m, uint16_t rows, uint16_t cols,
                           kf_real_t *data);

/** @brief Read element (r,c). No bounds checking on r/c. */
kf_real_t   kf_matrix_get(const kf_matrix_t *m, uint16_t r, uint16_t c);

/** @brief Write element (r,c). No bounds checking on r/c. */
void        kf_matrix_set(kf_matrix_t *m, uint16_t r, uint16_t c, kf_real_t v);

/* --------------------------------------------------------------------------
 * Basic element-wise operations
 * ------------------------------------------------------------------------ */

kf_status_t kf_matrix_zero(kf_matrix_t *m);
kf_status_t kf_matrix_identity(kf_matrix_t *m);                    /* square only */
kf_status_t kf_matrix_copy(kf_matrix_t *dst, const kf_matrix_t *src);
kf_status_t kf_matrix_scale(kf_matrix_t *m, kf_real_t s);
kf_status_t kf_matrix_add(kf_matrix_t *out, const kf_matrix_t *a,
                          const kf_matrix_t *b);                   /* out = a + b */
kf_status_t kf_matrix_add_scaled(kf_matrix_t *out, const kf_matrix_t *a,
                                 const kf_matrix_t *b, kf_real_t s); /* out = a + s*b */
kf_status_t kf_matrix_sub(kf_matrix_t *out, const kf_matrix_t *a,
                          const kf_matrix_t *b);                   /* out = a - b */
kf_status_t kf_matrix_symmetrize(kf_matrix_t *m);                  /* m = (m+m^T)/2 */

/* --------------------------------------------------------------------------
 * Multiplication (output must NOT alias inputs)
 * ------------------------------------------------------------------------ */

/** @brief out = a * b.  a:(p x q), b:(q x r), out:(p x r). */
kf_status_t kf_matrix_mul(kf_matrix_t *out, const kf_matrix_t *a,
                          const kf_matrix_t *b);

/** @brief out = a * b^T.  a:(p x q), b:(r x q), out:(p x r). */
kf_status_t kf_matrix_mul_transpose_b(kf_matrix_t *out, const kf_matrix_t *a,
                                      const kf_matrix_t *b);

/** @brief out = a^T * b.  a:(q x p), b:(q x r), out:(p x r). */
kf_status_t kf_matrix_mul_transpose_a(kf_matrix_t *out, const kf_matrix_t *a,
                                      const kf_matrix_t *b);

/* --------------------------------------------------------------------------
 * Transpose
 * ------------------------------------------------------------------------ */

/** @brief out = a^T. May be called with out == a (in-place). */
kf_status_t kf_matrix_transpose(kf_matrix_t *out, const kf_matrix_t *a);

/* --------------------------------------------------------------------------
 * Linear-system solves (numerically stable, no explicit inverse)
 * ------------------------------------------------------------------------ */

/**
 * @brief In-place Cholesky factorisation of a symmetric positive-definite
 *        matrix: on return the lower triangle holds L such that A = L L^T.
 *
 * Used for computing sqrt(P) for UKF sigma points and as the first step of a
 * stable linear solve. Returns KF_ERROR_NOT_POSITIVE_DEFINITE if a pivot is
 * non-positive or non-finite.
 */
kf_status_t kf_matrix_cholesky(kf_matrix_t *m);

/**
 * @brief Solves L L^T X = B for X, where L is lower-triangular (from
 *        kf_matrix_cholesky) and B is a (k x nrhs) matrix. B is overwritten
 *        with X in place.
 *
 * This is the preferred way to compute products of the form "B * A^-1" for a
 * symmetric A, e.g. the Kalman gain K = (P H^T) S^-1, without ever forming
 * the explicit inverse.
 */
kf_status_t kf_matrix_cholesky_solve(const kf_matrix_t *L, kf_matrix_t *B);

/* --------------------------------------------------------------------------
 * Vector operations (vectors are plain kf_real_t arrays)
 * ------------------------------------------------------------------------ */

void        kf_vec_copy(kf_real_t *dst, const kf_real_t *src, uint16_t n);
void        kf_vec_zero(kf_real_t *v, uint16_t n);
void        kf_vec_scale(kf_real_t *v, kf_real_t s, uint16_t n);
void        kf_vec_add(kf_real_t *out, const kf_real_t *a, const kf_real_t *b, uint16_t n);
void        kf_vec_sub(kf_real_t *out, const kf_real_t *a, const kf_real_t *b, uint16_t n);
void        kf_vec_axpy(kf_real_t *y, kf_real_t a, const kf_real_t *x, uint16_t n);
kf_real_t   kf_vec_dot(const kf_real_t *a, const kf_real_t *b, uint16_t n);

/** @brief out = A * x.  A:(m x n), x:(n), out:(m). */
kf_status_t kf_mat_vec_mul(kf_real_t *out, const kf_matrix_t *A, const kf_real_t *x);

/**
 * @brief Validates that all n elements of v are finite. Returns KF_OK if so,
 * KF_ERROR_NON_FINITE if any element is NaN/Inf, KF_ERROR_NULL_POINTER if v is
 * NULL. Used by the setters when KF_ENABLE_VALIDATION is on.
 */
kf_status_t kf_vec_validate_finite(const kf_real_t *v, uint16_t n);

/* --------------------------------------------------------------------------
 * Diagnostics (only meaningful when KF_ENABLE_DIAGNOSTICS or
 * KF_ENABLE_VALIDATION is enabled; otherwise thin wrappers returning success)
 * ------------------------------------------------------------------------ */

/** @brief Returns 1 if every element is finite, 0 otherwise. */
int kf_matrix_is_finite(const kf_matrix_t *m);

/** @brief Returns 1 if m is square and symmetric (within tolerance), else 0. */
int kf_matrix_is_symmetric(const kf_matrix_t *m);

#ifdef __cplusplus
}
#endif

#endif /* KALMAN_MATRIX_H */
