/**
 * @file    kalman_types.h
 * @brief   Common types and status codes shared by the whole library.
 */

#ifndef KALMAN_TYPES_H
#define KALMAN_TYPES_H

#include "kalman_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * STATUS CODES
 * ============================================================================
 *
 * Every public function that can fail returns a `kf_status_t`. The value
 * `KF_OK` (0) indicates success; all other values are negative error codes.
 * The library never "silently fails": callers should check the return value.
 */

typedef enum {
    KF_OK                       =  0,  /**< Success                          */
    KF_ERROR_NULL_POINTER       = -1,  /**< A required pointer was NULL      */
    KF_ERROR_INVALID_DIMENSION  = -2,  /**< Matrix/vector dimension mismatch */
    KF_ERROR_INVALID_PARAMETER  = -3,  /**< Parameter out of range           */
    KF_ERROR_NOT_INITIALIZED    = -4,  /**< Filter instance not initialised  */
    KF_ERROR_FEATURE_DISABLED   = -5,  /**< Requested feature compiled out   */
    KF_ERROR_SINGULAR_MATRIX    = -6,  /**< Matrix is singular               */
    KF_ERROR_NOT_POSITIVE_DEFINITE = -7, /**< Matrix is not positive definite */
    KF_ERROR_NUMERICAL          = -8,  /**< Numerical failure (see note)     */
    KF_ERROR_NON_FINITE         = -9,  /**< NaN or Inf detected              */
    KF_ERROR_DIVERGED           = -10  /**< Filter diverged (diagnostics)    */
} kf_status_t;

/**
 * @brief  Returns a human-readable, static string for a status code.
 *
 * The returned pointer is a constant string literal; it must not be freed and
 * is safe to use for logging. The function has no side effects and is safe to
 * call from any context.
 */
const char *kf_status_str(kf_status_t status);

/* ============================================================================
 * Internal guard helper (used by the implementations, not part of the public
 * contract). With KF_ENABLE_RUNTIME_CHECKS=0 the null-pointer guards compile
 * out for the fastest possible hot path.
 * ========================================================================== */

#if KF_ENABLE_RUNTIME_CHECKS
#define KF_NULL_CHECK(p) \
    do { if ((p) == NULL) { return KF_ERROR_NULL_POINTER; } } while (0)
#else
#define KF_NULL_CHECK(p) ((void)0)
#endif

#ifdef __cplusplus
}
#endif

#endif /* KALMAN_TYPES_H */
