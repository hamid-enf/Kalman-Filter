/**
 * @file    kalman_config.h
 * @brief   Compile-time configuration for the Kalman filtering library.
 *
 * This file is the single source of truth for all compile-time options. The
 * user is expected to edit this file (or provide these macros via the build
 * system, e.g. -DKF_ENABLE_UKF=0) to tailor the library to the application.
 *
 * Every non-essential feature is gated behind a configuration switch so that a
 * minimal build contains only the functionality that is actually required and
 * therefore has the smallest possible flash/RAM footprint.
 *
 * The library is written in C11 and is MISRA-C:2012 oriented. It performs no
 * dynamic memory allocation and contains no hidden global mutable state.
 */

#ifndef KALMAN_CONFIG_H
#define KALMAN_CONFIG_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * 1. NUMERIC PRECISION
 * ============================================================================
 *
 * The scalar type used by the whole library. The default is single precision
 * (float), which is the natural and fastest choice on STM32 Cortex-M cores
 * that feature a single-precision FPU (STM32F4/G4, etc.).
 *
 * Options (in order of precedence):
 *   - Define KF_REAL_TYPE to any arithmetic type (e.g. -DKF_REAL_TYPE=double)
 *     to fully control the scalar type.
 *   - Otherwise define KF_ENABLE_FLOAT64 to 1 to use double precision.
 *   - Otherwise single-precision float is used.
 */

#ifdef KF_REAL_TYPE
typedef KF_REAL_TYPE kf_real_t;
#elif defined(KF_ENABLE_FLOAT64) && (KF_ENABLE_FLOAT64 == 1)
typedef double kf_real_t;
#else
typedef float  kf_real_t;
#endif

/* ============================================================================
 * 2. MAXIMUM DIMENSIONS
 * ============================================================================
 *
 * The library uses statically-sized storage supplied by the caller, so the
 * maximum state and measurement dimensions must be known at compile time.
 * These values bound every internal buffer; keep them as small as the
 * application allows to minimise RAM usage.
 */

#ifndef KF_MAX_STATE_DIM
#define KF_MAX_STATE_DIM        6u   /**< Maximum state vector dimension n   */
#endif

#ifndef KF_MAX_MEASUREMENT_DIM
#define KF_MAX_MEASUREMENT_DIM  4u   /**< Maximum measurement dimension m    */
#endif

/* ============================================================================
 * 3. FILTER TYPE SWITCHES
 * ============================================================================
 *
 * Each filter type can be enabled/disabled independently. Disabling a filter
 * removes its public API and its implementation from the binary.
 */

#ifndef KF_ENABLE_KF
#define KF_ENABLE_KF            1     /**< Linear Kalman filter              */
#endif

#ifndef KF_ENABLE_EKF
#define KF_ENABLE_EKF           1     /**< Extended Kalman filter            */
#endif

#ifndef KF_ENABLE_UKF
#define KF_ENABLE_UKF           1     /**< Unscented Kalman filter           */
#endif

/* ============================================================================
 * 4. OPTIONAL FEATURES
 * ============================================================================
 */

#ifndef KF_ENABLE_FLOAT64
#define KF_ENABLE_FLOAT64       0     /**< Use double precision (see sec. 1) */
#endif

#ifndef KF_ENABLE_DIAGNOSTICS
#define KF_ENABLE_DIAGNOSTICS   0     /**< NaN/Inf/divergence diagnostics    */
#endif

#ifndef KF_ENABLE_VALIDATION
#define KF_ENABLE_VALIDATION    1     /**< Configuration/matrix validation   */
#endif

#ifndef KF_ENABLE_RUNTIME_CHECKS
#define KF_ENABLE_RUNTIME_CHECKS 1    /**< Cheap NULL/dimension guards       */
#endif

#ifndef KF_ENABLE_ADVANCED_API
#define KF_ENABLE_ADVANCED_API  1     /**< Low-level matrix/model access     */
#endif

/* ============================================================================
 * 5. NUMERICAL POLICIES
 * ============================================================================
 */

#ifndef KF_USE_JOSEPH_FORM
#define KF_USE_JOSEPH_FORM      1     /**< Use the Joseph-form covariance
                                           update (P = (I-KH)P(I-KH)^T + KRK^T)
                                           for guaranteed symmetry / PSD-ness.
                                           Slightly more expensive but more
                                           numerically robust. */
#endif

/*
 * Pivots smaller than this value are treated as (near-)singular when running
 * a Cholesky decomposition. Default is safe for single precision; raise it
 * only if you deliberately filter with very ill-conditioned covariances.
 */
#ifndef KF_MIN_PIVOT
#define KF_MIN_PIVOT            ((kf_real_t)1e-12)
#endif

/* ============================================================================
 * 6. LANGUAGE / TOOLING
 * ============================================================================
 */

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ < 201112L)
#error "kalman: C11 or later is required"
#endif

#ifdef __cplusplus
}
#endif

#endif /* KALMAN_CONFIG_H */
