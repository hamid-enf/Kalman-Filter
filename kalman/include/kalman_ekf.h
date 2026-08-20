/**
 * @file    kalman_ekf.h
 * @brief   Extended Kalman filter (EKF) for nonlinear systems.
 *
 * The EKF extends the linear KF to models where the state transition and/or
 * the measurement are **nonlinear**:
 *
 *      x(k+1) = f(x(k), u, dt) + w          w ~ N(0, Q)
 *      z(k)   = h(x(k))          + v          v ~ N(0, R)
 *
 * It linearises the models around the current estimate using the Jacobians
 *      F = df/dx   and   H = dh/dx
 * and then applies the standard linear KF recursion. It is a first-order
 * approximation: it performs best when the nonlinearity is mild over one step
 * or the estimate uncertainty is small.
 *
 * The user supplies four callback functions:
 *   - f:      the state-transition function
 *   - h:      the measurement function
 *   - F_jac:  the Jacobian of f w.r.t. the state
 *   - H_jac:  the Jacobian of h w.r.t. the state
 *
 * Complexity per call:
 *   - Predict: O(n^3) + O(n^2) + cost of f/F_jac
 *   - Update : O(n^2 m + n m^2 + m^3) + cost of h/H_jac
 */

#ifndef KALMAN_EKF_H
#define KALMAN_EKF_H

#include "kalman_config.h"
#include "kalman_types.h"
#include "kalman_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

#if KF_ENABLE_EKF

/* Size of the internal scratch buffer (in kf_real_t). */
#define KF_EKF_SCRATCH_FLOATS \
    ((2u * KF_MAX_STATE_DIM * KF_MAX_STATE_DIM) \
   + (3u * KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM) \
   + (KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM) \
   + (2u * KF_MAX_STATE_DIM) + (2u * KF_MAX_MEASUREMENT_DIM))

/* ==========================================================================
 * Model callbacks
 * ========================================================================== */

/**
 * @brief State-transition function. Computes x_out = f(x, u, dt).
 *
 * @param x      Current state (n values).
 * @param u      Control input (n values, may be NULL if unused).
 * @param dt     Time step.
 * @param x_out  Output predicted state (n values, caller-provided buffer).
 * @param ctx    Opaque user context.
 */
typedef void (*kf_ekf_f_fn)(const kf_real_t *x, const kf_real_t *u,
                            kf_real_t dt, kf_real_t *x_out, void *ctx);

/**
 * @brief Measurement function. Computes z_out = h(x).
 *
 * @param x      Current state (n values).
 * @param z_out  Output predicted measurement (m values, caller-provided).
 * @param ctx    Opaque user context.
 */
typedef void (*kf_ekf_h_fn)(const kf_real_t *x, kf_real_t *z_out, void *ctx);

/**
 * @brief Jacobian of the state-transition function. Fills F = df/dx evaluated
 * at x (an n x n matrix, pre-sized and bound to internal storage).
 *
 * Fill it with kf_matrix_set(&F, i, j, value). Do NOT resize it or change its
 * data pointer.
 */
typedef void (*kf_ekf_jacobian_f_fn)(const kf_real_t *x, const kf_real_t *u,
                                     kf_real_t dt, kf_matrix_t *F, void *ctx);

/**
 * @brief Jacobian of the measurement function. Fills H = dh/dx evaluated at x
 * (an m x n matrix, pre-sized and bound to internal storage).
 */
typedef void (*kf_ekf_jacobian_h_fn)(const kf_real_t *x, kf_matrix_t *H, void *ctx);

/**
 * @brief EKF instance. All storage is embedded; no dynamic allocation.
 */
typedef struct {
    uint16_t n;                /**< State dimension                              */
    uint16_t m;                /**< Measurement dimension                        */
    uint16_t initialized;      /**< Non-zero after a successful kf_ekf_init()   */

    kf_real_t x[KF_MAX_STATE_DIM];                                   /* n x 1   */
    kf_real_t P[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];                /* n x n   */
    kf_real_t Q[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];                /* n x n   */
    kf_real_t F[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];                /* n x n   */
    kf_real_t R[KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM];    /* m x m   */
    kf_real_t H[KF_MAX_MEASUREMENT_DIM * KF_MAX_STATE_DIM];          /* m x n   */

    kf_ekf_f_fn           f;         /**< state-transition function             */
    kf_ekf_h_fn           h;         /**< measurement function                  */
    kf_ekf_jacobian_f_fn  F_jac;     /**< df/dx                                 */
    kf_ekf_jacobian_h_fn  H_jac;     /**< dh/dx                                 */
    void                 *ctx_f;     /**< context for f / F_jac                 */
    void                 *ctx_h;     /**< context for h / H_jac                 */

#if KF_ENABLE_ADVANCED_API
    kf_real_t y[KF_MAX_MEASUREMENT_DIM];                  /* last innovation   */
    kf_real_t K[KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM];/* last gain n x m */
#endif

    kf_real_t scratch[KF_EKF_SCRATCH_FLOATS];   /**< internal temporaries        */
} kf_ekf_t;

/* ==========================================================================
 * Initialisation / reset
 * ========================================================================== */

/**
 * @brief Initialise an EKF for an n-state / m-measurement system.
 *
 * Applies defaults x = 0, F = I, H = 0, P = I, Q = 0, R = 0 and clears all
 * model callbacks. The caller must then register the four model callbacks
 * (kf_ekf_set_models) and set H/Q/R/P before filtering.
 */
kf_status_t kf_ekf_init(kf_ekf_t *ekf, uint16_t n, uint16_t m);

/** @brief Reset state/covariances to their initial (identity) values. */
kf_status_t kf_ekf_reset(kf_ekf_t *ekf);

/* ==========================================================================
 * Model registration
 * ========================================================================== */

/**
 * @brief Register the nonlinear model functions and their Jacobians, plus
 * optional per-callback contexts. Any callback may be NULL to unset it.
 */
kf_status_t kf_ekf_set_models(kf_ekf_t *ekf,
                              kf_ekf_f_fn f, kf_ekf_jacobian_f_fn F_jac,
                              kf_ekf_h_fn h, kf_ekf_jacobian_h_fn H_jac);

/** @brief Set the user context pointer passed to the transition callbacks. */
void        kf_ekf_set_context_f(kf_ekf_t *ekf, void *ctx);

/** @brief Set the user context pointer passed to the measurement callbacks. */
void        kf_ekf_set_context_h(kf_ekf_t *ekf, void *ctx);

/* ==========================================================================
 * State / covariance / noise (same semantics as the linear KF)
 * ========================================================================== */

kf_status_t kf_ekf_set_state(kf_ekf_t *ekf, const kf_real_t *x);
const kf_real_t *kf_ekf_get_state(const kf_ekf_t *ekf);

kf_status_t kf_ekf_set_covariance(kf_ekf_t *ekf, const kf_real_t *P);
kf_status_t kf_ekf_set_covariance_diagonal(kf_ekf_t *ekf, const kf_real_t *diag);
kf_status_t kf_ekf_set_covariance_scalar(kf_ekf_t *ekf, kf_real_t p);
const kf_real_t *kf_ekf_get_covariance(const kf_ekf_t *ekf);

kf_status_t kf_ekf_set_process_noise(kf_ekf_t *ekf, const kf_real_t *Q);
kf_status_t kf_ekf_set_process_noise_diagonal(kf_ekf_t *ekf, const kf_real_t *diag);
kf_status_t kf_ekf_set_process_noise_scalar(kf_ekf_t *ekf, kf_real_t q);

kf_status_t kf_ekf_set_measurement_noise(kf_ekf_t *ekf, const kf_real_t *R);
kf_status_t kf_ekf_set_measurement_noise_diagonal(kf_ekf_t *ekf, const kf_real_t *diag);
kf_status_t kf_ekf_set_measurement_noise_scalar(kf_ekf_t *ekf, kf_real_t r);

/* ==========================================================================
 * Filtering (real-time path)
 * ========================================================================== */

/**
 * @brief Predict step: linearise with F_jac, propagate the state with f, and
 * propagate the covariance:  P = F P F^T + Q.
 *
 * @param u   Optional control input (n values, may be NULL).
 * @param dt  Time step, passed to the model callbacks.
 */
kf_status_t kf_ekf_predict(kf_ekf_t *ekf, const kf_real_t *u, kf_real_t dt);

/**
 * @brief Update step: linearise with H_jac, predict the measurement with h,
 * then perform the standard Kalman correction (with Joseph-form covariance
 * update when KF_USE_JOSEPH_FORM is enabled).
 *
 * @param z  Measurement vector (m values).
 */
kf_status_t kf_ekf_update(kf_ekf_t *ekf, const kf_real_t *z);

#if KF_ENABLE_ADVANCED_API
kf_status_t kf_ekf_get_innovation(const kf_ekf_t *ekf, kf_real_t *y_out);
kf_status_t kf_ekf_get_gain(const kf_ekf_t *ekf, kf_real_t *K_out);
#endif

#if KF_ENABLE_DIAGNOSTICS
kf_status_t kf_ekf_check(const kf_ekf_t *ekf);
#endif

#endif /* KF_ENABLE_EKF */

#ifdef __cplusplus
}
#endif

#endif /* KALMAN_EKF_H */
