/**
 * @file    kalman_ukf.h
 * @brief   Unscented Kalman filter (UKF) for nonlinear systems.
 *
 * The UKF handles nonlinear state-transition and measurement models
 *
 *      x(k+1) = f(x(k), u, dt) + w          w ~ N(0, Q)
 *      z(k)   = h(x(k))          + v          v ~ N(0, R)
 *
 * without requiring Jacobians. Instead it deterministically samples the state
 * distribution at a small set of *sigma points* (2n + 1 points), propagates
 * them through the nonlinear functions, and reconstructs the mean and
 * covariance from the results. This captures the mean and covariance
 * accurately to 2nd order (and 3rd order for Gaussian inputs) without
 * linearisation error, at the cost of 2n+1 evaluations of f and h per step.
 *
 * The sigma-point spread is controlled by three scalar parameters:
 *   - alpha: spread of the sigma points (default 1.0; smaller = tighter,
 *            commonly 1e-3..1). Very small alpha introduces a large negative
 *            Wm[0] weight, which is numerically delicate.
 *   - beta : incorporates prior knowledge of the distribution (2 is optimal
 *            for Gaussians).
 *   - kappa: secondary scaling parameter (0 is a safe default).
 *
 * The sigma points are generated from the Cholesky factor of P (P = L L^T),
 * which is numerically stable and avoids eigen-decomposition.
 *
 * Complexity per call:
 *   - Predict: (2n+1) evaluations of f + O(n^3) + O(n^2 * (2n+1))
 *   - Update : (2n+1) evaluations of h + O(n * m * (2n+1)) + O(m^3)
 */

#ifndef KALMAN_UKF_H
#define KALMAN_UKF_H

#include "kalman_config.h"
#include "kalman_types.h"
#include "kalman_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

#if KF_ENABLE_UKF

/* Number of sigma points for the maximum state dimension. */
#define KF_UKF_MAX_SIGMA   (2u * KF_MAX_STATE_DIM + 1u)

/* Size of the internal scratch buffer (in kf_real_t). */
#define KF_UKF_SCRATCH_FLOATS \
    (2u * KF_MAX_STATE_DIM * KF_MAX_STATE_DIM \
   + 3u * KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM \
   +       KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM \
   + 2u * KF_MAX_STATE_DIM + 2u * KF_MAX_MEASUREMENT_DIM)

/* ==========================================================================
 * Model callbacks (same signatures as the EKF process/measurement functions;
 * the types are compatible, so the same function may be reused).
 * ========================================================================== */

/** @brief State-transition function: x_out = f(x, u, dt). */
typedef void (*kf_ukf_f_fn)(const kf_real_t *x, const kf_real_t *u,
                            kf_real_t dt, kf_real_t *x_out, void *ctx);

/** @brief Measurement function: z_out = h(x). */
typedef void (*kf_ukf_h_fn)(const kf_real_t *x, kf_real_t *z_out, void *ctx);

/**
 * @brief UKF instance. All storage is embedded; no dynamic allocation.
 */
typedef struct {
    uint16_t n;                /**< State dimension                              */
    uint16_t m;                /**< Measurement dimension                        */
    uint16_t initialized;      /**< Non-zero after a successful kf_ukf_init()   */

    kf_real_t x[KF_MAX_STATE_DIM];                                   /* n x 1   */
    kf_real_t P[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];                /* n x n   */
    kf_real_t Q[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];                /* n x n   */
    kf_real_t R[KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM];    /* m x m   */

    /* Sigma-point configuration. */
    kf_real_t alpha;           /**< sigma-point spread                          */
    kf_real_t beta;            /**< distribution prior (2 for Gaussian)         */
    kf_real_t kappa;           /**< secondary scaling                           */
    kf_real_t lambda;          /**< = alpha^2 (n + kappa) - n (derived)         */

    kf_real_t Wm[KF_UKF_MAX_SIGMA];   /**< mean weights (2n+1)                  */
    kf_real_t Wc[KF_UKF_MAX_SIGMA];   /**< covariance weights (2n+1)            */

    /* Sigma-point storage (row-major: point i at row i). */
    kf_real_t sigma[KF_UKF_MAX_SIGMA * KF_MAX_STATE_DIM];            /* (2n+1) x n */
    kf_real_t sigma_prop[KF_UKF_MAX_SIGMA * KF_MAX_STATE_DIM];       /* (2n+1) x n */
    kf_real_t sigma_z[KF_UKF_MAX_SIGMA * KF_MAX_MEASUREMENT_DIM];    /* (2n+1) x m */

    kf_ukf_f_fn f;             /**< state-transition function                   */
    kf_ukf_h_fn h;             /**< measurement function                        */
    void       *ctx_f;         /**< context for f                               */
    void       *ctx_h;         /**< context for h                               */

#if KF_ENABLE_ADVANCED_API
    kf_real_t y[KF_MAX_MEASUREMENT_DIM];                  /* last innovation   */
    kf_real_t K[KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM];/* last gain n x m */
#endif

    kf_real_t scratch[KF_UKF_SCRATCH_FLOATS];   /**< internal temporaries        */
} kf_ukf_t;

/* ==========================================================================
 * Initialisation / reset
 * ========================================================================== */

/**
 * @brief Initialise a UKF for an n-state / m-measurement system.
 *
 * Defaults: x = 0, P = I, Q = 0, R = 0, sigma-point parameters
 * alpha = 1, beta = 2, kappa = 0 (a robust, all-positive-weight configuration).
 * The caller must then set the model functions and Q/R/P.
 */
kf_status_t kf_ukf_init(kf_ukf_t *ukf, uint16_t n, uint16_t m);

/** @brief Reset state/covariance to their initial values (keeps parameters). */
kf_status_t kf_ukf_reset(kf_ukf_t *ukf);

/* ==========================================================================
 * Model / parameters
 * ========================================================================== */

/** @brief Register the nonlinear process and measurement functions. */
kf_status_t kf_ukf_set_models(kf_ukf_t *ukf, kf_ukf_f_fn f, kf_ukf_h_fn h);

/** @brief Set the context pointer passed to the transition function. */
void        kf_ukf_set_context_f(kf_ukf_t *ukf, void *ctx);

/** @brief Set the context pointer passed to the measurement function. */
void        kf_ukf_set_context_h(kf_ukf_t *ukf, void *ctx);

/**
 * @brief Set the sigma-point scaling parameters (alpha, beta, kappa) and
 * recompute lambda and the weights. Call before filtering when tuning.
 *
 * Typical values: alpha = 1.0 (robust) or 1e-3 (tight), beta = 2.0, kappa = 0.
 */
kf_status_t kf_ukf_set_parameters(kf_ukf_t *ukf, kf_real_t alpha,
                                  kf_real_t beta, kf_real_t kappa);

/* ==========================================================================
 * State / covariance / noise (same semantics as the linear KF)
 * ========================================================================== */

kf_status_t kf_ukf_set_state(kf_ukf_t *ukf, const kf_real_t *x);
const kf_real_t *kf_ukf_get_state(const kf_ukf_t *ukf);

kf_status_t kf_ukf_set_covariance(kf_ukf_t *ukf, const kf_real_t *P);
kf_status_t kf_ukf_set_covariance_diagonal(kf_ukf_t *ukf, const kf_real_t *diag);
kf_status_t kf_ukf_set_covariance_scalar(kf_ukf_t *ukf, kf_real_t p);
const kf_real_t *kf_ukf_get_covariance(const kf_ukf_t *ukf);

kf_status_t kf_ukf_set_process_noise(kf_ukf_t *ukf, const kf_real_t *Q);
kf_status_t kf_ukf_set_process_noise_diagonal(kf_ukf_t *ukf, const kf_real_t *diag);
kf_status_t kf_ukf_set_process_noise_scalar(kf_ukf_t *ukf, kf_real_t q);

kf_status_t kf_ukf_set_measurement_noise(kf_ukf_t *ukf, const kf_real_t *R);
kf_status_t kf_ukf_set_measurement_noise_diagonal(kf_ukf_t *ukf, const kf_real_t *diag);
kf_status_t kf_ukf_set_measurement_noise_scalar(kf_ukf_t *ukf, kf_real_t r);

/* ==========================================================================
 * Filtering (real-time path)
 * ========================================================================== */

/**
 * @brief Predict step: generate sigma points, propagate them through f, and
 * reconstruct the predicted mean and covariance (adding Q).
 *
 * @param u   Optional control input (n values, may be NULL).
 * @param dt  Time step, passed to the model callbacks.
 */
kf_status_t kf_ukf_predict(kf_ukf_t *ukf, const kf_real_t *u, kf_real_t dt);

/**
 * @brief Update step: propagate sigma points through h, reconstruct the
 * predicted measurement, and correct the state/covariance with measurement z.
 *
 * @param z  Measurement vector (m values).
 */
kf_status_t kf_ukf_update(kf_ukf_t *ukf, const kf_real_t *z);

#if KF_ENABLE_ADVANCED_API
kf_status_t kf_ukf_get_innovation(const kf_ukf_t *ukf, kf_real_t *y_out);
kf_status_t kf_ukf_get_gain(const kf_ukf_t *ukf, kf_real_t *K_out);
#endif

#if KF_ENABLE_DIAGNOSTICS
kf_status_t kf_ukf_check(const kf_ukf_t *ukf);
#endif

#endif /* KF_ENABLE_UKF */

#ifdef __cplusplus
}
#endif

#endif /* KALMAN_UKF_H */
