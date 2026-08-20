/**
 * @file    kalman_kf.h
 * @brief   Linear Kalman filter (KF).
 *
 * This is the simplest and most common estimator. It assumes both the state
 * transition and the measurement models are **linear**:
 *
 *      x(k+1) = F * x(k) + u          (state transition, u optional input)
 *      z(k)   = H * x(k) + v          (measurement, v ~ N(0, R))
 *
 * The filter maintains a Gaussian belief about the state, summarised by the
 * mean `x` and covariance `P`. `Q` is the process-noise covariance and `R` the
 * measurement-noise covariance.
 *
 * The implementation is generic in the state dimension `n` and the measurement
 * dimension `m`, which are set at initialisation time. It supports sequential
 * updates, i.e. calling kf_kf_update() multiple times with different sensors
 * (see examples 05 and 06).
 *
 * Complexity (big-O per call):
 *   - Predict: O(n^3) + O(n^2)
 *   - Update : O(n^2 * m) + O(n * m^2) + O(m^3)
 *
 * Thread-safety: a single instance must not be used concurrently from multiple
 * contexts (ISR + main loop + RTOS task) without external synchronisation. The
 * functions are reentrant and deterministic and touch only the instance passed
 * in; distinct instances may be used from distinct contexts freely.
 */

#ifndef KALMAN_KF_H
#define KALMAN_KF_H

#include "kalman_config.h"
#include "kalman_types.h"
#include "kalman_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Size of the internal scratch buffer (in kf_real_t) owned by each instance. */
#define KF_KF_SCRATCH_FLOATS \
    (2u * KF_MAX_STATE_DIM * KF_MAX_STATE_DIM \
   + 3u * KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM \
   +       KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM \
   + KF_MAX_STATE_DIM + KF_MAX_MEASUREMENT_DIM)

#if KF_ENABLE_KF

/**
 * @brief Linear Kalman filter instance.
 *
 * All storage is embedded in the instance: no heap allocation is performed.
 * The instance may be placed in static storage, on the stack, or in a
 * user-managed pool. Its size is fixed and known at compile time.
 */
typedef struct {
    uint16_t n;                /**< State dimension                              */
    uint16_t m;                /**< Measurement dimension                        */
    uint16_t initialized;      /**< Non-zero after a successful kf_kf_init()    */

    kf_real_t x[KF_MAX_STATE_DIM];                                   /* n x 1   */
    kf_real_t P[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];                /* n x n   */
    kf_real_t Q[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];                /* n x n   */
    kf_real_t F[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];                /* n x n   */
    kf_real_t R[KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM];    /* m x m   */
    kf_real_t H[KF_MAX_MEASUREMENT_DIM * KF_MAX_STATE_DIM];          /* m x n   */

#if KF_ENABLE_ADVANCED_API
    kf_real_t y[KF_MAX_MEASUREMENT_DIM];                  /* last innovation   */
    kf_real_t K[KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM];/* last gain n x m */
#endif

    kf_real_t scratch[KF_KF_SCRATCH_FLOATS];   /**< internal temporaries        */
} kf_kf_t;

/* ==========================================================================
 * Initialisation / reset
 * ========================================================================== */

/**
 * @brief Initialise a linear KF for an n-state / m-measurement system.
 *
 * Applies mathematically safe defaults:
 *   - x = 0
 *   - F = I (identity)
 *   - H = 0
 *   - P = I (identity)
 *   - Q = 0, R = 0
 *
 * The caller MUST subsequently set at least H, Q and R (and P if the identity
 * is inappropriate) before filtering. Defaults are kept explicit rather than
 * "magic" so no estimation assumption is hidden.
 */
kf_status_t kf_kf_init(kf_kf_t *kf, uint16_t n, uint16_t m);

/** @brief Reset the filter to its initial (identity) configuration. */
kf_status_t kf_kf_reset(kf_kf_t *kf);

/* ==========================================================================
 * State access
 * ========================================================================== */

/** @brief Set the full state vector (n values). */
kf_status_t kf_kf_set_state(kf_kf_t *kf, const kf_real_t *x);

/** @brief Return a pointer to the internal state vector (n values). */
const kf_real_t *kf_kf_get_state(const kf_kf_t *kf);

/** @brief Set a single state element. */
kf_status_t kf_kf_set_state_element(kf_kf_t *kf, uint16_t index, kf_real_t value);

/* ==========================================================================
 * Covariance access
 * ========================================================================== */

/** @brief Set the full covariance matrix P (n x n, row-major). */
kf_status_t kf_kf_set_covariance(kf_kf_t *kf, const kf_real_t *P);

/** @brief Set P from a diagonal vector (off-diagonals set to zero). */
kf_status_t kf_kf_set_covariance_diagonal(kf_kf_t *kf, const kf_real_t *diag);

/** @brief Set P to a scalar multiple of the identity (P = p * I). */
kf_status_t kf_kf_set_covariance_scalar(kf_kf_t *kf, kf_real_t p);

/** @brief Return a pointer to the internal covariance P (n x n, row-major). */
const kf_real_t *kf_kf_get_covariance(const kf_kf_t *kf);

/* ==========================================================================
 * Noise configuration
 * ========================================================================== */

/** @brief Set the full process-noise covariance Q (n x n, row-major). */
kf_status_t kf_kf_set_process_noise(kf_kf_t *kf, const kf_real_t *Q);

/** @brief Set Q from a diagonal vector. */
kf_status_t kf_kf_set_process_noise_diagonal(kf_kf_t *kf, const kf_real_t *diag);

/** @brief Set Q = q * I (scalar process noise). */
kf_status_t kf_kf_set_process_noise_scalar(kf_kf_t *kf, kf_real_t q);

/** @brief Set the full measurement-noise covariance R (m x m, row-major). */
kf_status_t kf_kf_set_measurement_noise(kf_kf_t *kf, const kf_real_t *R);

/** @brief Set R from a diagonal vector. */
kf_status_t kf_kf_set_measurement_noise_diagonal(kf_kf_t *kf, const kf_real_t *diag);

/** @brief Set R = r * I (scalar measurement noise). */
kf_status_t kf_kf_set_measurement_noise_scalar(kf_kf_t *kf, kf_real_t r);

/* ==========================================================================
 * Model configuration
 * ========================================================================== */

/** @brief Set the state-transition matrix F (n x n, row-major). */
kf_status_t kf_kf_set_transition_matrix(kf_kf_t *kf, const kf_real_t *F);

/** @brief Set the measurement matrix H (m x n, row-major). */
kf_status_t kf_kf_set_measurement_matrix(kf_kf_t *kf, const kf_real_t *H);

/* ==========================================================================
 * Filtering (real-time path)
 * ========================================================================== */

/**
 * @brief Perform a prediction (time update) step.
 *
 *      x = F x + u        (u may be NULL to skip the input term)
 *      P = F P F^T + Q
 *
 * @param kf  Filter instance.
 * @param u   Optional n-vector control input, already expressed in state
 *            space (i.e. pre-multiplied by any B matrix and dt). Pass NULL if
 *            there is no input.
 * @param dt  Time step. For a linear KF the step size is encoded in F (and Q);
 *            dt is accepted for API symmetry with EKF/UKF and is otherwise
 *            unused. For variable-dt models, update F (and Q) via the
 *            kf_kf_set_transition_matrix()/kf_kf_set_process_noise() setters
 *            before calling predict.
 */
kf_status_t kf_kf_predict(kf_kf_t *kf, const kf_real_t *u, kf_real_t dt);

/**
 * @brief Perform a measurement update (correction) step with measurement z.
 *
 *      y = z - H x              (innovation)
 *      S = H P H^T + R          (innovation covariance)
 *      K = P H^T S^-1           (Kalman gain, solved, not inverted)
 *      x = x + K y
 *      P = (I-KH) P (I-KH)^T + K R K^T   (Joseph form, if enabled)
 *
 * @param kf  Filter instance.
 * @param z   Measurement vector (m values). May be NULL to skip the update
 *            (prediction only).
 */
kf_status_t kf_kf_update(kf_kf_t *kf, const kf_real_t *z);

/* ==========================================================================
 * Advanced / diagnostics (optional feature sets)
 * ========================================================================== */

#if KF_ENABLE_ADVANCED_API
/** @brief Innovation y = z - H x, from the most recent update (m values). */
kf_status_t kf_kf_get_innovation(const kf_kf_t *kf, kf_real_t *y_out);

/** @brief Kalman gain K (n x m, row-major) from the most recent update. */
kf_status_t kf_kf_get_gain(const kf_kf_t *kf, kf_real_t *K_out);
#endif /* KF_ENABLE_ADVANCED_API */

#if KF_ENABLE_DIAGNOSTICS
/**
 * @brief Non-invasive health check. Returns KF_OK if the state/covariance are
 * finite and P is symmetric; KF_ERROR_NON_FINITE or KF_ERROR_NUMERICAL
 * otherwise. Not required in the real-time path.
 */
kf_status_t kf_kf_check(const kf_kf_t *kf);
#endif /* KF_ENABLE_DIAGNOSTICS */

#endif /* KF_ENABLE_KF */

#ifdef __cplusplus
}
#endif

#endif /* KALMAN_KF_H */
