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
    ((2u * KF_MAX_STATE_DIM * KF_MAX_STATE_DIM) \
   + (3u * KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM) \
   + (KF_MAX_MEASUREMENT_DIM * KF_MAX_MEASUREMENT_DIM) \
   + KF_MAX_STATE_DIM + (2u * KF_MAX_MEASUREMENT_DIM))

/* Scratch for the RTS smoother: 4 n x n temporaries + 1 n-vector. */
#define KF_KF_SMOOTHER_SCRATCH_FLOATS \
    ((4u * KF_MAX_STATE_DIM * KF_MAX_STATE_DIM) + KF_MAX_STATE_DIM)

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

#if KF_ENABLE_ADVANCED_API || KF_ENABLE_GATING || KF_ENABLE_ADAPTIVE_R
    kf_real_t y[KF_MAX_MEASUREMENT_DIM];  /**< innovation from last update    */
#endif

#if KF_ENABLE_ADVANCED_API
    kf_real_t K[KF_MAX_STATE_DIM * KF_MAX_MEASUREMENT_DIM];/* last gain n x m */
#endif

#if KF_ENABLE_GATING
    kf_real_t gate_threshold;             /**< chi-square gate; 0 = disabled  */
    kf_real_t last_nis;                   /**< last normalized innovation sq  */
#endif

#if KF_ENABLE_ADAPTIVE_R
    kf_real_t resid[KF_MAX_MEASUREMENT_DIM]; /**< residual r = z - H x^+ from
                                                  the last update               */
#endif

#if KF_ENABLE_SMOOTHER
    kf_real_t smoother_scratch[KF_KF_SMOOTHER_SCRATCH_FLOATS];
#endif

    kf_real_t scratch[KF_KF_SCRATCH_FLOATS];   /**< internal temporaries        */
} kf_kf_t;

/* ==========================================================================
 * Initialisation / reset
 * ==========================================================================
 *
 * kf_kf_init() is the general entry point. For the most common models, the
 * convenience constructors below configure the entire filter (F, H, Q, R, P)
 * in a single call, so a beginner never has to build matrices by hand.
 */

/**
 * @brief Convenience: 1-D constant-signal filter.
 *
 * Sets F = [1], H = [1], P = 1, Q = q, R = r and x = 0. Ideal for smoothing a
 * single noisy scalar (temperature, voltage, ...).
 *
 *      kf_kf_init_1d(&kf, q, r);
 */
kf_status_t kf_kf_init_1d(kf_kf_t *kf, kf_real_t q, kf_real_t r);

/**
 * @brief Convenience: N-D constant signal (F = I, H = I).
 *
 * Sets P = p0 * I, Q = q * I, R = r * I. Each state is measured directly and
 * drifts independently.
 */
kf_status_t kf_kf_init_constant(kf_kf_t *kf, uint16_t n,
                                kf_real_t q, kf_real_t r, kf_real_t p0);

/**
 * @brief Convenience: constant-velocity model (position + velocity).
 *
 *      x = [pos, vel],  F = [[1, dt],[0, 1]],  H = [1, 0]
 *      Q = q_accel * [[dt^4/4, dt^3/2],[dt^3/2, dt^2]]   (white acceleration)
 *      P = diag(p0_pos, p0_vel),  R = r
 *
 * @param q_accel  Variance of the driving white acceleration.
 */
kf_status_t kf_kf_init_constant_velocity(kf_kf_t *kf, kf_real_t dt,
                                         kf_real_t q_accel, kf_real_t r,
                                         kf_real_t p0_pos, kf_real_t p0_vel);

/**
 * @brief Convenience: constant-acceleration model (position + velocity +
 * acceleration).
 *
 *      x = [pos, vel, acc],  F has dt and dt^2/2 terms,  H = [1, 0, 0]
 *      Q = q_jerk * G G^T with G = [dt^3/6, dt^2/2, dt]^T   (white jerk)
 *      P = p0 * I,  R = r
 */
kf_status_t kf_kf_init_constant_acceleration(kf_kf_t *kf, kf_real_t dt,
                                             kf_real_t q_jerk, kf_real_t r,
                                             kf_real_t p0);

/* ==========================================================================
 * Reset
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

/** @brief Return a pointer to the process-noise covariance Q (n x n, row-major). */
const kf_real_t *kf_kf_get_process_noise(const kf_kf_t *kf);

/** @brief Set the full measurement-noise covariance R (m x m, row-major). */
kf_status_t kf_kf_set_measurement_noise(kf_kf_t *kf, const kf_real_t *R);

/** @brief Set R from a diagonal vector. */
kf_status_t kf_kf_set_measurement_noise_diagonal(kf_kf_t *kf, const kf_real_t *diag);

/** @brief Set R = r * I (scalar measurement noise). */
kf_status_t kf_kf_set_measurement_noise_scalar(kf_kf_t *kf, kf_real_t r);

/** @brief Return a pointer to the measurement-noise covariance R (m x m, row-major). */
const kf_real_t *kf_kf_get_measurement_noise(const kf_kf_t *kf);

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
 * Extensions (optional feature sets)
 * ========================================================================== */

#if KF_ENABLE_GATING
/**
 * @brief Set the innovation-gating threshold (chi-square value).
 *
 * With a non-zero threshold, kf_kf_update_gated() rejects measurements whose
 * normalized innovation squared (NIS) exceeds the threshold, protecting the
 * estimate from outliers (spikes). Typical 1-D chi-square values: 3.84 (95%),
 * 6.63 (99%), 10.83 (99.9%). Set 0 to disable gating.
 */
kf_status_t kf_kf_set_gate_threshold(kf_kf_t *kf, kf_real_t chi2);

/**
 * @brief Normalized innovation squared (NIS = y^T S^-1 y) from the most recent
 * update. Useful for monitoring filter health and for adaptive tuning.
 */
kf_real_t kf_kf_nis(const kf_kf_t *kf);

/**
 * @brief Measurement update with innovation gating.
 *
 * Computes the NIS first; if it exceeds the configured gate threshold the
 * measurement is treated as an outlier and the state/covariance are left
 * unchanged (returns KF_WARN_GATED). Otherwise behaves exactly like
 * kf_kf_update() and returns KF_OK.
 */
kf_status_t kf_kf_update_gated(kf_kf_t *kf, const kf_real_t *z);
#endif /* KF_ENABLE_GATING */

#if KF_ENABLE_ADAPTIVE_R
/**
 * @brief Adapt the measurement-noise covariance R online from the most recent
 * update (residual-based covariance matching).
 *
 *      R <- gamma * R + (1 - gamma) * (r r^T + H P H^T)
 *
 * where r = z - H x^+ is the post-update residual and P = P^+ the post-update
 * covariance. Because the residual has covariance R - H P H^T, the term
 * (r r^T + H P H^T) is an unbiased estimate of R and is positive
 * semi-definite by construction, so R stays well-conditioned.
 *
 * Call after an update. `gamma` in (0, 1) is the forgetting factor (closer to
 * 1 = slower adaptation); `r_min` is a floor applied to the diagonal to keep R
 * positive definite.
 */
kf_status_t kf_kf_adapt_r(kf_kf_t *kf, kf_real_t gamma, kf_real_t r_min);
#endif /* KF_ENABLE_ADAPTIVE_R */

#if KF_ENABLE_SMOOTHER
/**
 * @brief One backward step of the Rauch-Tung-Striebel (RTS) fixed-interval
 * smoother.
 *
 * Given the forward-pass filtered/predicted trajectory and the (already
 * smoothed) estimate at step k+1, computes the smoothed estimate at step k:
 *
 *      C   = P_k F^T (P_{k+1|k})^-1
 *      x_s = x_k + C (x_{s,k+1} - x_{k+1|k})
 *      P_s = P_k + C (P_{s,k+1} - P_{k+1|k}) C^T
 *
 * Usage: run the filter forward storing x_k, P_k, x_{k+1|k}, P_{k+1|k} and F_k
 * per step; then sweep backward from k = N-1 down to 0, seeding the sweep with
 * x_{s,N} = x_N, P_{s,N} = P_N. The result is a smoother, lower-variance
 * trajectory than the forward filter alone (see example 13).
 *
 * All matrices are n x n row-major; vectors are n elements.
 *
 * @param kf             Filter instance (provides n and scratch space).
 * @param x_filt         Forward filtered state at step k.
 * @param P_filt         Forward filtered covariance at step k.
 * @param F              State-transition matrix for step k -> k+1.
 * @param x_pred         Forward predicted state at step k+1 (x_{k+1|k}).
 * @param P_pred         Forward predicted covariance at step k+1.
 * @param x_smooth_next  Smoothed state at step k+1.
 * @param P_smooth_next  Smoothed covariance at step k+1.
 * @param x_smooth       [out] Smoothed state at step k.
 * @param P_smooth       [out] Smoothed covariance at step k.
 */
kf_status_t kf_kf_smooth_step(kf_kf_t *kf,
                              const kf_real_t *x_filt,
                              const kf_real_t *P_filt,
                              const kf_real_t *F,
                              const kf_real_t *x_pred,
                              const kf_real_t *P_pred,
                              const kf_real_t *x_smooth_next,
                              const kf_real_t *P_smooth_next,
                              kf_real_t *x_smooth,
                              kf_real_t *P_smooth);
#endif /* KF_ENABLE_SMOOTHER */

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
