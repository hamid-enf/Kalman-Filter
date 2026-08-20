/**
 * @file    kalman.h
 * @brief   Umbrella header for the Kalman filtering library.
 *
 * Include this single header in application code to pull in the whole library:
 *
 *      #include "kalman.h"
 *
 * The library provides three filters (each independently compile-time
 * switchable via kalman_config.h):
 *   - kf_kf_*  : linear Kalman filter
 *   - kf_ekf_* : extended Kalman filter
 *   - kf_ukf_* : unscented Kalman filter
 *
 * plus a small dependency-free matrix engine (kf_matrix_*) used internally and
 * exposed for advanced users.
 */

#ifndef KALMAN_H
#define KALMAN_H

#include "kalman_config.h"
#include "kalman_types.h"
#include "kalman_matrix.h"

#if KF_ENABLE_KF
#include "kalman_kf.h"
#endif

#if KF_ENABLE_EKF
#include "kalman_ekf.h"
#endif

#if KF_ENABLE_UKF
#include "kalman_ukf.h"
#endif

#endif /* KALMAN_H */
