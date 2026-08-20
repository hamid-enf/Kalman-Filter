/**
 * @file    kalman_common.c
 * @brief   Common helpers (status-code strings).
 */

#include "kalman_types.h"

const char *kf_status_str(kf_status_t status)
{
    switch (status) {
    case KF_OK:                        return "OK";
    case KF_WARN_GATED:                return "measurement gated (outlier)";
    case KF_ERROR_NULL_POINTER:        return "null pointer";
    case KF_ERROR_INVALID_DIMENSION:   return "invalid dimension";
    case KF_ERROR_INVALID_PARAMETER:   return "invalid parameter";
    case KF_ERROR_NOT_INITIALIZED:     return "not initialized";
    case KF_ERROR_FEATURE_DISABLED:    return "feature disabled";
    case KF_ERROR_SINGULAR_MATRIX:     return "singular matrix";
    case KF_ERROR_NOT_POSITIVE_DEFINITE: return "not positive definite";
    case KF_ERROR_NUMERICAL:           return "numerical failure";
    case KF_ERROR_NON_FINITE:          return "non-finite value";
    case KF_ERROR_DIVERGED:            return "filter diverged";
    default:                           return "unknown error";
    }
}
