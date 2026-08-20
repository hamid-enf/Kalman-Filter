/**
 * @file    util.h
 * @brief   Tiny shared helpers for the host examples: a deterministic PRNG
 *          (so example output is reproducible) and uniform/Gaussian-ish noise.
 *
 * The library core itself has no such helpers and no printf; these exist only
 * to make the examples self-contained and easy to read.
 */

#ifndef EXAMPLE_UTIL_H
#define EXAMPLE_UTIL_H

#include "kalman.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------------
 * xorshift32 PRNG — deterministic, so every example run prints the same result.
 * ------------------------------------------------------------------------ */
static inline uint32_t ex_rand(uint32_t *state)
{
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

/* Uniform noise in [-1, 1]. */
static inline kf_real_t ex_noise_unit(uint32_t *state)
{
    return (kf_real_t)((int32_t)(ex_rand(state) % 2000001u) - 1000000) / 1000000.0f;
}

/* Approximately Gaussian noise (sum of 3 uniforms -> bell-shaped). */
static inline kf_real_t ex_noise_gauss(uint32_t *state)
{
    kf_real_t s = (kf_real_t)0;
    int i;
    for (i = 0; i < 3; i++) {
        s += ex_noise_unit(state);
    }
    return s / (kf_real_t)3.0;
}

#ifdef __cplusplus
}
#endif

#endif /* EXAMPLE_UTIL_H */
