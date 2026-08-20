/**
 * @file    common.h
 * @brief   Shared helpers for the stm32_h743 examples:
 *            - a deterministic pseudo-random generator (so every run prints
 *              the same numbers — great for learning),
 *            - uniform / gaussian-ish noise,
 *            - printf retargeting for the UART.
 *
 * The random generator plays the role of the real sensors. In each example a
 * clearly-labelled "SENSOR (input)" function returns a *noisy measurement of a
 * known true value*, so you can see exactly what goes in and how the filter
 * recovers the truth on the output.
 */

#ifndef H743_COMMON_H
#define H743_COMMON_H

/* ---- select the environment: real HAL on target, mock on host ----------- */
#ifdef H743_HOST
#include "host_hal.h"
#else
#include "main.h"   /* CubeMX: provides stm32h7xx_hal.h, huart3, HAL_GetTick */
#endif

#include "kalman.h"
#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------------
 * Deterministic PRNG (xorshift32). Same numbers on every run.
 * ------------------------------------------------------------------------ */
static inline uint32_t h7_rand(uint32_t *state)
{
    uint32_t x = *state;
    x ^= x << 13u;
    x ^= x >> 17u;
    x ^= x << 5u;
    *state = x;
    return x;
}

/* Uniform noise in [-1, 1]. */
static inline float h7_noise(uint32_t *state)
{
    return (float)((int32_t)(h7_rand(state) % 2000001u) - 1000000) / 1000000.0f;
}

/* Approximately Gaussian noise (sum of 3 uniforms), std ~ 1/sqrt(3) ~ 0.577. */
static inline float h7_noise_gauss(uint32_t *state)
{
    return (h7_noise(state) + h7_noise(state) + h7_noise(state)) / 3.0f;
}

/* --------------------------------------------------------------------------
 * printf retarget: route printf() to the UART.
 *
 * On the host this is already handled by stdio (host_hal.h), so this function
 * is only needed on real hardware. It is kept here, guarded, so the examples
 * read identically on both platforms.
 * ------------------------------------------------------------------------ */
#ifdef H743_HOST
#define H7_UART_PUTS(s) fputs((s), stdout)
#else
static inline void h7_uart_puts(const char *s)
{
    HAL_UART_Transmit(&huart3, (uint8_t *)(uintptr_t)s,
                      (uint16_t)strlen(s), 1000u);
}
#define H7_UART_PUTS(s) h7_uart_puts(s)
#endif

/*
 * On real hardware, ALSO retarget the C library printf via _write() in main.c:
 *
 *     int _write(int fd, char *ptr, int len) {
 *         HAL_UART_Transmit(&huart3, (uint8_t*)ptr, len, 1000);
 *         return len;
 *     }
 *
 * The examples call printf() directly; that _write() hook sends it to UART3.
 */

#ifdef __cplusplus
}
#endif

#endif /* H743_COMMON_H */
