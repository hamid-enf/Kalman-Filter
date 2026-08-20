/**
 * @file    host_hal.h
 * @brief   Host-only mock of the STM32 HAL, so the stm32_h743 examples can be
 *          built and run on a PC (for learning/verification) without hardware.
 *
 * On real hardware this header is NOT used — the examples include <main.h>
 * (CubeMX) instead, which pulls in the genuine stm32h7xx_hal.h.
 *
 * Build the examples on the host with -DH743_HOST:
 *     gcc -DH743_HOST -I ../../kalman/include -I. 01_simple_temperature.c ../../kalman/src/kalman_*.c -lm
 */

#ifndef H743_HOST_HAL_H
#define H743_HOST_HAL_H

#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* Minimal HAL types used by the examples. */
typedef uint32_t HAL_StatusTypeDef;
#define HAL_OK ((HAL_StatusTypeDef)0u)

typedef struct { int unused; } UART_HandleTypeDef;

/* Fake "UART" handle used by the retargeted printf in the examples. */
extern UART_HandleTypeDef huart3;

/* --------------------------------------------------------------------------
 * HAL_GetTick(): on the host we advance 1 ms per call, which mimics the
 * SysTick-based millisecond counter on the real board well enough for the
 * examples (they only use it for dt / timing display).
 * ------------------------------------------------------------------------ */
static inline uint32_t HAL_GetTick(void)
{
    static uint32_t tick = 0u;
    tick += 1u;
    return tick;
}

/* printf retarget on the host simply writes to stdout. */
static inline HAL_StatusTypeDef HAL_UART_Transmit(UART_HandleTypeDef *huart,
                                                  uint8_t *pData,
                                                  uint16_t Size,
                                                  uint32_t Timeout)
{
    (void)huart;
    (void)Timeout;
    fwrite(pData, 1u, (size_t)Size, stdout);
    return HAL_OK;
}

static inline void HAL_Delay(uint32_t ms)
{
    (void)ms;
}

#endif /* H743_HOST_HAL_H */
