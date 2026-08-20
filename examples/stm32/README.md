# STM32 integration examples

These files show how to integrate the library into a HAL-based STM32 project.
They are **reference skeletons**: they compile as part of a CubeMX / STM32CubeIDE
project that provides the STM32 HAL (`stm32xxxx_hal.h`, `main.h`, `HAL_*`
functions), so they are not built by the host `make` targets.

The library core is hardware-independent — it does **not** include any STM32
header. Integration is simply:

1. Add `kalman/src/*.c` to the project sources.
2. Add `kalman/include` to the include paths.
3. `#include "kalman.h"` in your application file.
4. Place the filter instance in static storage (never on a tiny stack, and no
   heap is used):

       static kf_kf_t g_kf;   /* zero-initialised; call kf_kf_init() at boot */

## Files

| File | Target | Demonstrates |
|------|--------|--------------|
| `stm32f4_1d_sensor.c`  | STM32F4 (Cortex-M4F) | 1D sensor filtering in the main loop, `HAL_GetTick()` for dt |
| `stm32g4_pos_vel.c`    | STM32G4 (Cortex-M4F) | position+velocity from an ADC reading, periodic (timer/SysTick) execution |
| `stm32h7_fusion.c`     | STM32H7 (Cortex-M7)  | two-sensor fusion at different rates (fast predict + slow update) |

## Timing on target

For cycle-accurate timing use the DWT counter (see `docs/benchmarks.md`):

```c
CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
DWT->CYCCNT = 0;
DWT->CTRL |= DWT_CTRL_CYCCTENA_Msk;
uint32_t t0 = DWT->CYCCNT;
kf_kf_predict(&kf, u, dt);
uint32_t cycles = DWT->CYCCNT - t0;
```

## Real-time notes

- Run predict() from the periodic timer/SysTick callback and update() when a
  sensor sample is ready, or run everything in the main loop (the examples do
  this for clarity).
- If a filter instance is touched from both an ISR and the main loop, protect
  the shared calls with a critical section (`__disable_irq()`/`__enable_irq()`)
  or an RTOS mutex — the library itself performs no locking.
