# STM32 integration guide

The core library is hardware-independent: it never includes `stm32xxxx_hal.h`.
Integration is just adding source files and include paths (see
[Installation](installation.md)). This page covers the *application* patterns
specific to STM32.

## Where to place the instance

Use static storage so the instance is zero-initialised at boot and lives in
`.bss`:

```c
static kf_kf_t g_kf;   /* in .bss; call kf_kf_init() before first use */
```

Never place a large filter instance on a tiny stack (a default
`kf_ukf_t` is ~2 KB). The library itself uses no heap.

## Timing

`kf_kf_predict()` takes `dt` as a parameter; the library never assumes a fixed
sample rate. Obtain `dt` from the hardware:

```c
uint32_t now = HAL_GetTick();
kf_real_t dt = (kf_real_t)(now - last) / 1000.0f;
last = now;
```

For a fixed-rate control loop, run the filter from a timer callback with a
constant `dt` (see `examples/stm32/stm32g4_pos_vel.c`). For exact cycle counts
use the DWT counter (see [Benchmarks](benchmarks.md)).

## Execution contexts

The library performs no locking and creates no threads. The rules:

- **Single context** (all calls from the main loop, or all from one timer IRQ):
  no action needed — this is the common case.
- **Split predict/update across contexts** (e.g. `predict` in a timer ISR,
  `update` in the main loop): protect the calls with a critical section, since
  each call is a read-modify-write of the whole instance:

  ```c
  __disable_irq();
  kf_kf_predict(&g_kf, NULL, dt);
  __enable_irq();
  ```

  Or, with an RTOS, use a mutex around every call.
- **Multiple instances**: independent instances may be used from different
  contexts without any synchronisation.

## Typical patterns

### 1. Simple sensor filter (main loop, variable dt)

See `examples/stm32/stm32f4_1d_sensor.c`: read the ADC, compute `dt` from
`HAL_GetTick()`, `predict`, `update`, read the state.

### 2. Periodic estimation from a timer

See `examples/stm32/stm32g4_pos_vel.c`: `HAL_TIM_PeriodElapsedCallback`
performs the fixed-rate predict/update.

### 3. Sensor fusion at different rates

See `examples/stm32/stm32h7_fusion.c`: a fast loop predicts many times, and a
slow sensor updates occasionally — the classic IMU/GNSS split.

## ADC data

```c
HAL_ADC_Start(&hadc1);
if (HAL_ADC_PollForConversion(&hadc1, 10) == HAL_OK) {
    uint32_t raw = HAL_ADC_GetValue(&hadc1);
    kf_real_t volts = (kf_real_t)raw * (3.3f / 4095.0f);
    kf_kf_update(&g_kf, &volts);
}
```

If the ADC is DMA-driven, run the filter in the DMA-complete callback.

## Simulated sensor data (bring-up)

Before the real sensor is wired up, you can validate the filter with a
simulated signal (this is exactly what the host examples do):

```c
kf_real_t truth = 5.0f;
kf_real_t z = truth + noise();   /* your own noise generator */
```

## FPU notes

- STM32F4/G4: single-precision FPU — the default `float` build is the right
  choice.
- STM32H7: has a double-precision FPU; `float` is still faster and usually
  sufficient, but `KF_ENABLE_FLOAT64=1` is available if you need it.
- Enable the FPU and hard-float ABI in your toolchain settings; the library
  needs no special flags beyond C11.

## Enabling hardware floating point

Make sure your project uses the FPU (CubeIDE: `-mfloat-abi=hard -mfpu=fpv4-sp-d16`
for F4/G4, `-mfpu=fpv5-d16` for H7). Without it, float operations are
emulated in software and much slower.
