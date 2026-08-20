# STM32H743IIT6 examples — «ملموس» (concrete, easy to follow)

Three examples for your **STM32H743IIT6** board, written so the input/output
of the filter is completely obvious. Every example uses **simulated sensor
data** (a deterministic random generator) instead of real wiring, so you can
run and understand it immediately — even with nothing connected but a UART.

| Example | Level | What you filter | Input (sensor) | Output (estimate) |
|---------|-------|-----------------|----------------|-------------------|
| `01_simple_temperature.c` | Simple | temperature | noisy reading + spikes | cleaned temperature |
| `02_position_velocity.c` | Medium | position + velocity | noisy position | position + *inferred* velocity |
| `03_sensor_fusion.c` | Advanced | position + velocity | odometry (fast) + GNSS (slow) | fused position + velocity |

All three print a readable table over UART:

```
 step |  raw sensor  |  filtered  | note
------+--------------+------------+--------------------
   20 |        44.98 |      25.56 | <-- spike rejected
```

## Your board (STM32H743IIT6)

- **Core**: Cortex-M7 @ up to 480 MHz.
- **FPU**: *double*-precision (FPv5-D16). The library's default `float` is the
  fastest choice, but `double` (`-DKF_ENABLE_FLOAT64=1`) is also cheap here.
- **Flash**: 2 MB, **RAM**: 1 MB — filter instances (a few KB) are negligible.
- Compile flags: `-mfloat-abi=hard -mfpu=fpv5-d16 -mcpu=cortex-m7`.

## Run on the host first (no hardware needed)

The examples build and run on your PC via a small HAL mock, so you can see the
exact output before flashing:

```sh
cd examples/stm32_h743
for ex in 01_simple_temperature 02_position_velocity 03_sensor_fusion; do
    gcc -std=c11 -O2 -Wall -DH743_HOST -I. -I../../kalman/include \
        $ex.c ../../kalman/src/kalman_*.c -lm -o $ex
    ./$ex
done
```

`-DH743_HOST` switches `common.h` to the host HAL mock (`host_hal.h`).

## Run on the STM32H743IIT6

1. **Add the library** (once per project):
   - Copy `kalman/` into your CubeIDE project (e.g. `Core/`).
   - Add `kalman/src/*.c` to the build and `kalman/include` to the include
     paths (see `docs/installation.md`).
2. **Copy one example** into `Core/Src/main.c` (or add the `.c` file and call
   its `demo_xxx()` function from your main loop).
3. **Configure UART3** (or any USART) in CubeMX at **115200 8N1**, TX on the
   USART TX pin.
4. **Retarget printf** in `main.c`:

   ```c
   int _write(int fd, char *ptr, int len) {
       HAL_UART_Transmit(&huart3, (uint8_t*)ptr, (uint16_t)len, 1000);
       return len;
   }
   ```

   (Also check `Use MicroLIB` or implement `_write`/`_read`/`__io_putchar` per
   your toolchain.)
5. **Boot order** in `main()` (CubeMX generates most of this):

   ```c
   HAL_Init();
   SystemClock_Config();
   MX_GPIO_Init();
   MX_USART3_UART_Init();

   demo_temperature();   /* the example's function */
   while (1) {}
   ```

6. Open a serial terminal (PuTTY/TeraTerm/`minicom`) at 115200 and watch the
   table print.

## How to read the examples (the important part)

Each file has the same four clearly-marked sections:

1. **Model** — what we are estimating and the matrices `F`/`H` (for simple
   cases the one-call constructors build them for you).
2. **INPUT** — a `read_xxx()` function returning a *noisy measurement of a
   known truth*. This is the only place you'd touch to plug in a real sensor
   (ADC, encoder, UART GNSS…).
3. **Filter loop** — `kf_kf_predict()` then `kf_kf_update()`.
4. **OUTPUT** — `kf_kf_get_state()` printed to the UART.

The `TRUE_xxx` constants are the values the filter is *trying to recover* — in
a real system you don't know them, but here they let you see convergence
happening step by step.

## Example 1 — Simple: temperature

- Truth is a constant 25.0 °C; the "sensor" adds noise (σ = 0.5) and a +20 °C
  spike every 15th sample.
- `kf_kf_init_1d()` builds the whole 1-D filter.
- After a short warm-up, `kf_kf_update_gated()` **rejects the spikes** (the
  innovation gate), so the filtered column stays smooth.

## Example 2 — Medium: position + velocity

- A cart moves at 1.5 m/s; we measure only *position* (σ = 0.3 m).
- `kf_kf_init_constant_velocity()` builds the 2-state model.
- Watch the **velocity column converge to 1.5 m/s even though no velocity is
  ever measured** — the filter infers it from the position history.

## Example 3 — Advanced: sensor fusion

- Two sensors, different rates and qualities:
  - **Odometry** (100 Hz): noisy *velocity* (σ = 0.2 m/s) — precise but drifts
    if integrated alone.
  - **GNSS** (1 Hz): noisy *position* (σ = 2 m) — noisy but drift-free.
- The filter predicts fast, corrects velocity with odometry, and corrects
  position with GNSS (two `H` matrices, two `R` values, sequential updates).
- The fused position stays smooth **and** drift-free.

## Measuring execution time (optional)

On the M7, time the loop with the DWT cycle counter:

```c
CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
DWT->CYCCNT = 0;
DWT->CTRL |= DWT_CTRL_CYCCTENA_Msk;

uint32_t t0 = DWT->CYCCNT;
kf_kf_predict(&kf, NULL, dt);
kf_kf_update(&kf, &z);
uint32_t cycles = DWT->CYCCNT - t0;
```

At 480 MHz, a small KF step costs a few thousand cycles — microseconds of CPU.

---

### راهنمای سریع (فارسی)

این سه مثال روی برد STM32H743IIT6 شما کار می‌کنند و چون از داده‌های **تصادفی
شبیه‌سازی‌شده** به جای سنسور واقعی استفاده می‌کنند، بدون هیچ سیم‌کشی هم قابل
اجرا و فهم‌اند. هر مثال چهار بخش دارد: مدل، ورودی (سنسور)، حلقه‌ی فیلتر، و
خروجی (چاپ روی UART). اول روی کامپیوتر با `-DH743_HOST` اجرا کنید تا خروجی را
ببینید، بعد روی برد فلش کنید (UART3 با 115200 و retarget کردن printf).
