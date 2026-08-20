# Benchmarks

## What is provided

`benchmarks/benchmark.c` measures the execution time of `predict`+`update`
pairs for several state/measurement sizes across the KF, EKF, and UKF, and
reports the static instance sizes. Build and run with:

```sh
make bench
```

> **Important**: the numbers printed are measured on the **host (x86)** and are
> a *reference for relative scaling only*. They are **not STM32 measurements**
> and are not claimed to be. Absolute times will differ on target. What carries
> over to STM32 is the *shape* — e.g. that the cost grows ~`O(n³)`, and that the
> UKF costs ~`(2n+1)` model evaluations per step.

## Methodology (host)

1. Warm up each measured path (`WARMUP` iterations) to avoid cold-cache bias.
2. Time `REPEAT` (20 000) iterations with `clock_gettime(CLOCK_MONOTONIC)`.
3. Report the per-pair median.

## Reproducing on STM32

The DWT cycle counter is the correct instrument on Cortex-M. In a HAL project:

```c
/* Enable the cycle counter once, at start-up. */
CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
DWT->CYCCNT = 0;
DWT->CTRL  |= DWT_CTRL_CYCCTENA_Msk;

uint32_t t0 = DWT->CYCCNT;
kf_kf_predict(&kf, NULL, dt);
kf_kf_update(&kf, &z);
uint32_t cycles = DWT->CYCCNT - t0;
```

- Run at a fixed CPU clock and disable interrupts around the measurement (or
  subtract the ISR overhead) for reproducibility.
- Average many runs; report *cycles*, which are architecture-meaningful, in
  addition to nanoseconds.
- Measure each configuration (1/2/4/6-state KF, EKF, UKF) exactly as the host
  benchmark does.

## Interpreting the numbers

- **KF/EKF** cost is dominated by `O(n³)` matrix products in `F P Fᵀ` and the
  `O(n²m + nm² + m³)` innovation path; the EKF adds the user's callbacks.
- **UKF** additionally evaluates `f`/`h` at `2n+1` points, so its cost scales
  with `n` even at small `m` — which is why the UKF is preferred mainly when
  its accuracy benefits outweigh the extra work.
- On an FPU-equipped Cortex-M, a small KF step is typically a few thousand
  cycles; the host `ns` figures are not directly comparable to cycles.

## Footprint methodology

Flash/ROM footprint is best measured per-target with `arm-none-eabi-size` on
the linked ELF. The host benchmark reports `sizeof(kf_*_t)` for RAM budgeting
and the build produces per-object sizes (`size build/kalman/src/*.o`) for a
coarse flash estimate. With `-O2`, the host object sizes were approximately:

| Configuration | Approx. code size |
|---------------|-------------------|
| All filters | ~20 KB |
| KF only (`KF_ENABLE_EKF=0 -DKF_ENABLE_UKF=0`) | ~10 KB |

Disabling EKF/UKF removes their objects entirely (they compile to zero bytes).
