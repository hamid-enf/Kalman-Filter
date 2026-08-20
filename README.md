# Kalman-Filter

A **dependency-free, allocation-free Kalman filtering library for STM32**
microcontrollers — written in C11, with a tiny high-level API and a
sophisticated, numerically robust core.

It provides three estimators, each independently switchable at compile time:

| Filter | When to use it |
|--------|----------------|
| **KF** — Linear Kalman filter | Linear models (position, velocity, temperature, …) — the default choice |
| **EKF** — Extended Kalman filter | Mildly nonlinear models with computable Jacobians |
| **UKF** — Unscented Kalman filter | Strongly nonlinear models, or when Jacobians are hard/impossible |

## Highlights

- **Zero dependencies** — no HAL, CMSIS-DSP, FreeRTOS, `malloc`, or `printf`.
- **Static memory only** — every filter instance is a plain C struct with a
  compile-time-known size; place it in `.bss` or the stack, never on the heap.
- **Deterministic** — no hidden global state, no recursion, reentrant.
- **Numerically robust** — Joseph-form covariance update, Cholesky-based
  solves instead of naive matrix inversions, covariance symmetrisation,
  optional NaN/Inf diagnostics.
- **Simple outside, sophisticated inside** — a beginner needs only
  `init → predict → update → get_state`; experts get full matrix access.
- **MISRA-C:2012-oriented**, C11, explicit integer types, `const`-correct.
- **Tested** — a host-runnable unit suite (KF/EKF/UKF/matrix, numerical
  edge cases) plus examples and a benchmark framework.

## Quick start (1-D signal smoothing)

```c
#include "kalman.h"

kf_kf_t kf;                         /* 1 state, 1 measurement        */
kf_real_t F[1] = {1.0f};            /* x(k+1) = x(k)                 */
kf_real_t H[1] = {1.0f};            /* z(k)   = x(k) + noise         */

kf_kf_init(&kf, 1, 1);
kf_kf_set_transition_matrix(&kf, F);
kf_kf_set_measurement_matrix(&kf, H);
kf_kf_set_process_noise_scalar(&kf, 1e-3f);
kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
kf_kf_set_covariance_scalar(&kf, 1.0f);

for (;;) {
    kf_real_t z = read_sensor();   /* your noisy measurement         */

    kf_kf_predict(&kf, NULL, dt);  /* time update (dt varies freely) */
    kf_kf_update(&kf, &z);         /* measurement update             */

    kf_real_t x = kf_kf_get_state(&kf)[0];   /* filtered estimate    */
}
```

The same pattern, with different prefixes, applies to `kf_ekf_*` and `kf_ukf_*`
(the EKF/UKF additionally register model callback functions — see the
examples).

## Building

```sh
make            # build library, examples, tests, benchmark
make test       # run the unit test suite
make run-examples
make bench
```

For a **minimal** float-only linear-KF build:

```sh
make CFLAGS="-std=c11 -O2 -DKF_ENABLE_EKF=0 -DKF_ENABLE_UKF=0 \
             -DKF_ENABLE_DIAGNOSTICS=0 -DKF_ENABLE_ADVANCED_API=0" lib
```

A CMake build is also provided (`CMakeLists.txt`).

## Project layout

```
kalman/
├── include/           public headers (kalman.h is the umbrella)
│   ├── kalman_config.h   compile-time configuration
│   ├── kalman_types.h    status codes + scalar type
│   ├── kalman_matrix.h   internal matrix engine (also for advanced use)
│   ├── kalman_kf.h       linear KF
│   ├── kalman_ekf.h      EKF
│   └── kalman_ukf.h      UKF
└── src/               implementations (pure C, no HAL)
examples/              host examples (01..10) + stm32/ integration skeletons
tests/                 host unit-test suite
benchmarks/            execution-time benchmark framework
docs/                  full documentation (see below)
```

## Documentation

| Document | Contents |
|----------|----------|
| [Quick start](docs/quickstart.md) | Step-by-step for the three filters |
| [Installation](docs/installation.md) | Adding the library to a project |
| [Configuration](docs/configuration.md) | Every compile-time option |
| [API reference](docs/api_reference.md) | Complete function reference |
| [Architecture](docs/architecture.md) | Layers, memory model, design decisions |
| [Math concepts](docs/math_concepts.md) | State, P, Q, R, H, Kalman gain — in plain language |
| [Tuning guide](docs/tuning_guide.md) | How to choose Q, R, P |
| [Numerical stability](docs/numerical_stability.md) | What we do and why |
| [Performance](docs/performance.md) | Complexity and optimisation notes |
| [STM32 integration](docs/stm32_integration.md) | HAL wiring, ISR/RTOS, timing |
| [Troubleshooting](docs/troubleshooting.md) | Common problems and fixes |
| [FAQ](docs/faq.md) | Frequently asked questions |
| [Examples](docs/examples.md) | Walkthrough of examples 01–10 |
| [Testing](docs/testing.md) | What the test suite covers |
| [Benchmarks](docs/benchmarks.md) | Methodology and measured results |

## License

MIT — see [LICENSE](LICENSE).
