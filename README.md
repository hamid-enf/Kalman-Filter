# Kalman-Filter

[![CI](https://github.com/hamid-enf/Kalman-Filter/actions/workflows/ci.yml/badge.svg)](https://github.com/hamid-enf/Kalman-Filter/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Language: C11](https://img.shields.io/badge/language-C11-555555.svg)](kalman/)

A **dependency-free, allocation-free Kalman filtering library for STM32**
microcontrollers — written in C11, with a tiny high-level API and a
sophisticated, numerically robust core.

It provides three estimators, each independently switchable at compile time,
plus a set of practical extensions in the Kalman family:

| Component | When to use it |
|-----------|----------------|
| **KF** — Linear Kalman filter | Linear models (position, velocity, temperature, …) — the default choice |
| **EKF** — Extended Kalman filter | Mildly nonlinear models with computable Jacobians |
| **UKF** — Unscented Kalman filter | Strongly nonlinear models, or when Jacobians are hard/impossible |
| **Innovation gating (NIS)** | Rejecting sensor outliers/spikes |
| **Adaptive R** | Unknown or drifting sensor noise |
| **RTS smoother** | Offline post-processing of a recorded trajectory |

## Highlights

- **Zero dependencies** — no HAL, CMSIS-DSP, FreeRTOS, `malloc`, or `printf`.
- **Static memory only** — every filter instance is a plain C struct with a
  compile-time-known size; place it in `.bss` or the stack, never on the heap.
- **Deterministic** — no hidden global state, no recursion, reentrant.
- **Numerically robust** — Joseph-form covariance update, Cholesky-based
  solves instead of naive matrix inversions, covariance symmetrisation,
  optional NaN/Inf diagnostics.
- **Simple outside, sophisticated inside** — one-call constructors
  (`kf_kf_init_1d`, `kf_kf_init_constant_velocity`, …) for the common cases;
  `init → predict → update → get_state` is all a beginner needs, while experts
  get full matrix access.
- **MISRA-C:2012** — verified with cppcheck + the official `misra.py` addon:
  every Mandatory/Required rule cppcheck checks is clean; advisory findings are
  documented deviations (see `docs/misra_compliance.md`). C11, explicit integer
  types, `const`-correct.
- **Tested** — a host-runnable unit suite plus examples and a benchmark
  framework. Coverage includes: exact comparison against double-precision
  reference implementations, property-based stress tests (covariance stays
  symmetric/positive-definite/finite), boundary tests at maximum dimensions,
  the control-input path, and a matrix-engine fuzz; the suite also runs clean
  under AddressSanitizer + UndefinedBehaviorSanitizer.

## Quick start (1-D signal smoothing)

```c
#include "kalman.h"

kf_kf_t kf;
kf_kf_init_1d(&kf, /*q=*/1e-3f, /*r=*/1.0f);   /* whole model in one call */

for (;;) {
    kf_real_t z = read_sensor();   /* your noisy measurement         */

    kf_kf_predict(&kf, NULL, dt);  /* time update (dt varies freely) */
    kf_kf_update(&kf, &z);         /* measurement update             */

    kf_real_t x = kf_kf_get_state(&kf)[0];   /* filtered estimate    */
}
```

For position+velocity use `kf_kf_init_constant_velocity(&kf, dt, q, r, p0_pos,
p0_vel)`; for position+velocity+acceleration use
`kf_kf_init_constant_acceleration(...)`. These build `F`, `H`, `Q`, `R` and `P`
for you — no matrix algebra required. For anything custom, the granular setters
(`kf_kf_set_transition_matrix`, …) are still available.

The same pattern, with different prefixes, applies to `kf_ekf_*` and `kf_ukf_*`
(the EKF/UKF additionally register model callback functions — see the
examples).

## Building

```sh
make            # build library, examples, tests, benchmark
make test       # run the unit test suite
make run-examples
make bench
make sanitize   # run the suite under ASan + UBSan (memory safety / UB)
make check      # cppcheck + MISRA-C:2012 addon (if available)
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
| [MISRA compliance](docs/misra_compliance.md) | Rule-by-rule compliance matrix |

## Language ports

The same algorithms (identical formulas and numerical strategy) are also
available for offline prototyping and validation:

| Port | Directory | Notes |
|------|-----------|-------|
| **Python** (numpy) | [`python/`](python/) | fully tested here |
| **MATLAB** | [`matlab/`](matlab/) | reviewed, run `tests/test_kalman.m` to verify |
| **C# (.NET)** | [`csharp/`](csharp/) | reviewed, run `dotnet run --project Tests` to verify |

A [cross-language calibration](calibration/README.md) scenario implements the
same sensor-fusion filter in all four languages with an identical deterministic
RNG, so their numerical outputs can be diffed directly.

## License

MIT — see [LICENSE](LICENSE).
