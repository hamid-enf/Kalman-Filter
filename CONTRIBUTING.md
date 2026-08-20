# Contributing

Thanks for your interest! This library targets production STM32 firmware, so
the bar for changes is: **no hidden dependencies, no dynamic allocation, and
numerically-tested behaviour.**

## Ground rules (non-negotiable)

The core (`kalman/src/`, `kalman/include/`) must stay:

- **Dependency-free** — no HAL, CMSIS-DSP, FreeRTOS, or libc beyond `math.h`.
- **Allocation-free** — no `malloc`/`free`; all memory is static or
  caller-supplied.
- **C11, MISRA-C:2012-oriented** — explicit integer types, `const`-correct,
  no recursion, no hidden global mutable state, no magic numbers in the core.
- **Numerically robust** — prefer solves over inversions, keep the Joseph-form
  covariance update, symmetrise covariances.

If a change can't meet these, it belongs in an example or a language port, not
the C core.

## How to contribute

1. Fork the repo and create a branch.
2. Make your change and **add a test** in `tests/`.
3. Verify locally:

   ```sh
   make test                 # full C suite (default config)
   ./tools/sanitize.sh       # ASan + UBSan
   # plus the alternate configs:
   make CFLAGS="-std=c11 -O2 -DKF_ENABLE_FLOAT64=1" test
   make CFLAGS="-std=c11 -O2 -DKF_ENABLE_EKF=0 -DKF_ENABLE_UKF=0" test
   ```

   If you touched a language port, run its tests too (`python/tests`,
   `matlab/tests`, `csharp/Tests`).
4. Update the docs if the public API or configuration changed.
5. Open a pull request with a clear description.

## Test philosophy

See `docs/testing.md`. The suite is designed to catch *numerical* bugs, not
just compile errors: it includes exact comparisons against double-precision
reference implementations, property-based stress tests, boundary tests at the
maximum dimensions, and the UKF≡KF equivalence check. Please extend those
rather than only adding "smoke" tests.

## Style

- Follow the existing naming (`kf_<filter>_<verb>`, `KF_*` macros).
- Every public function returns `kf_status_t` (or is documented `void`).
- Document complexity/thread-safety/aliasing in the header, as the existing
  headers do.

## Reporting issues

Please include: the filter type, state/measurement dimensions, the model
matrices, `Q`/`R`/`P`, a minimal reproduction, and the observed vs expected
behaviour. A short C (or Python) reproducer is ideal.

## Code of conduct

Be kind and constructive. Keep discussion technical.
