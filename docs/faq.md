# FAQ

### Does the library use dynamic memory?
No. Every instance is a fixed-size struct with all storage embedded. There are
no calls to `malloc`/`calloc`/`realloc`/`free` anywhere in the code.

### Is it thread-safe / ISR-safe?
The library is reentrant and deterministic and creates no threads/locks. A
*single* instance must not be used concurrently from multiple contexts without
external synchronisation (each call is a read-modify-write of the instance).
Distinct instances are fully independent. See
[STM32 integration](stm32_integration.md).

### Can I use it without STM32 HAL?
Yes. The core has no HAL/RTOS/CMSIS dependency and compiles as plain C. The
examples and tests build and run on a host PC.

### Do I need CMSIS-DSP?
No. The library implements its own minimal matrix engine. The architecture
would allow an optional accelerated backend later, but none is required.

### Which filter should I use?
KF for linear models, EKF for mildly nonlinear models with computable
Jacobians, UKF for strongly nonlinear models or when Jacobians are
inconvenient. See [Math concepts](math_concepts.md).

### How do I choose Q, R, P?
See the [Tuning guide](tuning_guide.md). In short: `R` from sensor variance,
`P₀` from initial uncertainty, then tune `Q` for responsiveness.

### float or double?
`float` by default (fastest on STM32 FPUs, usually accurate enough given the
stable formulations). Enable `KF_ENABLE_FLOAT64=1` for `double`.

### Why is there no matrix inverse function?
Because the filters never need one — they solve `S X = B` with Cholesky, which
is more numerically stable. See [Architecture](architecture.md) and
[Numerical stability](numerical_stability.md).

### Does predict() need dt?
`predict` accepts `dt` for API symmetry with the EKF/UKF (whose model callbacks
receive it). For the linear KF the step size is encoded in `F`/`Q`; for
variable `dt`, re-set `F`/`Q` before each predict (example 08). It is safe to
pass any value.

### Can I use multiple sensors?
Yes — either as one vector measurement, or as sequential `update` calls (each
with its own `R`). See examples 05 and 06.

### What is the Joseph form, and can I turn it off?
It is the covariance update `P = (I−KH)P(I−KH)ᵀ + KRKᵀ`, symmetric and PSD by
construction. It is on by default; disable with `KF_USE_JOSEPH_FORM=0` if you
prefer the cheaper `P = (I−KH)P` and have verified stability.

### What C standard is required?
C11 or later. `kalman_config.h` asserts this.

### Is it MISRA-C compliant?
The code is written to MISRA-C:2012 and has been run through a static-analysis
pass (cppcheck + the official `misra.py` addon). Every Mandatory/Required rule
cppcheck can check is clean; the remaining findings are advisory and are
documented as explicit deviations (early-return style, precedence parentheses)
or false positives. See [MISRA compliance](misra_compliance.md) for the full
matrix and the reproducible check (`tools/misra_check.sh`). Formal
certification still requires a certified tool (LDRA, QAC, PC-lint, Parasoft).

### How do I report bugs or contribute?
Open an issue or pull request on the repository. Contributions should keep the
core dependency-free and allocation-free and extend the test suite.
