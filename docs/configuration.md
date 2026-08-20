# Configuration

All compile-time options live in `kalman/include/kalman_config.h`. You can edit
that file directly, or override any macro from the build system (e.g.
`-DKF_ENABLE_UKF=0`). Every non-essential feature is optional so a minimal
build contains only what you use.

## Scalar precision

| Macro | Effect |
|-------|--------|
| *(default)* | `kf_real_t` = `float` (single precision — fastest on M4F/G4) |
| `KF_ENABLE_FLOAT64=1` | `kf_real_t` = `double` |
| `KF_REAL_TYPE=<type>` | Fully custom scalar (highest precedence) |

## Dimensions

| Macro | Default | Meaning |
|-------|---------|---------|
| `KF_MAX_STATE_DIM` | 6 | Maximum state-vector dimension `n` |
| `KF_MAX_MEASUREMENT_DIM` | 4 | Maximum measurement dimension `m` |

These bound every internal buffer, so keep them as small as your application
allows. The per-instance RAM is documented in [Performance](performance.md).

## Filter types

| Macro | Default | Effect |
|-------|---------|--------|
| `KF_ENABLE_KF` | 1 | Linear Kalman filter |
| `KF_ENABLE_EKF` | 1 | Extended Kalman filter |
| `KF_ENABLE_UKF` | 1 | Unscented Kalman filter |

## Optional features

| Macro | Default | Effect |
|-------|---------|--------|
| `KF_ENABLE_DIAGNOSTICS` | 0 | NaN/Inf/symmetry checks (`kf_*_check()`) |
| `KF_ENABLE_VALIDATION` | 1 | Reject non-finite (NaN/Inf) inputs in the setters |
| `KF_ENABLE_RUNTIME_CHECKS` | 1 | NULL-pointer guards on every call |
| `KF_ENABLE_ADVANCED_API` | 1 | `kf_*_get_gain`, `get_innovation`, matrix access |
| `KF_ENABLE_GATING` | 1 | Innovation gating / NIS outlier rejection (linear KF) |
| `KF_ENABLE_ADAPTIVE_R` | 1 | Online measurement-noise (R) adaptation (linear KF) |
| `KF_ENABLE_SMOOTHER` | 1 | Rauch-Tung-Striebel fixed-interval smoother (linear KF) |

> Disable `KF_ENABLE_RUNTIME_CHECKS` only after you are confident all callers
> pass valid pointers — passing `NULL` with the checks off is undefined
> behaviour. Dimension checks are always performed (they are cheap and prevent
> buffer overruns).

## Numerical policy

| Macro | Default | Effect |
|-------|---------|--------|
| `KF_USE_JOSEPH_FORM` | 1 | Use the Joseph-form covariance update (`P = (I−KH)P(I−KH)ᵀ + KRKᵀ`), which guarantees symmetry/PSD-ness at a small extra cost |
| `KF_MIN_PIVOT` | `1e-12` | Cholesky pivots below this are treated as (near-)singular |

## Example configurations

**Minimal — 1-state KF, float, no extras**

```c
#define KF_MAX_STATE_DIM       1
#define KF_MAX_MEASUREMENT_DIM 1
#define KF_ENABLE_KF           1
#define KF_ENABLE_EKF          0
#define KF_ENABLE_UKF          0
#define KF_ENABLE_DIAGNOSTICS  0
#define KF_ENABLE_VALIDATION   0
#define KF_ENABLE_ADVANCED_API 0
```

**Medium — 4-state KF, float, diagnostics on**

```c
#define KF_MAX_STATE_DIM       4
#define KF_MAX_MEASUREMENT_DIM 2
#define KF_ENABLE_DIAGNOSTICS  1
```

**Advanced — EKF with a custom nonlinear model**

```c
#define KF_MAX_STATE_DIM       6
#define KF_MAX_MEASUREMENT_DIM 3
#define KF_ENABLE_KF           0   /* only need the EKF */
#define KF_ENABLE_EKF          1
#define KF_ENABLE_UKF          0
```

**UKF — configurable sigma points, float**

```c
#define KF_MAX_STATE_DIM       6
#define KF_MAX_MEASUREMENT_DIM 4
#define KF_ENABLE_KF           0
#define KF_ENABLE_EKF          0
#define KF_ENABLE_UKF          1
```

> `KF_ENABLE_FLOAT64` is independent: define it to switch the whole library to
> `double`. It is off by default because single precision is faster on most
> STM32 FPUs and is usually accurate enough with the library's stable
> formulations.
