# Performance

## Computational complexity

Let `n` = state dimension, `m` = measurement dimension.

| Operation | Complexity |
|-----------|------------|
| KF predict (`F P Fᵀ + Q`) | `O(n³)` (two `n×n×n` products) + `O(n²)` |
| KF update (Joseph form) | `O(n²m + nm² + m³)` (dominated by `S = HPHᵀ + R` and the Joseph products) |
| EKF predict / update | same as KF, plus the cost of the user's `f`/`h`/Jacobians |
| UKF predict | `O(n³)` (Cholesky) + `(2n+1)` calls to `f` + `O(n²(2n+1))` |
| UKF update | `O(m³)` (Cholesky) + `(2n+1)` calls to `h` + `O(nm(2n+1))` |
| Matrix multiply (general) | `O(p q r)` |
| Cholesky factorisation | `O(k³/3)` for a `k×k` matrix |
| Cholesky solve | `O(k²)` per right-hand side |

In practice the EKF and UKF are dominated by the *user's model callbacks*, so
keep those lightweight (no `printf`, no trig if avoidable) in the real-time
path.

## Measured results

Measured on the **host (x86)** with `-O2` as a *reference* only — **these are
not STM32 numbers**. They show relative scaling, which is what matters for
budgeting. Median ns per `predict`+`update` pair (see `benchmarks/benchmark.c`
for the full methodology):

| Filter | n | m | ns/pair |
|--------|---|---|---------|
| KF | 1 | 1 | ~106 |
| KF | 2 | 1 | ~161 |
| KF | 4 | 1 | ~382 |
| KF | 6 | 1 | ~945 |
| KF | 6 | 4 | ~1480 |
| EKF | 2 | 1 | ~159 |
| EKF | 4 | 1 | ~382 |
| EKF | 6 | 4 | ~1462 |
| UKF | 2 | 1 | ~258 |
| UKF | 4 | 1 | ~589 |
| UKF | 6 | 4 | ~1756 |

For STM32, measure cycles directly with the DWT counter (see
[Benchmarks](benchmarks.md)) — a Cortex-M4F typically executes a small KF step
in the low single-digit thousands of cycles.

## Memory footprint

All memory is embedded in the instance structs (plus whatever the caller
allocates for the instance itself). With `KF_MAX_STATE_DIM=6`,
`KF_MAX_MEASUREMENT_DIM=4`, `float`:

| Instance | Bytes (host build) | Composition |
|----------|--------------------|-------------|
| `kf_kf_t` | ~1.4 KB | `x` + `P,Q,F` (n² each) + `R,H` + scratch |
| `kf_ekf_t` | ~1.5 KB | KF storage + model callbacks + extra scratch |
| `kf_ukf_t` | ~2.2 KB | KF storage + sigma points `(2n+1)×(n + n + m)` + weights |

Scale linearly/quadratically with `n`,`m`. To shrink, lower
`KF_MAX_STATE_DIM`/`KF_MAX_MEASUREMENT_DIM` (see
[Configuration](configuration.md)). The memory for features you disable is not
compiled in at all: with `KF_ENABLE_EKF=0 -DKF_ENABLE_UKF=0` the EKF/UKF
objects are empty and only the KF + matrix code remains (~10 KB of flash
total, vs ~20 KB with all filters — measured on the host with `-O2`).

## What is deliberately *not* done for speed

- No `-ffast-math` assumptions; arithmetic follows IEEE 754 semantics.
- No hand-rolled SIMD/assembly. The row-major loops are written to
  auto-vectorise; if a target needs more, an optional CMSIS-DSP backend could
  be slotted in behind the matrix engine (the architecture allows this, but it
  is not required and not included).
- No unrolling of the generic loops (they are dimension-generic).

## Optimisation checklist for STM32

1. Build with `-O2` (or `-Os` if flash-limited) and enable the FPU.
2. Use `float` (the M4F/M7 single-precision FPU path is the fastest).
3. Keep `KF_MAX_*_DIM` at the exact sizes you use.
4. Turn off `KF_ENABLE_DIAGNOSTICS` and (once validated) the optional features.
5. Keep model callbacks minimal; hoist constant subexpressions out of the
   per-step path.
6. If the simplified covariance form is acceptable for your problem,
   `KF_USE_JOSEPH_FORM=0` trims a little work from the update.
