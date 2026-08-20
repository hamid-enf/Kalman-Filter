# Numerical stability

Kalman filters are prone to a small set of well-known numerical failure modes.
This library is designed to avoid them, and to *detect* them when they cannot
be avoided.

## The failure modes

1. **Covariance loses symmetry** — floating-point rounding makes `P` slightly
   asymmetric over many updates; an asymmetric "covariance" is meaningless and
   can drive divergence.
2. **Covariance loses positive semi-definiteness** — the simplified update
   `P = (I−KH)P` can produce a matrix with tiny negative eigenvalues under
   rounding, which later makes the Cholesky factorisation fail.
3. **Naive inversion of near-singular matrices** — computing `S⁻¹` explicitly
   is numerically much worse than solving `S X = B`.
4. **Singular / indefinite innovation covariance** — when `R` and `P` are both
   (near) zero, `S = HPHᵀ + R` can be singular; the standard formulas then
   produce NaNs.

## What the library does about them

### 1. Joseph-form covariance update (default)

Instead of `P = (I−KH)P`, the default is the **Joseph form**:

```
P = (I − KH) P (I − KH)ᵀ + K R Kᵀ
```

This form is symmetric and positive semi-definite *by construction*, so it is
far more robust to rounding. It costs a little extra work (a few more matrix
products) but is recommended for production. Disable it with
`KF_USE_JOSEPH_FORM=0` only for the absolute smallest/hottest path and only if
you have verified stability on your target.

### 2. Symmetrisation

After every covariance update the library symmetrises `P` (`P ← (P + Pᵀ)/2`)
to cancel the tiny asymmetry that accumulates. This is cheap and prevents
drift.

### 3. Cholesky solves instead of inversions

The Kalman gain is obtained by solving `S Kᵀ = (P Hᵀ)ᵀ` with a Cholesky
factorisation (`S = L Lᵀ`), never by forming `S⁻¹`. Cholesky is the standard
stable method for symmetric positive-definite systems. The UKF likewise builds
its sigma points from the Cholesky factor of `P` rather than an
eigendecomposition.

### 4. Explicit singular/indefinite detection

Cholesky returns `KF_ERROR_NOT_POSITIVE_DEFINITE` when a pivot is non-positive
(or below `KF_MIN_PIVOT`, default `1e-12`), and non-finite values are detected
when diagnostics are enabled. A failed factorisation is reported as an error
rather than silently producing NaN/Inf and corrupting the state.

## Practical guidance

- **Use `float` with care.** Single precision has ~7 significant digits. With
  the Joseph form and symmetrisation, single precision is adequate for most
  embedded filters (state dimension up to ~10). Switch to `double`
  (`KF_ENABLE_FLOAT64=1`) if you observe divergence, have an ill-conditioned
  system, or a very high-order model.
- **Keep `R` bounded away from zero** unless you truly have a noiseless sensor.
  A zero `R` combined with a zero initial covariance yields a singular `S`.
- **Prefer modest state dimensions.** High-order models are more prone to
  ill-conditioning; they rarely buy accuracy in practice.
- **Scale your state.** If your state mixes units of hugely different
  magnitude (e.g. nanovolts and megavolts), rescale so entries are comparable;
  this improves conditioning dramatically.
- **Use the diagnostics in development.** Compile with
  `KF_ENABLE_DIAGNOSTICS=1` and call `kf_kf_check()` (etc.) periodically; it
  flags non-finite values and asymmetric `P`. Turn it off for production.

## Known trade-offs

- The Joseph form is ~2–3× the floating-point work of the simplified form in
  the covariance step (the dominant cost is already `O(n³)` from `F P Fᵀ`, so
  the overall impact is modest). This is a deliberate robustness-over-speed
  choice, and it can be turned off if needed.
- Symmetrisation adds one `n×n` pass per update — negligible.
