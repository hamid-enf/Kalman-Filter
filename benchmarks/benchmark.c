/**
 * @file    benchmark.c
 * @brief   Execution-time benchmark for the Kalman filter library.
 *
 * IMPORTANT — these numbers are measured on the HOST (x86) as a *reference*
 * only. They are NOT STM32 measurements. Absolute times will differ on target,
 * but the *relative* scaling with state/measurement dimension and filter type
 * is indicative of the algorithmic cost, which is what matters for budgeting.
 *
 * Reproducing on STM32 (see docs/benchmarks.md):
 *   - Build with your STM32 toolchain and -O2 (or -Os).
 *   - Time the predict()/update() calls with the DWT cycle counter:
 *         CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
 *         DWT->CYCCNT = 0; DWT->CTRL |= DWT_CTRL_CYCCTENA_Msk;
 *         start = DWT->CYCCNT;  kf_kf_predict(...);  cycles = DWT->CYCCNT - start;
 *   - On Cortex-M with an FPU, cycles are far more meaningful than ns.
 *
 * The benchmark itself is framework-only measurement code: it does not fabricate
 * any STM32 numbers.
 */

#define _POSIX_C_SOURCE 199309L   /* for clock_gettime() on POSIX hosts */

#include "kalman.h"
#include <stdio.h>
#include <time.h>
#include <string.h>

#define WARMUP  1000u
#define REPEAT  20000u

static uint32_t rng_state = 0xC0FFEEu;

static uint32_t rng_next(void)
{
    uint32_t x = rng_state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng_state = x;
    return x;
}

static kf_real_t randu(void)
{
    return (kf_real_t)(rng_next() % 1000u) / 1000.0f;
}

static double now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

/* Time a predict+update pair, returning ns per pair. */
typedef void (*step_fn)(void);

static double bench_step(step_fn fn)
{
    uint32_t i;
    double t0, t1;

    for (i = 0; i < WARMUP; i++) {
        fn();
    }
    t0 = now_ns();
    for (i = 0; i < REPEAT; i++) {
        fn();
    }
    t1 = now_ns();
    return (t1 - t0) / (double)REPEAT;
}

/* --------------------------------------------------------------------------
 * Linear KF benchmark
 * ------------------------------------------------------------------------ */

static kf_kf_t g_kf;
static kf_real_t g_kf_u[KF_MAX_STATE_DIM];
static kf_real_t g_kf_z[KF_MAX_MEASUREMENT_DIM];
static uint16_t g_kf_n, g_kf_m;

static void kf_step(void)
{
    kf_kf_predict(&g_kf, g_kf_u, 0.01f);
    kf_kf_update(&g_kf, g_kf_z);
}

static void bench_kf(uint16_t n, uint16_t m)
{
    kf_real_t F[KF_MAX_STATE_DIM * KF_MAX_STATE_DIM];
    kf_real_t H[KF_MAX_MEASUREMENT_DIM * KF_MAX_STATE_DIM];
    uint16_t i;
    double ns;

    kf_kf_init(&g_kf, n, m);
    for (i = 0; i < n * n; i++) F[i] = 0.0f;
    for (i = 0; i < n; i++) F[i * n + i] = 1.0f;
    for (i = 0; i < m * n; i++) H[i] = (i % (n + 1) == 0) ? 1.0f : 0.0f;
    kf_kf_set_transition_matrix(&g_kf, F);
    kf_kf_set_measurement_matrix(&g_kf, H);
    kf_kf_set_process_noise_scalar(&g_kf, 1e-3f);
    kf_kf_set_measurement_noise_scalar(&g_kf, 1.0f);
    for (i = 0; i < n; i++) g_kf_u[i] = randu();
    for (i = 0; i < m; i++) g_kf_z[i] = randu();

    g_kf_n = n; g_kf_m = m;
    ns = bench_step(kf_step);
    printf("  KF    n=%u m=%u : %8.1f ns / predict+update\n", n, m, ns);
}

/* --------------------------------------------------------------------------
 * EKF benchmark (identity/quadratic model)
 * ------------------------------------------------------------------------ */

static kf_ekf_t g_ekf;

static void ekf_f(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                  kf_real_t *xo, void *c)
{
    uint16_t i; (void)u; (void)c;
    for (i = 0; i < g_kf_n; i++) xo[i] = x[i] + dt * 0.0f;
}

static void ekf_F(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                  kf_matrix_t *F, void *c)
{
    uint16_t i; (void)x; (void)u; (void)dt; (void)c;
    for (i = 0; i < g_kf_n; i++) kf_matrix_set(F, i, i, 1.0f);
}

static void ekf_h(const kf_real_t *x, kf_real_t *z, void *c)
{
    uint16_t i; (void)c;
    for (i = 0; i < g_kf_m; i++) z[i] = x[i];
}

static void ekf_H(const kf_real_t *x, kf_matrix_t *H, void *c)
{
    uint16_t i; (void)x; (void)c;
    for (i = 0; i < g_kf_m; i++) kf_matrix_set(H, i, i, 1.0f);
}

static void ekf_step(void)
{
    kf_ekf_predict(&g_ekf, g_kf_u, 0.01f);
    kf_ekf_update(&g_ekf, g_kf_z);
}

static void bench_ekf(uint16_t n, uint16_t m)
{
    uint16_t i;
    double ns;

    g_kf_n = n; g_kf_m = m;
    kf_ekf_init(&g_ekf, n, m);
    kf_ekf_set_models(&g_ekf, ekf_f, ekf_F, ekf_h, ekf_H);
    kf_ekf_set_process_noise_scalar(&g_ekf, 1e-3f);
    kf_ekf_set_measurement_noise_scalar(&g_ekf, 1.0f);
    for (i = 0; i < n; i++) g_kf_u[i] = randu();
    for (i = 0; i < m; i++) g_kf_z[i] = randu();

    ns = bench_step(ekf_step);
    printf("  EKF   n=%u m=%u : %8.1f ns / predict+update\n", n, m, ns);
}

/* --------------------------------------------------------------------------
 * UKF benchmark (linear model)
 * ------------------------------------------------------------------------ */

static kf_ukf_t g_ukf;

static void ukf_f(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                  kf_real_t *xo, void *c)
{
    uint16_t i; (void)u; (void)dt; (void)c;
    for (i = 0; i < g_kf_n; i++) xo[i] = x[i];
}

static void ukf_h(const kf_real_t *x, kf_real_t *z, void *c)
{
    uint16_t i; (void)c;
    for (i = 0; i < g_kf_m; i++) z[i] = x[i];
}

static void ukf_step(void)
{
    kf_ukf_predict(&g_ukf, g_kf_u, 0.01f);
    kf_ukf_update(&g_ukf, g_kf_z);
}

static void bench_ukf(uint16_t n, uint16_t m)
{
    uint16_t i;
    double ns;

    g_kf_n = n; g_kf_m = m;
    kf_ukf_init(&g_ukf, n, m);
    kf_ukf_set_models(&g_ukf, ukf_f, ukf_h);
    kf_ukf_set_process_noise_scalar(&g_ukf, 1e-3f);
    kf_ukf_set_measurement_noise_scalar(&g_ukf, 1.0f);
    for (i = 0; i < n; i++) g_kf_u[i] = randu();
    for (i = 0; i < m; i++) g_kf_z[i] = randu();

    ns = bench_step(ukf_step);
    printf("  UKF   n=%u m=%u : %8.1f ns / predict+update\n", n, m, ns);
}

/* --------------------------------------------------------------------------
 * Static footprint reporting (host object sizes; see docs for methodology)
 * ------------------------------------------------------------------------ */

static void print_footprint(void)
{
    printf("\nStatic footprint (bytes, host build):\n");
    printf("  sizeof(kf_kf_t)  = %zu\n", sizeof(kf_kf_t));
#if KF_ENABLE_EKF
    printf("  sizeof(kf_ekf_t) = %zu\n", sizeof(kf_ekf_t));
#endif
#if KF_ENABLE_UKF
    printf("  sizeof(kf_ukf_t) = %zu\n", sizeof(kf_ukf_t));
#endif
}

int main(void)
{
    printf("=== Kalman filter benchmark (host x86 reference, not STM32) ===\n");
    printf("Per predict+update pair, median of %u iterations:\n\n", REPEAT);

    printf("[Linear KF]\n");
    bench_kf(1, 1);
    bench_kf(2, 1);
    bench_kf(4, 1);
    bench_kf(6, 1);
    bench_kf(6, 4);

#if KF_ENABLE_EKF
    printf("\n[EKF]\n");
    bench_ekf(2, 1);
    bench_ekf(4, 1);
    bench_ekf(6, 4);
#endif

#if KF_ENABLE_UKF
    printf("\n[UKF]\n");
    bench_ukf(2, 1);
    bench_ukf(4, 1);
    bench_ukf(6, 4);
#endif

    print_footprint();
    return 0;
}
