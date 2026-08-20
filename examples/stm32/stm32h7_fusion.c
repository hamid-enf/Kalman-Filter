/**
 * @file    stm32h7_fusion.c
 * @brief   STM32H7 example: sensor fusion at different rates (fast predict,
 *          slow update) — the classic IMU/GNSS pattern.
 *
 * Target : any STM32H7 (Cortex-M7, single- or double-precision FPU) with HAL.
 * Model  : 2-state KF, x = [position, velocity].
 *
 * The "fast" sensor (e.g. an IMU/odometry) runs at 100 Hz and drives
 * prediction; the "slow" sensor (e.g. GNSS/vision) arrives at 1 Hz and drives
 * an update. Because predict() and update() are separate calls, the filter
 * simply predicts many times between updates.
 */

#include "main.h"
#include "kalman.h"

static kf_kf_t g_kf;

#define FAST_HZ   100u
#define FAST_DT   (1.0f / (kf_real_t)FAST_HZ)
#define SLOW_EVERY 100u               /* one slow update per 100 fast steps */

static void filter_init(void)
{
    kf_real_t F[4] = {1.0f, FAST_DT, 0.0f, 1.0f};
    kf_real_t H[2] = {1.0f, 0.0f};
    kf_real_t Q[4] = {0.0f, 0.0f, 0.0f, 0.01f};
    kf_real_t P0[4] = {10.0f, 0.0f, 0.0f, 10.0f};

    kf_kf_init(&g_kf, 2, 1);
    kf_kf_set_transition_matrix(&g_kf, F);
    kf_kf_set_measurement_matrix(&g_kf, H);
    kf_kf_set_process_noise(&g_kf, Q);
    kf_kf_set_measurement_noise_scalar(&g_kf, 4.0f);  /* GNSS is coarse    */
    kf_kf_set_covariance(&g_kf, P0);
}

/* Fast odometry input: position increment since last step (velocity * dt). */
static kf_real_t read_odometry_increment(void)
{
    /* Return the control/odometry contribution in state space, or 0 if none. */
    return 0.0f;
}

/* Slow absolute-position sensor. */
static kf_real_t read_gnss(void)
{
    return 0.0f;   /* replace with the actual GNSS position reading */
}

int main(void)
{
    uint32_t fast_count = 0;

    HAL_Init();
    /* SystemClock_Config(); peripheral init ... */

    filter_init();

    while (1) {
        /* Fast loop: run at FAST_HZ (e.g. scheduled by a timer). */
        kf_real_t u[2];
        u[0] = read_odometry_increment();   /* optional control input      */
        u[1] = 0.0f;

        kf_kf_predict(&g_kf, u, FAST_DT);

        /* Slow sensor: correct every SLOW_EVERY fast steps. */
        if (++fast_count >= SLOW_EVERY) {
            kf_real_t z = read_gnss();
            fast_count = 0;
            kf_kf_update(&g_kf, &z);
        }

        {
            const kf_real_t *x = kf_kf_get_state(&g_kf);
            (void)x;   /* x[0] = fused position, x[1] = fused velocity     */
        }

        /* Wait for the next 100 Hz tick (e.g. HAL_Delay(10) for illustration,
           or better: a timer-driven event/flag). */
        HAL_Delay(10);
    }
}
