/**
 * @file    03_sensor_fusion.c
 * @brief   ADVANCED example — fuse two sensors with different rates and
 *          qualities (the classic wheel-odometry + GNSS pattern).
 *
 * Two inputs, each measuring a DIFFERENT part of the state:
 *   - ODOMETRY (fast, 100 Hz): a noisy VELOCITY reading (wheel encoder).
 *     Precise short-term, but integrating it alone would drift in position.
 *   - GNSS (slow, 1 Hz): a noisy absolute POSITION (sigma 2 m).
 *     Noisy, but drift-free.
 *
 * The Kalman filter fuses them with a 2-state model x = [position, velocity]:
 *   - every 10 ms: predict, then correct VELOCITY with the odometry
 *     (H = [0, 1], small R).
 *   - every 1 s:   correct POSITION with the GNSS fix (H = [1, 0], large R).
 *
 * Result: a position estimate that is smooth (from odometry) and drift-free
 * (from GNSS). The library supports this naturally via sequential updates with
 * different measurement matrices and noise — no giant combined matrix needed.
 */

#include "common.h"

#define TRUE_VEL   2.0f                  /* m/s, true constant speed */
#define FAST_HZ    100u                  /* odometry rate */
#define FAST_DT    (1.0f / (float)FAST_HZ)
#define GNSS_EVERY 100u                  /* one GNSS fix per 100 fast steps */

/* ---------------------------------------------------------------------------
 * INPUT 1 — fast odometry: noisy velocity reading.
 * ------------------------------------------------------------------------- */
static float read_odometry(uint32_t *rng)
{
    return TRUE_VEL + 0.2f * h7_noise_gauss(rng);   /* sigma 0.2 m/s */
}

/* ---------------------------------------------------------------------------
 * INPUT 2 — slow GNSS: noisy absolute position.
 * ------------------------------------------------------------------------- */
static float read_gnss(uint32_t *rng, float true_pos)
{
    return true_pos + 2.0f * h7_noise_gauss(rng);   /* sigma 2 m */
}

/* ---------------------------------------------------------------------------
 * The filter demo.
 * ------------------------------------------------------------------------- */
static void demo_fusion(void)
{
    kf_kf_t kf;
    uint32_t rng = 0xAB0BA9Cu;
    float F[4] = {1.0f, FAST_DT, 0.0f, 1.0f};   /* pos += vel*dt; vel const */
    float H_vel[2] = {0.0f, 1.0f};              /* odometry measures velocity */
    float H_pos[2] = {1.0f, 0.0f};              /* GNSS measures position     */
    float Q[4] = {0.0f, 0.0f, 0.0f, 0.05f};     /* small velocity noise       */
    float P0[4] = {10.0f, 0.0f, 0.0f, 10.0f};
    float true_pos = 0.0f;
    uint32_t fast_count = 0u;
    uint32_t step;

    kf_kf_init(&kf, 2, 1);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_process_noise(&kf, Q);
    kf_kf_set_covariance(&kf, P0);

    printf("\r\n=== Example 3: sensor fusion (odometry 100 Hz + GNSS 1 Hz) ===\r\n");
    printf("  t(s) | odom vel | gnss pos | fused pos | fused vel\r\n");
    printf("-------+----------+----------+-----------+----------\r\n");

    for (step = 0; step < 400; step++) {
        float odom = read_odometry(&rng);             /* INPUT 1 (fast)  */
        float t;

        true_pos += TRUE_VEL * FAST_DT;

        /* 1) predict (integrate velocity into position). */
        kf_kf_predict(&kf, NULL, FAST_DT);

        /* 2) fast update: correct VELOCITY with the odometry. */
        kf_kf_set_measurement_matrix(&kf, H_vel);
        kf_kf_set_measurement_noise_scalar(&kf, 0.04f);   /* 0.2^2 */
        kf_kf_update(&kf, &odom);

        /* 3) slow update (once per second): correct POSITION with GNSS. */
        if (++fast_count >= GNSS_EVERY) {
            float gnss = read_gnss(&rng, true_pos);       /* INPUT 2 (slow) */
            fast_count = 0u;

            kf_kf_set_measurement_matrix(&kf, H_pos);
            kf_kf_set_measurement_noise_scalar(&kf, 4.0f); /* 2^2 */
            kf_kf_update(&kf, &gnss);

            t = (float)step * FAST_DT;
            printf(" %5.2f | %8.2f | %8.2f | %9.2f | %9.2f\r\n",
                   (double)t,
                   (double)odom,
                   (double)gnss,
                   (double)kf_kf_get_state(&kf)[0],
                   (double)kf_kf_get_state(&kf)[1]);
        }
    }
    printf("\r\nFinal: fused pos %.2f m (true %.2f m), fused vel %.2f m/s (true %.1f)\r\n",
           (double)kf_kf_get_state(&kf)[0], (double)true_pos,
           (double)kf_kf_get_state(&kf)[1], TRUE_VEL);
}

int main(void)
{
    /* On real hardware: HAL_Init(); SystemClock_Config(); MX_USART3_UART_Init(); */
    demo_fusion();
    return 0;
}
