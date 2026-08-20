/**
 * @file    02_position_velocity.c
 * @brief   MEDIUM example — estimate position AND velocity from position-only
 *          measurements.
 *
 * A cart moves at a constant speed of 1.5 m/s. We can only measure its
 * position (e.g. an encoder), and the measurement is noisy (sigma 0.3 m).
 *
 *   INPUT  : noisy position measurement (z = true position + noise).
 *   OUTPUT : filtered position AND inferred velocity (which is never measured
 *            directly — the filter derives it from the position history).
 *
 *   Model  : constant-velocity, x = [position, velocity], H = [1, 0].
 *
 * Watch the velocity column: it starts at 0 and converges to the true 1.5 m/s
 * even though the filter never sees a velocity reading.
 */

#include "common.h"

#define TRUE_VEL 1.5f                    /* m/s */
#define SAMPLE_HZ 100u                   /* 100 Hz measurement rate */

/* ---------------------------------------------------------------------------
 * INPUT — the (simulated) position sensor. Noisy reading of a moving cart.
 * ------------------------------------------------------------------------- */
static float read_position(uint32_t *rng, float true_pos)
{
    return true_pos + 0.3f * h7_noise_gauss(rng);      /* sigma 0.3 m */
}

/* ---------------------------------------------------------------------------
 * The filter demo.
 * ------------------------------------------------------------------------- */
static void demo_position_velocity(void)
{
    kf_kf_t kf;
    uint32_t rng = 0x5EED5EEDu;
    float dt = 1.0f / (float)SAMPLE_HZ;   /* 0.01 s */
    float true_pos = 0.0f;
    int step;

    /* Builds F, H, Q, R, P for the constant-velocity model in one call. */
    kf_kf_init_constant_velocity(&kf, dt,
                                 /*q_accel=*/0.1f,   /* unknown accelerations */
                                 /*r=*/0.09f,        /* 0.3^2                */
                                 /*p0_pos=*/1.0f, /*p0_vel=*/1.0f);

    printf("\r\n=== Example 2: position + velocity (true vel = %.1f m/s) ===\r\n", TRUE_VEL);
    printf("  t(s) | true pos | measured | filt pos | filt vel\r\n");
    printf("-------+----------+----------+----------+----------\r\n");

    for (step = 0; step < 100; step++) {
        float z;

        true_pos += TRUE_VEL * dt;                      /* cart keeps moving */
        z = read_position(&rng, true_pos);              /* INPUT             */

        kf_kf_predict(&kf, NULL, dt);
        kf_kf_update(&kf, &z);

        /* OUTPUT: print every 10th step (100 Hz is too fast to watch). */
        if (step % 10 == 0) {
            float t = (float)step * dt;
            printf(" %5.2f | %8.2f | %8.2f | %8.2f | %8.2f\r\n",
                   (double)t,
                   (double)true_pos,
                   (double)z,
                   (double)kf_kf_get_state(&kf)[0],
                   (double)kf_kf_get_state(&kf)[1]);
        }
    }
    printf("\r\nFinal: pos %.2f m (true %.2f m), vel %.2f m/s (true %.1f m/s)\r\n",
           (double)kf_kf_get_state(&kf)[0], (double)true_pos,
           (double)kf_kf_get_state(&kf)[1], TRUE_VEL);
}

int main(void)
{
    /* On real hardware: HAL_Init(); SystemClock_Config(); MX_USART3_UART_Init(); */
    demo_position_velocity();
    return 0;
}
