/**
 * Example 12 — Adaptive measurement noise (R).
 *
 * Suppose the sensor noise is not known in advance, or changes over time
 * (temperature drift, different operating conditions). Instead of hand-tuning
 * R, kf_kf_adapt_r() adjusts it online from the observed innovation:
 *
 *      R <- gamma * R + (1 - gamma) * (y y^T + H P H^T)
 *
 * Here the filter starts with an over-optimistic R = 0.01 while the sensor is
 * actually noisy (sigma ~ 1). The adaptation drives R upward toward the true
 * noise level, so the filter stops over-trusting the noisy sensor.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    kf_real_t truth = 20.0f;
    kf_real_t z;
    uint32_t rng = 12u;
    int i;

    kf_kf_init_1d(&kf, /*q=*/1e-4f, /*r=*/0.01f);   /* R starts too small */

    printf("step   est R     estimate\n");
    for (i = 0; i < 80; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);
        z = truth + 1.0f * ex_noise_gauss(&rng);    /* sensor noise ~ 1.0 */
        kf_kf_update(&kf, &z);
        kf_kf_adapt_r(&kf, /*gamma=*/0.95f, /*r_min=*/0.01f);

        if (i % 10 == 0) {
            printf("%4d   %6.3f   %8.2f\n", i,
                   (double)kf_kf_get_measurement_noise(&kf)[0],
                   (double)kf_kf_get_state(&kf)[0]);
        }
    }
    printf("\nAdapted R: %.3f (true sensor variance ~ 1.0)\n",
           (double)kf_kf_get_measurement_noise(&kf)[0]);
    return 0;
}
