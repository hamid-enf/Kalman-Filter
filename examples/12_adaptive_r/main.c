/**
 * Example 12 — Adaptive measurement noise (R).
 *
 * Suppose the sensor noise is not known in advance, or changes over time
 * (temperature drift, different operating conditions). Instead of hand-tuning
 * R, kf_kf_adapt_r() adjusts it online from the post-update residual
 * (residual-based covariance matching):
 *
 *      R <- gamma * R + (1 - gamma) * (r r^T + H P H^T)
 *
 * where r = z - H x^+ is the residual. Here the filter starts with an
 * over-optimistic R = 0.01 while the sensor actually has unit variance
 * (sigma = 1). The adaptation drives R upward toward the true value (1.0), so
 * the filter stops over-trusting the noisy sensor.
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
    for (i = 0; i < 200; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);
        /* ex_noise_gauss has std ~ 1/3, so scale by 3 -> unit variance. */
        z = truth + 3.0f * ex_noise_gauss(&rng);
        kf_kf_update(&kf, &z);
        kf_kf_adapt_r(&kf, /*gamma=*/0.98f, /*r_min=*/0.01f);

        if (i % 20 == 0) {
            printf("%4d   %6.3f   %8.2f\n", i,
                   (double)kf_kf_get_measurement_noise(&kf)[0],
                   (double)kf_kf_get_state(&kf)[0]);
        }
    }
    printf("\nAdapted R: %.2f (true sensor variance ~ 1.0)\n",
           (double)kf_kf_get_measurement_noise(&kf)[0]);
    return 0;
}
