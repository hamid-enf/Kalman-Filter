/**
 * Example 6 — Two sensors with different noise levels (sequential fusion).
 *
 * State       : x = the (unknown) true value            (n = 1)
 * Measurements: two independent sensors of the SAME quantity (m = 1 each)
 * Model       : F = [1], H = [1] for both sensors
 * Q           : small drift
 * R1 (sensor A): 4.0  (coarse, noisy)
 * R2 (sensor B): 0.25 (precise)
 *
 * Instead of building one giant combined measurement vector, the two sensors
 * are applied as *sequential* updates: each kf_kf_update() call corrects the
 * state using one sensor, with that sensor's own R. Mathematically, for
 * uncorrelated sensors this is equivalent to a single combined update, but it
 * is simpler to set up and lets sensors arrive at different rates.
 *
 * Expected: the fused estimate is better (lower variance) than either sensor
 * alone, and naturally weights the precise sensor more heavily.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    kf_real_t F[1] = {1.0f};
    kf_real_t H[1] = {1.0f};
    kf_real_t x0[1] = {0.0f};
    kf_real_t truth = 10.0f;
    uint32_t rng = 6u;
    int i;

    kf_kf_init(&kf, 1, 1);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise_scalar(&kf, 1e-4f);
    kf_kf_set_covariance_scalar(&kf, 1.0f);
    kf_kf_set_state(&kf, x0);

    printf("step   sensor A   sensor B   fused\n");
    for (i = 0; i < 40; i++) {
        kf_real_t za, zb;
        kf_kf_predict(&kf, NULL, 1.0f);

        /* Coarse sensor (R = 4). */
        kf_kf_set_measurement_noise_scalar(&kf, 4.0f);
        za = truth + 2.0f * ex_noise_gauss(&rng);
        kf_kf_update(&kf, &za);

        /* Precise sensor (R = 0.25). */
        kf_kf_set_measurement_noise_scalar(&kf, 0.25f);
        zb = truth + 0.5f * ex_noise_gauss(&rng);
        kf_kf_update(&kf, &zb);

        if (i % 5 == 0) {
            printf("%4d   %8.2f   %8.2f   %6.2f\n", i, (double)za,
                   (double)zb, (double)kf_kf_get_state(&kf)[0]);
        }
    }
    printf("\nFinal fused estimate: %.2f (true %.2f)\n",
           (double)kf_kf_get_state(&kf)[0], (double)truth);
    return 0;
}
