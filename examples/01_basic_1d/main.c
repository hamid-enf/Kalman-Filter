/**
 * Example 1 — 1D noisy signal filtering (the simplest possible use).
 *
 * State       : x = the (unknown) true signal value        (n = 1)
 * Measurement : z = x + noise                              (m = 1)
 * Model       : x(k+1) = x(k)   (constant signal), F = [1], H = [1]
 * Q           : 1e-3  (small process noise; lets the filter follow drift)
 * R           : 1.0   (measurement noise variance)
 * P0          : 1.0   (initial uncertainty)
 *
 * Predict : x = x, P = P + Q
 * Update  : y = z - x; K = P/(P+R); x += K y; P = (1-K) P
 *
 * Expected: the estimate converges to the true value (5.0) and stays close
 * even though each individual measurement is noisy (sigma = 1.0).
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
    kf_real_t truth = 5.0f;
    kf_real_t z;
    uint32_t rng = 1u;
    int i;

    /* Initialise: 1 state, 1 measurement. */
    kf_kf_init(&kf, 1, 1);

    /* Configure the model and noise. */
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise_scalar(&kf, 1e-3f);
    kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
    kf_kf_set_covariance_scalar(&kf, 1.0f);
    kf_kf_set_state(&kf, x0);

    printf("step   measurement   estimate\n");
    for (i = 0; i < 40; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);

        z = truth + ex_noise_gauss(&rng);          /* simulate the sensor   */
        kf_kf_update(&kf, &z);

        if (i % 5 == 0) {
            printf("%4d   %10.3f   %9.3f\n", i, (double)z,
                   (double)kf_kf_get_state(&kf)[0]);
        }
    }

    printf("\nFinal estimate: %.3f (true value %.3f)\n",
           (double)kf_kf_get_state(&kf)[0], (double)truth);
    return 0;
}
