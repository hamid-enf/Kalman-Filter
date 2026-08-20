/**
 * Example 1 — 1D noisy signal filtering (the simplest possible use).
 *
 * State       : x = the (unknown) true signal value        (n = 1)
 * Measurement : z = x + noise                              (m = 1)
 * Model       : x(k+1) = x(k)   (constant signal)
 * Q           : 1e-3  (small process noise; lets the filter follow drift)
 * R           : 1.0   (measurement noise variance)
 *
 * The one-line constructor kf_kf_init_1d() configures the whole filter. No
 * matrix is built by hand.
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
    kf_real_t truth = 5.0f;
    kf_real_t z;
    uint32_t rng = 1u;
    int i;

    /* One call: 1 state, 1 measurement, process noise q, sensor noise r. */
    kf_kf_init_1d(&kf, /*q=*/1e-3f, /*r=*/1.0f);

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
