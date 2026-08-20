/**
 * Example 2 — Position estimation (1D random-walk position).
 *
 * State       : x = position                               (n = 1)
 * Measurement : z = position + noise                       (m = 1)
 * Model       : x(k+1) = x(k) + w   (position drifts), F = [1], H = [1]
 * Q           : 0.01   (random-walk intensity)
 * R           : 1.0    (sensor noise variance)
 * P0          : 10.0   (large initial uncertainty)
 *
 * Expected: the estimate tracks the wandering true position, smoothing out
 * the sensor noise. A larger Q makes it follow faster but noisier; a smaller
 * Q smooths more but lags real motion. This is the classic tuning trade-off.
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
    kf_real_t truth = 0.0f;
    kf_real_t z;
    uint32_t rng = 2u;
    int i;

    kf_kf_init(&kf, 1, 1);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise_scalar(&kf, 0.01f);
    kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
    kf_kf_set_covariance_scalar(&kf, 10.0f);
    kf_kf_set_state(&kf, x0);

    printf("step   true pos   measurement   estimate\n");
    for (i = 0; i < 40; i++) {
        truth += 0.5f * ex_noise_gauss(&rng);      /* position random-walks */
        kf_kf_predict(&kf, NULL, 1.0f);
        z = truth + ex_noise_gauss(&rng);          /* noisy measurement     */
        kf_kf_update(&kf, &z);

        if (i % 5 == 0) {
            printf("%4d   %8.3f   %12.3f   %8.3f\n", i, (double)truth,
                   (double)z, (double)kf_kf_get_state(&kf)[0]);
        }
    }
    return 0;
}
