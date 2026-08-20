/**
 * Example 2 — Position estimation (1D random-walk position).
 *
 * State       : x = position                               (n = 1)
 * Measurement : z = position + noise                       (m = 1)
 * Model       : x(k+1) = x(k) + w   (position drifts)
 * Q           : 0.01   (random-walk intensity)
 * R           : 1.0    (sensor noise variance)
 * P0          : 10.0   (large initial uncertainty)
 *
 * Shows the convenience constructor followed by one tweak: the default
 * kf_kf_init_1d() sets P = 1, and we override the initial covariance to 10 to
 * reflect a genuinely unknown starting position.
 *
 * Expected: the estimate tracks the wandering true position, smoothing out the
 * sensor noise. A larger Q follows faster but noisier; a smaller Q smooths
 * more but lags. This is the classic tuning trade-off.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    kf_real_t truth = 0.0f;
    kf_real_t z;
    uint32_t rng = 2u;
    int i;

    kf_kf_init_1d(&kf, /*q=*/0.01f, /*r=*/1.0f);
    kf_kf_set_covariance_scalar(&kf, 10.0f);   /* unknown start position */

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
