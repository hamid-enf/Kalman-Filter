/**
 * Example 5 — Multiple measurements (vector measurement update).
 *
 * State       : x = [x_position, y_position]            (n = 2)
 * Measurement : z = [x_pos, y_pos] + noise              (m = 2)
 * Model       : F = I (positions drift), H = I (both axes measured together)
 * Q           : small diagonal random-walk noise
 * R           : diagonal measurement noise (different per axis)
 * P0          : diagonal
 *
 * Both axes are updated in a single kf_kf_update() call with a 2-element
 * measurement vector. Note R can encode different noise per axis.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    kf_real_t F[4] = {1.0f, 0.0f, 0.0f, 1.0f};
    kf_real_t H[4] = {1.0f, 0.0f, 0.0f, 1.0f};
    kf_real_t Q[4] = {0.01f, 0.0f, 0.0f, 0.01f};
    kf_real_t R[4] = {1.0f, 0.0f, 0.0f, 4.0f};   /* y-axis is noisier */
    kf_real_t P0[4] = {10.0f, 0.0f, 0.0f, 10.0f};
    kf_real_t x0[2] = {0.0f, 0.0f};
    kf_real_t z[2];
    uint32_t rng = 5u;
    int i;

    kf_kf_init(&kf, 2, 2);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise(&kf, Q);
    kf_kf_set_measurement_noise(&kf, R);
    kf_kf_set_covariance(&kf, P0);
    kf_kf_set_state(&kf, x0);

    printf("step   true(x,y)    estimate(x,y)\n");
    for (i = 0; i < 60; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);
        z[0] = 3.0f + ex_noise_gauss(&rng);        /* x target = 3  */
        z[1] = -2.0f + 2.0f * ex_noise_gauss(&rng); /* y target = -2, noisier */
        kf_kf_update(&kf, z);

        if (i % 10 == 0) {
            printf("%4d   (%4.1f,%4.1f)    (%5.2f,%5.2f)\n", i,
                   (double)z[0], (double)z[1],
                   (double)kf_kf_get_state(&kf)[0],
                   (double)kf_kf_get_state(&kf)[1]);
        }
    }
    return 0;
}
