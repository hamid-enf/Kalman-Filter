/**
 * Example 3 — Position + velocity (constant-velocity model).
 *
 * State       : x = [position, velocity]                (n = 2)
 * Measurement : z = position + noise                    (m = 1)
 * Model       : x(k+1) = [pos + vel*dt; vel], dt = 1
 *               F = [[1, dt],[0, 1]]
 *               H = [1, 0]   (only position is measured)
 * Q           : process noise on the velocity (unknown accelerations)
 * R           : 1.0   (position sensor noise variance)
 * P0          : large diagonal
 *
 * Even though only position is measured, the filter *infers* the velocity from
 * the sequence of position measurements, and it can predict position between
 * measurements. This is the workhorse model for many tracking applications.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    kf_real_t F[4] = {1.0f, 1.0f, 0.0f, 1.0f};   /* dt = 1 */
    kf_real_t H[2] = {1.0f, 0.0f};
    kf_real_t Q[4] = {0.0f, 0.0f, 0.0f, 0.1f};   /* velocity noise only */
    kf_real_t P0[4] = {10.0f, 0.0f, 0.0f, 10.0f};
    kf_real_t x0[2] = {0.0f, 0.0f};
    kf_real_t truth_pos = 0.0f, truth_vel = 2.0f; /* true velocity = 2 */
    kf_real_t z;
    uint32_t rng = 3u;
    int i;

    kf_kf_init(&kf, 2, 1);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise(&kf, Q);
    kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
    kf_kf_set_covariance(&kf, P0);
    kf_kf_set_state(&kf, x0);

    printf("step   true pos   est pos   est vel\n");
    for (i = 0; i < 60; i++) {
        truth_pos += truth_vel;                       /* pos += vel * dt */
        kf_kf_predict(&kf, NULL, 1.0f);
        z = truth_pos + ex_noise_gauss(&rng);
        kf_kf_update(&kf, &z);

        if (i % 10 == 0) {
            printf("%4d   %8.2f   %7.2f   %7.2f\n", i, (double)truth_pos,
                   (double)kf_kf_get_state(&kf)[0],
                   (double)kf_kf_get_state(&kf)[1]);
        }
    }
    printf("\nEstimated velocity: %.2f (true %.2f)\n",
           (double)kf_kf_get_state(&kf)[1], (double)truth_vel);
    return 0;
}
