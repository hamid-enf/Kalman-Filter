/**
 * Example 3 — Position + velocity (constant-velocity model).
 *
 * State       : x = [position, velocity]                (n = 2)
 * Measurement : z = position + noise                    (m = 1)
 * Model       : x(k+1) = [pos + vel*dt; vel], dt = 1, H = [1, 0]
 *
 * Even though only position is measured, the filter *infers* the velocity from
 * the sequence of position measurements, and it can predict position between
 * measurements. The constructor kf_kf_init_constant_velocity() builds F, H, Q,
 * R and P in one call.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    kf_real_t truth_pos = 0.0f, truth_vel = 2.0f; /* true velocity = 2 */
    kf_real_t z;
    uint32_t rng = 3u;
    int i;

    /* dt=1, acceleration noise q=0.1, sensor noise r=1, initial P diag(10,10). */
    kf_kf_init_constant_velocity(&kf, /*dt=*/1.0f, /*q_accel=*/0.1f,
                                 /*r=*/1.0f, /*p0_pos=*/10.0f, /*p0_vel=*/10.0f);

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
