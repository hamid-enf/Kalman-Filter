/**
 * Example 4 — Position + velocity + acceleration (constant-acceleration model).
 *
 * State       : x = [position, velocity, acceleration]   (n = 3)
 * Measurement : z = position + noise                     (m = 1)
 * Model       : constant acceleration (a_k+1 = a_k); F has dt and dt^2/2 terms.
 *
 * The constructor kf_kf_init_constant_acceleration() builds F, H, Q (white
 * jerk), R and P in one call. The filter recovers position, velocity AND
 * acceleration from position measurements alone, converging on the true
 * constant acceleration (1.0).
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    kf_real_t pos = 0.0f, vel = 1.0f, acc = 1.0f; /* true: v grows by 1/step */
    kf_real_t z;
    uint32_t rng = 4u;
    int i;

    kf_kf_init_constant_acceleration(&kf, /*dt=*/1.0f, /*q_jerk=*/0.01f,
                                     /*r=*/1.0f, /*p0=*/10.0f);

    printf("step   est pos   est vel   est acc\n");
    for (i = 0; i < 150; i++) {
        vel += acc;              /* v += a * dt        */
        pos += vel;              /* p += v * dt        */
        kf_kf_predict(&kf, NULL, 1.0f);
        z = pos + ex_noise_gauss(&rng);
        kf_kf_update(&kf, &z);

        if (i % 25 == 0) {
            printf("%4d   %7.2f   %7.2f   %7.2f\n", i,
                   (double)kf_kf_get_state(&kf)[0],
                   (double)kf_kf_get_state(&kf)[1],
                   (double)kf_kf_get_state(&kf)[2]);
        }
    }
    printf("\nEstimated acceleration: %.2f (true %.2f)\n",
           (double)kf_kf_get_state(&kf)[2], (double)acc);
    return 0;
}
