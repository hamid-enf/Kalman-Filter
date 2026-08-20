/**
 * Example 4 — Position + velocity + acceleration (constant-acceleration model).
 *
 * State       : x = [position, velocity, acceleration]   (n = 3)
 * Measurement : z = position + noise                     (m = 1)
 * Model       : constant acceleration (a_k+1 = a_k):
 *               F = [[1, dt, dt^2/2],[0, 1, dt],[0, 0, 1]],  dt = 1
 *               H = [1, 0, 0]
 * Q           : small process noise on the acceleration (jerk)
 * R           : 1.0
 * P0          : large diagonal
 *
 * Expected: the filter recovers position, velocity AND acceleration from
 * position measurements alone, converging on the true constant acceleration
 * (1.0). This generalises the constant-velocity model of example 3.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    /* dt = 1 -> dt^2/2 = 0.5 */
    kf_real_t F[9] = {1.0f, 1.0f, 0.5f,
                      0.0f, 1.0f, 1.0f,
                      0.0f, 0.0f, 1.0f};
    kf_real_t H[3] = {1.0f, 0.0f, 0.0f};
    kf_real_t Q[9] = {0.0f, 0.0f, 0.0f,
                      0.0f, 0.0f, 0.0f,
                      0.0f, 0.0f, 0.01f};   /* jerk noise on acceleration */
    kf_real_t P0[9] = {10.0f, 0.0f, 0.0f,
                       0.0f, 10.0f, 0.0f,
                       0.0f, 0.0f, 10.0f};
    kf_real_t x0[3] = {0.0f, 0.0f, 0.0f};
    kf_real_t pos = 0.0f, vel = 1.0f, acc = 1.0f; /* true: v grows by 1/step */
    kf_real_t z;
    uint32_t rng = 4u;
    int i;

    kf_kf_init(&kf, 3, 1);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise(&kf, Q);
    kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
    kf_kf_set_covariance(&kf, P0);
    kf_kf_set_state(&kf, x0);

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
