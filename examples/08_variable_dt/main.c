/**
 * Example 8 — Variable dt.
 *
 * In real systems the sampling interval is rarely perfectly constant (timer
 * jitter, varying loop duration). The library does not assume a fixed dt: the
 * step size is encoded in the transition matrix F, so the application simply
 * rebuilds F with the current dt before each predict().
 *
 * State       : x = [position, velocity]   (n = 2)
 * Measurement : z = position + noise       (m = 1)
 * Model       : F = [[1, dt],[0,1]] rebuilt each step, H = [1,0]
 *
 * Expected: despite irregular dt values the filter tracks the true trajectory
 * (position grows quadratically because velocity also grows linearly here).
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    kf_real_t H[2] = {1.0f, 0.0f};
    kf_real_t Q[4] = {0.0f, 0.0f, 0.0f, 0.05f};
    kf_real_t P0[4] = {10.0f, 0.0f, 0.0f, 10.0f};
    kf_real_t x0[2] = {0.0f, 0.0f};
    kf_real_t truth_pos = 0.0f, truth_vel = 1.0f;
    kf_real_t t = 0.0f, z;
    uint32_t rng = 8u;
    int i;

    kf_kf_init(&kf, 2, 1);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise(&kf, Q);
    kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
    kf_kf_set_covariance(&kf, P0);
    kf_kf_set_state(&kf, x0);

    printf("step    dt     true pos   est pos   est vel\n");
    for (i = 0; i < 30; i++) {
        kf_real_t dt = (i % 3 == 0) ? 0.5f : 1.0f;   /* irregular steps */
        kf_real_t F[4] = {1.0f, dt, 0.0f, 1.0f};

        truth_vel += 0.05f;                          /* slight acceleration */
        truth_pos += truth_vel * dt;
        t += dt;

        kf_kf_set_transition_matrix(&kf, F);         /* dt-dependent model */
        kf_kf_predict(&kf, NULL, dt);
        z = truth_pos + ex_noise_gauss(&rng);
        kf_kf_update(&kf, &z);

        printf("%4d   %.2f   %8.2f   %7.2f   %7.2f\n", i, (double)dt,
               (double)truth_pos, (double)kf_kf_get_state(&kf)[0],
               (double)kf_kf_get_state(&kf)[1]);
    }
    return 0;
}
