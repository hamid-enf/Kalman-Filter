/**
 * Example 7 — Different predict and update rates.
 *
 * A very common embedded pattern: a fast sensor drives prediction (e.g. an IMU
 * at 100 Hz) while a slow sensor drives correction (e.g. GNSS at 1 Hz). The
 * library's separate predict()/update() API supports this directly: you simply
 * call predict() many times between update() calls.
 *
 * Here we simulate a position+velocity system where prediction runs 10x per
 * measurement.
 *
 * State       : x = [position, velocity]   (n = 2)
 * Measurement : z = position + noise       (m = 1), arrives every 10 predicts
 * Model       : F = [[1, dt],[0,1]] with dt = 0.1, H = [1,0]
 *
 * Expected: the filter keeps propagating the state (and growing P) between
 * measurements, then corrects sharply when each measurement arrives.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

#define PREDICT_PER_UPDATE 10

int main(void)
{
    kf_kf_t kf;
    kf_real_t dt = 0.1f;   /* fast-loop step */
    kf_real_t F[4] = {1.0f, dt, 0.0f, 1.0f};
    kf_real_t H[2] = {1.0f, 0.0f};
    kf_real_t Q[4] = {0.0f, 0.0f, 0.0f, 0.001f};
    kf_real_t P0[4] = {1.0f, 0.0f, 0.0f, 1.0f};
    kf_real_t x0[2] = {0.0f, 0.0f};
    kf_real_t truth_pos = 0.0f, truth_vel = 2.0f;
    kf_real_t z;
    uint32_t rng = 7u;
    int update = 0, i;

    kf_kf_init(&kf, 2, 1);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_measurement_matrix(&kf, H);
    kf_kf_set_process_noise(&kf, Q);
    kf_kf_set_measurement_noise_scalar(&kf, 1.0f);
    kf_kf_set_covariance(&kf, P0);
    kf_kf_set_state(&kf, x0);

    printf("update   true pos   est pos   est vel\n");
    for (i = 0; i < 200; i++) {
        truth_pos += truth_vel * dt;             /* fast true dynamics */
        kf_kf_predict(&kf, NULL, dt);            /* fast prediction     */

        if ((i + 1) % PREDICT_PER_UPDATE == 0) { /* slow measurement    */
            z = truth_pos + ex_noise_gauss(&rng);
            kf_kf_update(&kf, &z);
            update++;
            if (update % 2 == 0) {
                printf("%6d   %8.2f   %7.2f   %7.2f\n", update,
                       (double)truth_pos,
                       (double)kf_kf_get_state(&kf)[0],
                       (double)kf_kf_get_state(&kf)[1]);
            }
        }
    }
    return 0;
}
