/* Cross-language calibration scenario — C reference implementation.
 * See README.md for the full spec. Build with the library:
 *   gcc -O2 fusion.c ../kalman/src/kalman_*.c -I../kalman/include -lm
 */

#include "kalman.h"
#include <stdio.h>
#include <stdint.h>

#define TRUE_VEL    2.0f
#define DT          0.01f
#define GNSS_EVERY  100u

static uint32_t rng = 0xC0FFEEu;

static uint32_t xorshift32(void)
{
    uint32_t x = rng;
    x ^= x << 13u;
    x ^= x >> 17u;
    x ^= x << 5u;
    rng = x;
    return x;
}

static float uniform(void)
{
    return (float)((int32_t)(xorshift32() % 2000001u) - 1000000) / 1000000.0f;
}

static float gauss(void)
{
    return (uniform() + uniform() + uniform()) / 3.0f;
}

int main(void)
{
    kf_kf_t kf;
    kf_real_t F[4] = {1.0f, DT, 0.0f, 1.0f};
    kf_real_t Hvel[2] = {0.0f, 1.0f};
    kf_real_t Hpos[2] = {1.0f, 0.0f};
    kf_real_t Q[4] = {0.0f, 0.0f, 0.0f, 0.05f};
    kf_real_t P0[4] = {10.0f, 0.0f, 0.0f, 10.0f};
    float true_pos = 0.0f;
    uint32_t fast = 0u;
    uint32_t step;

    kf_kf_init(&kf, 2u, 1u);
    kf_kf_set_transition_matrix(&kf, F);
    kf_kf_set_process_noise(&kf, Q);
    kf_kf_set_covariance(&kf, P0);

    printf("  t(s) | fused pos | fused vel\n");
    for (step = 0u; step < 400u; step++) {
        float odom = TRUE_VEL + 0.2f * gauss();
        true_pos += TRUE_VEL * DT;

        kf_kf_predict(&kf, NULL, DT);

        kf_kf_set_measurement_matrix(&kf, Hvel);
        kf_kf_set_measurement_noise_scalar(&kf, 0.04f);
        kf_kf_update(&kf, &odom);

        if (++fast >= GNSS_EVERY) {
            float gnss = true_pos + 2.0f * gauss();
            fast = 0u;
            kf_kf_set_measurement_matrix(&kf, Hpos);
            kf_kf_set_measurement_noise_scalar(&kf, 4.0f);
            kf_kf_update(&kf, &gnss);

            printf(" %5.2f | %9.5f | %9.5f\n",
                   (double)((float)step * DT),
                   (double)kf_kf_get_state(&kf)[0],
                   (double)kf_kf_get_state(&kf)[1]);
        }
    }
    return 0;
}
