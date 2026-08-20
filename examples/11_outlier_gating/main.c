/**
 * Example 11 — Outlier rejection with innovation gating (NIS).
 *
 * Real sensors occasionally emit spikes (a loose wire, a glitch). A plain
 * Kalman filter will happily follow a spike, corrupting the estimate. Gating
 * rejects measurements whose normalized innovation squared (NIS) exceeds a
 * chi-square threshold before they can corrupt the state.
 *
 * Here the true value is 5.0; every 10th measurement is corrupted to ~500.
 * With gating enabled the filter ignores those spikes and stays locked on 5.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

int main(void)
{
    kf_kf_t kf;
    kf_real_t truth = 5.0f;
    kf_real_t z;
    uint32_t rng = 11u;
    int i;

    kf_kf_init_1d(&kf, /*q=*/1e-3f, /*r=*/1.0f);

    /* Converge first WITHOUT gating (a gate rejects the far-off first reads). */
    for (i = 0; i < 30; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);
        z = truth + ex_noise_gauss(&rng);
        kf_kf_update(&kf, &z);
    }

    /* Now enable the chi-square gate: 6.63 ~ 99% for 1 degree of freedom. */
    kf_kf_set_gate_threshold(&kf, 6.63f);

    printf("step   measurement   NIS    estimate\n");
    for (i = 0; i < 60; i++) {
        kf_kf_predict(&kf, NULL, 1.0f);

        z = truth + ex_noise_gauss(&rng);
        if (i % 10 == 5) {
            z = 500.0f;                     /* a spike */
        }

        if (kf_kf_update_gated(&kf, &z) == KF_WARN_GATED) {
            printf("%4d   %10.1f   ----   %7.2f   (spike rejected)\n",
                   i, (double)z, (double)kf_kf_get_state(&kf)[0]);
        } else if (i % 10 == 0) {
            printf("%4d   %10.2f   %5.2f   %7.2f\n", i, (double)z,
                   (double)kf_kf_nis(&kf), (double)kf_kf_get_state(&kf)[0]);
        }
    }

    printf("\nFinal estimate: %.2f (true %.2f) -- spikes did not corrupt it\n",
           (double)kf_kf_get_state(&kf)[0], (double)truth);
    return 0;
}
