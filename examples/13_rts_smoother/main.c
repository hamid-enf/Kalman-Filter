/**
 * Example 13 — RTS fixed-interval smoother.
 *
 * A Kalman filter estimates each state using only past + current measurements.
 * A *smoother* revisits the whole trajectory using all measurements (past AND
 * future), producing a smoother, lower-variance estimate. It is an offline
 * (post-processing) step and needs no new measurements — just the stored
 * forward-pass data.
 *
 *   forward  : predict/update, store x_k, P_k, x_{k+1|k}, P_{k+1|k}, F_k
 *   backward : seed x_{s,N}=x_N, P_{s,N}=P_N, then smooth k = N-1 .. 0
 *
 * This example filters a random walk forward, then smooths it backward and
 * prints both trajectories alongside the truth.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>

#define N 40

int main(void)
{
    kf_kf_t kf;
    kf_real_t x_filt[N + 1], P_filt[N + 1];
    kf_real_t x_pred[N + 1], P_pred[N + 1];
    kf_real_t x_smooth[N + 1], P_smooth[N + 1];
    kf_real_t F[1] = {1.0f};
    kf_real_t truth[N + 1];
    uint32_t rng = 13u;
    int k;

    kf_kf_init_1d(&kf, /*q=*/0.01f, /*r=*/1.0f);

    /* Forward pass. */
    truth[0] = 0.0f;
    x_filt[0] = 0.0f;
    P_filt[0] = 1.0f;
    for (k = 0; k < N; k++) {
        kf_real_t z;
        truth[k + 1] = truth[k] + 0.5f * ex_noise_gauss(&rng);
        kf_kf_predict(&kf, NULL, 1.0f);
        x_pred[k + 1] = kf_kf_get_state(&kf)[0];
        P_pred[k + 1] = kf_kf_get_covariance(&kf)[0];
        z = truth[k + 1] + ex_noise_gauss(&rng);
        kf_kf_update(&kf, &z);
        x_filt[k + 1] = kf_kf_get_state(&kf)[0];
        P_filt[k + 1] = kf_kf_get_covariance(&kf)[0];
    }

    /* Backward smoothing pass. */
    x_smooth[N] = x_filt[N];
    P_smooth[N] = P_filt[N];
    for (k = N - 1; k >= 0; k--) {
        kf_kf_smooth_step(&kf,
                          &x_filt[k], &P_filt[k], F,
                          &x_pred[k + 1], &P_pred[k + 1],
                          &x_smooth[k + 1], &P_smooth[k + 1],
                          &x_smooth[k], &P_smooth[k]);
    }

    printf("step    truth   filtered   smoothed\n");
    for (k = 0; k <= N; k += 4) {
        printf("%4d   %6.2f   %8.2f   %8.2f\n", k, (double)truth[k],
               (double)x_filt[k], (double)x_smooth[k]);
    }
    return 0;
}
