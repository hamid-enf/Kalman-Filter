/**
 * Example 10 — Basic UKF with a simple nonlinear model.
 *
 * Same problem as example 9 (estimate an angle from sin(theta) measurements),
 * but solved with the Unscented Kalman filter. The key difference: the UKF
 * does NOT need Jacobians. You only provide the raw nonlinear functions f and
 * h; the filter samples the distribution at sigma points and reconstructs the
 * mean/covariance from the propagated points.
 *
 * State       : x = theta                      (n = 1)
 * Measurement : z = sin(theta) + noise         (m = 1)
 * Transition  : f(x) = x  (angle constant)
 * Measurement : h(x) = sin(x)
 * Sigma params: alpha = 1, beta = 2, kappa = 0 (robust default)
 * Q / R / P0  : same as the EKF example
 *
 * Expected: same convergence to the true angle 0.8 rad, without any Jacobian
 * code. Use the UKF when the model is strongly nonlinear or hard to
 * differentiate; it costs 2n+1 evaluations of f and h per step.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>
#include <math.h>

/* --- nonlinear model callbacks (no Jacobians needed) ------------------- */

static void f_const(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                    kf_real_t *x_out, void *ctx)
{
    (void)u; (void)dt; (void)ctx;
    x_out[0] = x[0];
}

static void h_sin(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = (kf_real_t)sin((double)x[0]);
}

int main(void)
{
    kf_ukf_t ukf;
    kf_real_t x0[1] = {0.3f};
    kf_real_t truth = 0.8f;
    kf_real_t z;
    uint32_t rng = 10u;
    int i;

    kf_ukf_init(&ukf, 1, 1);
    kf_ukf_set_models(&ukf, f_const, h_sin);
    /* Optional: tune the sigma-point spread (defaults are alpha=1,beta=2,kappa=0). */
    kf_ukf_set_parameters(&ukf, 1.0f, 2.0f, 0.0f);
    kf_ukf_set_process_noise_scalar(&ukf, 1e-4f);
    kf_ukf_set_measurement_noise_scalar(&ukf, 0.01f);
    kf_ukf_set_covariance_scalar(&ukf, 0.5f);
    kf_ukf_set_state(&ukf, x0);

    printf("step   measurement   estimate\n");
    for (i = 0; i < 40; i++) {
        kf_ukf_predict(&ukf, NULL, 1.0f);
        z = (kf_real_t)sin((double)truth) + 0.1f * ex_noise_gauss(&rng);
        kf_ukf_update(&ukf, &z);

        if (i % 5 == 0) {
            printf("%4d   %11.3f   %8.3f\n", i, (double)z,
                   (double)kf_ukf_get_state(&ukf)[0]);
        }
    }
    printf("\nEstimated angle: %.3f rad (true %.3f rad)\n",
           (double)kf_ukf_get_state(&ukf)[0], (double)truth);
    return 0;
}
