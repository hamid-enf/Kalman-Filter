/**
 * Example 9 — Basic EKF with a simple nonlinear model.
 *
 * Problem: estimate an angle theta from a sensor that measures sin(theta)
 * (e.g. a resolver or a pendulum). The measurement model is nonlinear:
 *
 *      h(theta) = sin(theta),   H = dh/dtheta = cos(theta)
 *
 * A linear KF cannot handle this because the relationship between the state
 * (angle) and the measurement (sine) is nonlinear. The EKF linearises it
 * around the current estimate using the Jacobian H.
 *
 * State       : x = theta                      (n = 1)
 * Measurement : z = sin(theta) + noise         (m = 1)
 * Transition  : f(x) = x  (angle is constant)  F = 1
 * Measurement : h(x) = sin(x)                  H = cos(x)
 * Q           : 1e-4 (small drift)
 * R           : 0.01 (sensor noise variance)
 * P0          : 0.5
 *
 * Expected: from an initial guess of 0.3 the estimate converges to the true
 * angle 0.8 rad, even though each measurement is a (noisy) sine, not an angle.
 */

#include "kalman.h"
#include "util.h"
#include <stdio.h>
#include <math.h>

/* --- nonlinear model callbacks ---------------------------------------- */

static void f_const(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                    kf_real_t *x_out, void *ctx)
{
    (void)u; (void)dt; (void)ctx;
    x_out[0] = x[0];                       /* angle is constant */
}

static void F_const(const kf_real_t *x, const kf_real_t *u, kf_real_t dt,
                    kf_matrix_t *F, void *ctx)
{
    (void)x; (void)u; (void)dt; (void)ctx;
    kf_matrix_set(F, 0, 0, 1.0f);         /* df/dtheta = 1 */
}

static void h_sin(const kf_real_t *x, kf_real_t *z, void *ctx)
{
    (void)ctx;
    z[0] = (kf_real_t)sin((double)x[0]);  /* h(theta) = sin(theta) */
}

static void H_sin(const kf_real_t *x, kf_matrix_t *H, void *ctx)
{
    (void)ctx;
    kf_matrix_set(H, 0, 0, (kf_real_t)cos((double)x[0])); /* dh/dtheta = cos */
}

int main(void)
{
    kf_ekf_t ekf;
    kf_real_t x0[1] = {0.3f};
    kf_real_t truth = 0.8f;               /* true angle (rad) */
    kf_real_t z;
    uint32_t rng = 9u;
    int i;

    kf_ekf_init(&ekf, 1, 1);
    kf_ekf_set_models(&ekf, f_const, F_const, h_sin, H_sin);
    kf_ekf_set_process_noise_scalar(&ekf, 1e-4f);
    kf_ekf_set_measurement_noise_scalar(&ekf, 0.01f);
    kf_ekf_set_covariance_scalar(&ekf, 0.5f);
    kf_ekf_set_state(&ekf, x0);

    printf("step   measurement   estimate\n");
    for (i = 0; i < 40; i++) {
        kf_ekf_predict(&ekf, NULL, 1.0f);
        z = (kf_real_t)sin((double)truth) + 0.1f * ex_noise_gauss(&rng);
        kf_ekf_update(&ekf, &z);

        if (i % 5 == 0) {
            printf("%4d   %11.3f   %8.3f\n", i, (double)z,
                   (double)kf_ekf_get_state(&ekf)[0]);
        }
    }
    printf("\nEstimated angle: %.3f rad (true %.3f rad)\n",
           (double)kf_ekf_get_state(&ekf)[0], (double)truth);
    return 0;
}
