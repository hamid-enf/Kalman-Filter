/**
 * @file    01_simple_temperature.c
 * @brief   SIMPLE example — filter one noisy temperature reading.
 *
 * The most tangible possible case:
 *
 *   INPUT  : a temperature sensor that reads ~25.0 C but is noisy (sigma 0.5)
 *            and occasionally glitches (a spike of +20 C every 15th sample).
 *   OUTPUT : the filtered temperature, printed over UART every step.
 *
 *   Model  : 1-D constant signal, F = [1], H = [1].
 *
 * Because the "sensor" is simulated with a deterministic random generator, you
 * can see — step by step — the noisy reading go IN and the cleaned estimate
 * come OUT, and watch the spike get rejected by the innovation gate.
 *
 * On real hardware, replace read_temperature() with your ADC/one-wire read;
 * the filtering loop stays identical.
 */

#include "common.h"

/* The value we are trying to estimate (in a real system this is unknown). */
#define TRUE_TEMP 25.0f

/* ---------------------------------------------------------------------------
 * INPUT — the (simulated) sensor.
 * Returns the "raw" reading: true temperature + noise (+ a spike sometimes).
 * ------------------------------------------------------------------------- */
static float read_temperature(uint32_t *rng, int step)
{
    float reading = TRUE_TEMP + 0.5f * h7_noise_gauss(rng);   /* noisy sensor */

    if (step % 15 == 5) {
        reading += 20.0f;               /* occasional wiring glitch (spike) */
    }
    return reading;
}

/* ---------------------------------------------------------------------------
 * The filter demo. Call this from main() (or from the CubeMX main loop).
 * ------------------------------------------------------------------------- */
static void demo_temperature(void)
{
    kf_kf_t kf;
    uint32_t rng = 0x1A2B3C4Du;         /* fixed seed -> reproducible output */
    int step;

    /* One call sets up the whole 1-D filter (Q = slow drift, R = sensor var). */
    kf_kf_init_1d(&kf, /*q=*/1e-3f, /*r=*/0.25f);

    printf("\r\n=== Example 1: temperature filtering (true = %.1f C) ===\r\n", TRUE_TEMP);

    /* Warm-up: converge WITHOUT gating first. The estimate starts at 0 C, so
       the first readings (near 25 C) look like huge "outliers" — gating them
       from the start would reject everything. */
    for (step = 0; step < 20; step++) {
        float z = read_temperature(&rng, step);
        kf_kf_predict(&kf, NULL, 0.001f);
        kf_kf_update(&kf, &z);          /* plain update, no gate yet */
    }

    /* Now enable the gate: 6.63 ~ 99% chi-square for 1 degree of freedom. */
    kf_kf_set_gate_threshold(&kf, 6.63f);

    printf(" step |  raw sensor  |  filtered  | note\r\n");
    printf("------+--------------+------------+--------------------\r\n");

    for (; step < 45; step++) {
        float z = read_temperature(&rng, step);       /* INPUT  */
        kf_status_t st;

        kf_kf_predict(&kf, NULL, 0.001f);             /* time update (1 ms) */
        st = kf_kf_update_gated(&kf, &z);             /* measurement update */

        /* OUTPUT */
        if (st == KF_WARN_GATED) {
            printf(" %4d | %12.2f | %10.2f | <-- spike rejected\r\n",
                   step, (double)z, (double)kf_kf_get_state(&kf)[0]);
        } else {
            printf(" %4d | %12.2f | %10.2f |\r\n",
                   step, (double)z, (double)kf_kf_get_state(&kf)[0]);
        }
    }
    printf("\r\nFinal estimate: %.2f C  (true %.1f C)\r\n",
           (double)kf_kf_get_state(&kf)[0], TRUE_TEMP);
}

int main(void)
{
    /* On real hardware: HAL_Init(); SystemClock_Config(); MX_USART3_UART_Init(); */
    demo_temperature();
    return 0;
}
