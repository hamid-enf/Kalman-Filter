/**
 * @file    stm32f4_1d_sensor.c
 * @brief   STM32F4 example: filter a single noisy sensor in the main loop.
 *
 * Target : any STM32F4 (Cortex-M4F, single-precision FPU) with HAL.
 * Model  : 1-state KF, x = sensor value, F = [1], H = [1].
 *
 * Build  : drop this file + kalman/src/*.c into a CubeMX project, add
 *          kalman/include to the include path. Provide `main.h`, the HAL and a
 *          `read_sensor()` returning the raw ADC-converted value.
 *
 * This is a complete main-loop integration: HAL_GetTick() supplies the real
 * elapsed time for the (variable) dt.
 */

#include "main.h"        /* CubeMX: provides HAL handles, main() decl, etc. */
#include "kalman.h"

static kf_kf_t g_kf;                     /* static: no heap, no large stack   */

/* ---------------------------------------------------------------------------
 * Application-specific sensor read (replace with your ADC/SPI/I2C code).
 * ------------------------------------------------------------------------- */
static float read_sensor(void)
{
    /* Example: a 12-bit ADC channel scaled to a physical value. */
    extern ADC_HandleTypeDef hadc1;
    HAL_ADC_Start(&hadc1);
    if (HAL_ADC_PollForConversion(&hadc1, 10) == HAL_OK) {
        uint32_t raw = HAL_ADC_GetValue(&hadc1);      /* 0..4095 */
        return (float)raw * (3.3f / 4095.0f);         /* volts   */
    }
    return 0.0f;
}

/* ---------------------------------------------------------------------------
 * Filter initialisation (call once, e.g. at the top of main()).
 * ------------------------------------------------------------------------- */
static void filter_init(void)
{
    kf_real_t F[1] = {1.0f};
    kf_real_t H[1] = {1.0f};

    kf_kf_init(&g_kf, 1, 1);
    kf_kf_set_transition_matrix(&g_kf, F);
    kf_kf_set_measurement_matrix(&g_kf, H);
    kf_kf_set_process_noise_scalar(&g_kf, 1e-4f);   /* tune for your signal  */
    kf_kf_set_measurement_noise_scalar(&g_kf, 0.05f);/* sensor variance      */
    kf_kf_set_covariance_scalar(&g_kf, 1.0f);        /* initial uncertainty  */
}

/*
 * Example main() (in a real CubeMX project `main()` already exists in main.c;
 * call filter_init() and the loop body from there instead).
 */
int main(void)
{
    uint32_t last_tick = 0;
    float measurement;
    float estimate;

    HAL_Init();
    /* SystemClock_Config(); MX_GPIO_Init(); MX_ADC1_Init(); ... */

    filter_init();

    while (1) {
        uint32_t now = HAL_GetTick();
        kf_real_t dt = (kf_real_t)(now - last_tick) / 1000.0f;  /* seconds */
        last_tick = now;
        if (dt <= 0.0f) {
            dt = 0.001f;
        }

        measurement = read_sensor();

        kf_kf_predict(&g_kf, NULL, dt);      /* time update (dt varies)     */
        kf_kf_update(&g_kf, &measurement);   /* measurement update          */

        estimate = kf_kf_get_state(&g_kf)[0];
        (void)estimate;                       /* use estimate (e.g. control) */
    }
}
