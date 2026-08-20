/**
 * @file    stm32g4_pos_vel.c
 * @brief   STM32G4 example: position + velocity estimation from an encoder/
 *          ADC, executed periodically from a timer (or SysTick) callback.
 *
 * Target : any STM32G4 (Cortex-M4F, single-precision FPU) with HAL.
 * Model  : 2-state KF, x = [position, velocity], measure position only.
 *
 * The STM32G4's FPU makes the float path very fast, and its ADC/timers suit
 * motor-control style sampling. The filter runs at a fixed rate from a timer
 * callback, so dt is constant here (which also lets you use a constant F).
 */

#include "main.h"
#include "kalman.h"

static kf_kf_t g_kf;

/* ---------------------------------------------------------------------------
 * Initialise the constant-velocity filter (fixed sample period T_s).
 * ------------------------------------------------------------------------- */
static void filter_init(float Ts)
{
    kf_real_t F[4] = {1.0f, Ts, 0.0f, 1.0f};     /* dt = Ts                 */
    kf_real_t H[2] = {1.0f, 0.0f};               /* measure position only   */
    kf_real_t Q[4] = {0.0f, 0.0f, 0.0f, 0.1f};   /* velocity process noise  */
    kf_real_t P0[4] = {1.0f, 0.0f, 0.0f, 1.0f};

    kf_kf_init(&g_kf, 2, 1);
    kf_kf_set_transition_matrix(&g_kf, F);
    kf_kf_set_measurement_matrix(&g_kf, H);
    kf_kf_set_process_noise(&g_kf, Q);
    kf_kf_set_measurement_noise_scalar(&g_kf, 0.01f);
    kf_kf_set_covariance(&g_kf, P0);
}

/* ---------------------------------------------------------------------------
 * Read position from an encoder counter (TIMx) or an ADC; project-specific.
 * ------------------------------------------------------------------------- */
static float read_position(void)
{
    extern TIM_HandleTypeDef htim2;
    int32_t counts = (int32_t)__HAL_TIM_GET_COUNTER(&htim2);
    return (float)counts * 0.001f;               /* counts -> position      */
}

/*
 * Periodic filter step. In a real project call this from the timer IRQ
 * (HAL_TIM_PeriodElapsedCallback) or from the SysTick handler, keeping the
 * period fixed. Alternatively drive it from the main loop with HAL_GetTick()
 * as in stm32f4_1d_sensor.c if exact timing is not required.
 */
void HAL_TIM_PeriodElapsedCallback(TIM_HandleTypeDef *htim)
{
    const float Ts = 0.001f;                     /* 1 kHz control loop      */
    float z;
    float *x;

    if (htim->Instance != TIM2) {
        return;
    }

    z = read_position();
    kf_kf_predict(&g_kf, NULL, Ts);              /* x = F x                */
    kf_kf_update(&g_kf, &z);                     /* correct with position  */

    x = (float *)kf_kf_get_state(&g_kf);         /* x[0]=pos, x[1]=vel     */
    (void)x;                                     /* feed into your control */
}

/*
 * main(): configure clocks/timer/ADC, then call filter_init(Ts) and start the
 * timer. The actual filtering happens in the callback above.
 */
int main(void)
{
    HAL_Init();
    /* SystemClock_Config(); MX_GPIO_Init(); MX_TIM2_Init(); MX_ADC1_Init(); */

    filter_init(0.001f);

    /* Start the 1 kHz timer. */
    HAL_TIM_Base_Start_IT(&htim2);

    while (1) {
        /* Low-priority housekeeping only; filtering is ISR-driven. */
    }
}
