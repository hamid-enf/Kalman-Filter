/**
 * @file    main.c
 * @brief   Test-suite entry point.
 */

#include "test_framework.h"
#include "kalman.h"

#include <stdio.h>

/* Test suites (each lives in its own translation unit). */
void test_matrix(void);
void test_kf(void);
void test_ekf(void);
void test_ukf(void);

int main(void)
{
    printf("=== Kalman filter test suite ===\n\n");

    test_matrix();
    test_kf();
    test_ekf();
    test_ukf();

    printf("\n=== %s (%d total failures) ===\n",
           t_total() == 0 ? "ALL TESTS PASSED" : "TESTS FAILED", t_total());
    return (t_total() == 0) ? 0 : 1;
}
