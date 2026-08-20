/**
 * @file    test_framework.c
 * @brief   Minimal assertion framework implementation.
 */

#include "test_framework.h"

#include <stdio.h>
#include <math.h>

static int      g_checks   = 0;
static int      g_failures = 0;
static int      g_total    = 0;
static const char *g_name  = "";

void t_begin(const char *name)
{
    g_name = name;
    g_checks = 0;
    g_failures = 0;
    printf("  [TEST] %s\n", name);
}

void t_check(int cond, const char *expr, const char *file, int line)
{
    g_checks++;
    if (!cond) {
        g_failures++;
        printf("    FAIL %s:%d: %s\n", file, line, expr);
    }
}

void t_check_near(double a, double b, double tol,
                  const char *ae, const char *be, const char *file, int line)
{
    g_checks++;
    {
        double diff = fabs(a - b);
        if (diff > tol) {
            g_failures++;
            printf("    FAIL %s:%d: %s = %g vs %s = %g (diff %g, tol %g)\n",
                   file, line, ae, a, be, b, diff, tol);
        }
    }
}

int t_end(void)
{
    g_total += g_failures;
    printf("    %d/%d checks passed\n", g_checks - g_failures, g_checks);
    return g_failures;
}

int t_total(void)
{
    return g_total;
}
