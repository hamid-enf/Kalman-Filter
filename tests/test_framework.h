/**
 * @file    test_framework.h
 * @brief   Minimal assertion framework for the Kalman library test suite.
 *
 * Kept intentionally tiny: no external dependencies, no dynamic allocation.
 */

#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

/* Begin a named test case. */
void t_begin(const char *name);

/* Record a pass/fail for a boolean condition. */
void t_check(int cond, const char *expr, const char *file, int line);

/* Record a pass/fail for |a - b| <= tol (floating point comparison). */
void t_check_near(double a, double b, double tol,
                  const char *ae, const char *be, const char *file, int line);

/* Finish the current test case; returns the number of failures. */
int t_end(void);

/* Total failures across all test cases so far. */
int t_total(void);

#define CHECK(cond)             t_check((cond) ? 1 : 0, #cond, __FILE__, __LINE__)
#define CHECK_NEAR(a, b, tol)   t_check_near((double)(a), (double)(b), (tol), \
                                             #a, #b, __FILE__, __LINE__)
#define TEST(name)              void name(void)

#endif /* TEST_FRAMEWORK_H */
