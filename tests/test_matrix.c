/**
 * @file    test_matrix.c
 * @brief   Matrix-engine unit tests (correctness against hand-computed values,
 *          numerical stability of Cholesky/solve, error handling).
 */

#include "test_framework.h"
#include "kalman.h"

#include <stddef.h>

/* Convenience stack buffers sized for these tests. */
static kf_real_t bufA[3 * 3];
static kf_real_t bufB[3 * 3];
static kf_real_t bufC[3 * 3];

static void test_construction(void)
{
    kf_matrix_t m;

    t_begin("matrix: construction & element access");
    CHECK(kf_matrix_init(&m, 2, 2, bufA) == KF_OK);
    CHECK(kf_matrix_init(&m, 0, 2, bufA) == KF_ERROR_INVALID_DIMENSION);
#if KF_ENABLE_RUNTIME_CHECKS
    CHECK(kf_matrix_init(&m, 2, 2, NULL) == KF_ERROR_NULL_POINTER);
    CHECK(kf_matrix_init(NULL, 2, 2, bufA) == KF_ERROR_NULL_POINTER);
#endif

    kf_matrix_init(&m, 2, 2, bufA);
    kf_matrix_set(&m, 0, 1, 3.0);
    CHECK_NEAR(kf_matrix_get(&m, 0, 1), 3.0, 1e-6);
    CHECK(kf_matrix_zero(&m) == KF_OK);
    CHECK_NEAR(kf_matrix_get(&m, 0, 1), 0.0, 1e-9);
    t_end();
}

static void test_identity_add_sub(void)
{
    kf_matrix_t A, B, C;

    t_begin("matrix: identity / add / sub / scale");
    kf_matrix_init(&A, 3, 3, bufA);
    kf_matrix_init(&B, 3, 3, bufB);
    kf_matrix_init(&C, 3, 3, bufC);

    CHECK(kf_matrix_identity(&A) == KF_OK);
    CHECK_NEAR(kf_matrix_get(&A, 0, 0), 1.0, 1e-9);
    CHECK_NEAR(kf_matrix_get(&A, 0, 1), 0.0, 1e-9);
    CHECK_NEAR(kf_matrix_get(&A, 2, 2), 1.0, 1e-9);

    /* 2x identity matrix is not 3x3 -> dimension error for non-square ops. */
    {
        kf_matrix_t small;
        kf_real_t sb[4];
        kf_matrix_init(&small, 2, 2, sb);
        CHECK(kf_matrix_identity(&small) == KF_OK);
        CHECK(kf_matrix_add(&C, &A, &small) == KF_ERROR_INVALID_DIMENSION);
    }

    kf_matrix_identity(&B);
    kf_matrix_scale(&B, 2.0);
    CHECK_NEAR(kf_matrix_get(&B, 1, 1), 2.0, 1e-9);

    kf_matrix_add(&C, &A, &B);
    CHECK_NEAR(kf_matrix_get(&C, 1, 1), 3.0, 1e-9);

    kf_matrix_sub(&C, &B, &A);
    CHECK_NEAR(kf_matrix_get(&C, 1, 1), 1.0, 1e-9);
    CHECK_NEAR(kf_matrix_get(&C, 0, 1), 0.0, 1e-9);

    kf_matrix_add_scaled(&C, &A, &A, -1.0);
    CHECK_NEAR(kf_matrix_get(&C, 0, 0), 0.0, 1e-9);
    t_end();
}

static void test_multiply(void)
{
    kf_matrix_t A, B, C;

    t_begin("matrix: multiplication (incl. transpose helpers)");
    kf_matrix_init(&A, 2, 2, bufA);
    kf_matrix_init(&B, 2, 2, bufB);
    kf_matrix_init(&C, 2, 2, bufC);

    /* A = [[1,2],[3,4]], B = [[5,6],[7,8]]
       A*B = [[19,22],[43,50]] */
    kf_matrix_set(&A, 0, 0, 1); kf_matrix_set(&A, 0, 1, 2);
    kf_matrix_set(&A, 1, 0, 3); kf_matrix_set(&A, 1, 1, 4);
    kf_matrix_set(&B, 0, 0, 5); kf_matrix_set(&B, 0, 1, 6);
    kf_matrix_set(&B, 1, 0, 7); kf_matrix_set(&B, 1, 1, 8);

    CHECK(kf_matrix_mul(&C, &A, &B) == KF_OK);
    CHECK_NEAR(kf_matrix_get(&C, 0, 0), 19.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&C, 0, 1), 22.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&C, 1, 0), 43.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&C, 1, 1), 50.0, 1e-6);

    /* A * B^T = [[1*5+2*6, 1*7+2*8],[3*5+4*6, 3*7+4*8]]
              = [[17, 23],[39, 53]] */
    CHECK(kf_matrix_mul_transpose_b(&C, &A, &B) == KF_OK);
    CHECK_NEAR(kf_matrix_get(&C, 0, 0), 17.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&C, 0, 1), 23.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&C, 1, 0), 39.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&C, 1, 1), 53.0, 1e-6);

    /* A^T * B = [[1*5+3*7, 1*6+3*8],[2*5+4*7, 2*6+4*8]]
              = [[26, 30],[38, 44]] */
    CHECK(kf_matrix_mul_transpose_a(&C, &A, &B) == KF_OK);
    CHECK_NEAR(kf_matrix_get(&C, 0, 0), 26.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&C, 0, 1), 30.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&C, 1, 0), 38.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&C, 1, 1), 44.0, 1e-6);

    /* Dimension mismatch. */
    {
        kf_matrix_t rect;
        kf_real_t rb[6];
        kf_matrix_init(&rect, 2, 3, rb);
        CHECK(kf_matrix_mul(&rect, &A, &B) == KF_ERROR_INVALID_DIMENSION);
    }
    t_end();
}

static void test_transpose_symmetrize(void)
{
    kf_matrix_t A, C;

    t_begin("matrix: transpose & symmetrize");
    kf_matrix_init(&A, 2, 3, bufA);
    kf_matrix_init(&C, 3, 2, bufC);

    kf_matrix_set(&A, 0, 0, 1); kf_matrix_set(&A, 0, 1, 2); kf_matrix_set(&A, 0, 2, 3);
    kf_matrix_set(&A, 1, 0, 4); kf_matrix_set(&A, 1, 1, 5); kf_matrix_set(&A, 1, 2, 6);

    CHECK(kf_matrix_transpose(&C, &A) == KF_OK);
    CHECK_NEAR(kf_matrix_get(&C, 0, 0), 1.0, 1e-9);
    CHECK_NEAR(kf_matrix_get(&C, 1, 0), 2.0, 1e-9);
    CHECK_NEAR(kf_matrix_get(&C, 2, 1), 6.0, 1e-9);

    /* Symmetrize a non-symmetric 2x2. */
    {
        kf_matrix_t S;
        kf_real_t sb[4];
        kf_matrix_init(&S, 2, 2, sb);
        kf_matrix_set(&S, 0, 0, 1); kf_matrix_set(&S, 0, 1, 5);
        kf_matrix_set(&S, 1, 0, 3); kf_matrix_set(&S, 1, 1, 2);
        CHECK(kf_matrix_symmetrize(&S) == KF_OK);
        CHECK_NEAR(kf_matrix_get(&S, 0, 1), 4.0, 1e-9);
        CHECK_NEAR(kf_matrix_get(&S, 1, 0), 4.0, 1e-9);
    }
    t_end();
}

static void test_cholesky_and_solve(void)
{
    kf_matrix_t A, L, B;

    t_begin("matrix: Cholesky factorisation");
    /* A = [[4,2],[2,3]]  ->  L = [[2,0],[1,sqrt(2)]] */
    kf_matrix_init(&A, 2, 2, bufA);
    kf_matrix_set(&A, 0, 0, 4); kf_matrix_set(&A, 0, 1, 2);
    kf_matrix_set(&A, 1, 0, 2); kf_matrix_set(&A, 1, 1, 3);

    CHECK(kf_matrix_cholesky(&A) == KF_OK);
    CHECK_NEAR(kf_matrix_get(&A, 0, 0), 2.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&A, 1, 0), 1.0, 1e-6);
    CHECK_NEAR(kf_matrix_get(&A, 1, 1), 1.41421356, 1e-5);

    /* Solve A x = b with b = [6, 5] -> x = [1, 1]. */
    kf_matrix_init(&A, 2, 2, bufA);
    kf_matrix_set(&A, 0, 0, 4); kf_matrix_set(&A, 0, 1, 2);
    kf_matrix_set(&A, 1, 0, 2); kf_matrix_set(&A, 1, 1, 3);
    kf_matrix_cholesky(&A);

    kf_matrix_init(&B, 2, 1, bufB);
    kf_matrix_set(&B, 0, 0, 6);
    kf_matrix_set(&B, 1, 0, 5);
    CHECK(kf_matrix_cholesky_solve(&A, &B) == KF_OK);
    CHECK_NEAR(kf_matrix_get(&B, 0, 0), 1.0, 1e-5);
    CHECK_NEAR(kf_matrix_get(&B, 1, 0), 1.0, 1e-5);

    /* A singular matrix -> NOT positive definite. */
    kf_matrix_init(&L, 2, 2, bufC);
    kf_matrix_set(&L, 0, 0, 0); kf_matrix_set(&L, 0, 1, 0);
    kf_matrix_set(&L, 1, 0, 0); kf_matrix_set(&L, 1, 1, 0);
    CHECK(kf_matrix_cholesky(&L) == KF_ERROR_NOT_POSITIVE_DEFINITE);

    /* Indefinite matrix ([[1,2],[2,1]]) -> NOT positive definite. */
    kf_matrix_set(&L, 0, 0, 1); kf_matrix_set(&L, 0, 1, 2);
    kf_matrix_set(&L, 1, 0, 2); kf_matrix_set(&L, 1, 1, 1);
    CHECK(kf_matrix_cholesky(&L) == KF_ERROR_NOT_POSITIVE_DEFINITE);
    t_end();
}

static void test_vec_and_matvec(void)
{
    kf_matrix_t A;
    kf_real_t x[3], y[3];
    int i;

    t_begin("matrix: vector ops & matrix-vector multiply");
    for (i = 0; i < 3; i++) {
        x[i] = (kf_real_t)(i + 1);   /* [1,2,3] */
    }

    kf_vec_axpy(x, 2.0, x, 3);       /* x = 3*[1,2,3] = [3,6,9] */
    CHECK_NEAR(x[2], 9.0, 1e-9);
    CHECK_NEAR(kf_vec_dot(x, x, 3), 3 * 3 + 6 * 6 + 9 * 9, 1e-6);

    /* A = [[1,0,0],[0,1,0],[0,0,1]] * [1,2,3] = [1,2,3] */
    kf_matrix_init(&A, 3, 3, bufA);
    kf_matrix_identity(&A);
    x[0] = 1; x[1] = 2; x[2] = 3;
    CHECK(kf_mat_vec_mul(y, &A, x) == KF_OK);
    CHECK_NEAR(y[0], 1.0, 1e-9);
    CHECK_NEAR(y[1], 2.0, 1e-9);
    CHECK_NEAR(y[2], 3.0, 1e-9);
    t_end();
}

void test_matrix(void)
{
    test_construction();
    test_identity_add_sub();
    test_multiply();
    test_transpose_symmetrize();
    test_cholesky_and_solve();
    test_vec_and_matvec();
}
