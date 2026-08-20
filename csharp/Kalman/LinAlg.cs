using System;

namespace Kalman
{
    /// <summary>
    /// Minimal dense linear-algebra helpers used by the KF/EKF/UKF ports.
    /// Matrices are <c>double[,]</c> (row, column). This mirrors the C matrix
    /// engine: Cholesky factorisation + solve instead of explicit inversions.
    /// </summary>
    public static class LinAlg
    {
        public static double[,] Identity(int n)
        {
            var a = new double[n, n];
            for (int i = 0; i < n; i++) a[i, i] = 1.0;
            return a;
        }

        public static double[,] Clone(double[,] a)
        {
            return (double[,])a.Clone();
        }

        public static double[,] Transpose(double[,] a)
        {
            int r = a.GetLength(0), c = a.GetLength(1);
            var t = new double[c, r];
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    t[j, i] = a[i, j];
            return t;
        }

        public static double[,] Multiply(double[,] a, double[,] b)
        {
            int p = a.GetLength(0), q = a.GetLength(1), r = b.GetLength(1);
            if (q != b.GetLength(0)) throw new ArgumentException("inner dimensions mismatch");
            var c = new double[p, r];
            for (int i = 0; i < p; i++)
                for (int j = 0; j < r; j++)
                {
                    double s = 0.0;
                    for (int k = 0; k < q; k++) s += a[i, k] * b[k, j];
                    c[i, j] = s;
                }
            return c;
        }

        public static double[,] Add(double[,] a, double[,] b)
        {
            int r = a.GetLength(0), c = a.GetLength(1);
            var o = new double[r, c];
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++) o[i, j] = a[i, j] + b[i, j];
            return o;
        }

        public static double[,] Subtract(double[,] a, double[,] b)
        {
            int r = a.GetLength(0), c = a.GetLength(1);
            var o = new double[r, c];
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++) o[i, j] = a[i, j] - b[i, j];
            return o;
        }

        public static double[,] Scale(double[,] a, double s)
        {
            int r = a.GetLength(0), c = a.GetLength(1);
            var o = new double[r, c];
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++) o[i, j] = s * a[i, j];
            return o;
        }

        public static void SymmetrizeInPlace(double[,] a)
        {
            int n = a.GetLength(0);
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                {
                    double v = 0.5 * (a[i, j] + a[j, i]);
                    a[i, j] = v;
                    a[j, i] = v;
                }
        }

        public static double[] MatVec(double[,] a, double[] x)
        {
            int r = a.GetLength(0), c = a.GetLength(1);
            var y = new double[r];
            for (int i = 0; i < r; i++)
            {
                double s = 0.0;
                for (int k = 0; k < c; k++) s += a[i, k] * x[k];
                y[i] = s;
            }
            return y;
        }

        public static double[] VecSub(double[] a, double[] b)
        {
            var o = new double[a.Length];
            for (int i = 0; i < a.Length; i++) o[i] = a[i] - b[i];
            return o;
        }

        public static double[] VecAdd(double[] a, double[] b)
        {
            var o = new double[a.Length];
            for (int i = 0; i < a.Length; i++) o[i] = a[i] + b[i];
            return o;
        }

        public static double Dot(double[] a, double[] b)
        {
            double s = 0.0;
            for (int i = 0; i < a.Length; i++) s += a[i] * b[i];
            return s;
        }

        /// <summary>
        /// In-place Cholesky factorisation (lower): on return <c>a</c> holds L
        /// with a = L Lᵀ. Throws if the matrix is not positive definite.
        /// </summary>
        public static double[,] Cholesky(double[,] a)
        {
            int n = a.GetLength(0);
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= i; j++)
                {
                    double sum = a[i, j];
                    for (int k = 0; k < j; k++) sum -= a[i, k] * a[j, k];
                    if (i > j)
                    {
                        a[i, j] = sum / a[j, j];
                    }
                    else
                    {
                        if (sum <= 1e-12 || double.IsNaN(sum))
                            throw new InvalidOperationException("matrix is not positive definite");
                        a[i, j] = Math.Sqrt(sum);
                    }
                }
            return a;
        }

        /// <summary>
        /// Solves L Lᵀ X = B for X, where L is lower-triangular (from
        /// <see cref="Cholesky"/>). B is overwritten with X in place.
        /// </summary>
        public static double[,] Solve(double[,] L, double[,] B)
        {
            int n = L.GetLength(0), nb = B.GetLength(1);
            for (int c = 0; c < nb; c++)
            {
                for (int i = 0; i < n; i++)               // forward: L y = b
                {
                    double s = B[i, c];
                    for (int k = 0; k < i; k++) s -= L[i, k] * B[k, c];
                    B[i, c] = s / L[i, i];
                }
                for (int i = n - 1; i >= 0; i--)           // backward: Lᵀ x = y
                {
                    double s = B[i, c];
                    for (int k = i + 1; k < n; k++) s -= L[k, i] * B[k, c];
                    B[i, c] = s / L[i, i];
                }
            }
            return B;
        }
    }
}
