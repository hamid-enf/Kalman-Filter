using System;

namespace Kalman
{
    /// <summary>
    /// Unscented Kalman filter (C# port of the C <c>kf_ukf_*</c> API).
    /// Nonlinear models without Jacobians, via 2n+1 sigma points.
    /// </summary>
    public class UnscentedKalmanFilter
    {
        public int N { get; }
        public int M { get; }

        public double[] x;
        public double[,] P, Q, R;

        public Func<double[], double[], double, double[]> F;   // f(x, u, dt)
        public Func<double[], double[]> H;                     // h(x)

        public double Alpha { get; private set; } = 1.0;
        public double Beta { get; private set; } = 2.0;
        public double Kappa { get; private set; } = 0.0;
        private double lambda;
        private double[] Wm, Wc;

        public UnscentedKalmanFilter(int n, int m,
            Func<double[], double[], double, double[]> f = null,
            Func<double[], double[]> h = null)
        {
            N = n; M = m;
            x = new double[n];
            P = LinAlg.Identity(n);
            Q = new double[n, n];
            R = new double[m, m];
            F = f; H = h;
            SetParameters(1.0, 2.0, 0.0);
        }

        public void SetModels(Func<double[], double[], double, double[]> f,
                              Func<double[], double[]> h)
        {
            F = f; H = h;
        }

        public void SetParameters(double alpha, double beta, double kappa)
        {
            if (alpha <= 0 || beta < 0) throw new ArgumentOutOfRangeException("alpha > 0 and beta >= 0");
            Alpha = alpha; Beta = beta; Kappa = kappa;
            lambda = alpha * alpha * (N + kappa) - N;
            var denom = N + lambda;
            int nsig = 2 * N + 1;
            Wm = new double[nsig]; Wc = new double[nsig];
            for (int i = 0; i < nsig; i++) { Wm[i] = 0.5 / denom; Wc[i] = 0.5 / denom; }
            Wm[0] = lambda / denom;
            Wc[0] = lambda / denom + (1.0 - alpha * alpha + beta);
        }

        public void SetState(double[] v) { Array.Copy(v, x, N); }
        public double[] GetState() { return (double[])x.Clone(); }
        public void SetCovariance(double[,] a) { P = (double[,])a.Clone(); }
        public void SetCovarianceScalar(double p) { P = LinAlg.Scale(LinAlg.Identity(N), p); }
        public double[,] GetCovariance() { return (double[,])P.Clone(); }
        public void SetProcessNoise(double[,] a) { Q = (double[,])a.Clone(); }
        public void SetProcessNoiseScalar(double q) { Q = LinAlg.Scale(LinAlg.Identity(N), q); }
        public void SetMeasurementNoise(double[,] a) { R = (double[,])a.Clone(); }
        public void SetMeasurementNoiseScalar(double r) { R = LinAlg.Scale(LinAlg.Identity(M), r); }

        private double[][] SigmaPoints()
        {
            var L = LinAlg.Cholesky(LinAlg.Clone(P));   // lower; P = L Lᵀ
            double c = Math.Sqrt(N + lambda);
            int nsig = 2 * N + 1;
            var sig = new double[nsig][];
            sig[0] = (double[])x.Clone();
            for (int i = 0; i < N; i++)
            {
                var plus = new double[N];
                var minus = new double[N];
                for (int j = 0; j < N; j++)
                {
                    plus[j] = x[j] + c * L[j, i];
                    minus[j] = x[j] - c * L[j, i];
                }
                sig[i + 1] = plus;
                sig[i + 1 + N] = minus;
            }
            return sig;
        }

        public void Predict(double[] u = null, double dt = 0.0)
        {
            if (F == null) throw new InvalidOperationException("call SetModels() first");
            var sig = SigmaPoints();
            int nsig = 2 * N + 1;

            var X = new double[nsig][];
            for (int i = 0; i < nsig; i++) X[i] = F(sig[i], u, dt);

            var xm = new double[N];
            for (int i = 0; i < nsig; i++)
                for (int j = 0; j < N; j++) xm[j] += Wm[i] * X[i][j];

            var Pn = (double[,])Q.Clone();
            for (int i = 0; i < nsig; i++)
            {
                var d = new double[N];
                for (int j = 0; j < N; j++) d[j] = X[i][j] - xm[j];
                for (int r = 0; r < N; r++)
                    for (int c = 0; c < N; c++)
                        Pn[r, c] += Wc[i] * d[r] * d[c];
            }
            LinAlg.SymmetrizeInPlace(Pn);
            P = Pn;
            x = xm;
        }

        public void Update(double[] z)
        {
            if (H == null) throw new InvalidOperationException("call SetModels() first");
            var sig = SigmaPoints();
            int nsig = 2 * N + 1;

            var Z = new double[nsig][];
            for (int i = 0; i < nsig; i++) Z[i] = H(sig[i]);

            var zm = new double[M];
            for (int i = 0; i < nsig; i++)
                for (int j = 0; j < M; j++) zm[j] += Wm[i] * Z[i][j];

            var S = (double[,])R.Clone();
            for (int i = 0; i < nsig; i++)
            {
                var dz = new double[M];
                for (int j = 0; j < M; j++) dz[j] = Z[i][j] - zm[j];
                for (int r = 0; r < M; r++)
                    for (int c = 0; c < M; c++)
                        S[r, c] += Wc[i] * dz[r] * dz[c];
            }

            var Pxz = new double[N, M];
            for (int i = 0; i < nsig; i++)
            {
                var dx = new double[N];
                for (int j = 0; j < N; j++) dx[j] = sig[i][j] - x[j];
                var dz = new double[M];
                for (int j = 0; j < M; j++) dz[j] = Z[i][j] - zm[j];
                for (int r = 0; r < N; r++)
                    for (int c = 0; c < M; c++)
                        Pxz[r, c] += Wc[i] * dx[r] * dz[c];
            }

            var K = LinAlg.Transpose(LinAlg.Solve(LinAlg.Cholesky(LinAlg.Clone(S)), LinAlg.Transpose(Pxz)));

            var innov = LinAlg.VecSub(z, zm);
            x = LinAlg.VecAdd(x, LinAlg.MatVec(K, innov));
            P = LinAlg.Subtract(P, LinAlg.Multiply(K, LinAlg.Transpose(Pxz)));
            LinAlg.SymmetrizeInPlace(P);
        }
    }
}
