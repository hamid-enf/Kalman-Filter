using System;

namespace Kalman
{
    /// <summary>
    /// Extended Kalman filter (C# port of the C <c>kf_ekf_*</c> API).
    /// Nonlinear models linearised with user-provided Jacobians.
    /// </summary>
    public class ExtendedKalmanFilter
    {
        public int N { get; }
        public int M { get; }

        public double[] x;
        public double[,] P, Q, R;

        public Func<double[], double[], double, double[]> F;      // f(x, u, dt)
        public Func<double[], double[], double, double[,]> Fjac;   // df/dx
        public Func<double[], double[]> H;                         // h(x)
        public Func<double[], double[,]> Hjac;                     // dh/dx

        public ExtendedKalmanFilter(int n, int m,
            Func<double[], double[], double, double[]> f = null,
            Func<double[], double[], double, double[,]> fjac = null,
            Func<double[], double[]> h = null,
            Func<double[], double[,]> hjac = null)
        {
            N = n; M = m;
            x = new double[n];
            P = LinAlg.Identity(n);
            Q = new double[n, n];
            R = new double[m, m];
            F = f; Fjac = fjac; H = h; Hjac = hjac;
        }

        public void SetModels(Func<double[], double[], double, double[]> f,
                              Func<double[], double[], double, double[,]> fjac,
                              Func<double[], double[]> h,
                              Func<double[], double[,]> hjac)
        {
            F = f; Fjac = fjac; H = h; Hjac = hjac;
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

        public void Predict(double[] u = null, double dt = 0.0)
        {
            if (F == null || Fjac == null) throw new InvalidOperationException("call SetModels() first");
            var Fm = Fjac(x, u, dt);
            x = F(x, u, dt);
            P = LinAlg.Add(LinAlg.Multiply(LinAlg.Multiply(Fm, P), LinAlg.Transpose(Fm)), Q);
            LinAlg.SymmetrizeInPlace(P);
        }

        public void Update(double[] z)
        {
            if (H == null || Hjac == null) throw new InvalidOperationException("call SetModels() first");
            var Hm = Hjac(x);
            var y = LinAlg.VecSub(z, H(x));                       // innovation

            var S = LinAlg.Add(LinAlg.Multiply(LinAlg.Multiply(Hm, P), LinAlg.Transpose(Hm)), R);
            var B = LinAlg.Multiply(P, LinAlg.Transpose(Hm));     // n x m
            var rhs = LinAlg.Transpose(B);
            var K = LinAlg.Transpose(LinAlg.Solve(LinAlg.Cholesky(LinAlg.Clone(S)), rhs));

            x = LinAlg.VecAdd(x, LinAlg.MatVec(K, y));
            var IKH = LinAlg.Subtract(LinAlg.Identity(N), LinAlg.Multiply(K, Hm));
            P = LinAlg.Add(
                    LinAlg.Multiply(LinAlg.Multiply(IKH, P), LinAlg.Transpose(IKH)),
                    LinAlg.Multiply(LinAlg.Multiply(K, R), LinAlg.Transpose(K)));
            LinAlg.SymmetrizeInPlace(P);
        }
    }
}
