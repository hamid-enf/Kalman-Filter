using System;

namespace Kalman
{
    /// <summary>
    /// Linear Kalman filter (C# port of the C <c>kf_kf_*</c> API).
    ///
    /// <code>
    ///   predict : x = F x + u ;  P = F P Fᵀ + Q
    ///   update  : y = z - H x ;  S = H P Hᵀ + R ;  K = (P Hᵀ) S⁻¹ (solved)
    ///             x = x + K y ;  P = (I-KH) P (I-KH)ᵀ + K R Kᵀ   (Joseph form)
    /// </code>
    ///
    /// Also provides innovation gating (NIS), adaptive R, and the RTS smoother,
    /// matching the C library.
    /// </summary>
    public class KalmanFilter
    {
        public int N { get; }
        public int M { get; }

        public double[] x;                     // state (N)
        public double[,] P, Q, F;              // (N x N)
        public double[,] R;                    // (M x M)
        public double[,] H;                    // (M x N)

        public double GateThreshold { get; set; }
        public double LastNis { get; private set; }

        private double[] y;                    // innovation
        private double[] resid;                // residual r = z - H x+

        public KalmanFilter(int n, int m)
        {
            if (n <= 0 || m <= 0) throw new ArgumentOutOfRangeException("dimensions must be positive");
            N = n; M = m;
            x = new double[n];
            P = LinAlg.Identity(n);
            Q = new double[n, n];
            F = LinAlg.Identity(n);
            R = new double[m, m];
            H = new double[m, n];
            y = new double[m];
            resid = new double[m];
        }

        // ---- convenience constructors ----
        public static KalmanFilter ConstantSignal(double q, double r)
        {
            var f = new KalmanFilter(1, 1);
            f.F[0, 0] = 1.0; f.H[0, 0] = 1.0;
            f.Q[0, 0] = q; f.R[0, 0] = r; f.P[0, 0] = 1.0;
            return f;
        }

        public static KalmanFilter ConstantVelocity(double dt, double qAccel,
                                                    double r, double p0Pos, double p0Vel)
        {
            var f = new KalmanFilter(2, 1);
            f.F = new double[,] { { 1.0, dt }, { 0.0, 1.0 } };
            f.H = new double[,] { { 1.0, 0.0 } };
            double d2 = dt * dt, d3 = d2 * dt, d4 = d3 * dt;
            f.Q = new double[,] { { qAccel * d4 / 4, qAccel * d3 / 2 },
                                  { qAccel * d3 / 2, qAccel * d2 } };
            f.R[0, 0] = r;
            f.P = new double[,] { { p0Pos, 0.0 }, { 0.0, p0Vel } };
            return f;
        }

        // ---- state / covariance / noise access ----
        public void SetState(double[] v)
        {
            if (v.Length != N) throw new ArgumentException("state length mismatch");
            Array.Copy(v, x, N);
        }

        public double[] GetState() { return (double[])x.Clone(); }

        public void SetCovariance(double[,] a)
        {
            CheckShape(a, N, N, "covariance");
            foreach (var v in a) if (double.IsNaN(v) || double.IsInfinity(v))
                    throw new ArgumentException("covariance has non-finite values");
            P = (double[,])a.Clone();
        }

        public void SetCovarianceDiagonal(double[] d) { P = Diag(d); }
        public void SetCovarianceScalar(double p) { P = LinAlg.Scale(LinAlg.Identity(N), p); }
        public double[,] GetCovariance() { return (double[,])P.Clone(); }

        public void SetProcessNoise(double[,] a) { CheckShape(a, N, N, "Q"); Q = (double[,])a.Clone(); }
        public void SetProcessNoiseDiagonal(double[] d) { Q = Diag(d); }
        public void SetProcessNoiseScalar(double q) { Q = LinAlg.Scale(LinAlg.Identity(N), q); }
        public double[,] GetProcessNoise() { return (double[,])Q.Clone(); }

        public void SetMeasurementNoise(double[,] a) { CheckShape(a, M, M, "R"); R = (double[,])a.Clone(); }
        public void SetMeasurementNoiseDiagonal(double[] d) { R = Diag(d); }
        public void SetMeasurementNoiseScalar(double r) { R = LinAlg.Scale(LinAlg.Identity(M), r); }
        public double[,] GetMeasurementNoise() { return (double[,])R.Clone(); }

        public void SetTransitionMatrix(double[,] a) { CheckShape(a, N, N, "F"); F = (double[,])a.Clone(); }
        public void SetMeasurementMatrix(double[,] a) { CheckShape(a, M, N, "H"); H = (double[,])a.Clone(); }

        private static double[,] Diag(double[] d)
        {
            var a = new double[d.Length, d.Length];
            for (int i = 0; i < d.Length; i++) a[i, i] = d[i];
            return a;
        }

        private static void CheckShape(double[,] a, int r, int c, string name)
        {
            if (a.GetLength(0) != r || a.GetLength(1) != c)
                throw new ArgumentException($"{name} must be {r}x{c}");
        }

        // ---- filtering ----
        public void Predict(double[] u = null, double dt = 0.0)
        {
            x = LinAlg.MatVec(F, x);
            if (u != null)
            {
                if (u.Length != N) throw new ArgumentException("control input length mismatch");
                x = LinAlg.VecAdd(x, u);
            }
            P = LinAlg.Add(LinAlg.Multiply(LinAlg.Multiply(F, P), LinAlg.Transpose(F)), Q);
            LinAlg.SymmetrizeInPlace(P);
        }

        private bool DoUpdate(double[] z, bool gate)
        {
            if (z.Length != M) throw new ArgumentException("measurement length mismatch");

            var Hx = LinAlg.MatVec(H, x);
            y = LinAlg.VecSub(z, Hx);                        // innovation
            var S = LinAlg.Add(LinAlg.Multiply(LinAlg.Multiply(H, P), LinAlg.Transpose(H)), R);

            if (gate && GateThreshold > 0.0)
            {
                // NIS = yᵀ S⁻¹ y  (solve S t = y, then y·t)
                var L = LinAlg.Cholesky(LinAlg.Clone(S));
                var t = LinAlg.Solve(L, ToColumn(y));
                LastNis = LinAlg.Dot(y, ColumnToVec(t));
                if (LastNis > GateThreshold) return false;   // rejected
            }

            // K = (P Hᵀ) S⁻¹   (solve S Kᵀ = (P Hᵀ)ᵀ, then transpose)
            var B = LinAlg.Multiply(P, LinAlg.Transpose(H));  // n x m
            var rhs = LinAlg.Transpose(B);                    // m x n
            var Ls = LinAlg.Cholesky(LinAlg.Clone(S));
            rhs = LinAlg.Solve(Ls, rhs);                      // = Kᵀ
            var K = LinAlg.Transpose(rhs);                    // n x m

            x = LinAlg.VecAdd(x, LinAlg.MatVec(K, y));
            resid = LinAlg.VecSub(z, LinAlg.MatVec(H, x));

            var KH = LinAlg.Multiply(K, H);
            var IKH = LinAlg.Subtract(LinAlg.Identity(N), KH);
            P = LinAlg.Add(
                    LinAlg.Multiply(LinAlg.Multiply(IKH, P), LinAlg.Transpose(IKH)),
                    LinAlg.Multiply(LinAlg.Multiply(K, R), LinAlg.Transpose(K)));
            LinAlg.SymmetrizeInPlace(P);
            return true;
        }

        public void Update(double[] z) { DoUpdate(z, gate: false); }

        public bool UpdateGated(double[] z) { return DoUpdate(z, gate: true); }

        public void SetGateThreshold(double chi2)
        {
            if (chi2 < 0) throw new ArgumentOutOfRangeException("threshold must be >= 0");
            GateThreshold = chi2;
        }

        public double Nis() { return LastNis; }

        public void AdaptR(double gamma, double rMin)
        {
            if (!(gamma > 0.0 && gamma < 1.0)) throw new ArgumentOutOfRangeException("gamma in (0,1)");
            var HPHt = LinAlg.Multiply(LinAlg.Multiply(H, P), LinAlg.Transpose(H));
            var rr = Outer(resid, resid);
            var s = LinAlg.Add(HPHt, rr);
            var newR = LinAlg.Add(LinAlg.Scale(R, gamma), LinAlg.Scale(s, 1.0 - gamma));
            for (int i = 0; i < M; i++)
                if (newR[i, i] < rMin) newR[i, i] = rMin;
            R = newR;
        }

        private static double[,] Outer(double[] a, double[] b)
        {
            var o = new double[a.Length, b.Length];
            for (int i = 0; i < a.Length; i++)
                for (int j = 0; j < b.Length; j++)
                    o[i, j] = a[i] * b[j];
            return o;
        }

        private static double[,] ToColumn(double[] v)
        {
            var c = new double[v.Length, 1];
            for (int i = 0; i < v.Length; i++) c[i, 0] = v[i];
            return c;
        }

        private static double[] ColumnToVec(double[,] c)
        {
            var v = new double[c.GetLength(0)];
            for (int i = 0; i < v.Length; i++) v[i] = c[i, 0];
            return v;
        }

        // ---- RTS smoother (one backward step) ----
        public static (double[] x, double[,] P) RtsSmoothStep(
            double[,] P_filt, double[,] F, double[,] P_pred,
            double[] x_filt, double[] x_pred,
            double[] x_smooth_next, double[,] P_smooth_next)
        {
            var A = LinAlg.Multiply(F, P_filt);          // = F P_k
            var C = LinAlg.Transpose(LinAlg.Solve(LinAlg.Cholesky(LinAlg.Clone(P_pred)), A));
            var d = LinAlg.VecSub(x_smooth_next, x_pred);
            var xs = LinAlg.VecAdd(x_filt, LinAlg.MatVec(C, d));
            var D = LinAlg.Subtract(P_smooth_next, P_pred);
            var Ps = LinAlg.Add(P_filt, LinAlg.Multiply(LinAlg.Multiply(C, D), LinAlg.Transpose(C)));
            LinAlg.SymmetrizeInPlace(Ps);
            return (xs, Ps);
        }
    }
}
