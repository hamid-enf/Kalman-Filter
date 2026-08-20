using System;
using Kalman;

// Numerical tests for the C# port — mirrors the C/Python suite.
// Run:  dotnet run --project Tests
class TestKalman
{
    static uint _rng;
    static uint Rand()
    {
        uint x = _rng;
        x ^= x << 13; x ^= x >> 17; x ^= x << 5;
        _rng = x; return x;
    }
    static double Noise() => ((int)(Rand() % 2000001u) - 1000000) / 1000000.0;
    static double Gauss() => (Noise() + Noise() + Noise()) / 3.0;

    static int _failed;

    static void Check(bool ok, string name)
    {
        Console.WriteLine((ok ? "  PASS  " : "  FAIL  ") + name);
        if (!ok) _failed++;
    }

    static void Close(double a, double b, double tol, string name)
        => Check(Math.Abs(a - b) <= tol, name);

    static int Main()
    {
        Check1DConvergence();
        CheckConstantVelocity();
        CheckUkfMatchesKf();
        CheckGating();
        CheckAdaptiveR();
        CheckRtsSmoother();
        CheckEkf();
        CheckUkfNonlinear();

        Console.WriteLine();
        Console.WriteLine(_failed == 0 ? "ALL TESTS PASSED" : $"{_failed} TEST(S) FAILED");
        return _failed;
    }

    static void Check1DConvergence()
    {
        _rng = 0x12345678u;
        var kf = KalmanFilter.ConstantSignal(1e-4, 1.0);
        for (int i = 0; i < 500; i++) { kf.Predict(); kf.Update(new[] { 10.0 + Gauss() }); }
        Close(kf.GetState()[0], 10.0, 0.3, "1D convergence");
    }

    static void CheckConstantVelocity()
    {
        _rng = 0x5EED5EEDu;
        var kf = KalmanFilter.ConstantVelocity(1.0, 0.05, 1.0, 100.0, 100.0);
        double pos = 0.0;
        for (int i = 0; i < 500; i++) { pos += 2.0; kf.Predict(); kf.Update(new[] { pos + 0.5 * Gauss() }); }
        var x = kf.GetState();
        Check(Math.Abs(x[0] - pos) < 5.0 && Math.Abs(x[1] - 2.0) < 0.2, "constant velocity");
    }

    static void CheckUkfMatchesKf()
    {
        _rng = 0x1234u;
        var F = new double[,] { { 1.0, 1.0 }, { 0.0, 1.0 } };
        var H = new double[,] { { 1.0, 0.0 } };
        var Q = new double[,] { { 0.0, 0.0 }, { 0.0, 0.05 } };
        var P0 = new double[,] { { 10.0, 0.0 }, { 0.0, 10.0 } };

        var kf = new KalmanFilter(2, 1);
        kf.SetTransitionMatrix(F); kf.SetMeasurementMatrix(H);
        kf.SetProcessNoise(Q); kf.SetMeasurementNoiseScalar(1.0); kf.SetCovariance(P0);

        var ukf = new UnscentedKalmanFilter(2, 1,
            (x, u, dt) => new[] { F[0,0]*x[0]+F[0,1]*x[1], F[1,0]*x[0]+F[1,1]*x[1] },
            x => new[] { H[0,0]*x[0]+H[0,1]*x[1] });
        ukf.SetProcessNoise(Q); ukf.SetMeasurementNoiseScalar(1.0); ukf.SetCovariance(P0);

        for (int i = 0; i < 200; i++)
        {
            double z = 2.0 * (i + 1) + 0.5 * Gauss();
            kf.Predict(); kf.Update(new[] { z });
            ukf.Predict(); ukf.Update(new[] { z });
        }
        var xk = kf.GetState(); var xu = ukf.GetState();
        Check(Math.Abs(xk[0]-xu[0]) < 0.05 && Math.Abs(xk[1]-xu[1]) < 0.02, "UKF == KF (linear)");
    }

    static void CheckGating()
    {
        _rng = 1u;
        var kf = KalmanFilter.ConstantSignal(1e-4, 1.0);
        for (int i = 0; i < 100; i++) { kf.Predict(); kf.Update(new[] { 10.0 + Gauss() }); }
        kf.SetGateThreshold(6.63);
        kf.Predict();
        bool a = kf.UpdateGated(new[] { 10.0 + Gauss() });
        double before = kf.GetState()[0];
        kf.Predict();
        bool b = kf.UpdateGated(new[] { 1000.0 });
        Check(a && !b && Math.Abs(kf.GetState()[0] - before) < 1e-4, "gating rejects spike");
    }

    static void CheckAdaptiveR()
    {
        _rng = 2u;
        var kf = KalmanFilter.ConstantSignal(1e-4, 0.05);
        for (int i = 0; i < 3000; i++)
        {
            kf.Predict();
            kf.Update(new[] { 10.0 + Math.Sqrt(2.0) * Gauss() });
            kf.AdaptR(0.995, 0.001);
        }
        double r = kf.GetMeasurementNoise()[0, 0];
        Check(r > 0.3 && r < 4.0, "adaptive R converges");
    }

    static void CheckRtsSmoother()
    {
        _rng = 3u;
        int N = 60;
        var kf = KalmanFilter.ConstantSignal(0.01, 1.0);
        var xf = new double[N + 1]; var Pf = new double[N + 1];
        var xp = new double[N + 1]; var Pp = new double[N + 1];
        xf[0] = 0.0; Pf[0] = 1.0;
        for (int k = 0; k < N; k++)
        {
            kf.Predict();
            xp[k + 1] = kf.GetState()[0]; Pp[k + 1] = kf.GetCovariance()[0, 0];
            kf.Update(new[] { Gauss() });
            xf[k + 1] = kf.GetState()[0]; Pf[k + 1] = kf.GetCovariance()[0, 0];
        }
        var xs = new double[N + 1]; var Ps = new double[N + 1];
        xs[N] = xf[N]; Ps[N] = Pf[N];
        bool ok = true;
        for (int k = N - 1; k >= 0; k--)
        {
            var r = KalmanFilter.RtsSmoothStep(
                new double[,] { { Pf[k] } }, new double[,] { { 1.0 } }, new double[,] { { Pp[k + 1] } },
                new[] { xf[k] }, new[] { xp[k + 1] },
                new[] { xs[k + 1] }, new double[,] { { Ps[k + 1] } });
            xs[k] = r.x[0]; Ps[k] = r.P[0, 0];

            double C = Pf[k] / Pp[k + 1];
            double xr = xf[k] + C * (xs[k + 1] - xp[k + 1]);
            double Pr = Pf[k] + C * C * (Ps[k + 1] - Pp[k + 1]);
            ok &= Math.Abs(xs[k] - xr) < 1e-4 && Math.Abs(Ps[k] - Pr) < 1e-4;
        }
        Check(ok, "RTS smoother vs reference");
    }

    static void CheckEkf()
    {
        _rng = 4u;
        var ekf = new ExtendedKalmanFilter(1, 1,
            (x, u, dt) => x,
            (x, u, dt) => new double[,] { { 1.0 } },
            x => new[] { x[0] * x[0] },
            x => new double[,] { { 2.0 * x[0] } });
        ekf.SetState(new[] { 2.0 });
        ekf.SetProcessNoiseScalar(1e-4);
        ekf.SetMeasurementNoiseScalar(1.0);
        ekf.SetCovarianceScalar(1.0);
        for (int i = 0; i < 300; i++) { ekf.Predict(); ekf.Update(new[] { 9.0 + 2.0 * Gauss() }); }
        Close(ekf.GetState()[0], 3.0, 0.15, "EKF h=x^2");
    }

    static void CheckUkfNonlinear()
    {
        _rng = 5u;
        var ukf = new UnscentedKalmanFilter(1, 1,
            (x, u, dt) => x,
            x => new[] { x[0] * x[0] });
        ukf.SetState(new[] { 2.0 });
        ukf.SetProcessNoiseScalar(1e-4);
        ukf.SetMeasurementNoiseScalar(1.0);
        ukf.SetCovarianceScalar(1.0);
        for (int i = 0; i < 300; i++) { ukf.Predict(); ukf.Update(new[] { 9.0 + 2.0 * Gauss() }); }
        Close(ukf.GetState()[0], 3.0, 0.15, "UKF h=x^2");
    }
}
