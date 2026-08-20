using System;
using Kalman;

// Cross-language calibration scenario — C# implementation.
// See calibration/README.md for the spec. Identical xorshift32 + noise to
// fusion.c so the printed numbers match the C reference.
class Program
{
    const double TrueVel = 2.0, Dt = 0.01;
    const int GnssEvery = 100;
    static uint _rng = 0xC0FFEE;

    static uint Xorshift32()
    {
        uint x = _rng;
        x ^= x << 13; x ^= x >> 17; x ^= x << 5;
        _rng = x; return x;
    }
    static double Uniform() => ((int)(Xorshift32() % 2000001u) - 1000000) / 1000000.0;
    static double Gauss() => (Uniform() + Uniform() + Uniform()) / 3.0;

    static void Main()
    {
        var F = new double[,] { { 1.0, Dt }, { 0.0, 1.0 } };
        var Hvel = new double[,] { { 0.0, 1.0 } };
        var Hpos = new double[,] { { 1.0, 0.0 } };
        var Q = new double[,] { { 0.0, 0.0 }, { 0.0, 0.05 } };
        var P0 = new double[,] { { 10.0, 0.0 }, { 0.0, 10.0 } };

        var kf = new KalmanFilter(2, 1);
        kf.SetTransitionMatrix(F);
        kf.SetProcessNoise(Q);
        kf.SetCovariance(P0);

        double truePos = 0.0;
        int fast = 0;

        Console.WriteLine("  t(s) | fused pos | fused vel");
        for (int step = 0; step < 400; step++)
        {
            double odom = TrueVel + 0.2 * Gauss();
            truePos += TrueVel * Dt;

            kf.Predict();

            kf.SetMeasurementMatrix(Hvel);
            kf.SetMeasurementNoiseScalar(0.04);
            kf.Update(new[] { odom });

            if (++fast >= GnssEvery)
            {
                double gnss = truePos + 2.0 * Gauss();
                fast = 0;
                kf.SetMeasurementMatrix(Hpos);
                kf.SetMeasurementNoiseScalar(4.0);
                kf.Update(new[] { gnss });
                var x = kf.GetState();
                Console.WriteLine($" {step * Dt,5:F2} | {x[0],9:F5} | {x[1],9:F5}");
            }
        }
    }
}
