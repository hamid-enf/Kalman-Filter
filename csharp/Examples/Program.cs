using System;
using Kalman;

// C# examples — mirrors python/examples and the C/STM32 examples.
// Run:  dotnet run --project Examples
class Program
{
    // deterministic xorshift32
    static uint _rng;
    static uint Rand()
    {
        uint x = _rng;
        x ^= x << 13; x ^= x >> 17; x ^= x << 5;
        _rng = x; return x;
    }
    static double Noise() => ((int)(Rand() % 2000001u) - 1000000) / 1000000.0;
    static double Gauss() => (Noise() + Noise() + Noise()) / 3.0;

    static void Main()
    {
        Example1D();
        Console.WriteLine();
        ExampleConstantVelocity();
        Console.WriteLine();
        ExampleFusion();
    }

    // (simple) 1-D temperature filter with spike rejection
    static void Example1D()
    {
        _rng = 0x1A2B3C4Du;
        const double truth = 25.0;
        var kf = KalmanFilter.ConstantSignal(1e-3, 0.25);

        for (int step = 0; step < 20; step++)          // warm-up, no gate
        {
            double z = truth + 0.5 * Gauss();
            if (step % 15 == 5) z += 20.0;
            kf.Predict();
            kf.Update(new[] { z });
        }
        kf.SetGateThreshold(6.63);

        Console.WriteLine("step |  raw  | filtered | note");
        for (int step = 20; step < 45; step++)
        {
            double z = truth + 0.5 * Gauss();
            if (step % 15 == 5) z += 20.0;
            kf.Predict();
            bool accepted = kf.UpdateGated(new[] { z });
            Console.WriteLine($"{step,4} | {z,5:F2} | {kf.GetState()[0],7:F2} |" +
                              (accepted ? "" : " <-- spike rejected"));
        }
        Console.WriteLine($"\nFinal estimate: {kf.GetState()[0]:F2} (true {truth:F1})");
    }

    // (medium) position + velocity from position-only measurements
    static void ExampleConstantVelocity()
    {
        _rng = 0x5EED5EEDu;
        const double vel = 1.5, dt = 0.01;
        var kf = KalmanFilter.ConstantVelocity(dt, 0.1, 0.09, 1.0, 1.0);
        double pos = 0.0;

        Console.WriteLine(" t(s) | true pos | measured | filt pos | filt vel");
        for (int step = 0; step < 100; step++)
        {
            pos += vel * dt;
            double z = pos + 0.3 * Gauss();
            kf.Predict();
            kf.Update(new[] { z });
            if (step % 10 == 0)
            {
                var x = kf.GetState();
                Console.WriteLine($"{step*dt,5:F2} | {pos,8:F2} | {z,8:F2} | {x[0],8:F2} | {x[1],8:F2}");
            }
        }
        var xf = kf.GetState();
        Console.WriteLine($"\nFinal: pos {xf[0]:F2} (true {pos:F2}), vel {xf[1]:F2} (true {vel:F1})");
    }

    // (advanced) odometry (velocity) + GNSS (position) fusion
    static void ExampleFusion()
    {
        _rng = 0xAB0BA9Cu;
        const double vel = 2.0, dt = 0.01;
        var kf = new KalmanFilter(2, 1);
        kf.SetTransitionMatrix(new double[,] { { 1.0, dt }, { 0.0, 1.0 } });
        kf.SetProcessNoise(new double[,] { { 0.0, 0.0 }, { 0.0, 0.05 } });
        kf.SetCovariance(new double[,] { { 10.0, 0.0 }, { 0.0, 10.0 } });
        var Hvel = new double[,] { { 0.0, 1.0 } };
        var Hpos = new double[,] { { 1.0, 0.0 } };

        double pos = 0.0;
        int fast = 0;

        Console.WriteLine(" t(s) | odom vel | gnss pos | fused pos | fused vel");
        for (int step = 0; step < 400; step++)
        {
            double odom = vel + 0.2 * Gauss();
            pos += vel * dt;
            kf.Predict();

            kf.SetMeasurementMatrix(Hvel);
            kf.SetMeasurementNoiseScalar(0.04);
            kf.Update(new[] { odom });

            if (++fast >= 100)
            {
                double gnss = pos + 2.0 * Gauss();
                fast = 0;
                kf.SetMeasurementMatrix(Hpos);
                kf.SetMeasurementNoiseScalar(4.0);
                kf.Update(new[] { gnss });
                var x = kf.GetState();
                Console.WriteLine($"{step*dt,5:F2} | {odom,8:F2} | {gnss,8:F2} | {x[0],9:F2} | {x[1],9:F2}");
            }
        }
        var xf = kf.GetState();
        Console.WriteLine($"\nFinal: fused pos {xf[0]:F2} (true {pos:F2}), vel {xf[1]:F2} (true {vel:F1})");
    }
}
