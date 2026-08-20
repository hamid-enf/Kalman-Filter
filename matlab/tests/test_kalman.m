function test_kalman()
%TEST_KALMAN Numerical tests for the MATLAB port (mirrors the C/Python suite).
%
%   Run:  addpath('..'); test_kalman();
%
% Uses the documented deviations: errors are MATLAB errors, and update_gated
% returns a logical. Checks: convergence, UKF==KF on a linear system, RTS
% smoother vs an exact scalar reference, gating, adaptive R, EKF/UKF nonlinear.

    failed = 0;
    failed = failed + check_1d_convergence();
    failed = failed + check_constant_velocity();
    failed = failed + check_ukf_matches_kf();
    failed = failed + check_gating();
    failed = failed + check_adaptive_r();
    failed = failed + check_rts_smoother();
    failed = failed + check_ekf();
    failed = failed + check_ukf_nonlinear();

    if failed == 0
        fprintf('\nALL TESTS PASSED\n');
    else
        fprintf('\n%d TEST(S) FAILED\n', failed);
    end
end

function f = check_1d_convergence()
    rng(0x12345678);
    kf = kalman.KF(1, 1);
    kf.set_constant_signal(1e-4, 1.0);
    for i = 1:500
        kf.predict();
        kf.update(10 + (randn + randn + randn) / 3);
    end
    f = report('1D convergence', abs(kf.get_state() - 10) < 0.3);
end

function f = check_constant_velocity()
    rng(0x5EED5EED);
    kf = kalman.KF(2, 1);
    kf.set_constant_velocity(1.0, 0.05, 1.0, 100, 100);
    pos = 0;
    for i = 1:500
        pos = pos + 2.0;
        kf.predict();
        kf.update(pos + 0.5 * (randn + randn + randn) / 3);
    end
    x = kf.get_state();
    ok = abs(x(1) - pos) < 5 && abs(x(2) - 2.0) < 0.2;
    f = report('constant velocity', ok);
end

function f = check_ukf_matches_kf()
    rng(0x1234);
    F = [1 1; 0 1]; H = [1 0]; Q = [0 0; 0 0.05]; P0 = diag([10 10]);

    kf = kalman.KF(2, 1);
    kf.set_transition_matrix(F); kf.set_measurement_matrix(H);
    kf.set_process_noise(Q); kf.set_measurement_noise_scalar(1); kf.set_covariance(P0);

    ukf = kalman.UKF(2, 1, @(x, u, dt) F * x, @(x) H * x);
    ukf.set_process_noise(Q); ukf.set_measurement_noise_scalar(1); ukf.set_covariance(P0);

    for i = 1:200
        z = 2 * i + 0.5 * (randn + randn + randn) / 3;
        kf.predict();  kf.update(z);
        ukf.predict(); ukf.update(z);
    end
    xk = kf.get_state(); xu = ukf.get_state();
    ok = abs(xk(1) - xu(1)) < 0.05 && abs(xk(2) - xu(2)) < 0.02;
    f = report('UKF == KF (linear)', ok);
end

function f = check_gating()
    rng(1);
    kf = kalman.KF(1, 1);
    kf.set_constant_signal(1e-4, 1.0);
    for i = 1:100
        kf.predict(); kf.update(10 + (randn + randn + randn) / 3);
    end
    kf.set_gate_threshold(6.63);
    kf.predict();
    a = kf.update_gated(10 + (randn + randn + randn) / 3);
    before = kf.get_state();
    kf.predict();
    b = kf.update_gated(1000);
    ok = a && ~b && abs(kf.get_state() - before) < 1e-4;
    f = report('gating rejects spike', ok);
end

function f = check_adaptive_r()
    rng(2);
    kf = kalman.KF(1, 1);
    kf.set_constant_signal(1e-4, 0.05);
    for i = 1:3000
        kf.predict();
        kf.update(10 + sqrt(2) * (randn + randn + randn) / 3);
        kf.adapt_r(0.995, 0.001);
    end
    r = kf.get_measurement_noise();
    f = report('adaptive R converges', r > 0.3 && r < 4.0);
end

function f = check_rts_smoother()
    rng(3);
    N = 60;
    kf = kalman.KF(1, 1);
    kf.set_constant_signal(0.01, 1.0);
    xf = zeros(N + 1, 1); Pf = zeros(N + 1, 1);
    xp = zeros(N + 1, 1); Pp = zeros(N + 1, 1);
    xf(1) = 0; Pf(1) = 1;
    for k = 1:N
        kf.predict();
        xp(k + 1) = kf.get_state(); Pp(k + 1) = kf.get_covariance();
        kf.update((randn + randn + randn) / 3);
        xf(k + 1) = kf.get_state(); Pf(k + 1) = kf.get_covariance();
    end
    xs = zeros(N + 1, 1); Ps = zeros(N + 1, 1);
    xs(N + 1) = xf(N + 1); Ps(N + 1) = Pf(N + 1);
    for k = N:-1:1
        [xs(k), Ps(k)] = kalman.KF.rts_smooth_step(1, xf(k), Pf(k), 1, ...
            xp(k + 1), Pp(k + 1), xs(k + 1), Ps(k + 1));
    end
    ok = true;
    for k = N:-1:1
        C = Pf(k) / Pp(k + 1);
        xr = xf(k) + C * (xs(k + 1) - xp(k + 1));
        Pr = Pf(k) + C * C * (Ps(k + 1) - Pp(k + 1));
        ok = ok && abs(xs(k) - xr) < 1e-4 && abs(Ps(k) - Pr) < 1e-4;
    end
    f = report('RTS smoother vs reference', ok);
end

function f = check_ekf()
    rng(4);
    ekf = kalman.EKF(1, 1, @(x, u, dt) x, @(x, u, dt) 1, ...
        @(x) x^2, @(x) 2 * x);
    ekf.set_state(2);
    ekf.set_process_noise_scalar(1e-4);
    ekf.set_measurement_noise_scalar(1.0);
    ekf.set_covariance_scalar(1.0);
    for i = 1:300
        ekf.predict();
        ekf.update(9 + 2 * (randn + randn + randn) / 3);
    end
    f = report('EKF h=x^2', abs(ekf.get_state() - 3) < 0.15);
end

function f = check_ukf_nonlinear()
    rng(5);
    ukf = kalman.UKF(1, 1, @(x, u, dt) x, @(x) x^2);
    ukf.set_state(2);
    ukf.set_process_noise_scalar(1e-4);
    ukf.set_measurement_noise_scalar(1.0);
    ukf.set_covariance_scalar(1.0);
    for i = 1:300
        ukf.predict();
        ukf.update(9 + 2 * (randn + randn + randn) / 3);
    end
    f = report('UKF h=x^2', abs(ukf.get_state() - 3) < 0.15);
end

function f = report(name, ok)
    if ok
        fprintf('  PASS  %s\n', name);
        f = 0;
    else
        fprintf('  FAIL  %s\n', name);
        f = 1;
    end
end
