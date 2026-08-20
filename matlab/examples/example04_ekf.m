%% Example 4 (nonlinear) — EKF: estimate an angle from sin(angle) readings.
% Mirrors python/examples/04_ekf.py.

rng(9);
TRUE = 0.8;

ekf = kalman.EKF(1, 1, ...
    @(x, u, dt) x, ...
    @(x, u, dt) 1, ...
    @(x) sin(x), ...
    @(x) cos(x));
ekf.set_state(0.3);
ekf.set_process_noise_scalar(1e-4);
ekf.set_measurement_noise_scalar(0.01);
ekf.set_covariance_scalar(0.5);

fprintf('step | measurement | estimate\n');
for i = 0:39
    ekf.predict();
    z = sin(TRUE) + 0.1 * (randn + randn + randn) / 3;
    ekf.update(z);
    if mod(i, 5) == 0
        fprintf('%4d | %11.3f | %8.3f\n', i, z, ekf.get_state());
    end
end
fprintf('\nEstimated angle: %.3f rad (true %.3f)\n', ekf.get_state(), TRUE);
