%% Example 5 (nonlinear) — UKF: same sin(angle) problem, no Jacobians.
% Mirrors python/examples/05_ukf.py.

rng(10);
TRUE = 0.8;

ukf = kalman.UKF(1, 1, ...
    @(x, u, dt) x, ...
    @(x) sin(x));
ukf.set_state(0.3);
ukf.set_process_noise_scalar(1e-4);
ukf.set_measurement_noise_scalar(0.01);
ukf.set_covariance_scalar(0.5);

fprintf('step | measurement | estimate\n');
for i = 0:39
    ukf.predict();
    z = sin(TRUE) + 0.1 * (randn + randn + randn) / 3;
    ukf.update(z);
    if mod(i, 5) == 0
        fprintf('%4d | %11.3f | %8.3f\n', i, z, ukf.get_state());
    end
end
fprintf('\nEstimated angle: %.3f rad (true %.3f)\n', ukf.get_state(), TRUE);
