%% Example 3 (advanced) — fuse odometry (velocity) with GNSS (position).
% Two sensors at different rates/qualities, fused via sequential updates.

rng(0xAB0BA9C);
TRUE_VEL = 2.0;
FAST_DT = 0.01;
GNSS_EVERY = 100;

F     = [1 FAST_DT; 0 1];
H_vel = [0 1];
H_pos = [1 0];
Q     = [0 0; 0 0.05];
P0    = diag([10 10]);

kf = kalman.KF(2, 1);
kf.set_transition_matrix(F);
kf.set_process_noise(Q);
kf.set_covariance(P0);

pos = 0;
fast_count = 0;
fprintf(' t(s) | odom vel | gnss pos | fused pos | fused vel\n');
for step = 0:399
    odom = TRUE_VEL + 0.2 * (randn + randn + randn) / 3;
    pos = pos + TRUE_VEL * FAST_DT;

    kf.predict();

    % fast update: correct velocity with odometry
    kf.set_measurement_matrix(H_vel);
    kf.set_measurement_noise_scalar(0.04);
    kf.update(odom);

    fast_count = fast_count + 1;
    if fast_count >= GNSS_EVERY
        gnss = pos + 2.0 * (randn + randn + randn) / 3;
        fast_count = 0;
        kf.set_measurement_matrix(H_pos);
        kf.set_measurement_noise_scalar(4.0);
        kf.update(gnss);

        t = step * FAST_DT;
        x = kf.get_state();
        fprintf('%5.2f | %8.2f | %8.2f | %9.2f | %9.2f\n', ...
            t, odom, gnss, x(1), x(2));
    end
end
x = kf.get_state();
fprintf('\nFinal: fused pos %.2f (true %.2f), vel %.2f (true %.1f)\n', ...
    x(1), pos, x(2), TRUE_VEL);
