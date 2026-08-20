%% Example 2 (medium) — position + velocity from position-only measurements.

rng(0x5EED5EED);
TRUE_VEL = 1.5;
DT = 0.01;

kf = kalman.KF(2, 1);
kf.set_constant_velocity(DT, 0.1, 0.09, 1.0, 1.0);

pos = 0;
fprintf(' t(s) | true pos | measured | filt pos | filt vel\n');
for step = 0:99
    pos = pos + TRUE_VEL * DT;
    z = pos + 0.3 * (randn + randn + randn) / 3;
    kf.predict();
    kf.update(z);
    if mod(step, 10) == 0
        t = step * DT;
        x = kf.get_state();
        fprintf('%5.2f | %8.2f | %8.2f | %8.2f | %8.2f\n', ...
            t, pos, z, x(1), x(2));
    end
end
x = kf.get_state();
fprintf('\nFinal: pos %.2f (true %.2f), vel %.2f (true %.1f)\n', ...
    x(1), pos, x(2), TRUE_VEL);
