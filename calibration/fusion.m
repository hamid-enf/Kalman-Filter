function fusion()
%FUSION Cross-language calibration scenario — MATLAB implementation.
%   See README.md for the spec. Identical xorshift32 + noise to fusion.c, so
%   the printed numbers match the C reference. Run: addpath('../matlab'); fusion;

TRUE_VEL = 2.0;
DT = 0.01;
GNSS_EVERY = 100;

F    = [1 DT; 0 1];
Hvel = [0 1];
Hpos = [1 0];
Q    = [0 0; 0 0.05];
P0   = diag([10 10]);

kf = kalman.KF(2, 1);
kf.set_transition_matrix(F);
kf.set_process_noise(Q);
kf.set_covariance(P0);

true_pos = 0;
fast = 0;

fprintf('  t(s) | fused pos | fused vel\n');
for step = 0:399
    odom = TRUE_VEL + 0.2 * gauss();
    true_pos = true_pos + TRUE_VEL * DT;

    kf.predict();

    kf.set_measurement_matrix(Hvel);
    kf.set_measurement_noise_scalar(0.04);
    kf.update(odom);

    fast = fast + 1;
    if fast >= GNSS_EVERY
        gnss = true_pos + 2.0 * gauss();
        fast = 0;
        kf.set_measurement_matrix(Hpos);
        kf.set_measurement_noise_scalar(4.0);
        kf.update(gnss);
        x = kf.get_state();
        fprintf(' %5.2f | %9.5f | %9.5f\n', step * DT, x(1), x(2));
    end
end
end

% ---- identical deterministic RNG to the C/Python reference ----------------

function x = xorshift32()
persistent rng;
if isempty(rng), rng = uint32(hex2dec('C0FFEE')); end
s = rng;
s = bitxor(s, bitshift(s, 13));
s = bitxor(s, bitshift(s, -17));
s = bitxor(s, bitshift(s, 5));
rng = s;
x = s;
end

function u = uniform()
u = (double(mod(xorshift32(), 2000001)) - 1000000) / 1000000.0;
end

function g = gauss()
g = (uniform() + uniform() + uniform()) / 3.0;
end
