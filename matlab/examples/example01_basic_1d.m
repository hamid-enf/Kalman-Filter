%% Example 1 (simple) — filter a noisy 1-D signal with spike rejection.
% Mirrors python/examples/01_basic_1d.py and the C/STM32 temperature example.

rng(0x1A2B3C4D);
TRUE = 25.0;

kf = kalman.KF(1, 1);
kf.set_constant_signal(1e-3, 0.25);

% warm-up WITHOUT gating (estimate starts at 0, far from 25)
for step = 0:19
    z = TRUE + 0.5 * (randn + randn + randn) / 3;
    if mod(step, 15) == 5, z = z + 20; end
    kf.predict();
    kf.update(z);
end

kf.set_gate_threshold(6.63);   % ~99% chi-square, 1 dof

fprintf('step |  raw  | filtered | note\n');
for step = 20:44
    z = TRUE + 0.5 * (randn + randn + randn) / 3;
    if mod(step, 15) == 5, z = z + 20; end
    kf.predict();
    accepted = kf.update_gated(z);
    note = '';
    if ~accepted, note = ' <-- spike rejected'; end
    fprintf('%4d | %5.2f | %7.2f |%s\n', step, z, kf.get_state(), note);
end
fprintf('\nFinal estimate: %.2f (true %.1f)\n', kf.get_state(), TRUE);
