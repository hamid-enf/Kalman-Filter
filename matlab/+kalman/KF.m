classdef KF < handle
%KF Linear Kalman filter (MATLAB port of the C kf_kf_* API).
%
%   Predict:  x = F x + u ;  P = F P F' + Q
%   Update:   y = z - H x ;  S = H P H' + R ;  K = (P H') / S   (solved)
%             x = x + K y ;  P = (I-KH) P (I-KH)' + K R K'     (Joseph form)
%
% Also provides the same extensions as the C library: innovation gating (NIS),
% adaptive measurement noise, and the RTS smoother.

    properties
        n (1,1) double = 0          % state dimension
        m (1,1) double = 0          % measurement dimension
        x  double                   % state (n x 1)
        P  double                   % covariance (n x n)
        Q  double                   % process noise (n x n)
        F  double                   % transition (n x n)
        R  double                   % measurement noise (m x m)
        H  double                   % measurement matrix (m x n)
        gate_threshold (1,1) double = 0
        last_nis (1,1) double = 0
    end

    properties (Access = private)
        y_  double                  % innovation
        resid_ double               % residual r = z - H x+
    end

    methods
        function obj = KF(n, m)
            %KF Construct a linear Kalman filter for n states / m measurements.
            if n <= 0 || m <= 0
                error('kalman:KF:badDim', 'dimensions must be positive');
            end
            obj.n = n;
            obj.m = m;
            obj.x = zeros(n, 1);
            obj.P = eye(n);
            obj.Q = zeros(n);
            obj.F = eye(n);
            obj.R = zeros(m);
            obj.H = zeros(m, n);
            obj.y_ = zeros(m, 1);
            obj.resid_ = zeros(m, 1);
        end

        % ---- convenience constructors ----
        function obj = set_constant_signal(obj, q, r)
            %1-D constant-signal model: F=[1], H=[1], P=1.
            obj.x = 0; obj.P = 1; obj.F = 1; obj.H = 1; obj.Q = q; obj.R = r;
        end

        function obj = set_constant_velocity(obj, dt, q_accel, r, p0_pos, p0_vel)
            %Constant-velocity model x = [pos; vel], measuring position only.
            obj.F = [1 dt; 0 1];
            obj.H = [1 0];
            obj.Q = q_accel * [dt^4/4 dt^3/2; dt^3/2 dt^2];
            obj.R = r;
            obj.P = diag([p0_pos p0_vel]);
        end

        % ---- state / covariance / noise access ----
        function set_state(obj, x)
            obj.x = double(x(:));
            if numel(obj.x) ~= obj.n, error('kalman:KF:badState', 'state length mismatch'); end
        end
        function x = get_state(obj), x = obj.x; end

        function set_covariance(obj, P)
            P = double(P);
            if ~isequal(size(P), [obj.n obj.n]), error('kalman:KF:badCov', 'covariance size mismatch'); end
            if ~all(isfinite(P(:))), error('kalman:KF:nonFinite', 'covariance has non-finite values'); end
            obj.P = P;
        end
        function set_covariance_diagonal(obj, d), obj.P = diag(double(d(:))); end
        function set_covariance_scalar(obj, p), obj.P = p * eye(obj.n); end
        function P = get_covariance(obj), P = obj.P; end

        function set_process_noise(obj, Q), obj.Q = double(Q); end
        function set_process_noise_diagonal(obj, d), obj.Q = diag(double(d(:))); end
        function set_process_noise_scalar(obj, q), obj.Q = q * eye(obj.n); end
        function Q = get_process_noise(obj), Q = obj.Q; end

        function set_measurement_noise(obj, R), obj.R = double(R); end
        function set_measurement_noise_diagonal(obj, d), obj.R = diag(double(d(:))); end
        function set_measurement_noise_scalar(obj, r), obj.R = r * eye(obj.m); end
        function R = get_measurement_noise(obj), R = obj.R; end

        function set_transition_matrix(obj, F), obj.F = double(F); end
        function set_measurement_matrix(obj, H), obj.H = double(H); end

        % ---- filtering ----
        function predict(obj, u, dt)
            %predict Time update. u (optional) is a control input in state space.
            if nargin < 3, dt = 0; end
            if nargin < 2, u = []; end
            obj.x = obj.F * obj.x;
            if ~isempty(u)
                obj.x = obj.x + double(u(:));
            end
            obj.P = obj.F * obj.P * obj.F' + obj.Q;
            obj.P = (obj.P + obj.P') / 2;          % symmetrize
        end

        function update(obj, z)
            %update Measurement update (raises an error if S is singular).
            obj.do_update(z, false);
        end

        function accepted = update_gated(obj, z)
            %update_gated Gated update; returns false if rejected as outlier.
            accepted = obj.do_update(z, true);
        end

        function set_gate_threshold(obj, chi2)
            if chi2 < 0, error('kalman:KF:badGate', 'threshold must be >= 0'); end
            obj.gate_threshold = chi2;
        end

        function v = nis(obj)
            v = obj.last_nis;
        end

        function adapt_r(obj, gamma, r_min)
            %adapt_r Residual-based covariance matching:
            %   R <- gamma R + (1 - gamma) (r r' + H P H')
            if gamma <= 0 || gamma >= 1, error('kalman:KF:badGamma', 'gamma in (0,1)'); end
            s = obj.H * obj.P * obj.H' + obj.resid_ * obj.resid_';
            obj.R = gamma * obj.R + (1 - gamma) * s;
            d = diag(obj.R);
            obj.R(1:(obj.m+1):end) = max(d, r_min);   % floor the diagonal
        end
    end

    methods (Access = private)
        function accepted = do_update(obj, z, gate)
            z = double(z(:));
            if numel(z) ~= obj.m, error('kalman:KF:badZ', 'measurement length mismatch'); end
            y = z - obj.H * obj.x;
            obj.y_ = y;
            S = obj.H * obj.P * obj.H' + obj.R;

            if gate && obj.gate_threshold > 0
                nis = y' * (S \ y);
                obj.last_nis = nis;
                if nis > obj.gate_threshold
                    accepted = false;      % reject: state unchanged
                    return;
                end
            end

            K = (S \ (obj.P * obj.H')')';  % K = (P H') S^-1, solved not inverted
            obj.x = obj.x + K * y;
            obj.resid_ = z - obj.H * obj.x;
            IKH = eye(obj.n) - K * obj.H;
            obj.P = IKH * obj.P * IKH' + K * obj.R * K';   % Joseph form
            obj.P = (obj.P + obj.P') / 2;
            accepted = true;
        end
    end

    methods (Static)
        function [x_s, P_s] = rts_smooth_step(n, x_filt, P_filt, F, ...
                                              x_pred, P_pred, x_sn, P_sn)
            %rts_smooth_step One backward RTS smoother step (see C kf_kf_smooth_step).
            A = F * P_filt;              % = F P_k
            C = (P_pred \ A)';           % C = P_k F' P_pred^-1
            d = x_sn - x_pred;
            x_s = x_filt + C * d;
            D = P_sn - P_pred;
            P_s = P_filt + C * D * C';
            P_s = (P_s + P_s') / 2;
        end
    end
end
