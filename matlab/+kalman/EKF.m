classdef EKF < handle
%EKF Extended Kalman filter (MATLAB port of the C kf_ekf_* API).
%
%   Nonlinear models, linearised with the user's Jacobians:
%     x(k+1) = f(x, u, dt) + w,   z = h(x) + v
%     F = df/dx,  H = dh/dx  (at the current estimate)

    properties
        n (1,1) double = 0
        m (1,1) double = 0
        x  double
        P  double
        Q  double
        R  double
        f     function_handle   % f(x, u, dt) -> x_out
        F_jac function_handle   % F_jac(x, u, dt) -> F (n x n)
        h     function_handle   % h(x) -> z_out
        H_jac function_handle   % H_jac(x) -> H (m x n)
    end

    methods
        function obj = EKF(n, m, f, F_jac, h, H_jac)
            if nargin < 5, H_jac = []; end
            if nargin < 4, h = []; end
            if nargin < 3, F_jac = []; end
            if nargin < 2, f = []; end
            obj.n = n;
            obj.m = m;
            obj.x = zeros(n, 1);
            obj.P = eye(n);
            obj.Q = zeros(n);
            obj.R = zeros(m);
            obj.f = f; obj.F_jac = F_jac; obj.h = h; obj.H_jac = H_jac;
        end

        function set_models(obj, f, F_jac, h, H_jac)
            obj.f = f; obj.F_jac = F_jac; obj.h = h; obj.H_jac = H_jac;
        end

        function set_state(obj, x), obj.x = double(x(:)); end
        function x = get_state(obj), x = obj.x; end
        function set_covariance(obj, P), obj.P = double(P); end
        function set_covariance_diagonal(obj, d), obj.P = diag(double(d(:))); end
        function set_covariance_scalar(obj, p), obj.P = p * eye(obj.n); end
        function P = get_covariance(obj), P = obj.P; end
        function set_process_noise(obj, Q), obj.Q = double(Q); end
        function set_process_noise_scalar(obj, q), obj.Q = q * eye(obj.n); end
        function set_measurement_noise(obj, R), obj.R = double(R); end
        function set_measurement_noise_scalar(obj, r), obj.R = r * eye(obj.m); end

        function predict(obj, u, dt)
            if nargin < 3, dt = 0; end
            if nargin < 2, u = []; end
            if isempty(obj.f) || isempty(obj.F_jac)
                error('kalman:EKF:noModel', 'call set_models() first');
            end
            F = double(obj.F_jac(obj.x, u, dt));
            obj.x = double(obj.f(obj.x, u, dt));
            obj.P = F * obj.P * F' + obj.Q;
            obj.P = (obj.P + obj.P') / 2;
        end

        function update(obj, z)
            z = double(z(:));
            if isempty(obj.h) || isempty(obj.H_jac)
                error('kalman:EKF:noModel', 'call set_models() first');
            end
            H = double(obj.H_jac(obj.x));
            y = z - double(obj.h(obj.x));
            S = H * obj.P * H' + obj.R;
            K = (S \ (obj.P * H')')';
            obj.x = obj.x + K * y;
            IKH = eye(obj.n) - K * H;
            obj.P = IKH * obj.P * IKH' + K * obj.R * K';
            obj.P = (obj.P + obj.P') / 2;
        end
    end
end
