classdef UKF < handle
%UKF Unscented Kalman filter (MATLAB port of the C kf_ukf_* API).
%
%   Nonlinear models without Jacobians: 2n+1 sigma points are sampled from the
%   Cholesky factor of P, propagated through f/h, and the mean/covariance
%   reconstructed. Parameters: alpha (spread), beta (prior), kappa (scaling).

    properties
        n (1,1) double = 0
        m (1,1) double = 0
        x  double
        P  double
        Q  double
        R  double
        f  function_handle     % f(x, u, dt) -> x_out
        h  function_handle     % h(x) -> z_out
        alpha (1,1) double = 1
        beta  (1,1) double = 2
        kappa (1,1) double = 0
        lambda_ (1,1) double
        Wm double
        Wc double
    end

    methods
        function obj = UKF(n, m, f, h)
            if nargin < 4, h = []; end
            if nargin < 3, f = []; end
            obj.n = n;
            obj.m = m;
            obj.x = zeros(n, 1);
            obj.P = eye(n);
            obj.Q = zeros(n);
            obj.R = zeros(m);
            obj.f = f; obj.h = h;
            % lambda = alpha^2 (n + kappa) - n = 0 for the defaults (matches C)
            obj.lambda_ = obj.alpha^2 * (obj.n + obj.kappa) - obj.n;
            obj = recompute_weights(obj);
        end

        function set_models(obj, f, h)
            obj.f = f; obj.h = h;
        end

        function set_parameters(obj, alpha, beta, kappa)
            if alpha <= 0 || beta < 0
                error('kalman:UKF:badParams', 'alpha > 0 and beta >= 0');
            end
            obj.alpha = alpha; obj.beta = beta; obj.kappa = kappa;
            obj.lambda_ = alpha^2 * (obj.n + kappa) - obj.n;
            obj = recompute_weights(obj);
        end

        function obj = recompute_weights(obj)
            denom = obj.n + obj.lambda_;
            nsig = 2 * obj.n + 1;
            obj.Wm = 0.5 / denom * ones(1, nsig);
            obj.Wc = 0.5 / denom * ones(1, nsig);
            obj.Wm(1) = obj.lambda_ / denom;
            obj.Wc(1) = obj.lambda_ / denom + (1 - obj.alpha^2 + obj.beta);
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

        function sig = sigma_points(obj)
            L = chol(obj.P)';            % lower; P = L L'
            c = sqrt(obj.n + obj.lambda_);
            sig = zeros(2*obj.n + 1, obj.n);
            sig(1, :) = obj.x';
            for i = 1:obj.n
                sig(i+1, :)     = (obj.x + c * L(:, i))';
                sig(i+1+obj.n, :) = (obj.x - c * L(:, i))';
            end
        end

        function predict(obj, u, dt)
            if nargin < 3, dt = 0; end
            if nargin < 2, u = []; end
            if isempty(obj.f), error('kalman:UKF:noModel', 'call set_models() first'); end
            sig = sigma_points(obj);
            nsig = 2 * obj.n + 1;
            X = zeros(nsig, obj.n);
            for i = 1:nsig
                X(i, :) = obj.f(sig(i, :)', u, dt)';
            end
            xm = obj.Wm * X;                         % 1 x n
            d = X - repmat(xm, nsig, 1);
            obj.P = d' * (obj.Wc' .* d) + obj.Q;
            obj.P = (obj.P + obj.P') / 2;
            obj.x = xm';
        end

        function update(obj, z)
            z = double(z(:));
            if isempty(obj.h), error('kalman:UKF:noModel', 'call set_models() first'); end
            sig = sigma_points(obj);
            nsig = 2 * obj.n + 1;
            Z = zeros(nsig, obj.m);
            for i = 1:nsig
                Z(i, :) = obj.h(sig(i, :)')';
            end
            zm = obj.Wm * Z;                         % 1 x m
            dz = Z - repmat(zm, nsig, 1);
            S = dz' * (obj.Wc' .* dz) + obj.R;
            dx = sig - repmat(obj.x', nsig, 1);
            Pxz = dx' * (obj.Wc' .* dz);             % n x m
            K = (S \ Pxz')';                         % K = Pxz / S
            obj.x = obj.x + K * (z - zm');
            obj.P = obj.P - K * Pxz';
            obj.P = (obj.P + obj.P') / 2;
        end
    end
end
