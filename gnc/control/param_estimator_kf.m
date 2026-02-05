function [params, K] = param_estimator_kf(time, w, d, c)

    % d = encoder angle (rad)
    % c = rho * area * arm / inertia
    % dw/dt = c * Cl * d + C0 * c = phi_k' * p
    % Cl = canard lift coeff
    % C0 = rocket induced angular acceleration / (rho * area * arm)

    % tuning parameters
    Q = diag([1e-3 1e-8]);
    R = 10;

    persistent t p_k_minus P_k_minus w_old d_hat w_dot
    if isempty(t)
        t = -0.01;
    end
    if isempty(p_k_minus)
        p_k_minus = [2; 0];
    end
    if isempty(P_k_minus)
        P_k_minus = Q;
    end
    if isempty(w_old)
        w_old = w;
    end
    if isempty(d_hat)
        d_hat = 0;
    end
    if isempty(w_dot)
       w_dot = 0;
    end

    % lowpass command and measurement
    d_hat = 0.75 * d_hat + 0.25 * d;
    w_dot = 0.75 * w_dot + 0.25 * (w - w_old) / (time - t);

    % Kalman filter
    phi = c * [d_hat; 1];
    P_k_minus = P_k_minus + Q;
    K = P_k_minus * phi / (phi' * P_k_minus * phi + R); % the stuff inside in brackets is just a scalar so you can just divide
    params = p_k_minus + K * (w_dot - phi' * p_k_minus);
    P_k = (eye(2) - K * phi') * P_k_minus * (eye(2) - K * phi')' + K * R * K'; % Joseph form for numerical stability
    
    % update for next cycle
    t = time;
    p_k_minus = params;
    P_k_minus = P_k;
    w_old = w;

end