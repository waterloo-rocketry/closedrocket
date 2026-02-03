function [params, K] = param_estimator_kf(time, w, d, c)

    % wrote out entire KF math by hand to learn, will consolidate onto one
    % line as we move forward.

    % d = encoder angle (rad)
    % c = rho * area * arm / inertia
    % dw/dt = 2 * c * Cl * d + C0 * c = phi_k' * p
    % Cl = canard lift coeff
    % C0 = rocket induced angular acceleration / (rho * area * arm)

    persistent t p_k_minus1 P_k_minus1 w_old
    if isempty(t)
        t = -0.01;
    end
    if isempty(p_k_minus1)
        p_k_minus1 = [2; 0];
    end
    if isempty(P_k_minus1)
        P_k_minus1 = diag([1, 1]);
    end
    if isempty(w_old)
        w_old = w;
    end
    
    % KF constants
    dt = time - t;
    y_k = (w - w_old)/dt; % measurement
    phi_k = c * [2 * d; 1];
    F = eye(2); 
    Q = zeros(2, 2); % used for forgetting factor
    R = 0.1;   % scalar (rad/s^2)^2

    % prediction
    p_k_minus = F * p_k_minus1;
    P_k_minus = F * P_k_minus1 * F' + Q;
    y_k_hat = phi_k' * p_k_minus;
    
    % prediction error
    nu_k = y_k - y_k_hat;
    S_k = phi_k' * P_k_minus * phi_k + R;

    % kalman gain
    K_k = P_k_minus * phi_k / S_k; % bless matlab i can divide matrices
    p_k_hat = p_k_minus + K_k * nu_k;
    P_k = (eye(2) - K_k * phi_k') * P_k_minus;
    
    % update for next cycle
    t = time;
    p_k_minus1 = p_k_hat;
    P_k_minus1 = P_k;
    w_old = w;

    params = p_k_hat;
    K = K_k;

end