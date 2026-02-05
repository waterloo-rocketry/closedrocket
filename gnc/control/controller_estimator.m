function [C_l_delta, C_l_0] = controller_estimator(time, w, delta, pdyn_params)
    % estimates the canard aerodynamic coefficients from canard angle, roll rates, air data
    % coeffs : canard coefficients C_l_delta and C_l_0
    % time : current time stamp (s)
    % w : angular rate measurement (rad/s)
    % d : canard angle measurement or command (rad)
    % pdyn_params : dynamic pressure * constant parameters (pressure * area * arm / inertia)
 
    % dw/dt = c * Cl * d + C0 * c = phi_k' * p
    % C_l_delta = canard lift coeff
    % C_l_0 = rocket induced angular acceleration / (rho * area * arm)

    %% tuning parameters
    Q = diag([1e-5, 1e-9]);

    %% initialize
    persistent t c P_minus w_old d_old w_dot_old
    if isempty(t)
        t = -0.01; % for /(time - t)
    end
    if isempty(c)
        c = [2; 0]; % initial coefficient guess
    end
    if isempty(P_minus)
        P_minus = Q; 
    end
    if isempty(w_old)
        w_old = w;
    end
    if isempty(d_old)
        d_old = 0;
    end
    if isempty(w_dot_old)
       w_dot_old = 0;
    end

    %% lowpass command and measurement
    alpha = 0.25;
    delta = (1-alpha) * d_old + alpha * delta;
    w_dot = (1-alpha) * w_dot_old + alpha * (w - w_old) / (time - t);

    %% Kalman filter
    r = pdyn_params * [delta; 1]; % regression 
    P = P_minus + Q; % covariance prediction
    K = P * r / (r' * P * r + 1); % correction gain. the stuff inside in brackets is just a scalar so you can just divide
    coeffs = c + K * (w_dot - r' * c); % coefficient correction
    P_plus = (eye(2) - K * r') * P * (eye(2) - K * r')' + K * 1* K'; % covariance correction. Joseph form for numerical stability
    
    %% update for next cycle
    t = time;
    c = coeffs;
    P_minus = P_plus;
    w_old = w;
    d_old = delta;
    w_dot_old = w_dot; 
    C_l_delta = coeffs(1);
    C_l_0 = coeffs(2);
end