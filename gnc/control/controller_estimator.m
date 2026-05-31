function [coeffs_ret, w_old_ret, P_minus_ret, d_old_ret, w_dot_old_ret] = controller_estimator(dt_ctrl, w, delta, pdyn_params, w_old, coeffs, P_minus, d_old, w_dot_old)
    %#codegen
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
    %%% covariance
    Q = diag([1e-5, 1e-9]);
    
    %%% small angle cutoff
    if abs(delta) < 0.005 % prevents high noise density for small delta from affecting estimate
        delta = 0;  % probably should make more rigorous
    end

    %%% lowpass
    tau = 0.25; % time constant

    %% lowpass command and measurement
    delta = (1 - tau) * d_old + tau * delta;
    w_dot = (1 - tau) * w_dot_old + tau * (w - w_old) / dt_ctrl;

    %% Kalman filter
    r = pdyn_params * [delta; 1]; % regression 
    P = P_minus + Q; % covariance prediction
    K = P * r / (r' * P * r + 1); % correction gain. the stuff inside in brackets is just a scalar so you can just divide
    coeffs = coeffs + K * (w_dot - r' * coeffs); % coefficient correction
    P_plus = (eye(2) - K * r') * P * (eye(2) - K * r')' + K * 1* K'; % covariance correction. Joseph form for numerical stability
    
    %% update for next cycle
    coeffs_ret = coeffs;
    P_minus_ret = P_plus;
    w_old_ret = w;
    d_old_ret = delta;
    w_dot_old_ret = w_dot; 
end