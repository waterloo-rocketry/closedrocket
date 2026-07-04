function ctrl_mem = controller_estimator(dt_ctrl, w, delta, pdyn_params, ctrl_mem)
    %#codegen
    % estimates the canard aerodynamic coefficients from canard angle, roll rates, air data
    % w : (rad/s) angular rate measurement
    % delta : (rad) canard angle measurement or command
    % pdyn_params : (pressure * area * arm / inertia) dynamic pressure * constant parameters
    % ctrl_mem : estimator state, including coeffs, covariance, and filtered previous values

    % dw/dt = c * Cl * d + C0 * c = phi_k' * p
    % C_l_delta = canard lift coeff
    % C_l_0 = (angular acceleration / (rho * area * arm)) rocket induced coefficient

    %% tuning parameters
    %%% covariance
    Q = diag([5e-5, 5e-9]);

    %%% small angle cutoff
    if abs(delta) < 0.005 % prevents high noise density for small delta from affecting estimate
        delta = 0;  % probably should make more rigorous
    end

    %%% lowpass
    tau = 0.25; % (1/s) time constant

    %% lowpass command and measurement
    delta_lp = (1 - tau) * ctrl_mem.delta_lp + tau * delta;
    w_dot_lp = (1 - tau) * ctrl_mem.w_dot_lp + tau * (w - ctrl_mem.w) / dt_ctrl;

    %% Kalman filter
    r = pdyn_params * [delta_lp; 1]; % regression
    P = ctrl_mem.P + Q; % covariance prediction
    K = P * r / (r' * P * r + 1); % correction gain. the stuff inside in brackets is just a scalar so you can just divide
    coeffs = ctrl_mem.coeffs + K * (w_dot_lp - r' * ctrl_mem.coeffs); % coefficient correction
    P_plus = (eye(2) - K * r') * P * (eye(2) - K * r')' + K * 1* K'; % covariance correction. Joseph form for numerical stability

    %% update for next cycle
    ctrl_mem.coeffs = coeffs;
    ctrl_mem.P = P_plus;
    ctrl_mem.w = w;
    ctrl_mem.delta_lp = delta_lp;
    ctrl_mem.w_dot_lp = w_dot_lp;
end
