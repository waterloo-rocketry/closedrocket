function [x_new, P_new] = ekf_correct_gyro(x, P, y, b, R)
    % Computes EKF correction step.
    % Inputs: estimates x, P; measurement y; sensor bias b;
    % Input parameters: weighting R; 
    % Outputs: new estimates x, P

    %% Correction
    % computes a-posteriori state and covariance estimates
    % Uses discrete-time measurement model and analytical Jacobian

    %%% compute expected measurement and difference to measured values
    % y_expected = meas_gyro(0,x,b); % hardcoded this
    % decompose state vector: [q(4); w(3); v(3); alt]
    w = x(5:7); 
    % decompose bias matrix: [b_A(3,i); b_W(3, i); M_E(3, i); b_P(1, i)]
    b_W = b(4:6);
    y_expected = w + b_W; 
    innovation = y - y_expected;

    %%% compute Jacobian: H = dh/dx
    % H = meas_gyro_jacobian(0, x, b); % hardcoded this
    H = zeros(3, 11);
    H(:, 5:7) = eye(3);

    %%% compute Kalman gain (and helper matrices)
    L = P(5:7, 5:7) + R;
    K = P(:, 5:7) * inv(L);
    E = eye(11) - [zeros(11,4), K, zeros(11,4)];
    
    %%% correct covariance estimate
    P_corr = E * P * E' + K * R * K'; % joseph stabilized

    %%% correct state estimate
    x_corr = x + K * innovation;
    x_corr(1:4) = x_corr(1:4) / norm(x_corr(1:4)); % norm quaternions

    %%% return a-posteriori estimates
    x_new = x_corr; P_new = P_corr;
end