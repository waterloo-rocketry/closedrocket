function [x_new, P_new] = ekf_update(dt, x, P, a_meas, w_meas, Q, R)
    % Computes EKF prediction+correction step for IMU data.
    % Input variables: time step dt, old state x, old covariance P; acceleration a, angular rate w; 
    % Input parameters: dynamics weighting Q, gyroscope weighting R;  
    % Outputs: new state x, new covariance P

    %% Prediction
    % computes a-priori state and covariance estimates
    % Uses discrete-time dynamics and analytical Jacobian
    
    %%% discrete dynamics update
    x_pred = dynamics(dt, x, a_meas);

    %%% discrete Jacobian: F = df/dx
    F = dynamics_jacobian(dt, x);

    %%% discrete covariance
    P_pred = F * P * F' + Q;


    %% Correction
    % computes a-posteriori state and covariance estimates
    % Uses discrete-time measurement model and analytical Jacobian

    %%% compute expected measurement and difference to measured values
    w = x_pred(5:7); % hardcoded measurement model, state vector: [q(4); w(3); v(3); alt]
    innovation = w_meas - w;

    %%% compute Jacobian: H = dh/dx
    % H = zeros(3, 11); % is hardcoded
    % H(:, 5:7) = eye(3); % hardcoded measurement jacobian

    %%% compute Kalman gain (and helper matrices)
    L = P_pred(5:7, 5:7) + R; % hardcoded H*P*H' 
    K = P_pred(:, 5:7) * inv(L); % hardcoded P*H
    E = eye(11) - [zeros(11,4), K, zeros(11,4)]; % hardcoded K*H
    
    %%% correct covariance estimate
    P_new = E * P_pred * E' + K * R * K'; % joseph stabilized

    %%% correct state estimate
    x_new = x_pred + K * innovation;
    x_new(1:4) = quaternion_norm(x_new(1:4)); % norm quaternions

end