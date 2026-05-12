function [x_new, P_new] = ekf_correct(model_measurement, model_jacobian, x, P, y, b, R)
    % Computes EKF correction step for other sensor data.
    % Input function handles: measurement model, model jacobian; 
    % Input variables: old state x, old covariance P, measurement y, sensor bias b;
    % Input parameters: sensor noise weighting R; 
    % Outputs: new state x, new covariance P

    %% Correction
    % computes a-posteriori state and covariance estimates
    % Uses discrete-time measurement model and analytical Jacobian

    %%% compute expected measurement and difference to measured values
    y_expected = model_measurement(x,b);
    innovation = y - y_expected;

    %%% compute Jacobian: H = dh/dx
    % H = jacobian(@model_measurement, 0, x, b); 
    H = model_jacobian(x, b);

    %%% compute Kalman gain (and helper matrices)
    L = H * P * H' + R;
    K = P * H' * inv(L);
    E = eye(length(x)) - K * H;
    
    %%% correct covariance estimate
    P_new = E * P * E' + K * R * K'; % joseph stabilized

    %%% correct state estimate
    x_new = x + K * innovation;
    x_new(1:4) = quaternion_norm(x_new(1:4)); % norm quaternions
end