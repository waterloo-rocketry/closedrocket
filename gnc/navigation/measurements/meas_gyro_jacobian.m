function [J] = meas_gyro_jacobian(t, x, bias)
    % Computes gradient of gyroscope measurement prediction using current state and sensor biases

    %% Initialize
    % Jacobian is of size: length(measurement) rows, length(x) columns
    % fill with zeros at first
    J = zeros(3, length(x));

    %% rates
    % as W = w + b_W, Jacobian is unity
    W_w = eye(3);

    J(:, 5:7) = W_w;
end