function [y] = meas_gyro(t, x, bias)
    % Computes gyroscope measurement prediction using current state and sensor biases

    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt]
    w = x(5:7); 

    % decompose bias matrix: [b_A(3,i); b_W(3, i); M_E(3, i); b_P(1, i)]
    b_W = bias(4:6);

    %% rates
    W = w + b_W; 

    %% measurement prediction
    y = [W];
end