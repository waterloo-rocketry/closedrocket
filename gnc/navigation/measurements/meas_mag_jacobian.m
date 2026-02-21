function [J] = meas_mag_jacobian(x, bias)
    % Computes measurement prediction using current state and sensor biases

    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt]
    q = x(1:4);

    % decompose bias matrix: [b_A(3,i); b_W(3, i); M_E(3, i); b_P(1, i)]
    M_E = bias(7:9);

    %% Initialize
    % Jacobian is of size: length(measurement) rows, length(x) columns
    % fill with zeros at first
    J = zeros(3, length(x));

    %% magnetic field model
    % S = quaternion_rotmatrix(q);
    % M = S * M_E; % Earth magnetic field in body frame
    % TODO: add iron corrections
    M_q = quaternion_rotate_jacobian(q, M_E);

    %% measurement prediction
    J(:, 1:4) = M_q;
end