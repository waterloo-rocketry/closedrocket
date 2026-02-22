function [y] = meas_mag(x, M_E)
    % Computes measurement prediction using current state and sensor biases

    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt; Cl; delta]
    q = x(1:4);

    %% magnetic field model
    S = quaternion_rotmatrix(q);
    M = S * M_E; % Earth magnetic field in body frame
    % TODO: add iron corrections

    %% measurement prediction
    y = [M];
end