function [y] = meas_mag(t, x, bias)
    % Computes measurement prediction using current state and sensor biases

    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt; Cl; delta]
    q = x(1:4);

    % decompose bias matrix: [b_A(3,i); b_W(3, i); M_E(3, i); b_P(1, i)]
    M_E = bias(7:9);

    %% magnetic field model
    S = quaternion_rotmatrix(q);
    M = S * M_E; % Earth magnetic field in body frame
    % TODO: add iron corrections

    %% measurement prediction
    y = [M];
end